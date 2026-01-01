# Zones Implementation Plan

## Overview

Refactor the binning system from a single global binning scheme to a list of zones, where each zone has its own binning config, results, and parameters.

---

## Phase 1: Core Zone Data Structures

### 1.1 Define ZoneType enum (`src/zones.rs` - new file)
```rust
enum ZoneType {
    Full,      // 0-180 theta coverage, computes asymmetry, scat cross-section
    Forward,   // Forward scattering, computes ext cross (optical theorem)
    Backward,  // Backscatter, computes lidar ratio
    Custom,    // User-defined arbitrary range
}
```

### 1.2 Define Zone struct
```rust
struct Zone {
    label: Option<String>,        // Optional user label
    zone_type: ZoneType,
    scheme: Scheme,               // Reuse existing Simple/Interval/Custom
    bins: Vec<SolidAngleBin>,     // Generated from scheme
    field_2d: Vec<ScattResult2D>, // Results per bin
    field_1d: Option<Vec<ScattResult1D>>,
    params: Params,               // Zone-specific parameters
}
```

### 1.3 Define ZoneConfig for deserialization
```rust
struct ZoneConfig {
    label: Option<String>,
    scheme: Scheme,  // Existing binning scheme config
}
```
- ZoneType is inferred: if theta range is 0-180, it's Full; otherwise Custom
- Forward/Backward zones are added programmatically, not from config

### 1.4 Update Results struct
```rust
struct Results {
    zones: Vec<Zone>,
    powers: Powers,  // Stays global
}
```

---

## Phase 2: Configuration & Parsing

### 2.1 Update TOML structure
**default.toml / local.toml:**
```toml
[[zones]]
label = "main"  # optional
[zones.scheme.Interval]
thetas = [0, 180]
theta_spacings = [2]
phis = [0, 360]
phi_spacings = [4]

[[zones]]
label = "high_res_forward"
[zones.scheme.Interval]
thetas = [0, 5]
theta_spacings = [0.1]
phis = [0, 360]
phi_spacings = [2]
```

### 2.2 Update Settings struct (`src/settings.rs`)
- Replace `binning: BinningScheme` with `zones: Vec<ZoneConfig>`
- Add methods for runtime zone manipulation

### 2.3 Update CLI parsing (`src/settings/cli.rs`)
- Adapt existing `--simple`, `--interval`, `--custom` to create zone configs
- Multiple zone specifications allowed
- Print `info!` as each zone is processed

### 2.4 Zone processing logic
- For each ZoneConfig:
  - Infer ZoneType from theta range (0-180 = Full, else Custom)
  - Log: `info!("Processing zone '{}': {:?}", label, zone_type)`
  - Generate bins using existing scheme logic
- Forward/Backward zones added programmatically (not from config)

---

## Phase 3: Results & Zone Integration

### 3.1 Refactor Results initialization
- `Results::new(zones_config)` creates zone list
- Each zone initializes its own `field_2d`, empty `params`

### 3.2 Zone-specific parameter computation
- Full zone: asymmetry, scattering cross-section
- Forward zone: extinction cross-section (optical theorem)
- Backward zone: lidar ratio, backscatter cross-section
- Custom zone: TBD based on coverage

### 3.3 Update `Params` struct if needed
- May need zone-aware getters/setters
- Or each Zone just owns its own Params instance (simpler)

---

## Phase 4: Problem Solving (REQUIRES GUIDANCE)

### 4.1 Review `solve_far()` in `problem.rs`
**CHECKPOINT: Pause here for user guidance before any edits**

Current flow:
1. Process beams through ray tracing
2. Map beams to bins (n2f_go or aperture diffraction)
3. Accumulate amplitudes per bin
4. Convert to Mueller matrices
5. Special handling for forward/backward bins

New flow considerations:
- Each zone processed independently or in parallel?
- Beam mapping needs to check all zones' bins
- Forward/Backward become regular zones

### 4.2 `solve_far_queue()` 
**DO NOT EDIT WITHOUT APPROVAL**

### 4.3 Update beam-to-bin mapping
- `n2f_go()` and diffraction need zone awareness
- Return `Vec<(zone_idx, bin_idx, Ampl)>` instead of `Vec<(bin_idx, Ampl)>`?
- Or process each zone separately?

---

## Phase 5: Convergence Integration

### 5.1 Design decision point
**Option A**: Boolean flags on zones to enable/disable processing
- Simpler convergence logic
- More complexity in problem.rs

**Option B**: Pass subset of zones into problem solve
- Cleaner problem.rs (processes whatever zones it's given)
- Convergence manages which zones are active
- **Preferred approach** - keeps inner pieces simple

### 5.2 Update Convergence solver
- Zones can be dynamically added/removed between iterations
- Convergence targets are zone-specific: `(zone_label, Param, threshold)`

---

## Phase 6: Output Format

### 6.1 Consolidate output to single JSON
**New format: `results.json`**
```json
{
  "powers": {
    "absorbed": 0.123,
    "reflected": 0.456,
    ...
  },
  "zones": [
    {
      "label": "main",
      "type": "Full",
      "params": {
        "ScatCross": { "Total": 1.23, "Beam": 0.8, "Ext": 0.43 },
        "Asymmetry": { "Total": 0.76, ... }
      },
      "mueller_file": "mueller_main.dat"
    },
    {
      "label": "forward",
      "type": "Forward",
      "params": {
        "ExtCross": { "OpticalTheorem": 2.34 }
      },
      "mueller_file": "mueller_forward.dat"
    }
  ]
}
```

### 6.2 Mueller output per zone
- Each zone writes its own `mueller_{label}.dat` or similar
- Naming scheme: use label if provided, else `zone_0`, `zone_1`, etc.

### 6.3 Deprecate old files
- `results.dat`, `powers.json`, `params.json` → single `results.json`
- Update any scripts/tests that depend on old format

---

## Phase 7: Python Bindings

### 7.1 Update PyO3 bindings for Results
```python
results.zones  # List of Zone objects
results.zones[0].params  # Zone-specific params
results.zones[0].mueller  # Zone Mueller data
results.powers  # Global powers
```

### 7.2 Zone access methods
- `results.get_zone("main")` - by label
- `results.get_zone_by_type(ZoneType.Full)` - returns list
- `results.full_zone` - convenience for first Full zone (common case)

---

## Implementation Order

1. **Phase 1**: Core structs (Zone, ZoneType, ZoneConfig) - foundation
2. **Phase 2**: Config parsing - can test zone loading
3. **Phase 3**: Results refactor - integrate zones into results
4. **Phase 4**: Problem solving - **CHECKPOINT FOR GUIDANCE**
5. **Phase 5**: Convergence - after problem.rs is stable
6. **Phase 6**: Output format - can be done in parallel with Phase 5
7. **Phase 7**: Python bindings - final polish

---

## Open Questions (to resolve during implementation)

1. How should overlapping zones handle the same beam? Process beam for each zone independently?
2. Should Forward/Backward zones be auto-added always, or configurable?
3. 1D integration (`field_1d`) - does it make sense for Custom zones with partial coverage?
4. Convergence: can we converge on multiple zones simultaneously with different targets?
5. Memory: for large zone counts, should we lazy-load Mueller data?

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/zones.rs` | NEW - Zone, ZoneType, ZoneConfig structs |
| `src/bins.rs` | Minor - may move some logic to zones.rs |
| `src/settings.rs` | Replace binning with zones |
| `src/settings/cli.rs` | Update CLI parsing for zones |
| `src/settings/loading.rs` | Update TOML loading for zones |
| `src/result.rs` | Major refactor - zones container |
| `src/problem.rs` | **CAREFUL** - solve_far changes |
| `src/convergence.rs` | Zone-aware convergence |
| `src/multiproblem.rs` | Update initialization |
| `src/output.rs` | New JSON output format |
| `config/default.toml` | New zones format |
| Python bindings | Zone access methods |
