"""Task 1: complex-argument edge sum vs float64 quadrature.

Sweeps the ramp strength s = k |A_t| r_max across the physically
reachable range (bounded by the beam power threshold: a beam attenuated
by more than -ln(P_th) in amplitude is dead; default P_th = 0.005 gives
s_max ~ 5.3) crossed with the (Re p, Re q) cases from
fraunhofer_limit_check.py. For each case, reports the relative error of:

  - f64 plain edge sum       (discretisation-free check of the identity)
  - f32 plain edge sum       (float32 rounding, no difference trick)
  - f32 difference-form sum  (what diff2.rs will actually compute)

against the dense float64 quadrature reference. Also validates the
effective cross-section identity D_eff = A(i mu, i nu) (pure imaginary
argument, real positive result).
"""

import numpy as np

from lossy import hexagon, polygon_ft, polygon_moments, reference_quadrature

RAMP_DIRECTION = np.array([0.6, 0.8])  # unit vector for (mu, nu)
RAMP_SCALES = [0.0, 0.25, 1.0, 2.65, 5.3, 8.0]  # s = k |A_t| r_max

PQ_CASES = [
    ("generic", 1.3, 0.7),
    ("oscillatory", 10.0, -6.0),
    ("q zero", 1.3, 0.0),
    ("both small", 5e-4, 5e-4),
    ("backscatter pole", 2e-5, -8e-8),
    ("both zero", 0.0, 0.0),
]


def main():
    verts = hexagon()
    area, mx, my, r_max = polygon_moments(verts)
    print(f"hexagon: area = {area:.6f}, r_max = {r_max:.4f}")
    print(f"ramp scales s = k|A_t| r_max: {RAMP_SCALES}")
    print(f"(s = 2.65 is half of -ln(0.005); s = 5.3 is the full power-threshold bound)\n")

    hdr = (f"{'case':<18} {'s':>5} {'|A_ref|':>12} "
           f"{'f64 err':>10} {'f32 plain':>10} {'f32 diff':>10}")
    print(hdr)
    print("-" * len(hdr))

    worst = {"f64": 0.0, "f32p": 0.0, "f32d": 0.0}
    for name, pr, qr in PQ_CASES:
        for s in RAMP_SCALES:
            mu, nu = RAMP_DIRECTION * (s / r_max)
            p = pr + 1j * mu
            q = qr + 1j * nu
            ref = reference_quadrature(verts, p, q)
            e64 = abs(polygon_ft(verts, p, q, f32=False, difference_form=False) - ref) / abs(ref)
            e32p = abs(polygon_ft(verts, p, q, f32=True, difference_form=False) - ref) / abs(ref)
            e32d = abs(polygon_ft(verts, p, q, f32=True, difference_form=True) - ref) / abs(ref)
            worst["f64"] = max(worst["f64"], e64)
            worst["f32p"] = max(worst["f32p"], e32p)
            worst["f32d"] = max(worst["f32d"], e32d)
            print(f"{name:<18} {s:>5.2f} {abs(ref):>12.5e} "
                  f"{e64:>10.2e} {e32p:>10.2e} {e32d:>10.2e}")
        print()

    print("worst-case relative errors:")
    print(f"  f64 plain:      {worst['f64']:.2e}   (identity + quadrature floor)")
    print(f"  f32 plain:      {worst['f32p']:.2e}")
    print(f"  f32 difference: {worst['f32d']:.2e}   (the diff2.rs form)")
    print("note: rows where the f64 and f32 errors coincide are limited by the")
    print("quadrature reference (O(h^2) midpoint rule), not by the edge sum; the")
    print("oscillatory case error falls 4x per doubling of the quadrature grid.")

    # D_eff identity: pure imaginary argument must give the real, positive
    # ramp-weighted area, equal to direct quadrature of the real integrand.
    print("\nD_eff = A(i mu, i nu) identity (pure imaginary argument):")
    print(f"{'s':>5} {'A(i mu, i nu)':>16} {'quadrature':>16} {'rel err':>10} {'Im/Re':>10}")
    for s in RAMP_SCALES:
        mu, nu = RAMP_DIRECTION * (s / r_max)
        val = polygon_ft(verts, 1j * mu, 1j * nu, f32=True, difference_form=True)
        ref = reference_quadrature(verts, 1j * mu, 1j * nu)
        rel = abs(val - ref) / abs(ref)
        imre = abs(val.imag) / abs(val.real) if val.real != 0 else float("inf")
        print(f"{s:>5.2f} {val.real:>16.6f} {ref.real:>16.6f} {rel:>10.2e} {imre:>10.2e}")


if __name__ == "__main__":
    main()
