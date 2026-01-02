# API Reference

This page documents the Python API for GOAD.

## Simulation Classes

These are the main classes for running light scattering simulations.

### Problem

Single-orientation scattering simulation.

::: goad.Problem
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"
        - "!^py_"

### MultiProblem

Multi-orientation averaging simulation.

::: goad.MultiProblem
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"
        - "!^py_"

### Convergence

Adaptive convergence-based simulation.

::: goad.Convergence
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"
        - "!^py_"

---

## Configuration

Classes for configuring simulations.

### Settings

Main configuration object for simulations.

::: goad.Settings
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"
        - "!^set_"

### Orientation

Particle orientation configuration.

::: goad.Orientation
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

### Euler

Euler angle representation.

::: goad.Euler
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

### BinningScheme

Angular binning configuration.

::: goad.BinningScheme
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

### ZoneConfig

Zone configuration for multi-zone simulations.

::: goad.ZoneConfig
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"
        - "!^set_"

---

## Results

Classes for accessing simulation results.

### Results

Complete simulation results container.

::: goad.Results
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"
        - "!^get_"
        - "!^set_"

### Zone

A single angular zone with its results.

::: goad.Zone
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

### Zones

Collection of zones.

::: goad.Zones
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

---

## Geometry

Classes for defining particle geometry.

### Geom

Complete particle geometry (collection of shapes).

::: goad.Geom
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

### Shape

A single 3D surface mesh.

::: goad.Shape
    options:
      show_source: false
      members_order: source
      show_root_heading: true
      show_root_full_path: false
      filters:
        - "!^_"

---

## Enumerations

### Param

Available scattering parameters for convergence targets.

::: goad.Param
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false

### ZoneType

Types of angular zones.

::: goad.ZoneType
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false

### Mapping

Near-to-far field mapping methods.

::: goad.Mapping
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false

### EulerConvention

Euler angle conventions.

::: goad.EulerConvention
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false

### Scheme

Orientation sampling schemes.

::: goad.Scheme
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false

---

## Helper Functions

### create_uniform_orientation

::: goad.create_uniform_orientation
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false

### create_discrete_orientation

::: goad.create_discrete_orientation
    options:
      show_source: false
      show_root_heading: true
      show_root_full_path: false
