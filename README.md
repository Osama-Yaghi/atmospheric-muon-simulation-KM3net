# atmospheric-muon-simulation-KM3net

This project simulates the flux of atmospheric muons underwater using a Monte Carlo method, tailored for the geometry and detection environment of the KM3NeT neutrino telescope.


## 📌 Overview

Atmospheric muons are a key background and calibration source for deep-sea neutrino telescopes. This simulation models the flux of muons arriving at the detector site by:
- Sampling muon energies and angles using a physics-based probability distribution.
- Propagating muons through seawater while accounting for continuous energy loss.
- Calculating intersections with a cylindrical detector volume.
- Filtering muons that reach the detector with energy above a threshold.

The result is a realistic sample of muon events at the detector’s boundary, with spatial, directional, and energy distributions consistent with physical expectations.


## 🌊 About KM3NeT

**KM3NeT** (Cubic Kilometre Neutrino Telescope) is a multi-site, deep-sea research infrastructure hosting neutrino telescopes in the Mediterranean Sea. Its primary goal is to detect neutrinos from astrophysical sources and study fundamental neutrino properties.

Key features:
- **Detector Environment:** Located at depths of 2500–5000 meters, minimizing background from atmospheric particles.
- **Detector Geometry:** A 3D array of photomultiplier tubes (PMTs) housed in digital optical modules (DOMs), arranged along vertical detection units.
- **Scientific Goals:**
  - Identify cosmic neutrino sources.
  - Investigate the neutrino mass hierarchy.
  - Study atmospheric neutrino interactions and oscillations.

Understanding the **flux of atmospheric muons**, a dominant background source, is essential for:
- Background rejection and calibration.
- Optimizing detector design and reconstruction algorithms.
- Interpreting low-energy neutrino events that can be mimicked by muon-induced signatures.

This simulation models those muons in realistic geometry, helping assess and tune the detector response.

## 🧠 Physics Modeled

- Muon flux based on the Gaisser parameterization.
- Solid angle correction via Jacobian transformation.
- Continuous energy loss using tabulated data (`loss.txt`) and interpolation.
- Propagation through a cylindrical geometry with geometric rejection sampling.
- Thresholding by muon range in water to ensure realism.

## 🔧 Features

- Highly parallelized using `joblib` to maximize performance.
- 3D visualization of muon hits on the detector.
- Histogram plots of muon energy and vertical momentum.
- Output CSV files for statistical post-processing and visualization.

## 📂 Repository Structure

```bash
.
├── Monte-carlo-simulation-generate_v1.0.py  # Main simulation script
├── loss.txt                                 # Energy loss table for muons in seawater
├── result_<run_id>.csv                      # Output: simulated muon hits
├── efficiency_<run_id>.csv                  # Output: efficiency metrics (optional)
├── hists1.png                               # Histogram of energy and direction
├── 3d1.png                                  # 3D scatter plot of muon hits
├── README.md                                # This file
