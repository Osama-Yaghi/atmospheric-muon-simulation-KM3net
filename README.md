# Atmospheric Muon Flux Simulation for KM3NeT

This project simulates the flux of atmospheric muons underwater using a Monte Carlo method, tailored for the geometry and detection environment of the KM3NeT neutrino telescope.
Written in 11/2020.
## 🌊 About KM3NeT

**KM3NeT** (Cubic Kilometre Neutrino Telescope) is a deep-sea neutrino observatory being built in the Mediterranean Sea. It is designed to detect high-energy neutrinos from astrophysical sources and study fundamental properties of neutrinos.

Key features:
- Detector depth of 2500–5000 meters reduces background noise.
- Photomultiplier tubes arranged in large 3D arrays across vertical detection lines.
- Dual scientific goals: cosmic neutrino detection (ARCA) and neutrino oscillation physics (ORCA).

Atmospheric muons represent a dominant background signal and a calibration tool. Accurately simulating them is critical for:
- Background rejection algorithms.
- Trigger optimization.
- Detector calibration and energy reconstruction studies.

## 📌 Overview

This code simulates the atmospheric muon flux arriving at an underwater detector. It:
- Samples muon energy and angle distributions based on physics-driven probability densities.
- Uses geometric sampling to calculate muon paths toward a cylindrical detector.
- Applies continuous energy loss models through seawater using real tabulated stopping powers.
- Filters muons that survive to the detector boundary with sufficient energy.

## 🧠 Physics Modeled

- Gaisser-like parameterization of atmospheric muon flux.
- Angular correction using Jacobian solid angle factor.
- Continuous energy loss in water, interpolated from a tabulated file (`loss.txt`).
- Propagation and intersection with a cylindrical detector volume.
- Muon survival filtering based on energy threshold at arrival.

## 🔧 Features

- Parallelized with `joblib` for faster simulation across CPU cores.
- Generates histograms of muon energies and directional cosines.
- Produces 3D scatter plots of detected muon hit locations.
- Outputs all data in clean CSV format for analysis or plotting.

## 🚀 How to Run

### Requirements

- Python 3.x
- numpy
- pandas
- matplotlib
- seaborn
- scipy
- joblib

Install requirements using pip or your environment manager of choice.

### Execution
-Modify the input file in input/input.txt
-Run the main script directly:
python Monte-carlo-simulation-generate_v1.0.py

This will launch the simulation, generate events, and produce output files and plots in output/.

## 📊 Output Description

The main output file (`result_<run_id>.csv`) contains one line per detected muon event. Columns:

- `x, y, z`: Position where the muon intersects the detector volume.
- `vx, vy, vz`: Directional components of the muon.
- `E`: Energy of the muon at the detector.
- `num_calls`: Number of random sampling iterations used to generate the muon (for efficiency diagnostics).

Two plots are also saved:
- `hists1.png`: Histograms of energy and vertical direction component.
- `3d1.png`: 3D scatter plot of muon entry points on the detector.

Additional optional output includes efficiency data in `efficiency_<run_id>.csv`.

## ⚠ Notes & Assumptions

- Simulation assumes maximum muon energy of 100 GeV.
- Only downward-going muons (zenith < 85°) are considered.
- `loss.txt` must be present and correctly formatted (columns: energy, range, stopping power, etc.).
- Energy loss model uses continuous approximation; stochastic effects are not included.

## ✍ Author

**Osama Yaghi**
