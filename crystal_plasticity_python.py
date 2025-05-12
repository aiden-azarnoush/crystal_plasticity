# Crystal Plasticity Integration Tool

## Description

This Python tool implements explicit point integration for Crystal Plasticity using Busso (1990) flow rule as modified by Luscher et al. The program simulates crystal deformation processes by tracking slip systems in crystalline materials.

The tool features three different hardening models:
1. Asaro and Needleman
2. Bassani and Wu
3. Luscher et al.

Users can specify custom crystal orientations, with several preset options available. The code assumes that the total deformation gradient (F) is known for each step of the simulation.

## Features

- Implementation of three different hardening models for crystal plasticity
- User-definable crystal orientation (<100>, <110>, <310>, <580>, or custom)
- Calculation of stress-strain relationships
- Tracking of slip systems activation
- Visualization of results through matplotlib
- MATLAB version included for additional compatibility

## Requirements

For Python version:
- Python 2.x
- NumPy
- SciPy
- Matplotlib

For MATLAB version:
- MATLAB (compatible with most recent versions)

## Installation

1. Ensure Python/MATLAB and required dependencies are installed
2. Download the repository to your local machine
3. Make the Python script executable (Unix/Linux):
   chmod +x crystal_plasticity.py

## Usage

For Python version, run the script from the command line:
python crystal_plasticity.py

For MATLAB version, open MATLAB and run:
crystal_plasticity

You'll be prompted to:
1. Select a crystal orientation (or provide a custom one)
2. Choose which hardening model to use

The program will run the simulation and display the results as plots showing:
- Normalized stress vs. strain
- Sp (slip resistance) vs. time
- Gamma_dot (slip rate) vs. time
- Tau (resolved shear stress) vs. gamma (slip)

## Technical Information

The code implements crystal plasticity based on the following key components:
- Rotation matrices for crystal orientation
- Elasticity tensor calculation in the rotated frame
- Flow rule for slip system activation based on Busso (1990)
- Three different hardening models:
  - Asaro and Needleman hardening model
  - Bassani and Wu hardening model
  - Luscher et al. hardening model
- Numerical integration of plastic deformation

> [!NOTE]
> The simulation uses small time steps (default: dt = 1e-12) to accurately capture the deformation process over the specified total time.

## Output

The code generates four plots:
1. Normalized stress vs. strain
2. Slip resistance (Sp) vs. time for each slip system
3. Slip rate (gamma_dot) vs. time for each slip system
4. Resolved shear stress (tau) vs. slip (gamma) for each slip system

## Author

Aiden Azarnoush (mazarnou@purdue.edu)

---

*Dedicated to my mother, Simin Nematpour*
