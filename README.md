# Trajectory-Design-and-Optimization-for-a-Flyby-Mission-to-Asteroid-2024-YR4

![MATLAB](https://img.shields.io/badge/MATLAB-R2021b+-blue)
![SPICE](https://img.shields.io/badge/NAIF-SPICE-orange)
![Status](https://img.shields.io/badge/status-active-success)

MATLAB tools for interplanetary trajectory optimisation and flyby mission design for asteroid 2024 YR4. Includes MGADSM trajectory design, local/global optimisation, and orbital visualisation.

## Project Abstract
This repository provides a high-fidelity trajectory design and optimisation framework for a flyby mission to Asteroid 2024 YR4. The mission is designed to intercept the asteroid before its peak Earth-impact risk window in December 2032. Two primary mission architectures are explored: an Earth-Resonant Gravity Assist and a Venus Gravity Assist, both utilising a Deep Space Manoeuvre (DSM) to minimise total propellant requirements.

## Intended Audience
This repository is intended for astrodynamics students, researchers, and mission analysts interested in interplanetary trajectory design and optimisation using gravity assists and deep space manoeuvres.

## Key Features:
* **MGADSM Engine:** A modular engine for Multi-Gravity Assist trajectories with segmented Deep Space Manoeuvres.
* **Global-Local Optimisation:** Uses Monotonic Basin Hopping (MBH) to navigate non-convex search spaces and find global Delta-V minima.
* **High-Fidelity Ephemeris:** Fully integrated with the NASA SPICE (Mice) toolkit for precise J2000 planetary and asteroid state vectors.
* **B-Plane Targeting:** Advanced flyby modelling using Rodrigues’ rotation formula for 3D trajectory steering.

## Results Summary:

| Parameter | Ballistic Transfer | Earth-Resonant Assist + DSM | Venus Flyby Assist + DSM |
| :--- | :--- | :--- | :--- |
| **Departure Date** | 17 Feb 2032 | 14 Oct 2028 | 27 May 2031 |
| **Arrival Date** | 17 Nov 2032 | 07 Nov 2032 | 13 Oct 2032 |
| **Total $\Delta V$ (km/s)** | 1.7491 | 3.807 | 2.796 |
| **Arrival $v_{\infty}$ (km/s)** | 7.446 | 6.998 | 6.984 |
| **Arrival Angle (deg)** | 83.87 | 59.92 | 60.00 |
| **Total TOF (days)** | 274 | 1485 | 505 |

> **Notes:**  
> - Results correspond to the current optimisation configuration and SPICE kernel versions; values may vary with different bounds or solver settings.  
> - The arrival angle quantifies whether the spacecraft approaches the sun-illuminated side of the asteroid. Values below **60°** were adopted as a criterion for favourable illumination conditions during the flyby.

### Mission Analysis
The results highlight a clear trade-off between **Time of Flight (TOF)** and **ΔV efficiency**. While the ballistic transfer provides the shortest TOF, the **Venus Gravity Assist** architecture offers more favourable arrival conditions on the sun-illuminated side of the asteroid. In addition, gravity assist trajectories yield lower arrival \( v_{\infty} \), which is critical for maximising observation time and scientific return during the flyby phase.

## Repository Structure

* `startup.m`: Initialises the MATLAB environment and maps all project subfolders. *This script must be executed before running any other code.*
* `/data`: Contains all the necessary small data files. BIG FILES TO DOWNLOAD: de430.bsp from: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
* `/core_solvers`: Universal physics engines, including the Lambert solver 
* `/mgadsm_engine`: Architecture-specific cost functions, nonlinear constraints, and GA calculations.
* `/ballistic-optimisation`:  ballistic trajectory optimisation cost function and nonlinear constraints
* `/plotting_utils`: Specialised tools for generating standardised trajectory visualisations for Earth-asteroid 2024 YR4 pork-chop results.
* `/trajectory_plots`: Contains plots of the optimisation results from ballistic and MGADSM optimisers.

## Getting Started

1. Dependencies: Download the Mice Toolkit and required planetary/asteroid kernels from the NASA NAIF website: https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
2. Clone the Repository:
> git clone https://github.com/AyaBelkat/Trajectory-Design-and-Optimization-for-a-Flyby-Mission-to-Asteroid-2024-YR4.git
3. Setup: Place the Mice files and kernels in your project root or update the paths in load_kernels.m
4. Initialise: Run startup.m to configure the environment and verify that the NAIF data is correctly loaded.

