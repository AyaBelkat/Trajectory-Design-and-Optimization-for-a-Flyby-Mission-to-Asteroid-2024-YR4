# Trajectory-Design-and-Optimization-for-a-Flyby-Mission-to-Asteroid-2024-YR4
MATLAB tools for interplanetary trajectory optimisation and flyby mission design for asteroid 2024 YR4. Includes MGADSM trajectory design, local/global optimisation, and orbital visualisation.

Project Abstract
This repository provides a high-fidelity trajectory design and optimisation framework for a flyby mission to Asteroid 2024 YR4. The mission is designed to intercept the asteroid before its peak Earth-impact risk window in December 2032. Two primary mission architectures are explored: an Earth-Resonant Gravity Assist and a Venus Gravity Assist, both utilising a Deep Space Manoeuvre (DSM) to minimise total propellant requirements.

Key Features:
-> MGADSM Engine: A modular engine for Multi-Gravity Assist trajectories with segmented Deep Space Manoeuvres.
-> Global-Local Optimisation: Uses Monotonic Basin Hopping (MBH) to navigate non-convex search spaces and find global Delta-V minima.
-> High-Fidelity Ephemeris: Fully integrated with the NASA SPICE (Mice) toolkit for precise J2000 planetary and asteroid state vectors.
-> B-Plane Targeting: Advanced flyby modelling using Rodrigues’ rotation formula for 3D trajectory steering.

Results Summary (check code version):

Parameter           |  Ballistic Transfer  | Earth-Resonant Assist+DSM  |  Venus Flyby Assist+DSM
---------------------------------------------------------------------------------------------------
Departure Date      |  17 Feb 2032         |  14 Oct 2028               |  27 May 2031
Arrival Date        |  17 Nov 2032         |  07 Nov 2032               |  13 Oct 2032
Total ΔV (km/s)     |  1.7491              |  3.807                     |  2.796
Arrival v∞​ (km/s)   |  7.446               |  6.998                     |  6.984
Arrival angle (deg) |  83.87               |  59.92                     |  60.00
Total TOF (days)    |  274                 |  1485                      |  505

Repository Structure:

