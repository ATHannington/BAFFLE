# BAFFLE (https://github.com/ATHannington/BAFFLE)
**This project has concluded. Further work is unlikely. Please be aware that this project was for learning purposes, and the software is not optimised, fully debugged, or complete.**
BAFFLE - BAsic Fortan FLuids Engine: A basic fixed cartesian grid, flux based, second order hydrodynamics code in Fortran90. Developed in collaboration with Dr B. Vandenbroucke, University of St Andrews.

BAFFLE is, as the name suggests, a rather basic fluids engine written in Fortran90. 

BAFFLE is a direct derivative of work by Bert Vandenbroucke ().Find more work by Bert Vandenbroucke at: (https://github.com/bwvdnbro)

BAFFLE was written by Andrew Hannington () as part of Summer project work with Prof. Ian Bonnell, St Andrews University ().  

BAFFLE was a collaboration of the author Andrew Hannington () with Bert Vandenbroucke and William Lucas (), under the supervision of Prof. Ian Bonnell, St Andrews University ().  

Copyright (C) 2017 Andrew Hannington ()
!----------------------------------------------------------------------------------------------------------------------------------

To run BAFFLE is rather simple, as it is designed to be as modular as possible. To run:

1. Download the relevant files from this repository (https://github.com/ATHannington/BAFFLE). You will have to download **AT LEAST ONE FILE FROM EACH FOLDER** and **ADAPT TO RELEVANT NUMBER OF DIMENSIONS IF NECESSARY**. The BAFFLE code provided is **NOT** exhaustive, and was written by the author within a limited time-frame, and with limited experience.
2. Setup Initial Conditions subroutine file to choose the setup you wish to model
3. Configure the relevant constants module file for desired resolution, accuracy (first or second order), run time, isothermal/adiabatic setup, and number of outputs and output location etc. .
4. Compile your files using a your choice of compiler (the author used Gfortran on Linux)
5. Run the executable file generated by your compiler, and perform data analysis as you desire.

For further information on the theory behind BAFFLE, and the software it is a derivative of, please see the work of Bert Vandenbroucke at (https://github.com/bwvdnbro), specifically his python_finite_volume_solver @:(https://github.com/bwvdnbro/python_finite_volume_solver)
