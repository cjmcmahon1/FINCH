# FINCH
## Introduction
This is a MATLAB repository designed to simulate Fresnel Incoherent Correlation Holography (FINCH). The starting point for most of the code is a system with a bi-refringent lens.\
![image](Images/Presentation_Images/Figures/BRL_diagram.png)
Because of the self-interference structure of FINCH, 
## Code Structure
Relevant scripts for plotting simulation are in the top directory `./`. Most helper functions are in the subfolder `./MATLAB_functions`. Some scripts for processing actual data are in `./Data_Scripts/`, along with some data-specific helper functions. Scripts designed to test basic functionality of the Fresnel propagator, Fourier transform, etc., are in `./Test_Scripts`.
## Bench Setup
Images from the bench setup can be found in `./Images`. The actual setup uses a triangle interferometer, which produces the same effect as the simulated bi-refringent lens. The images taken are of a pinhole illuminated by 490nm incoherent (LED) light. Different pinhole sizes have been used, but 5 micron is most recent.
## References
Introduction to Fourier Optics - Joseph Goodman\
Computational Fourier Optics: A MATLAB Tutorial - David Voelz
