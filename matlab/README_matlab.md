# MagnetostaticDesign

1.0     This is the Readme file for the Matlab magnet array and coil design toolbox. It includes:

- libraries for designing magnet arrays/EM coils, generating other files for simulation (Biot-Savart files; COMSOL model files), and generating volumetric STL models for importing into CAD drawing programs
- example scripts to produce the "Head-Optimized MRI" designs described in my thesis: McDaniel, Patrick C. "Computational design and fabrication of portable MRI systems", Chapter 7. (Ph.D Thesis), MIT (Sept 2020)
- CAD models for the "Head-Optimized MRI" mechanical components, as described in my thesis
- Documentation on how to use the code and post-process the generated STL models

2.0     Using this code requires several other libraries, linked below. These are all portable Matlab toolboxes that you should put in your Matlab search path.
    
- Gael Bringout's Coil Design Toolbox: https://github.com/gBringout/CoilDesign/
    - This in turn requires other toolboxes, as described on the github page
        - regu: http://www.imm.dtu.dk/~pcha/Regutools/
        - SphericalHarmonics: https://github.com/gBringout/SphericalHarmonics
        - YALMIP: http://users.isy.liu.se/johanl/yalmip/
- HarmonicY: https://github.com/jmontalt/harmonicY
- stlwrite (if using pre-R2018b Matlab): https://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-write-ascii-or-binary-stl-files
- stlread (if using pre-2018b Matlab): https://www.mathworks.com/matlabcentral/fileexchange/6678-stlread

2.1     The permanent magnet array design code includes scripts to interface with COMSOL via the "Matlab live link" capability. This code requires having Comsol, the Comsol EM toolbox, and Matlab live link installed. You also need to start Matlab using "comsol server matlab" (ie to run it with the Comsol-Matlab link), and not just run vanilla Matlab.

3.0     Notes on use:

- Some of the examples will precompute large (~1Gb) matrices and save them as *.mat files. These files will all be saved in the ~/precomp/ directory. Clear this directory after running the example scripts if you no longer need these files.

4.0     Quick start guide

- Download and link the libraries listed in (2.0)
- Follow the included slides README_ppt.pptx
	- They provide instructions on how to quickly run the examples in ~/ExampleScripts/

