# MagnetostaticDesign

The MagnetostaticDesign toolbox allows for target-field design of permanent magnet arrays and quasi-static EM coils. This toolbox accompanies my PhD thesis:

   Patrick C McDaniel. "Computational design and fabrication of portable MRI systems". MIT dept of EECS (2020) [URL forthcoming]
   
This toolbox presently exists only in the Matlab implementation (found in /matlab/), the details of which can be found in the ~/matlab/README_matlab.md file. The python implementation will be found in the ~/python/ directory (eventually)
        
This toolbox allows for:
    - Field simulation of permanent magnet/electromagnet structures
    - Specification of arbitrary target volume+ magnet/coil geometry
    - Magnet/coil design optimization
    - Automatic generation/simulation of Comsol models (requires Comsol + Matlab/Comsol live link)
    - Conversion to STL model for subsequent CAD, processing, or manufacturing steps
