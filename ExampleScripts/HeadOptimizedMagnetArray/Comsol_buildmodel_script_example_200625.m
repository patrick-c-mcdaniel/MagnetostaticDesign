
%% Generates and simulates comsol model; specify simulation parameters here
%   Nround - number of discrete block sizes to "round" continuous solution
%               to
%   wRand  - width of the uniform random error distribution for Monte-Carlo
%               simulations of block size errors
%   Br     - magnet material remanent flux density
%   mur    - magnet material relative permeability

Nround = 26;
wRand  = 0e-3;
Br     = 1.42;
mur    = 1.05;


%% Run one iteration of the comsol simulation
%   NOTE: this function is set up to simulate the NON-ROUNDED block sizes;
%   you will need to comment out a line in this function in order to
%   simulate the rounded block sizes
Nmc = 1;

for imc = 1:Nmc
    Comsol_buildmodel_func_example_200625( Nround, wRand, Br, mur );
end
