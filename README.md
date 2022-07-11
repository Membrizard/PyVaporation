# PyVaporation

This solution is designed specifically to assist Researchers in the field of Pervaporation membranes development.
By means of the proposed instrument one can easily model a performance of a particular membrane with known permeance (Pi - GPU, SI, kg/(m2 * h * kPa)) and apparent energy of transport activiation (Ea - J/mol) values for each component of a considered binary mixture, if the transport is considered Ideal (Permeances are not dependent on the mixture composition)

Or, given that the diffusion curve set of a non-ideal process is measured one can model the non-ideal process in isothermal or non-isothermal (adiabatic) mode.
Non-isothermal modelling for both type of processes supports self-cooling mode, or temperature program mode.

# Following mixtures are Currently built in the solution:
(Current version supports only binary mixtures)


* H2O/MeOH
* H2O/EtOH
* H2O/IPOH
* H2O/Acetic acid
* MeOH/toluene
* MeOH/Methyl-tert-butyl ether
* MeOH/Dimethylcarbonate
* EtOH/Ethyl-tert-butyl ether


# Assumptions and applicability

* The activity coefficients of the binary mixture are calculated by means of NRTL model;
* Saturated vapour pressure could be assessed using Antoine or Frost equations
* Vaporisation/Condensation heat  values are calculated using Clapeyron-Clausius equation
* Specific heat capcities are calculated using polynomial approximation
* The ideal modelling process is applicable only for the modelling of processes, where permeance values does not depend significantly on mixture composition
* The non-ideal modelling is performed only based on the basis of specified diffusion curves (Fluxes/Permeances of each component as a function of first component concentration in feed)
* Non Ideal modelling supports non-linear dependencies of permeances and activation energies on feed composition 
* Non-Isothermal processes support pre-defined temperature program (feed temperature as a function of process time may be specified for process modelling)

# Installation

Rquirements:

python 3.8 or higher

To install:
```
pip install pyvaporation
```

# Code examples

# Hints for general usage

* Pre-configured default membranes are located in 
```
   ./tests/default_membranes
``` 
* To run automated tests for all the modules: 
```
   python -m pytest -sv tests/
```

