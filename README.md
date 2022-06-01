# PVmodelling

This solution is designed specifically to assist Researchers in the field of Pervaporation membranes development.
By means of the proposed instrument one can easily model a performance of a particular membrane with known permeance (Pi) and apparent energy of activiation  of transport (Ea) values for each component of a considered binary mixture, if the transport is considered Ideal (Permeances are not dependant on the mixture composition)

Or, given that the diffusion curve of a non-ideal process is measured one can model the non ideal process at a given temperature.

# The Following mixtures are Currently built in the solution:
(Current version supports only binary mixtures)


* H2O/MeOH
* H2O/EtOH
* H2O/IPOH
* MeOH/toluene
* MeOH/Methyl-tert-butyl ether
* MeOH/Dimethylcarbonate
* EtOH/Ethyl-tert-butyl ether


# Assumptions and applicability

* The activity coefficients of the binary mixture are calculated by means of NRTL model;
* Saturated vapour pressure could be assesd using assessed using Antoine or Frost equations
* Vaporisation/Condensation heat is calculated using Clapeyron-Clausius equations
* Specific heat capcities are calculated using polynomial approximation
* The ideal modelling process is applicable only for the modelling of processes, where permeance values does not depend significantly on mixture composition
* The non-ideal modelling is performed only based on the specified diffusion curves (Permeances of each component as a function of first component concentration in feed)
* Non Ideal modelling supports non-linear dependencies of permeances and activation energies on feed composition 
