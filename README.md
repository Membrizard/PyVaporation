# PVmodelling

This solution is designed specifically to assist Researchers in the field of Pervaporation membranes development.
By means of the proposed instrument one can easily model a performance of a particular membrane with known permeance (Pi) and appearing energy of activiation  of transport (Ea) values for each component of a considered binary mixture, if the transport is considered Ideal (Permeances are not dependant on the mixture composition)

Or, given that the diffusion curve of a non-ideal process is measured one can model the non ideal process at a given temperature.

# The Following mixtures are Currently built in the solution:

binary:

* H2O/MeOH
* H2O/EtOH
* H2O/IPOH
* MeOH/toluene
* MeOH/Methyl-tertbuthyl ether
* MeOH/DimethoxyEthane
* MeOH/Dimethylcarbonate
* MeOH/CycloHexane


# Assumptions and applicability

* The activity coefficients of the binary mixture are calculated by means of NRTL model;
* Saturated vapour pressure and Vaporisation/Condensation heat are assessed usnig Antoine and Clapeyron-Clausius equations
* The ideal modelling process is applicable only for the modelling of processes, where permeance values does depend on mixture composition
* The non-ideal modelling is performed only based on the specified diffusion curve (The overall and component fluxes as a function of separated mixture's composition)
