import numpy as np
R = 8.314462
class Component:
    def __init__(self,A,B,C):
        self.AntoineConstants = [A,B,C]
        # Calculate Pressure in kPa from Antoine's equation Component = [A,B,C,name] (by the base of 10)
    def getAntoinePressure (self,T):
           return 10**(self.AntoineConstants[0]-self.AntoineConstants[1]/(T+self.AntoineConstants[2]))*100
        #Calculate Vaporisation/Condensation enthalpy in J/mol from Clapeyron-Clausius equation
    def getVaporisationHeat (self,T):
        return (T/(T+self.AntoineConstants[2]))**2*R*self.AntoineConstants[1]*np.log(10)
        #Calculate molar heat Capacity from definition of isobaric heat capacity
class Mixture:
    def __init__(self,g12,g21,alpha,component1,component2):
        self.Components = [component1,component2]
        self.NRTLParameters = [g12,g21,alpha]
        # Calculate Activity Coefficients using NRTL model for binary mixture
        # g=[g12,g21], mixture = [g12,g21,alpha,Component1,Component2,name]
    def getNRTLPartialPressures(self,T,Composition):
        #Composition is in molar fractions
        g = self.NRTLParameters[0:2]
        t = np.multiply(1/(R*T),g)
        G = np.exp(-t*self.NRTLParameters[2])
        ActivityCoefficients=[]
        ComponentsPartialPressure=[]
        ActivityCoefficients.append(np.exp((Composition[1]**2)*(t[1]*(G[1]/(Composition[0]+Composition[1]*G[1]))**2+t[0]*G[0]/(Composition[1]+Composition[0]*G[0])**2)))
        ActivityCoefficients.append(np.exp((Composition[0]**2)*(t[0]*(G[0]/(Composition[1]+Composition[0]*G[0]))**2+t[1]*G[1]/(Composition[0]+Composition[1]*G[1])**2)))
        ComponentsPartialPressure.append(self.Components[0].getAntoinePressure(T)*ActivityCoefficients[0]*Composition[0])
        ComponentsPartialPressure.append(self.Components[1].getAntoinePressure(T)*ActivityCoefficients[1]*Composition[1])  
        return ComponentsPartialPressure

class Membrane:
    def __init__(self, Permeances, Energies) -> None:
        
        pass

