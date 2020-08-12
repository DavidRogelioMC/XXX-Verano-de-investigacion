# This is CDM cosmology with w, wa and Ok


import math as N
from simplemc.models.LCDMCosmology import LCDMCosmology
from simplemc.cosmo.paramDefs import w_par, wa_par, Ok_par, q0, j0, Omega0

class owa0CDMCosmology(LCDMCosmology):
    def __init__(self, varyw=False, varywa=False, varyOk=False, varyq0=True, varyj0=True, varyOmega0=True):
        # three parameters: w, wa, Ok

        self.varyw  = varyw
        self.varywa = varywa
        self.varyOk = varyOk
        
        self.varyq0 = varyq0
        self.varyj0 = varyj0
        self.varyOmega0 = varyOmega0

        self.Ok = Ok_par.value
        self.w0 = w_par.value
        self.wa = wa_par.value
        
        self.q0=q0.value
        self.j0=j0.value
        self.Omega0 = Omega0.value
        
        LCDMCosmology.__init__(self)


    # my free parameters. We add Ok on top of LCDM ones (we inherit LCDM)
    def freeParameters(self):
        l = LCDMCosmology.freeParameters(self)
        if (self.varyw):  l.append(w_par)
        if (self.varywa): l.append(wa_par)
        if (self.varyOk): l.append(Ok_par)
        if (self.varyq0): l.append(q0)
        if (self.varyj0): l.append(j0)
        if (self.varyOmega0): l.append(Omega0)
        return l


    def updateParams(self, pars):
        ok = LCDMCosmology.updateParams(self, pars)
        if not ok:
            return False
        for p in pars:
            if p.name == "w":
                self.w0 = p.value
            elif p.name == "wa":
                self.wa = p.value
            elif p.name == "Ok":
                self.Ok = p.value
                self.setCurvature(self.Ok)
                if (abs(self.Ok) > 1.0):
                    return False
            elif p.name == "q0":
                self.q0 = p.value
            elif p.name =="j0":
                self.j0 = p.value
            elif p.name == "Omega0":
                self.Omega0 = p.value
        return True


    # this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    #def RHSquared_a(self, a):
     #   NuContrib = self.NuDensity.rho(a)/self.h**2
      #  rhow = a**(-3*(1.0+self.w0+self.wa))*N.exp(-3*self.wa*(1-a))
       # return (self.Ocb/a**3+self.Ok/a**2+self.Omrad/a**4+NuContrib+(1.0-self.Om-self.Ok)*rhow)
     
    #Weak    
   # def RHSquared_a(self,a):
    #    c=3e8 # m/s
    #    beta = -2
     #   gamma = -3
      #  H0 = 67.4*1000/3.085e22
      #  kk = -2.3e-5 #Producto de la k y k^prime cuando gamma=-3 y beta = -2
       # a0 = 1.2e-10 # m/s^2
       # Z= 1 + (gamma-1)*((self.q0-1)/gamma + (self.j0-self.q0-2)/(1-self.q0))
       # return ((a0/H0*c)**2)*(((kk*(self.q0-1)**(1-gamma))/((6**gamma)*gamma*Z))*((3*self.Omega0/(8*N.pi))**beta))**(1/(gamma-beta))
    
    #Strong
    def RHSquared_a(self,a):
        c=3e8 # m/s
        beta = 3
        gamma = -3
        tau=-3.0 
        H0 = 70*1000/3.085e22
        kk = 9/(((N.pi)**2)*(4**5)) #Producto de la k y k^prime cuando gamma=-3 y beta = 3
        a0 = 1.2e-10 # m/s^2
        Z= 1 + (1-gamma)*((1-self.q0)/gamma - (self.j0-self.q0-2)/(1-self.q0))+tau*beta
        return ((a0/H0*c)**2)*(((6*(self.q0-1))**(1-gamma))*((8*N.pi*kk)/(3*gamma*Z))*((3*self.Omega0/(8*N.pi))**(1-beta)))**(1/(gamma+beta-1))

