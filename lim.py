import numpy as np
import math
from scipy.integrate import dblquad
import sys

if len(sys.argv) < 2:
        print("It is correct")
        sys.exit(1)

interaction = sys.argv[1].upper()

if interaction not in ["SI", "SD"]:
        print("Error, use SI or SD")
        sys.exit(1)

##Definition of units of measurement

s = 1.
kg = 1.
m = 1.
m3 = 1.
fm = 10**(-15)*m
gr = 10**(-3)
Joule = 1.*kg*m*m/s/s
u = 1.66054*10**(-27)*kg
eV = 1.60218*10**(-19)*Joule
GeV = 1*10**(9)*eV
MeV = 1*10**(6)*eV
keV = 1*10**(3)*eV
c = 3*10**(8)*m/s
ht = 1.055*10**(-34)*Joule*s
GeVc2 = 1*GeV/c/c
cm2 = 1*10**(-4)*m*m
cm3 = 1*10**(-6)*m*m*m
day = 86400.*s
y = 365*day
bar = 1.
atmpressure = 1.01325*bar

## Parameters that can be changed
LowEThre = 1.5*keV   		#energy threshold in ee		
percCF4 = 0.4			# percentage of CF4
percHe = 0.6			# percentage He
Daqtime = 1502711 # 25*day		# correct exposuretime
workingpressure = 0.907*bar     # working pressure
volume = 0.043*m3 #0.022*m3 #0.03*m3 # 0.027*m3		# sensitive volume (it should include geometrical cuts

## Constant values

rho0 = 0.3*GeVc2/cm3		# DM density at Earth
vc = 230*10**(3)*m/s		# average velocity of Halo at Earth
vesc = 544*10**(3)*m/s		# galactic escape velocity
ve = 242*10**(3)*m/s		# Tangent velocity of Earth
mp = 1.007*u			# proton mass
mn = 1.009*u			# neutron mass
NA = 6.022*10**(26)		# Avogadro's number
AHe = 4			# mass number He
AC = 12			# mass number C
AF = 19			# mass number F
AmolHe = 4			# Molar mass He
AmolC = AmolF = 88		# Molar mass CF4
NatomiHe = NatomiC = 1		# number of atoms of He in He molecule and C atoms in CF4 molecule
NatomiF = 4			# number of F atoms in CF4 molecule
rhoHe = 0.1662*kg/m3		# density of helium gas (all densities need to be used at atmospheric pressure
rhoC = rhoF = 3.72*kg/m3	# density of CF4
massHe = 4.0026*u		# Atomic mass He
massC = 12.011*u		# Atomic mass C
massF = 18.998*u		# Atomic mass F
zbar = vesc/vc
AN3 = 1/((np.pi*vc**2)**(3/2)*(math.erf(zbar) - 2/np.sqrt(np.pi)*zbar*np.exp(-zbar**2)))  # Integral normalisation factor

## Importing tables containing DM masses under test and the corresponding signal events found by the fit to be the 95% CI 

masses = np.loadtxt("masses.txt")
Nev = 650
TabMassEvt = np.column_stack((masses * GeVc2, np.full(len(masses), Nev)))

## Quenching factor part
## a b and c are parameters of the function (qf) used to fit Flaminia's simulation of QF for He, C, F
## ap bp cp and dp are parameters of the inverse function of qf
#print('mx', TabMassEvt[0,0])

aHe = 0.131
bHe = 3.76
cHe = 0.345
apHe = 0.250
bpHe = 0.768
cpHe = 3.90*keV
dpHe = 0.67
def qfHe(x):
	return aHe*(x + bHe*keV**(1 - cHe)*x**cHe)/(1*keV +  aHe*(x + bHe*keV**(1 - cHe)*x**cHe))
def qfpHe(x):
	return apHe + bpHe/(1 + (cpHe/x)**dpHe)

aC = 0.0203
bC = 14.8
cC = 0.306
apC = 0.157
bpC = 0.93
cpC = 21.52*keV
dpC = 0.506
def qfC(x):
	return aC*(x + bC*keV**(1 - cC)*x**cC)/(1*keV + aC*(x + bC*keV**(1 - cC)*x**cC))
def qfpC(x):
	return apC + bpC/(1 + (cpC/x)**dpC)

aF = 0.00857
bF = 27.7
cF = 0.284
apF = 0.127
bpF = 1.03
cpF = 68.3*keV
dpF = 0.446
def qfF(x):
	return aF*(x + bF*keV**(1 - cF)*x**cF)/(1*keV + aF*(x + bF*keV**(1 - cF)*x**cF))
def qfpF(x):
	return apF + bpF/(1 + (cpF/x)**dpF);
	
LowEThreHe = LowEThre/qfpHe(LowEThre)
#LowEThreHe = LowEThreHe/keV		# Energy threshold for He after QF
LowEThreC = LowEThre/qfpC(LowEThre)
#LowEThreC = LowEThreC/keV
LowEThreF = LowEThre/qfpF(LowEThre)
#LowEThreF = LowEThreF/keV

## Other factors required in the calculation of the limit
def muHe2(x):
	return massHe*TabMassEvt[x, 0]/(massHe + TabMassEvt[x, 0])		# reduced mass
def muC2(x):
	return massC*TabMassEvt[x, 0]/(massC + TabMassEvt[x, 0])
def muF2(x):
	return massF*TabMassEvt[x, 0]/(massF + TabMassEvt[x, 0])
def mup2(x):
	return mp*TabMassEvt[x, 0]/(mp + TabMassEvt[x, 0])
def rHe2(x):
	return 4*TabMassEvt[x, 0]*massHe/(TabMassEvt[x, 0] + massHe)**2		# kinematic ratio
def rC2(x):
	return 4*TabMassEvt[x, 0]*massC/(TabMassEvt[x, 0] + massC)**2
def rF2(x):
	return 4*TabMassEvt[x, 0]*massF/(TabMassEvt[x, 0] + massF)**2
def maxiHe2(x):
	return 0.5*TabMassEvt[x, 0]*(vesc + ve)*(vesc + ve)*rHe2(x)/keV
def maximinHe2(x):
	return 0.5*TabMassEvt[x, 0]*(vesc - ve)*(vesc - ve)*rHe2(x)/keV
def maxiC2(x):
	return 0.5*TabMassEvt[x, 0]*(vesc + ve)*(vesc + ve)*rC2(x)/keV
def maximinC2(x):
	return 0.5*TabMassEvt[x, 0]*(vesc - ve)*(vesc - ve)*rC2(x)/keV
def maxiF2(x):
	return 0.5*TabMassEvt[x, 0]*(vesc + ve)*(vesc + ve)*rF2(x)/keV
def maximinF2(x):
	return 0.5*TabMassEvt[x, 0]*(vesc - ve)*(vesc - ve)*rF2(x)/keV
	
exposureHe = rhoHe*workingpressure/atmpressure * volume*percHe*Daqtime		#exposure
exposureC = rhoC*workingpressure/atmpressure *volume*percCF4*Daqtime
exposureF = rhoF*workingpressure/atmpressure *volume*percCF4*Daqtime

RHe = (0.3 + 0.91*AHe**(1/3))						# Radius He (form factor calculation)
def SHe(En):
	qHe = np.sqrt(2*1000*massHe/u*En/MeV)
	return abs(9*(math.sin(qHe*RHe/(197)) - qHe*RHe*math.cos(qHe*RHe/(197))/197)**2/(qHe*RHe/(197))**6)**2

RC = (0.3 + 0.91*AC**(1/3))
def SC(En):
	qC = np.sqrt(2*1000*massC/u*En/MeV)
	return abs(9*(math.sin(qC*RC/(197)) - qC*RC*math.cos(qC*RC/(197))/197)**2/(qC*RC/(197))**6)**2

RF = (0.3 + 0.91*AF**(1/3))
def SF(En):
	qF = np.sqrt(2*1000*massF/u*En/MeV)
	return abs(9*(math.sin(qF*RF/(197)) - qF*RF*math.cos(qF*RF/(197))/197)**2/(qF*RF/(197))**6)**2

def compute_integrals(i):
        SHem , _ = dblquad(lambda v, En: SHe(En)*np.exp(-(v - ve)**2/vc**2), min(LowEThreHe, maxiHe2(i)*keV), maxiHe2(i)*keV, lambda En: np.sqrt(2*massHe*En)/(2*muHe2(i)), vesc + ve, epsrel=1.49e-50)
        SHep , _ = dblquad(lambda v, En: SHe(En)*np.exp(-(v + ve)**2/vc**2), min(LowEThreHe, maximinHe2(i)*keV), maximinHe2(i)*keV, lambda En: np.sqrt(2*massHe*En)/(2*muHe2(i)), vesc - ve, epsrel=1.49e-50)
        SHe0, _ = dblquad(lambda v, En: SHe(En), min(LowEThreHe, maxiHe2(i)*keV), maxiHe2(i)*keV, lambda En: max(vesc - ve, np.sqrt(2*massHe*En)/(2*muHe2(i))), vesc + ve, epsrel=1.49e-50)

        SCm , _ = dblquad(lambda v, En: SC(En)*np.exp(-(v - ve)**2/vc**2), min(LowEThreC, maxiC2(i)*keV), maxiC2(i)*keV, lambda En: np.sqrt(2*massC*En)/(2*muC2(i)), vesc + ve, epsrel=1.49e-50)
        SCp , _ = dblquad(lambda v, En: SC(En)*np.exp(-(v + ve)**2/vc**2), min(LowEThreC, maximinC2(i)*keV), maximinC2(i)*keV, lambda En: np.sqrt(2*massC*En)/(2*muC2(i)), vesc - ve, epsrel=1.49e-50)
        SC0, _ = dblquad(lambda v, En: SC(En), min(LowEThreC, maxiC2(i)*keV), maxiC2(i)*keV, lambda En: max(vesc - ve, np.sqrt(2*massC*En)/(2*muC2(i))), vesc + ve, epsrel=1.49e-50)

        SFm , _ = dblquad(lambda v, En: SF(En)*np.exp(-(v - ve)**2/vc**2), min(LowEThreF, maxiF2(i)*keV), maxiF2(i)*keV, lambda En: np.sqrt(2*massF*En)/(2*muF2(i)), vesc + ve, epsrel=1.49e-50)
        SFp , _ = dblquad(lambda v, En: SF(En)*np.exp(-(v + ve)**2/vc**2), min(LowEThreF, maximinF2(i)*keV), maximinF2(i)*keV, lambda En: np.sqrt(2*massF*En)/(2*muF2(i)), vesc - ve, epsrel=1.49e-50)
        SF0, _ = dblquad(lambda v, En: SF(En), min(LowEThreF, maxiF2(i)*keV), maxiF2(i)*keV, lambda En: max(vesc - ve, np.sqrt(2*massF*En)/(2*muF2(i))), vesc + ve, epsrel=1.49e-50)

        return SHem, SHep, SHe0, SCm, SCp, SC0, SFm, SFp, SF0

results = []

spinHe = 0
spinC = 0
spinF = 0.91

if interaction == "SI":
        pHe = AHe
        pC = AC
        pF = AF
if interaction == "SD":
        pHe = spinHe
        pC = spinC
        pF = spinF

for j in range(len(TabMassEvt)):
        SHem, SHep, SHe0, SCm, SCp, SC0, SFm, SFp, SF0 = compute_integrals(j)

        den =(exposureHe*NA/AmolHe*NatomiHe*2*
              rho0*(massHe + TabMassEvt[j, 0])**2/(4*TabMassEvt[j, 0]**3*
                                                   massHe)*pHe**2*(muHe2(j)/mup2(j))**2*AN3*np.pi*
              vc**2/ve*(SHem - SHep - np.exp(-vesc**2/vc**2)* SHe0) +
              exposureC*NA/AmolC*NatomiC*2*
              rho0*(massC + TabMassEvt[j, 0])**2/(4*TabMassEvt[j, 0]**3*
                                                  massC)*pC**2*(muC2(j)/mup2(j))**2*AN3*np.pi*
              vc**2/ve*(SCm - SCp - np.exp(-vesc**2/vc**2)* SC0) + 
              exposureF*NA/AmolF*NatomiF*2*
              rho0*(massF + TabMassEvt[j, 0])**2/(4*TabMassEvt[j, 0]**3*
                                                  massF)*pF**2*(muF2(j)/mup2(j))**2*AN3*np.pi*
              vc**2/ve*(SFm - SFp - np.exp(-vesc**2/vc**2)* SF0))

        if den <= 0:
                sigma_val = np.inf
        else:
                sigma_val = TabMassEvt[j,1]/cm2/den

        m = TabMassEvt[j,0] / GeVc2
        results.append([m,sigma_val])
        print(m, sigma_val)



np.savetxt("lim_output.txt",
           results,
#           header = "m_DM[GeV] sigma[cm2]"
           )




