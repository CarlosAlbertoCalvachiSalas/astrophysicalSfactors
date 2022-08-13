import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Resonant import * 
from sfactors.nonResonant import * 
from sfactors.Plots.reactionPlot import *

from sfactors.Utils.formulas import sfactorFromCrossSection

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "sans-serif"
})


r16O18O = ReactionPlot('16O', '18O', crossSection = True)
O16O18Yakovlev = Yakovlev(Z1=8, Z2=8, A1=r16O18O.r1.A, A2=r16O18O.r2.A)
plotFuncOverData(O16O18Yakovlev.sFactor, r16O18O.E, r16O18O.S, label = 'Yakovlev')
plt.yscale('log')
r16O18O.show()

r16O17O = ReactionPlot('16O', '17O', crossSection = True)
O16O17Yakovlev = Yakovlev(Z1=8, Z2=8, A1=15.99491461956, A2=17)
plotFuncOverData(O16O17Yakovlev.sFactor, r16O17O.E, r16O17O.S, label = 'Yakovlev')
plt.yscale('log')
r16O17O.show()

rp16O16 = ReactionPlot('16O', '16O')
O16Yakovlev = Yakovlev(Z1=8, Z2=8, A1=15.99491461956, A2=15.99491461956)
plotFuncOverData(O16Yakovlev.sFactor, rp16O16.E, rp16O16.S, label = 'Yakovlev')
plt.ylim(1e22, 1e26)
plt.xlim(6, 13)
plt.yscale('log')
#plt.title('Astrophysical S-factor energy curve for the 16O(16O) reaction')
rp16O16.show()

plt.show()


r12C16O = ReactionPlot('12C', '16O')
C12O16Yakovlev = Yakovlev(Z1=6, Z2=8, A1=12, A2=15.99491461956)
plotFuncOverData(C12O16Yakovlev.sFactor, r12C16O.E, r12C16O.S, label = 'Yakovlev')
plt.yscale('log')
r12C16O.show()


rp12C12C = ReactionPlot('12C', '12C')
C12Yakovlev = Yakovlev(Z1=6, Z2=6, A1=12, A2=12)
plotFuncOverData(C12Yakovlev.sFactor, rp12C12C.E, rp12C12C.S, label = 'Yakovlev')
plt.yscale('log')
rp12C12C.show()


rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(exponential1, p0 = [0.7E-1, 0.7], color = 'red', label = 'exp1')
rp2Hd.fit(exponential2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'exp2')
rp2Hd.fit(exponential3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'exp3')
rp2Hd.fit(exponential4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'exp4')
#plotFuncOverData(lambda E: exponential1(E,0.7E-1, 0.7), rp2Hp.E, rp2Hp.S)
plt.xscale('log')
plt.yscale('log')
rp2Hd.save()
rp2Hd.show()



rp7Bep = ReactionPlot('7Be','p','gamma','8B')
#rp7Bep.fit(breitWigner, p0 = [110, 0.65, 0.08], Emin = 0.5, Emax = 0.8, color = 'red')
rp7Bep.fit(breitWigner, p0 = [110, 0.65, 0.08], Emin = 0.5, Emax = 0.8, Smin = 30 , color = 'red', label = 'BW')
#plotFuncOverData(lambda E: breitWigner(E, 110, 0.65, 0.08), rp7Bep.E, rp7Bep.S)
#plotFuncOverData(lambda E:  breitWignerBackground(E, 110, 0.65, 0.08, 19.09067254, 0.30103139), rp7Bep.E, rp7Bep.S)
rp7Bep.fit(exponential1, color = 'green', label = 'exp')
rp7Bep.save()
rp7Bep.show()


rp7Bep = ReactionPlot('7Be','p','gamma','8B')
rp7Bep.fit(breitWignerBackground, p0 = [110, 0.65, 0.08, 19.09067254, 0.30103139], color = 'blue', label = 'exp+BW')
rp7Bep.show()

rp13Cp = ReactionPlot('13C','p','gamma','14N')
#plotFuncOverData(lambda E: breitWigner(E, 2.03956298, 0.51380586, 0.0325387), rp7Bep.E, rp7Bep.S)
rp13Cp.fit(breitWigner, p0 = [1, 0.5, 0.2], color = 'red')
rp13Cp.save()
rp13Cp.show()

plt.scatter(rp13Cp.E,rp13Cp.residuals)
plotFitFunc(exponential1, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00] )  
#rp13Cp.fit(exponential1, color = 'green', label = 'exp', Emax = 0.44)
plt.title('Residuals BW fit')
plt.ylabel('S(MeV-b)')
plt.xlabel('E(MeV)')
plt.yscale('log')
plt.show()


rp13Cp = ReactionPlot('13C','p','gamma','14N')
rp13Cp.fit(breitWigner, p0 = [1, 0.5, 0.2], color = 'red')
rp13Cp.fit(breitWignerBackground, p0 = [2.03956298,0.51380586,0.03253878,0.00862043,0.18103796], color = 'blue', label = 'exp+BW')
#plotFuncOverData(lambda E: breitWigner(E, 2.03956298, 0.51380586, 0.0325387), rp7Bep.E, rp7Bep.S)
#rp13Cp.fit(breitWigner, p0 = [1, 0.5, 0.2], color = 'red')
rp13Cp.show()



rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(exponential1, p0 = [0.7E-1, 0.7], color = 'red', label = 'exp1')
rp2Hd.fit(exponential2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'exp2')
rp2Hd.fit(exponential3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'exp3')
rp2Hd.fit(exponential4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'exp4')
#plotFuncOverData(lambda E: exponential1(E,0.7E-1, 0.7), rp2Hp.E, rp2Hp.S)
plt.xscale('log')
plt.yscale('log')
rp2Hd.save()
rp2Hd.show()
#plt.show()

rp2Hp = ReactionPlot('2H','p','gamma','3He')
#rp2Hp.fit(exponential1, p0 = [0.2, 3], color = 'red', label = 'exp1')
#rp2Hp.fit(exponential2, p0 = [0.2, 3, -0.2], color = 'green', label = 'exp2')
rp2Hp.fit(exponential3, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428], color = 'brown', label = 'exp3-lowE', Emax = 1)
rp2Hp.fit(exponential3, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428], color = 'blue', label = 'exp3-higE')
#rp2Hp.fit(exponential4, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428, 0.05], label = 'exp4', Emax = 0.5)
#plotFuncOverData(lambda E: exponential3(E,  0.275367003, 1.33359118, -0.14065138, 0.01865428), rp2Hp.E, rp2Hp.S)
plt.xscale('log')
plt.yscale('log')
rp2Hp.save()
rp2Hp.show()



[2.03956298, 0.51380586, 0.0325387]
[99.1531399, 0.631883871, 0.0535006075]



#H2pYakovlev = Yakovlev(Z1=1, Z2=1, A1=1, A2=1, S0=1, Ec=0.5, delta=0.6263, xi=10)
#plt.show()


#H2nYakovlev = Yakovlev(Z1=1, Z2=1, A1=1, A2=2, S0=0.5, Ec=1, delta=0.1, xi=10)
#plt.show()


"""
DdData = pd.read_csv('Databases/Dd.csv')
DpData = pd.read_csv('Databases/Dp.csv')
C13Data = pd.read_csv('Databases/C13p.csv')
C13DataSquare = pd.read_csv('Databases/C13pSquare.csv')
C13DataCircle = pd.read_csv('Databases/C13pCircle.csv')
C12C12StarData = pd.read_csv('Databases/C12C12Star.csv')

plt.scatter(C13DataSquare['E(MeV)'], C13DataSquare['S(MeV-b)'], color = 'black', s = 7, marker ='s', label = 'WO52') 
plt.scatter(C13DataCircle['E(MeV)'], C13DataCircle['S(MeV-b)'], color = 'black', s = 7, marker ='o', label = 'KI94')
params = fitBreitWigner(C13Data['E(MeV)'], C13Data['S(MeV-b)'], p0=[1, 0.5, 0.2])
print(params)
continumE = np.linspace(np.min(C13Data['E(MeV)']), np.max(C13Data['E(MeV)']), 1000)
plt.plot(continumE, breitWigner(continumE, *params[0]), color = 'red')
plt.xlabel('E(MeV)')
plt.ylabel('S(MeV-b)')
plt.yscale('log')
plt.title(r'Astrophysical S-factor energy curve for the $\mathrm{{}^{13}C(p,\gamma){}^{14}N}$ reaction')
plt.legend()
plt.savefig('ReconstructedPlots/C13.eps', format = 'eps')
plt.show()





plt.scatter(DpData['E(MeV)'], DpData['S(eV-b)'], color = 'blue', s = 7)
plt.xlabel('E(keV)')
plt.ylabel('S(eV-b)')
plt.xscale('log')
plt.yscale('log')
plt.title(r'Astrophysical S-factor energy curve for the $\mathrm{D(p,\gamma){}^{3}He}$ reaction')
plt.savefig('ReconstructedPlots/Dp.eps', format = 'eps')
plt.show()

plt.scatter(DdData['E(MeV)'], DdData['S(MeV-b)'], color = 'red', s = 7)
plt.xlabel('E(MeV)')
plt.ylabel('S(MeV-b)')
plt.xscale('log')
plt.yscale('log')
plt.title(r'Astrophysical S-factor energy curve for the $\mathrm{D(d,n){}^{3}He}$ reaction')
plt.savefig('ReconstructedPlots/Dd.eps', format = 'eps')
plt.show()

plt.scatter(C12C12StarData['E(MeV)'], C12C12StarData['S(MeV-b)'], color = 'purple', s = 7, label = r'$S^{\star}$')
plt.scatter(C12C12StarData['E(MeV)'], C12C12StarData['S(MeV-b)']*np.exp(-0.46*C12C12StarData['E(MeV)']), color = 'orange', s = 7, label = r'$S$')
plt.xlabel('E(MeV)')
plt.ylabel('S(MeV-b)')
plt.yscale('log')
plt.title(r'Astrophysical S-factor energy curve for the $\mathrm{{}^{12}C({}^{12}C,{}^{4}He){}^{20}Ne}$ reaction')
plt.legend()
plt.savefig('ReconstructedPlots/C12C12.eps', format = 'eps')
plt.show()

"""
