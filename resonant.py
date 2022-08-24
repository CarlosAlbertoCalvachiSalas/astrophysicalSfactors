import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Resonant import * 
from sfactors.nonResonant import * 
from sfactors.Plots.reactionPlot import *

rp7Bep = ReactionPlot('7Be','p','gamma','8B')
rp7Bep.fit(polynomial1, color = 'grey', label = 'poly1')
rp7Bep.fit(polynomial2, color = 'dodgerblue', label = 'poly2')
rp7Bep.fit(polynomial3, color = 'chocolate', label = 'poly3', p0 = [18.31552157,8.04220547,-0.21736566, 0])
rp7Bep.fit(polynomial4, color = 'seagreen', label = 'poly4', p0 = [18.31552157,8.04220547,-0.21736566, 0, 0])
rp7Bep.fit(polynomial(5), color = 'darkorange', label = 'poly5', p0 = [18.31552157,8.04220547,-0.21736566, 0, 0, 0])
rp7Bep.fit(polynomial(6), color = 'violet', label = 'poly6',  p0 = [18.31552157,8.04220547,-0.21736566, 0, 0, 0, 0])
rp7Bep.save()
rp7Bep.show()

rp7Bep = ReactionPlot('7Be','p','gamma','8B')
#rp7Bep.fit(breitWigner, p0 = [110, 0.65, 0.08], Emin = 0.5, Emax = 0.8, color = 'red')
rp7Bep.fit(breitWigner, p0 = [110, 0.65, 0.08], Emin = 0.5, Emax = 0.8, Smin = 30 , color = 'red', label = 'BW')
#plotFuncOverData(lambda E: breitWigner(E, 110, 0.65, 0.08), rp7Bep.E, rp7Bep.S)
#plotFuncOverData(lambda E:  breitWignerBackground(E, 110, 0.65, 0.08, 19.09067254, 0.30103139), rp7Bep.E, rp7Bep.S)
rp7Bep.fit(exponential1, color = 'green', label = 'exp1')
rp7Bep.fit(polynomial1, color = 'grey', label = 'poly1')
rp7Bep.save()
rp7Bep.show()

rp7Bep = ReactionPlot('7Be','p','gamma','8B')
rp7Bep.fit(breitWignerBackground1, p0 = [110, 0.65, 0.08, 19.09067254, 0.30103139], color = 'blue', label = 'exp+BW')
rp7Bep.fit(rmatrixPolynomial1, p0 = [110, 0.65, 0.08, (0.08**2/4), 19.09067254, 0.30103139], color = 'red', label = 'poly+BW')
rp7Bep.show()


rp13Cp = ReactionPlot('13C','p','gamma','14N')
#plotFuncOverData(lambda E: breitWigner(E, 2.03956298, 0.51380586, 0.0325387), rp7Bep.E, rp7Bep.S)
rp13Cp.fit(breitWigner, p0 = [1, 0.5, 0.2], color = 'red', label = 'BW')
rp13Cp.save()
rp13Cp.show()


plt.scatter(rp13Cp.E,rp13Cp.residuals)
plotFitFunc(exponential1, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00], label = 'exp1')  
plotFitFunc(polynomial1, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00], color = 'red', label = 'poly1')  
#rp13Cp.fit(exponential1, color = 'green', label = 'exp', Emax = 0.44)
plt.title('Residuals BW fit')
plt.ylabel('S(MeV-b)')
plt.legend()
plt.xlabel('E(MeV)')
plt.yscale('log')
plt.show()


plt.scatter(rp13Cp.E,rp13Cp.residuals)
plotFitFunc(polynomial1, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00], color = 'red', label = 'poly1')  
plotFitFunc(polynomial2, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0], color = 'black', label = 'poly2')  
plotFitFunc(polynomial3, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0], color = 'seagreen', label = 'poly3')  
plotFitFunc(polynomial4, rp13Cp.E, rp13Cp.residuals, color = 'grey', label = 'poly4')  
# p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0],
#plotFitFunc(polynomial(5), rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0, 0], color = 'navy', label = 'poly5')  
#plotFitFunc(polynomial(6), rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0, 0, 0], color = 'fuchsia', label = 'poly6')  
#rp13Cp.fit(exponential1, color = 'green', label = 'exp', Emax = 0.44)
plt.title('Residuals BW fit')
plt.ylabel('S(MeV-b)')
plt.xlabel('E(MeV)')
plt.legend()
plt.yscale('log')
plt.show()

plt.scatter(rp13Cp.E,rp13Cp.residuals)
plotFitFunc(exponential1, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00], color = 'red', label = 'exp1')  
plotFitFunc(exponential2, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0], color = 'black', label = 'exp2')  
plotFitFunc(exponential3, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0], color = 'seagreen', label = 'exp3')  
plotFitFunc(exponential4, rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0], color = 'grey', label = 'exp4')  
# p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0],
#plotFitFunc(exponential(5), rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0, 0], color = 'navy', label = 'exp5')  
#plotFitFunc(exponential(6), rp13Cp.E, rp13Cp.residuals, p0 = [2.49185090e-03,7.81050299e+00, 0, 0, 0, 0, 0], color = 'fuchsia', label = 'exp6')  
#rp13Cp.fit(exponential1, color = 'green', label = 'exp', Emax = 0.44)
plt.title('Residuals BW fit')
plt.ylabel('S(MeV-b)')
plt.xlabel('E(MeV)')
plt.legend()
plt.yscale('log')
plt.show()


rp13Cp = ReactionPlot('13C','p','gamma','14N')
rp13Cp.fit(breitWigner, p0 = [1, 0.5, 0.2], color = 'red', label = 'BW')
rp13Cp.fit(breitWignerBackground1, p0 = [2.03956298,0.51380586,0.03253878,0.00862043,0.18103796], color = 'blue', label = 'exp1+BW')
rp13Cp.fit(rmatrixPolynomial1, p0 = [2.03956298,0.51380586,0.03253878, (0.03253878**2/4), 0.0084133, 0.0021128], color = 'indianred', label = 'poly1+BW')
rp13Cp.fit(breitWignerBackground2, color = 'forestgreen', p0 = [2.03956298,0.51380586,0.03253878,-0.00203961,0.06281975,-0.06971001], bounds = ([2.0, 0.5, 0.03, -1, -1, -1], [2.5, 0.6, 0.04, 1, 1, 1]), label = 'exp2+BW')
rp13Cp.fit(rmatrixPolynomial2, color = 'darkorange', p0 = [2.03956298,0.51380586,0.03253878, (0.03253878**2/4),2.92356824e-03,6.04883663e+00,-6.52640185e+00],  bounds = ([2.0, 0.5, 0.03, -10, -10, -10, -10], [2.5, 0.6, 0.04, 10, 10, 10, 10]), label = 'poly2+BW')
#plotFuncOverData(lambda E: breitWigner(E, 2.03956298, 0.51380586, 0.0325387), rp7Bep.E, rp7Bep.S)
#rp13Cp.fit(breitWigner, p0 = [1, 0.5, 0.2], color = 'red')
#plotFitFunc(breitWignerBackground2, rp13Cp.E, rp13Cp.S, p, color = 'fuchsia', label = 'exp6')  
#plotFuncOverData(lambda E: breitWignerBackground2(E, *[2.03956298,0.51380586,0.03253878,-0.00203961,0.06281975,-0.06971001] ), rp13Cp.E, rp13Cp.S, color = 'forestgreen', label = 'exp2+BW')
#plotFuncOverData(lambda E: rmatrixPolynomial2(E, *[2.03956298,0.51380586,0.03253878, (0.03253878**2/4),2.92356824e-03,6.04883663e+00,-6.52640185e+00] ), rp13Cp.E, rp13Cp.S, color = 'darkorange', label = 'poly2+BW')

rp13Cp.show()





