import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Resonant import * 
from sfactors.nonResonant import * 
from sfactors.Plots.reactionPlot import *

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