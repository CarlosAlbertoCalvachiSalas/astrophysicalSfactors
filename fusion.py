import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Resonant import * 
from sfactors.nonResonant import * 
from sfactors.nonResonant.fusion import empiricalFusion34
from sfactors.Plots.reactionPlot import *

"""
yPrediction16O18O = lambda E: empiricalFusion34(E, 
						b1 = 66.9, 
						b2 = -0.710, 
						b3 = -0.01137, 
						c1 = 4, 
						c2 = -0.19, 
						c3 = -0.05928, 
						c4 = 0.002232, 
						Ec = 10.871, 
						D = 0.89)

"""
# 66 - 67 

"""
						b1 = 66.9, 
						b2 = -0.710, 
						b3 = -0.01137, 
						c1 = 4, 
						c2 = -0.19, 
						c3 = -0.05928, 
						c4 = 0.002232, 
						Ec = 10.871, 
						D = 0.89)

Checkpoint

"""



yPrediction16O18O = lambda E: empiricalFusion34(E, 
						b1 = 63.550, 
						b2 = -0.4610, 
						b3 = -0.01137, 
						c1 = 5.498, 
						c2 = -0.2623, 
						c3 = -0.05928, 
						c4 = 0.002232, 
						Ec = 10.871, 
						D = 0.89)





r16O18O = ReactionPlot('16O', '18O', crossSection = True)
r16O18O.noSigma()
O16O18Yakovlev = Yakovlev(Z1=8, Z2=8, A1=r16O18O.r1.A, A2=r16O18O.r2.A)
params, errors = r16O18O.fit(empiricalFusion34, p0 = 	[66.9, -0.710,  -0.01137,  4, -0.19, -0.05928,  0.002232, 10.871, 0.89], 
							bounds = (	[66, -0.8, -0.1, 3.5, -0.25, -0.1, 0, 10.86, 0.88 ],
										[67.5, -0.6, 0,  4.5, -0.10, 0, 0.1, 10.88, 0.90 ]),  color = 'red', label = 'empirical-fit')
plt.close()


r16O18O = ReactionPlot('16O', '18O', crossSection = True)
O16O18Yakovlev = Yakovlev(Z1=8, Z2=8, A1=r16O18O.r1.A, A2=r16O18O.r2.A)
plotFuncOverData(O16O18Yakovlev.sFactor, r16O18O.E, r16O18O.S, label = 'potential-Yakovlev')
plotFuncOverData(yPrediction16O18O,  r16O18O.E,  r16O18O.S, label = 'empirical-Yakovlev', color = 'green')
plotFuncOverData(lambda E: empiricalFusion34(E, *params),  r16O18O.E,  r16O18O.S, label = 'empirical-fit', color = 'red')
plt.yscale('log')
r16O18O.save()
r16O18O.show()

"""
yPrediction16O17O = lambda E: empiricalFusion34(E, 
						b1 = 70, 
						b2 = -1.15, 
						b3 =  -0.01137, 
						c1 =  4, 
						c2 = -0.155, 
						c3 = -0.046, 
						c4 = 0.002232, 
						Ec = 10.9, 
						D = 0.89)


"""

"""
yPrediction16O17O = lambda E: empiricalFusion34(E, 
						b1 = 70, 
						b2 = -1.15, 
						b3 =  -0.01137, 
						c1 =  4, 
						c2 = -0.09, 
						c3 = -0.045, 
						c4 = 0.002232, 
						Ec = 10.871, 
						D = 0.89)

"""


r16O17O = ReactionPlot('16O', '17O', crossSection = True)
r16O17O.noSigma()
params, errors = r16O17O.fit(empiricalFusion34, p0 = 		[70, -1.15, -0.01137, 4.0, -0.155, -0.046, 0.002232, 10.9, 0.89], 
								bounds = (	[69.5, -1.2,   -0.015, 3.9,  -0.16,  -0.05, 0.002,    10.8, 0.88],
											[70.5, -1.0,    -0.010, 4.1,  -0.14, -0.04, 0.004,    11.0, 0.90]),  color = 'red', label = 'empirical-fit')

plt.close()

r16O17O = ReactionPlot('16O', '17O', crossSection = True)
O16O17Yakovlev = Yakovlev(Z1=8, Z2=8, A1=15.99491461956, A2=17)
plotFuncOverData(O16O17Yakovlev.sFactor, r16O17O.E, r16O17O.S, label = 'potential-Yakovlev')
#plotFuncOverData(yPrediction16O17O,  r16O17O.E,  r16O17O.S, label = 'empirical-Yakovlev', color = 'green')
plotFuncOverData(lambda E : empiricalFusion34(E, *params),  r16O17O.E,  r16O17O.S, label = 'empirical-fit', color = 'red')
plt.yscale('log')
r16O17O.save()
r16O17O.show()

"""

yPrediction16O16O = lambda E: empiricalFusion34(E, 
						b1 = 65, 
						b2 = -0.6, 
						b3 = -0.03, 
						c1 = 2.485, 
						c2 = 0.3363, 
						c3 = -0.080, 
						c4 = 0.002830, 
						Ec = 11.079, 
						D = 0.88)

"""
# 66 - 70
# -1.5 - 


yPrediction16O16O = lambda E: empiricalFusion34(E, 
						b1 = 60.932, 
						b2 = -0.4236, 
						b3 = -0.01018, 
						c1 = 2.485, 
						c2 = 0.3363, 
						c3 = -0.09320, 
						c4 = 0.002830, 
						Ec = 11.079, 
						D = 0.88)


rp16O16O = ReactionPlot('16O', '16O')
rp16O16O.noSigma()
rp16O16O.excludeRef(['SaoPablo'])
params, errors = rp16O16O.fit(empiricalFusion34, p0 = 		[65, -0.6,  -0.03, 2.485,  0.3363, -0.080, 0.002830, 11.079, 0.88], 
								bounds = (	[64.9, -0.61, -0.032, 2.48,  0.32,  -0.081, 0.0028,    10.9, 0.87],
											[65.1, -0.59, -0.028, 2.49,  0.34, -0.079, 0.003,    11.1, 0.89]),  color = 'red', label = 'empirical-fit')

plt.close()

rp16O16O = ReactionPlot('16O', '16O')
O16Yakovlev = Yakovlev(Z1=8, Z2=8, A1=15.99491461956, A2=15.99491461956)
plotFuncOverData(O16Yakovlev.sFactor, rp16O16O.E, rp16O16O.S, label = 'potential-Yakovlev')
plotFuncOverData(yPrediction16O16O,  rp16O16O.E,  rp16O16O.S, label = 'empirical-Yakovlev', color = 'green')
plotFuncOverData(lambda E: empiricalFusion34(E, *params),  rp16O16O.E,  rp16O16O.S, label = 'empirical-fit', color = 'red')
plt.ylim(1e22, 1e26)
plt.xlim(6, 13)
plt.yscale('log')
rp16O16O.save()
rp16O16O.show()



"""
yPrediction12C16O = lambda E: empiricalFusion34(E, 
						b1 = 49, 
						b2 = -0.5158, 
						b3 = -0.01195, 
						c1 = 7.454, 
						c2 = -1.4683, 
						c3 = 0.04019, 
						c4 = -0.000005, 
						Ec = 8.998, 
						D = 1.00)


"""
yPrediction12C16O = lambda E: empiricalFusion34(E, 
						b1 = 47.541, 
						b2 = -0.4158, 
						b3 = -0.01195, 
						c1 = 7.454, 
						c2 = -1.3683, 
						c3 = 0.04019, 
						c4 = -0.000005, 
						Ec = 8.998, 
						D = 1.00)



r12C16O = ReactionPlot('12C', '16O')
C12O16Yakovlev = Yakovlev(Z1=6, Z2=8, A1=12, A2=15.99491461956)
plotFuncOverData(C12O16Yakovlev.sFactor, r12C16O.E, r12C16O.S, label = 'potential-Yakovlev')
plotFuncOverData(yPrediction12C16O,  r12C16O.E,  r12C16O.S, label = 'empirical-Yakovlev', color = 'green')
plt.yscale('log')
r12C16O.fit(empiricalFusion34, p0 = 		[47.541, -0.4158,  -0.01195, 7.454, -1.3683, 0.04019, -0.000005, 8.998, 1.00], 
								bounds = (	[45, -0.6,  -0.02, 7.4, -np.inf, 0, -1e-5, 8.9, 0.99],
											[49, -0.2,  -0.01, 7.6, 0, np.inf, 1e-5, 9.0, 1.01]),  color = 'red', label = 'empirical-fit')
#r12C16O.fit(empiricalFusion34, p0 = [47.541, -0.4158,  -0.01195, 7.454, -1.3683, 0.04019, -0.000005, 8.998, 1.00], color = 'red', label = 'empirical-fit')
r12C16O.save()
r12C16O.show()


yPrediction12C12C = lambda E: empiricalFusion34(E, 
						b1 = 37.333, 
						b2 = -0.4065, 
						b3 = -0.0137, 
						c1 = 4.881, 
						c2 = -1.1909, 
						c3 = 0.0418, 
						c4 = 0.00014, 
						Ec = 7.134, 
						D = 0.94)


"""
yPrediction12C12C = lambda E: empiricalFusion34(E, 
						b1 = 38.5, 
						b2 = -0.4065, 
						b3 = -0.0137, 
						c1 = 4.881, 
						c2 = -1.1909, 
						c3 = 0.0418, 
						c4 = 0.00014, 
						Ec = 7.134, 
						D = 0.94)

"""
rp12C12C = ReactionPlot('12C', '12C')
C12Yakovlev = Yakovlev(Z1=6, Z2=6, A1=12, A2=12)
plotFuncOverData(C12Yakovlev.sFactor, rp12C12C.E, rp12C12C.S, label = 'potential-Yakovlev')
plotFuncOverData(yPrediction12C12C,  rp12C12C.E,  rp12C12C.S, color = 'green', label = 'empirical-Yakovlev')
rp12C12C.fit(lambda E, b1, b2, b3: empiricalFusion34(E, b1, b2, b3, 4.881, -1.1909, 0.0418, 0.00014, 7.134, 0.94), p0 = [38.5, -0.4065, -0.0137], color = 'red', label = 'empirical-fit', Emin = 2.25)
plt.yscale('log')
rp12C12C.save()
rp12C12C.show()