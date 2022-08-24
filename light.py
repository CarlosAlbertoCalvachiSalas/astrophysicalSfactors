import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Resonant import * 
from sfactors.nonResonant import * 
from sfactors.Plots.reactionPlot import *


rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(exponential1, p0 = [0.7E-1, 0.7], color = 'red', label = 'exp1')
rp2Hd.fit(exponential2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'exp2')
rp2Hd.fit(exponential3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'exp3')
rp2Hd.fit(exponential4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'exp4')
rp2Hd.fit(exponential(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'exp5', color = 'orangered')
#plotFuncOverData(lambda E: exponential1(E,0.7E-1, 0.7), rp2Hp.E, rp2Hp.S)
#rp2Hd.excludeRef('RA02')
#rp2Hd.fit(exponential(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'exp5', color = 'dodgerblue')
#rp2Hd.excludeRef(['RA02', 'GR95'])
#rp2Hd.fit(exponential(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'exp5', color = 'mediumseagreen')
plt.xscale('log')
plt.yscale('log')
rp2Hd.save()
rp2Hd.show()

rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(polynomial1, p0 = [0.7E-1, 0.7], color = 'red', label = 'poly1')
rp2Hd.fit(polynomial2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'poly2')
rp2Hd.fit(polynomial3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'poly3')
rp2Hd.fit(polynomial4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'poly4', color = 'orangered')
rp2Hd.fit(polynomial(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'poly5', color = 'midnightblue')
#plotFuncOverData(lambda E: exponential1(E,0.7E-1, 0.7), rp2Hp.E, rp2Hp.S)
plt.xscale('log')
plt.yscale('log')
rp2Hd.save()
rp2Hd.show()


rp2Hd = ReactionPlot('2H','d','p','3H')
#rp2Hd.fit(exponential1, p0 = [0.7E-1, 0.7], color = 'red', label = 'exp1')
#rp2Hd.fit(exponential2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'exp2')
#rp2Hd.fit(exponential3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'exp3')
#rp2Hd.fit(exponential4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'exp4')
rp2Hd.fit(exponential(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'exp5', color = 'orangered')
params, conv = rp2Hd.fit(polynomial(5), label = 'poly5', color = 'navy')

plotFuncOverData(lambda E: screenedFunction(
						E=E, 
						f=lambda E: polynomial(5)(E, *params),
						S=1,
						Ue=309E-6, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), rp2Hd.E, rp2Hd.S,  color = 'mediumspringgreen')
rp2Hd.excludeRef(['RA02', 'GR95'])
params, conv = rp2Hd.fit(polynomial(5), label = 'poly5', color = 'magenta')
plotFuncOverData(lambda E: screenedFunction(
						E=E, 
						f=lambda E: polynomial(5)(E, *params),
						S=1,
						Ue=309E-6, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), rp2Hd.E, rp2Hd.S,  color = 'mediumvioletred')

#plotFuncOverData(lambda E: polynomial(5)(E, 5.5325E-02,0.18293,0.28256,0.62121,0.44865,0.61893),rp2Hd.E, rp2Hd.S, color = 'pink',  label = 'poly5')
#plotFuncOverData(lambda E: exponential1(E,0.7E-1, 0.7), rp2Hp.E, rp2Hp.S)
plt.xscale('log')
plt.yscale('log')
rp2Hd.save()
rp2Hd.show()

rp2Hp = ReactionPlot('2H','p','gamma','3He')
rp2Hp.fit(polynomial(1), p0 = [ 4.34462668e-01,1.38914847e+00], label = 'poly1', color = 'red')
rp2Hp.fit(polynomial(2), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00], label = 'poly2', color = 'green')
rp2Hp.fit(polynomial(3), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00], label = 'poly3', color = 'blue')
rp2Hp.fit(polynomial(4), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00, 1.68928099e-01], label = 'poly4', color = 'orangered')
rp2Hp.fit(polynomial(5), label = 'poly5', color = 'midnightblue')
plt.xscale('log')
plt.yscale('log')
#rp2Hd.save()
rp2Hp.show()

rp2Hp = ReactionPlot('2H','p','gamma','3He')
rp2Hp.fit(polynomial(1), p0 = [ 4.34462668e-01,1.38914847e+00], bounds = ([1e-5, -10], [1, 10]), label = 'poly1-c', color = 'red')
rp2Hp.fit(polynomial(2), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00], bounds = ([1e-5, -10, -10], [1, 10, 10]), label = 'poly2-c', color = 'green')
rp2Hp.fit(polynomial(3), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00], bounds = ([1e-5, -10, -10, -10], [1, 10, 10, 10]), label = 'poly3-c', color = 'blue')
rp2Hp.fit(polynomial(4), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00, 1.68928099e-01], bounds = ([1e-5, -10, -10, -10, -10], [1, 10, 10, 10, 10]), label = 'poly4-c', color = 'orangered')
rp2Hp.fit(polynomial(5),  bounds = ([1e-3, -10, -10, -10, -10, -10], [1, 10, 10, 10, 10, 10]), color = 'midnightblue',  label = 'poly5-c')
plt.xscale('log')
plt.yscale('log')
#rp2Hd.save()
rp2Hp.show()


rp2Hp = ReactionPlot('2H','p','gamma','3He')
rp2Hp.fit(exponential1, p0 = [0.2, 3], color = 'red', label = 'exp1')
rp2Hp.fit(exponential2, p0 = [0.2, 3, -0.2], color = 'green', label = 'exp2')
#rp2Hp.fit(exponential3, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428], color = 'brown', label = 'exp3-lowE', Emax = 1)
rp2Hp.fit(exponential3, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428], color = 'blue', label = 'exp3')
#rp2Hp.fit(exponential3, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428], color = 'blue', label = 'exp3-higE')
rp2Hp.fit(exponential(4), p0 = [2.75369389,1.33358604,-0.18065048,0.00865423, 0], label = 'exp4', color = 'dodgerblue')
rp2Hp.fit(exponential(5), label = 'exp5', p0 = [ 2.75369389,1.33358604,-0.18065048,0.00865423, 0, 0], color = 'darkviolet')
#rp2Hp.fit(exponential4, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428, 0.05], label = 'exp4', Emax = 0.5)
plt.xscale('log')
plt.yscale('log')
rp2Hp.save()
rp2Hp.show()

rp2Hp = ReactionPlot('2H','p','gamma','3He')

rp2Hp.fit(exponential(5), label = 'exp5', p0 = [ 2.75369389,1.33358604,-0.18065048,0.00865423, 0, 0], color = 'darkviolet')
rp2Hp.fit(polynomial(4), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00, 1.68928099e-01], bounds = ([1e-5, -10, -10, -10, -10], [1, 10, 10, 10, 10]), color = 'orangered', label = 'poly4')
rp2Hp.fit(polynomial(5), label = 'poly5', color = 'midnightblue')
#plotFuncOverData(lambda E: exponential3(E,  0.275367003, 1.33359118, -0.14065138, 0.01865428), rp2Hp.E, rp2Hp.S)

#rp2Hp.fit(exponential(5), label = 'exp5', p0 = [ 2.75369389,1.33358604,-0.18065048,0.00865423, 0, 0], color = 'darkorange', Emax = 1)
#rp2Hp.fit(polynomial(5), label = 'poly5', color = 'olivedrab',  Emax = 1)


plt.xscale('log')
plt.yscale('log')
rp2Hp.show()


