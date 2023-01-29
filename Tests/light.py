import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sfactors.Resonant import * 
from sfactors.nonResonant import * 
from sfactors.Plots.reactionPlot import *
from sfactors.Utils.fitFuncMerger import MergeFitFuncs

rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(exponential1, p0 = [0.7E-1, 0.7], color = 'red', label = 'exp1')
rp2Hd.fit(exponential2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'exp2')
rp2Hd.fit(exponential3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'exp3')
rp2Hd.fit(exponential4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'exp4')
rp2Hd.fit(exponential(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'exp5', color = 'orangered')
plt.xscale('log')
plt.yscale('log')
rp2Hd.setTags(['exp12345'])
rp2Hd.save()
rp2Hd.show()


rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(polynomial1, p0 = [0.7E-1, 0.7], color = 'red', label = 'poly1')
rp2Hd.fit(polynomial2, p0 = [0.7E-1, 0.7, 1], color = 'green', label = 'poly2')
rp2Hd.fit(polynomial3, p0 = [0.7E-1, 0.7, -0.5, 0.1], color = 'brown', label = 'poly3')
rp2Hd.fit(polynomial4, p0 = [0.7E-1, 1, -0.7, 0.1, 1], label = 'poly4', color = 'orangered')
params, conv = rp2Hd.fit(polynomial(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'poly5', color = 'midnightblue')
plt.xscale('log')
plt.yscale('log')
rp2Hd.setTags(['poly12345'])
rp2Hd.save()
rp2Hd.show()

rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(exponential(5), p0 = [0.7E-1, 1, -0.7, 0.1, 1, 0], label = 'exp5', color = 'orangered')
rp2Hd.fit(polynomial(5), label = 'poly5', color = 'navy')
rp2Hd.excludeRef(['RA02', 'GR95'])
params, conv = rp2Hd.fit(polynomial(5), label = 'poly5-exclude', color = 'magenta')
rp2Hd.setTags(['exp5', 'poly5', 'best'])
plt.xscale('log')
plt.yscale('log')
rp2Hd.save()
rp2Hd.show()


rp2Hd = ReactionPlot('2H','d','p','3H')
rp2Hd.fit(lambda E, Ue: screenedFunction(
						E=E, 
						f=lambda E: polynomial(5)(E, *params),
						Ue=Ue, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), p0 = [309E-6], bounds = ([100E-6],[400E-6]),
						label = 'screening-all-fit', 
						color = 'royalblue')

#rp2Hd.includeRef(['RA02'])
"""
plotFuncOverData(lambda E: screenedFunction(
						E=E, 
						f=lambda E: polynomial(5)(E, *params),
						Ue=309E-6, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), rp2Hd.E, rp2Hd.S,  label = 'screening', color = 'mediumspringgreen')
"""
rp2Hd.includeRef(['RA02'])
rp2Hd.fit(lambda E, Ue: screenedFunction(
						E=E, 
						f=lambda E: polynomial(5)(E, *params),
						Ue=Ue, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), p0 = [309E-6], bounds = ([100E-6],[400E-6]),
						label = 'screening-only-RA02-fit', 
						color = 'blue')
#rp2Hd.fit(MergeFitFuncs([exponential(5), lambda E, Ue: screeningFactor(E, Ue, Z1=1, Z2=1, mu=rp2Hd.mu)]).prod(), p0=[0.0554494,0.25410999,-0.08265845,-0.00711975,0.01966242,-0.00417225, 309E-6],  Emax = 1E-2, label = 'screening-exclude', color = 'mediumvioletred')
rp2Hd.allRef()
plotFuncOverData(lambda E: screenedFunction(
						E=E, 
						f=lambda E: polynomial(5)(E, *params),
						Ue=309E-6, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), rp2Hd.E, rp2Hd.S, label = 'screening-Ue-theory', color = 'mediumvioletred')

plotFuncOverData(lambda E: screenedFunction(
						E=E, 
						f=lambda E: polynomial1(E, 43E-3, 0.54),
						Ue=309E-6, 
						Z1=1,
						Z2=1,
						mu=rp2Hd.mu), rp2Hd.E, rp2Hd.S, label = 'RA02-paper', color = 'sandybrown')

#plotFuncOverData(lambda E: polynomial(5)(E, 5.5325E-02,0.18293,0.28256,0.62121,0.44865,0.61893),rp2Hd.E, rp2Hd.S, color = 'pink',  label = 'poly5')
#plotFuncOverData(lambda E: exponential1(E,0.7E-1, 0.7), rp2Hp.E, rp2Hp.S)
plt.xscale('log')
plt.yscale('log')
rp2Hd.setTags(['screening', 'exp5', 'poly5', 'exclude'])
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
rp2Hp.setTags(['poly12345'])
rp2Hp.save()
rp2Hp.show()

rp2Hp = ReactionPlot('2H','p','gamma','3He')
rp2Hp.fit(polynomial(1), p0 = [ 4.34462668e-01,1.38914847e+00], bounds = ([1e-5, -10], [1, 10]), label = 'poly1-c', color = 'red')
rp2Hp.fit(polynomial(2), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00], bounds = ([1e-5, -10, -10], [1, 10, 10]), label = 'poly2-c', color = 'green')
rp2Hp.fit(polynomial(3), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00], bounds = ([1e-5, -10, -10, -10], [1, 10, 10, 10]), label = 'poly3-c', color = 'blue')
rp2Hp.fit(polynomial(4), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00, 1.68928099e-01], bounds = ([1e-5, -10, -10, -10, -10], [1, 10, 10, 10, 10]), label = 'poly4-c', color = 'orangered')
rp2Hp.fit(polynomial(5),  bounds = ([1e-3, -10, -10, -10, -10, -10], [1, 10, 10, 10, 10, 10]), color = 'midnightblue',  label = 'poly5-c')
plt.xscale('log')
plt.yscale('log')
rp2Hp.setTags(['poly12345', 'bounds'])
rp2Hp.save()
rp2Hp.show()


rp2Hp = ReactionPlot('2H','p','gamma','3He')
rp2Hp.fit(exponential1, p0 = [0.2, 3], color = 'red', label = 'exp1')
rp2Hp.fit(exponential2, p0 = [0.2, 3, -0.2], color = 'green', label = 'exp2')
rp2Hp.fit(exponential3, p0 = [0.275367003, 1.33359118, -0.14065138, 0.01865428], color = 'blue', label = 'exp3')
rp2Hp.fit(exponential(4), p0 = [2.75369389,1.33358604,-0.18065048,0.00865423, 0], label = 'exp4', color = 'dodgerblue')
rp2Hp.fit(exponential(5), label = 'exp5', p0 = [ 2.75369389,1.33358604,-0.18065048,0.00865423, 0, 0], color = 'darkviolet')
plt.xscale('log')
plt.yscale('log')
rp2Hp.setTags(['exp12345'])
rp2Hp.save()
rp2Hp.show()

rp2Hp = ReactionPlot('2H','p','gamma','3He')
rp2Hp.fit(exponential(5), label = 'exp5', p0 = [ 2.75369389,1.33358604,-0.18065048,0.00865423, 0, 0], color = 'darkviolet')
rp2Hp.fit(polynomial(4), p0 = [ 4.34462668e-01,1.38914847e+00,7.98034243e+00,-1.82914588e+00, 1.68928099e-01], bounds = ([1e-5, -10, -10, -10, -10], [1, 10, 10, 10, 10]), color = 'orangered', label = 'poly4')
rp2Hp.fit(polynomial(5), label = 'poly5', color = 'midnightblue')
plt.xscale('log')
plt.yscale('log')
rp2Hp.setTags(['selected', 'exp5', 'poly45'])
rp2Hp.save()
rp2Hp.show()


