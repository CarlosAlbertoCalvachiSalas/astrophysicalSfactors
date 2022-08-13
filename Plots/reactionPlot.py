import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import regex
import sfactors.Resonant.single
from scipy.optimize import curve_fit

from sfactors.Utils.formulas import sfactorFromCrossSection
from sfactors.Utils.strings import decodeNucleusStr
from sfactors.Databases.Reader.nuclearData import Nucleus

from sfactors import FILE_PATH

import os

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "sans-serif"
})

def plotFromDataFrameVect(dataFrames, r1, r2, r3, r4, Eaxis = 'lin', Saxis = 'lin', Eunit = 'MeV', Sunit = 'MeV-b', color = 'black'):
	np.vectorize(plotFromDataFrame, otypes = [np.dtype('O')], cache = True, 
		excluded = ['r1', 'r2', 'r3', 'r4', 'Eaxis', 'Saxis', 'Eunit', 'Sunit', 'color'])(**{
		'dataFrame': dataFrames,
		'r1':r1, 'r2':r2, 'r3':r3, 'r4':r4, 
		'Eaxis':Eaxis, 'Saxis':Saxis, 
		'Eunit':Eunit, 'Sunit':Sunit, 
		'color':color})


def plotFromDataFrame(dataFrame, r1, r2, r3, r4, Eaxis = 'lin', Saxis = 'lin', Eunit = 'MeV', Sunit = 'MeV-b', color = 'black'):

	plotASeries(r1, r2, r3, r4, dataFrame.iloc[:, 0],dataFrame.iloc[:, 1], dataFrame.iloc[:, 2], dataFrame.iloc[0, 3], dataFrame.iloc[0, 4], Eaxis , Saxis , Eunit, Sunit, color)


def plotASeries(r1, r2, r3, r4, E, S, sigma, marker, reference, Eaxis = 'lin', Saxis = 'lin', Eunit = 'MeV', Sunit = 'MeV-b', color = 'black'):
	
	plt.errorbar(E, S, sigma, color = color, linestyle = 'None', capsize = 2, marker = marker, label = reference)
	plt.xlabel('E({})'.format(Eunit))
	plt.ylabel('S({})'.format(Sunit))

	if(Eaxis != 'lin'):
		plt.xscale('log')

	if(Eaxis != 'lin'):
		plt.yscale('log')	

	if(r3 and r4):
		plt.title('Astrophysical S-factor energy curve for the ${}({},{}){}$ reaction'.format(r1, r2, r3, r4))
	else:
		plt.title('Astrophysical S-factor energy curve for the ${}+{}$ reaction'.format(r1, r2))

def plotFuncOverData(f, E, S, N = 10000, color = 'blue', label = None):
	minE = np.min(E)
	maxE = np.max(E)

	Elin = np.linspace(minE, maxE, N)

	y = np.vectorize(f)(Elin)

	plt.plot(Elin, y, color = color, label = label)

def plotFitFunc(f, E, S, p0 = None, sigma = None, Emin = None, Emax = None, N = 10000, label = None, color = 'blue'):

	params, conv = curve_fit(f, E, S, p0, sigma)
	print(params)
	EminLin = np.min(E)
	EmaxLin = np.max(E)

	if(Emin):
		EminLin = Emin

	if(Emax):
		EmaxLin = Emax

	Elin = np.linspace(EminLin, EmaxLin, N)

	plt.plot(Elin, f(Elin, *params), color = color, label = label)
	
	return params, conv

class ReactionPlot():
	def __init__(self, r1, r2, r3 = None, r4 = None, crossSection = False):
		
		self.r1 = Nucleus(r1)
		self.r2 = Nucleus(r2)
		self.r3 = Nucleus(r3)
		self.r4 = Nucleus(r4)
	
		self.link = FILE_PATH + 'Databases/ReactionData/' + r1 + r2 + '.csv'
		#print(self.link)
		try:
			self.dataFrame = pd.read_csv(self.link)	
			self.E = self.dataFrame.iloc[:, 0]
			self.S = self.dataFrame.iloc[:, 1]
			self.sigma = self.dataFrame.iloc[:, 2]
			#print(self.E)
		except Exception as e:
			print('No database available')
		
		self.groupData()

		if(crossSection):
			self.S = sfactorFromCrossSection(crossSection=self.S, 
																			E = self.E,
																			Z1 = self.r1.Z, 	
																			Z2 = self.r2.Z, 
																			A1 = self.r1.A,
																			A2 = self.r2.A)

			self.groups[0].iloc[:, 1] = self.S

		
		self.plot()

	def groupData(self):
		self.groups = self.dataFrame.groupby(by=["Reference"])

		self.groups = np.array(self.groups, dtype = np.dtype('O'))[:,1]

	def plot(self, Eaxis = 'lin', Saxis = 'lin', Eunit = 'MeV', Sunit = 'MeV-b', color = 'black'):
		
		r1 = self.r1.__str__(latex = True)
		r2 = self.r2
		r3 = self.r3
		r4 = self.r4

		if(r2.synonym):
			r2 = r2.synonym
			if(len(r2) > 1):
				r2 = '\\' + r2
		else:
			r2 = r2.__str__(latex = True)

		if(r3 and r4):

			if(r3.synonym):
				r3 = r3.synonym
				if(len(r3) > 1):
					r3 = '\\' + r3
			else:
				r3 = r3.__str__(latex = True)

			r4 = r4.__str__(latex = True)


		plotFromDataFrameVect(**{
			'dataFrames': self.groups, 
			'r1':r1, 'r2':r2, 'r3':r3, 'r4':r4, 
			'Eaxis':Eaxis, 'Saxis':Saxis, 
			'Eunit':Eunit, 'Sunit':Sunit, 
			'color':color
		})



	def fit(self, f, p0 = None, Emin = None, Emax = None, Smin = None, Smax = None , N = 10000, color = 'blue', label = None):
		
		indexes = np.ones(len(self.E), dtype = int) == 1

		if(Emin):
			indexes = indexes & (self.E >= Emin)
			
		if(Emax):
			indexes = indexes & (self.E <= Emax)

		if(Smin):
			indexes = indexes & (self.S >= Smin)

		if(Smax):
			indexes = indexes & (self.S <= Smax)

		if(np.sum(self.sigma) == 0):
			sigmaFit = None
		else:
			sigmaFit = self.sigma[indexes]

		params, conv = plotFitFunc(f, self.E[indexes], self.S[indexes], p0 = p0, sigma = sigmaFit, Emin = Emin, Emax = Emax, N = N, color = color, label = label)

		logE = np.log10(np.sort(self.E))
		logS = np.log10(np.sort(self.S))

		if(np.max(logE) - np.min(logE) > 1.5):
			plt.xscale('log')

		if(np.max(logS) - np.min(logS) > 1.5):
			plt.yscale('log')

		self.residuals = self.S[indexes] - f(self.E[indexes], *params)

		return self.residuals


	def show(self):
		plt.legend()
		plt.show()

	def save(self, close = False):
		plt.legend()
		if(self.r2.synonym):
				r2 = self.r2
		else:
				r2 = self.r2.__str__()

		if(self.r3 and self.r4):
			if(self.r3.synonym):
				r3 = self.r3
			else:
				r3 = self.r3.__str__()

			plt.savefig(FILE_PATH + 'ReconstructedPlots/{}({},{}){}.eps'.format(self.r1, r2, r3, self.r4), format = 'eps')
		else:
			plt.savefig(FILE_PATH + 'ReconstructedPlots/{}{}.eps'.format(self.r1, r2), format = 'eps')

		if(close):
			plt.close()

