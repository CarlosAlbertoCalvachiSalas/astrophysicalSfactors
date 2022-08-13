import numpy as np
import pandas as pd 
from sfactors import FILE_PATH
from sfactors.Utils.strings import decodeNucleusStr, getSynonym

nuclearData = pd.read_csv(FILE_PATH + 'Databases/NuclearData/nuclearData.csv')

atomicNames 	= nuclearData.iloc[:, 0]
atomicNumbers 	= nuclearData.iloc[:, 1]
massNumbers   	= nuclearData.iloc[:, 2]

def getNuclearDataByNumbers(Z, A):

	index = np.arange(len(nuclearData))[(atomicNumbers == Z) & (massNumbers == A)]

	if(len(index) > 0):

		index = index[0]

		return nuclearData.iloc[index]

	else:
		return None

def getNuclearDataByStr(str):

	decoder = decodeNucleusStr(str)

	if(type(decoder) != type(None)):

		A, strZ = decoder

		index = np.arange(len(nuclearData))[(atomicNames == strZ) & (massNumbers == A)]

		if(len(index) > 0):
			
			index = index[0]

			return nuclearData.iloc[index]
		else:
			
			return None
	else:
		
		return None

class Nucleus():
	def __init__(self, *args):

		if(len(args) == 1):
			query = getNuclearDataByStr(args[0])
			if(type(query) != type(None)):
				self.name, self.Z, self.A, self.m, self.errorM = query.values
				self.data = query
				self.synonym = getSynonym(args[0])
			elif(type(decodeNucleusStr(args[0])) != type(None)):
				decoded = decodeNucleusStr(args[0])
				self.name, self.Z, self.A, self.m, self.errorM = decoded[1], decoded[0], decoded[0], decoded[0], decoded[0]
				self.data = pd.Series([self.name, self.Z, self.A, self.m, self.errorM], index = nuclearData.columns)
				self.synonym = getSynonym('{}'.format(self.name))
			else:
				return None
		elif(len(args) == 2):
			query = getNuclearDataByNumbers(args[0], args[1])
			if(type(query) != type(None)):
				self.name, self.Z, self.A, self.m, self.errorM = query.values
				self.data = query
				self.synonym = getSynonym('{}{}'.format(self.A, self.name))
			elif(args[0] == 0 and args[1] == 0):
				decoded = decodeNucleusStr('gamma')
				self.name, self.Z, self.A, self.m, self.errorM = decoded[1], decoded[0], decoded[0], decoded[0], decoded[0]
				self.data = pd.Series([self.name, self.Z, self.A, self.m, self.errorM], index = nuclearData.columns)
				self.synonym = getSynonym('{}'.format(self.name))
			else:
				return None
		else:
			return None
	
	def __str__(self, latex = False):

		if(hasattr(self, 'A') and hasattr(self, 'name')):
			if(latex):
				return  r'\mathrm{{}^' + '{' + str(self.A) + '}' + self.name + '}'
			else:
				if(self.A > 0):
					return str(self.A) + self.name 
				else:
					return self.name
		else:
			return 'None'

	def __eq__(self, test):

		if(type(test) == type(None)):
			return not(hasattr(self, 'Z') or hasattr(self, 'A') or hasattr(self, 'name'))
		elif(isinstance(test, bool)):
			return (hasattr(self, 'Z') or hasattr(self, 'A') or hasattr(self, 'name'))
		elif(isinstance(test, Nucleus)):
			return (self.Z == test.Z) and (self.A == test.A) and (self.name == test.name)
		elif(isinstance(test, str)):
			return self.__eq__(Nucleus(test))
		elif(hasattr(test, 'dtype')):
			if('str' in test.dtype.name and test.ndim == 0):
				test = str(test)
				return self.__eq__(Nucleus(test))
			else:
				return False
		elif(hasattr(test, '__len__')):	
			if(np.array(test).ndim == 1 and len(test) >= 2):
				return self.__eq__(Nucleus(*np.array(test)[:2]))
		else:
			return False

	def __and__(self, test):

		if(hasattr(self, 'Z') or hasattr(self, 'A') or hasattr(self, 'name')):
			if(type(test) == type(None)):
				return False
			elif(isinstance(test, Nucleus)):
				return not(hasattr(test, 'Z') or hasattr(test, 'A') or hasattr(test, 'name'))
			else:
				return None
		else:
			return False

	def __bool__(self):
		return (hasattr(self, 'Z') or hasattr(self, 'A') or hasattr(self, 'name'))
#print(Nucleus('gamma'))


"""

if(Nucleus(1, 1)):
	print('Y')
else:
	print('N')


if(Nucleus(12, 1)):
	print('Y')
else:
	print('N')

"""

"""
print(Nucleus(1, 1 ) == np.array('2H'),
Nucleus(1, 2) == np.array('2H'),
Nucleus(0, 0) == np.array('2H'),
Nucleus(2, 4) == np.array('2H'),
Nucleus(92, 238) == np.array('2H'),
 sep = '\n')
"""

"""
print(Nucleus(1, 1 ),
Nucleus(1, 2),
Nucleus(0, 0),
Nucleus(2, 4),
Nucleus(92, 238), sep = '\n')
"""

"""
print(Nucleus(' 2 N '),
Nucleus(' 7 7 Yb'),
Nucleus(' Li 7 '),
Nucleus(' Yb  77 '),
Nucleus(' 13 C *'),
Nucleus('p'),
Nucleus('alpha'),
Nucleus(' gamma'), sep = '\n')
"""

"""
print(getNuclearDataByStr(' 2 N '),
getNuclearDataByStr(' 7 7 Yb'),
getNuclearDataByStr(' Yb  77 '),
getNuclearDataByStr(' 13 C *'),
getNuclearDataByStr('p'),
getNuclearDataByStr('alpha'),
getNuclearDataByStr(' gamma'), sep = '\n')
"""


"""
print(getNuclearDataByNumbers(1, 0),
	getNuclearDataByNumbers(1, 1),
	getNuclearDataByNumbers(8, 9),
	getNuclearDataByNumbers(1000, 1000), sep = '\n')
"""
