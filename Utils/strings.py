import numpy as np
import re
import pandas as pd

synonyms = {
	'p': [1, 'H'],
	'd': [2, 'H'],
	't': [3, 'H'],
	'alpha': [4, 'He'],
	'α': [4, 'He'],
	'gamma': [0, 'gamma'],
	'γ': [0, 'gamma'],
}


dataFrameSynonyms = pd.DataFrame(synonyms, index = ['A', 'name']).transpose()

def decodeNucleusStr(str):

	if(type(str) == type(None)):
		return None 
		
	regex = r'\s*(\d+)\s*([a-zA-Z]+)\s*|\s*([a-zA-Z]+)\s*(\d+)\s*'

	query = re.match(regex, str)

	if(type(query) == type(None)):

		if(type(dictionaryCheck(str, synonyms)) != type(None)):
			return dictionaryCheck(str, synonyms)
		else:
			return None

		return None
	else:
		g1, g2, g3, g4 = query.groups()

		if(g1 and g2):
			return int(g1), g2
		else:
			return int(g4), g3


def getSynonym(str):

	if(type(str) == type(None)):
		return None
	else:
		query = decodeNucleusStr(str)

		if(type(query) != type(None)):
			A, name = query

			selection = dataFrameSynonyms[(dataFrameSynonyms.A == A) & (dataFrameSynonyms.name == name)]

			if(len(selection) > 0):
				choices = pd.Series(selection.index)
				synonym = sorted(choices.values, key = lambda x: -len(x))[0]

				return synonym		

			else:
				return None
		else:
			return None

def dictionaryCheck(str, dictionary):

	#regex =  r'\s*(' + '|'.join(list(dictionary.keys())) + ')' + r'\s*'

	regex = r'^[^a-zA-Z0-9]*(' +  '|'.join(list(dictionary.keys())) + ')' + r'[^a-zA-Z0-9]*$'
	#print(regex)
	query = re.match(regex, str)

	if(type(query) == type(None)):
		return None
	else:

		character = query.groups()[0]
		
		return tuple(dictionary[character])

"""
print(getSynonym('gamma'),
	getSynonym('4He'),
	getSynonym('1H'),
	getSynonym('2H'),
	getSynonym('3H'), 
	getSynonym('6Li'),
 sep = '\n')
"""


"""
print(dictionaryCheck('p', synonyms),
		dictionaryCheck(' dd ', synonyms),
		dictionaryCheck(' t', synonyms),
		dictionaryCheck(' alpha', synonyms),
		dictionaryCheck('α', synonyms),
		dictionaryCheck(' gamma', synonyms),
		dictionaryCheck('γ ', synonyms), sep = '\n')

"""

"""
print(decodeNucleusStr(' 2 N '),
decodeNucleusStr(' 7 7 Yb'),
decodeNucleusStr(' Yb  77 '),
decodeNucleusStr(' 13 C *'),
decodeNucleusStr('p'),
decodeNucleusStr('alpha'),
decodeNucleusStr(' gamma'), sep = '\n')
"""
