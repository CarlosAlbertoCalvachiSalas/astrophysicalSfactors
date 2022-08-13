import numpy as np
import pandas as pd
#from bs4 import BeautifulSoup
import re

import requests

baseUrl = 'https://wwwndc.jaea.go.jp/cgi-bin/nucltab14?'

massRegexString = r'(\d\.\d+)\s*\+-\s*(\d\.\d+)'
isotopeRegexString = r'([a-zA-Z]+)\-\s*(\d+)\s*<\/a>'
regexString = r'([a-zA-Z]+)\-\s*(\d+)\s*<\/a>\s*(\d+\.\d+)\s*\+-\s*(\d+\.\d+)'


def getContent(url):
	try:
		response = requests.get(url)
	except Exception as e:
		return 'error'

	return response.text

def htmlParser(text):

	searchData = re.findall(regexString, text)
		
	names = np.array([])
	massNumbers = np.array([])
	masses = np.array([])
	errors = np.array([])

	if(len(searchData) > 0):
		names, massNumbers, masses, errors = np.array(searchData, dtype = np.dtype('O')).T 

		names = np.array(names, dtype = np.dtype('U'))
		massNumbers = np.array(massNumbers, dtype = int)
		masses = np.array(masses, dtype = float)
		errors = np.array(errors, dtype = float)

	return names, massNumbers, masses, errors

	"""
	searchMassData = re.findall(massRegexString, text)
	searchIsotopeData = re.findall(isotopeRegexString, text)

	names, massNumbers = np.array(searchIsotopeData, dtype = np.dtype('O')).T 
	masses, errors = np.array(searchMassData, dtype = float).T

	massNumbers = np.unique(np.array(massNumbers, dtype = int))
	names = np.array(names, dtype = np.dtype('U'))[:len(massNumbers)]
	
	if(len(errors) >= len(names)):
		masses = masses[len(errors)-len(names):]
		errors = errors[len(errors)-len(names):]
	else:


	minimumLength = np.min([len(masses), len(errors), len(massNumbers), len(names)])

	return names[:minimumLength], massNumbers[:minimumLength], masses[:minimumLength], errors[:minimumLength]

	"""

def extractNucleusData(Z):

	url = baseUrl + str(Z)

	names, massNumbers, masses, errors = htmlParser(getContent(url))
	print(names, massNumbers, masses, errors, sep = '\n')
	data = pd.DataFrame({'name': names, 'Z': Z*np.ones(len(massNumbers), dtype = int),'A': massNumbers, 'm': masses, 'error': errors})

	return data

def extractNucleusDataVect(Zarray):

	return np.vectorize(extractNucleusData, otypes = [np.dtype('O')])(Zarray)

allData = pd.DataFrame(columns = ['name', 'Z', 'A', 'm', 'error'])

allData = allData.append(list(extractNucleusDataVect(np.arange(1, 119))))

allData.to_csv('nuclearData.csv', index = False)