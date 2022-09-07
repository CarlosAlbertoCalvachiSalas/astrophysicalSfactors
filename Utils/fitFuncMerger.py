import numpy as np
import inspect

def removeFirstElement(array):

	return array[1:]

def replaceStrVect(string, toReplaceAr, replacementAr):

	if(len(toReplaceAr) == 0):
		return string 
	else:
		string = string.replace(toReplaceAr[0], replacementAr[0])
		return replaceStrVect(string, toReplaceAr[1:], replacementAr[1:])

#print(removeFirstElement([1, 2, 3]))

class MergeFitFuncs():

	def __init__(self, fs):

		functions = np.array(fs, dtype = np.dtype('O'))

		if(functions.ndim > 0):
			self.fs = functions
		else:
			self.fs = np.array([fs], dtype = np.dtype('O'))

		self.fNum = len(self.fs)

		self.fullArgs = np.vectorize(inspect.getfullargspec, 
									otypes = [np.dtype('O')])(self.fs)

		self.argsObj  = np.vectorize(getattr, 
									otypes = [np.dtype('O')])(self.fullArgs, 'args')

		self.argsComma = np.vectorize(str.join, otypes = [str])(',', self.argsObj)

		self.arg0 	  = np.vectorize(np.take, otypes = [str])(self.argsObj, [0])[0]
		
		self.othersArgs = np.concatenate(np.vectorize(removeFirstElement, otypes = [np.dtype('O')])(self.argsObj))

		self.args = np.append(self.arg0, self.othersArgs)
		argUnique, indexes = np.unique(self.args, return_index = True)
		self.args = self.args[np.sort(indexes)]

		self.setSignatures()

	def setSignatures(self):

		toReplaceAr = np.arange(self.fNum).astype(str)
		self.toReplaceAr = np.char.add('#', toReplaceAr)

		replacementAr = np.arange(self.fNum).astype(str)
		replacementAr = np.char.add(replacementAr, '](')
		replacementAr = np.char.add('fs[', replacementAr)
		replacementAr = np.char.add(replacementAr, self.argsComma)
		self.replacementAr = np.char.add(replacementAr, ')')
		
		self.paramsStr = ', '.join(self.args)

	def funcStr(self, string):

		expression = 'lambda ' + self.paramsStr + ': ' + replaceStrVect(string, self.toReplaceAr, self.replacementAr)
		
		global fs 
		fs = self.fs

		return eval(expression)

	def sum(self):
		return self.funcStr(' + '.join(self.toReplaceAr))

	def prod(self):
		return self.funcStr(' * '.join(self.toReplaceAr))
"""
f1 = lambda E, x: E + x + 1
f2 = lambda E, y, a, b: E * y + 2 + a + b
f3 = lambda E, ab, abb: E + ab * abb

m1 = MergeFitFuncs(f1)
m2 = MergeFitFuncs([f1, f2])
m3 = MergeFitFuncs([f1, f2, f3])
"""

"""
print(m1.arg0, m2.arg0, m3.arg0, sep = '\n')
print(m1.argsObj, m2.argsObj, m3.argsObj, sep = '\n')
print(m1.argsComma, m2.argsComma, m3.argsComma, sep = '\n')
print(m1.othersArgs, m2.othersArgs, m3.othersArgs, sep = '\n')
print(m1.args, m2.args, m3.args, sep = '\n')
"""

#print(m1.funcStr('2*#0'), m2.funcStr('#0+#1'), m3.funcStr('(#0+#1)*(#2)') )

#m1.sum(), m2.sum(), m3.sum()
#m1.prod(), m2.prod(), m3.prod()
