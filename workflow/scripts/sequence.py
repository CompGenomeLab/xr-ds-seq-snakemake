import os
import re
import itertools

kDNAcomplementaryDictionary = {
	'A':'T', 'T':'A', 'G':'C', 'C':'G',
	'a':'t', 't':'a', 'g':'c', 'c':'g',
	'n':'n', 'N':'N', 'x':'x', 'X':'X',
	'B':'V', 'V':'B', 'D':'H', 'H':'D',
	'b':'v', 'v':'b', 'd':'h', 'h':'d'
}

kDNAmotifDictionary = {
	'A': '(A)',
	'C': '(C)',
	'G': '(G)',
	'T': '(T)',
	'U': '(U)',
	'W': '(A|T)',
	'S': '(C|G)',
	'M': '(A|C)',
	'K': '(G|T)',
	'R': '(A|G)',
	'Y': '(C|T)',
	'B': '(C|G|T)',
	'D': '(A|G|T)',
	'H': '(A|C|T)',
	'V': '(A|C|G)',
	'N': '.'
	}

class Sequence:
	def __init__(self, input, up=True):
		self.sequence = input
		if up:
			self.sequence = self.sequence.upper()
		self.header = None

	def subseq(self, start, end):
		# Zero-based subsequence function
		return DNA(self.sequence[start : end])

	def motifCount(self, motif):
		count = start = 0
		while True:
			start = self.sequence.find(motif, start) + 1
			if start > 0:
				count+=1
			else:
				return count

	def assignHeader(self, header):
		self.header = header

	def getHeader(self):
		return self.header

	def getLength(self):
		return len(self.getSequence())

	def getSequence(self):
		return self.sequence


class DNA(Sequence):
	def reverseComplement(self, reverse=True):
		seq_dict = kDNAcomplementaryDictionary
		if reverse:
			string = reversed(self.sequence)
		else:
			string = self.sequence
		return DNA("".join([seq_dict[base] for base in string]), False)

class reMotif(object):
	def __init__(self, input, ignorecaseFlag=True):
		self.patternAsString = input
		if ignorecaseFlag:
			self.pattern = re.compile(self.patternAsString, re.IGNORECASE)
		else:
			self.pattern = re.compile(self.patternAsString)

	def getPatternAsString(self):
		return self.patternAsString

	def DNA_complement(self, reverseFlag = False):
		def complementIfPossibe(key, dictionary):
			if key in dictionary.keys():
				return dictionary[key]
			return key
		if reverseFlag: # Use with caution. It doesn't work for non-symmetric patterns like: .{3}ATG.{4} will give .{3}GTA.{4}
			string = reverse(self.getPatternAsString())
		else:
			string = self.getPatternAsString()
		return reMotif("".join([complementIfPossibe(base, kDNAcomplementaryDictionary) for base in string]))

	def reverse(self):
		replacements = {"[": "]", "]": "[", "(": ")", ")": "(", "{": "}", "}": "{"}
		reversedString = reversed(self.getPatternAsString())
		replacedString = "".join(replacements.get(c, c) for c in reversedString)
		return reMotif(replacedString)

	def getIndexList(self, subject):
		queryPattern = self.pattern
		matchObjects = re.finditer(queryPattern, subject)
		startList = []
		endList = []
		for matchObject in matchObjects:
			start = matchObject.start()
			end = matchObject.end()
			startList.append(start)
			endList.append(end)
		return startList, endList

class consensus(reMotif):
	def __init__(self, input):
		convertedInput = self.toRegex(input)
		reMotif.__init__(self, convertedInput)

	def toRegex(self, string):
		newString = ''
		for letter in string:
			upperLetter = letter.upper()
			if upperLetter in kDNAmotifDictionary.keys():
				newString += kDNAmotifDictionary[upperLetter]
			else:
				newString += upperLetter
		return newString



