import os
import re
import sys
import itertools
import random
from sequence import Sequence

def nucleotideAbundanceDict2percentageDict(theDict):
	newDict = {}
	for key in theDict.keys():
		newDict[key] = {}
		positionTotalCount = 0
		for letter in theDict[key].keys():
			positionTotalCount += theDict[key][letter]
		for letter in theDict[key].keys():
			newDict[key][letter] = float(theDict[key][letter])/positionTotalCount
	return newDict


class fasta:
	def __init__(self, input):
		self.file = input

	def stream(self, bufsize=4096):
		def chunk2seqDict(chunk):
			lines = chunk.split('\n')
			header = lines[0]
			del lines[0]
			sequence = ''.join(lines)
			seqObject = {'h': header, 's': sequence}
			return seqObject

		filein = open(self.file, 'r')
		delimiter = '\n>'
		buf = ''
		justStarted = True
		while True:
			newbuf = filein.read(bufsize)
			if not newbuf:
				yield chunk2seqDict(buf)
				return
			buf += newbuf
			sequenceChunks = buf.split(delimiter)
			for chunk in sequenceChunks[0:-1]:
				if justStarted and chunk.startswith('>'):
					chunk = chunk[1:]
					justStarted = False
				yield chunk2seqDict(chunk)
			buf = sequenceChunks[-1]

	def stream2(self, bufsize=4096):
		def chunk2seqDict(chunk):
			lines = chunk.split('\n')
			header = lines[0]
			del lines[0]
			sequence = ''.join(lines)
			seqObject = Sequence(sequence)
			seqObject.assignHeader(header)
			return seqObject

		filein = open(self.file, 'r')
		delimiter = '\n>'
		buf = ''
		justStarted = True
		while True:
			newbuf = filein.read(bufsize)
			if not newbuf:
				yield chunk2seqDict(buf)
				return
			buf += newbuf
			sequenceChunks = buf.split(delimiter)
			for chunk in sequenceChunks[0:-1]:
				if justStarted and chunk.startswith('>'):
					chunk = chunk[1:]
					justStarted = False
				yield chunk2seqDict(chunk)
			buf = sequenceChunks[-1]

	def getGenomeLength(self):
		totalLength = 0
		seqDicts = self.stream()
		for seqDict in seqDicts:
			sequenceLength = len(seqDict['s'])
			totalLength += sequenceLength
		return totalLength

	def getRandomIndexList(self, numberOfReads):
		genomeLength = self.getGenomeLength
		indexList = []
		for i in numberOfReads:
			value = random.random * genomeLength 
			

	def getSequenceCount(self):
		count = 0
		filein = open(self.file, 'r')
		for line in filein:
			if line.startswith('>'):
				count += 1
		return count

	def singleEntry2chromSize(self):
		d = self.read()
		if len(d.keys()) != 1:
			raise ValueError('Fasta file does not have a single entry')
		return 'chr\t' + str(len(d[d.keys()[0]]))

	def getFirstSeqLength(self):
		gtCount = 0
		sequence = ''
		filein = open(self.file, 'r')
		for line in filein:
			if line.startswith('>'):
				gtCount += 1
			else:
				sequence += line.strip()
			if gtCount == 2:
				break
		return len(sequence)

	def getMaxSeqLength(self):

		sequence = ''
		maxLength = 0
		filein = open(self.file, 'r')
		for line in filein:
			if line.startswith('>'):
				maxLength = max(maxLength, len(sequence))
				sequence = ''
			else:
				sequence += line.strip()
		else:
			maxLength = max(maxLength,len(sequence))
		return maxLength

	def getKmerAbundance(self, kmer, firstNletters=None):
		nucleotides = 'ATGC'
		subseqTupleList = itertools.product(nucleotides, repeat=int(kmer))
		subseqList = sorted(list(subseqTupleList))

		def subseqHasValidNucleotides(subseq):
			allowedLetters = 'ATCGatcg'
			for letter in subseq:
				if letter not in allowedLetters:
					return False
			return True

		def addSeqToPosDict(sequence, theDict):
			for i in range(min(firstNletters,len(sequence)) - (kmer - 1)):
				subseq = sequence[i : i + kmer]
				if subseqHasValidNucleotides(subseq):
						 theDict[i + 1][subseq.upper()] += 1

		if firstNletters == None:
			firstNletters = self.getMaxSeqLength()
		else:
			firstNletters = int(firstNletters)

		positionDict = {}
		for i in range(1, firstNletters + 1 - (kmer - 1)):
			positionDict[i] = {}
			for subseq in subseqList:
				positionDict[i][''.join(subseq)] = 0


		filein = open(self.file, 'r')
		seq = header = ''
		for line in filein:
			if line.startswith('>'):
				header = line.strip()
				addSeqToPosDict(seq, positionDict)
				seq = ''
			else:
				seq += line.strip()
		else:
			addSeqToPosDict(seq,positionDict)

		return positionDict

	def headerDict2kmerAbundance(self, headerDict, kmer, firstNletters=None):
		nucleotides = 'ATGC'
		print(kmer)
		subseqTupleList = itertools.product(nucleotides, repeat=int(kmer))
		subseqList = sorted(list(subseqTupleList))

		def subseqHasValidNucleotides(subseq):
			allowedLetters = 'ATCGatcg'
			for letter in subseq:
				if letter not in allowedLetters:
					return False
			return True

		def addSeqToPosDict(sequence, theDict):
			for i in range(min(firstNletters,len(sequence)) - (kmer - 1)):
				subseq = sequence[i : i + kmer]
				if subseqHasValidNucleotides(subseq):
						 theDict[i + 1][subseq.upper()] += 1

		if firstNletters == None:
			firstNletters = self.getMaxSeqLength()
		else:
			firstNletters = int(firstNletters)

		positionDict = {}
		for i in range(1, firstNletters + 1 - (kmer - 1)):
			positionDict[i] = {}
			for subseq in subseqList:
				positionDict[i][''.join(subseq)] = 0

		for header in headerDict:
			addSeqToPosDict(headerDict[header], positionDict)

		return positionDict


	def getNucleotideAbundance(self, firstNletters=None):

		if firstNletters == None:
			firstNletters = self.getMaxSeqLength()
		else:
			firstNletters = int(firstNletters)

		def letterIsValidNucleotide(letter):
			if re.match('A|T|G|C|a|t|g|c', letter):
				return True
			else:
				return False

		def addSeqToPosDict(sequence, theDict):
			for i in range(min(firstNletters,len(sequence))):
				letter = sequence[i]
				if letterIsValidNucleotide(letter):
						 theDict[i + 1][letter.upper()] += 1

		positionDict = {}
		for i in range(1, firstNletters + 1):
			positionDict[i] = {
				'A': 0,
				'T': 0,
				'G': 0,
				'C': 0
				}
		filein = open(self.file, 'r')
		seq = header = ''
		for line in filein:
			if line.startswith('>'):
				header = line.strip()
				addSeqToPosDict(seq, positionDict)
				seq = ''
			else:
				seq += line.strip()
		else:
			addSeqToPosDict(seq,positionDict)
		return positionDict

	def getNucleotidePercentages(self, abuDict):
		return nucleotideAbundanceDict2percentageDict(abuDict)

	def getKmerPercentages(self, abuDict):
		return nucleotideAbundanceDict2percentageDict(abuDict)


	def printSeqLength(self):
		filein = open(self.file)
		header = False
		for line in filein:
			if line.startswith('>'):
				if header:
					print(header + '\t' + str(sequenceLength))
				header = line.strip()[1:]
				sequenceLength = 0
			else:
				sequenceLength += len(line.strip())
		else:
			print(header + '\t' + str(sequenceLength))

	def read(self):
		filein = open(self.file)
		header = False
		myDict = {}
		sequence = ''		
		for line in filein:
			if line.startswith('>'):
				if header:
					myDict[header] = sequence
				header = line.strip()[1:]
				sequence = ''
			else:
				sequence += line.strip()
		else:
			myDict[header] = sequence
		return myDict

	def separateBySequenceLength(self):
		headerDict = self.read()
		lengthDict = {}
		for header in headerDict.keys():
			sequence = headerDict[header]
			sequenceLength = len(sequence)
			if not sequenceLength in lengthDict.keys():
				lengthDict[sequenceLength] = {}
				lengthDict[sequenceLength][header] = sequence
			else:
				lengthDict[sequenceLength][header] = sequence
		return lengthDict

	def separateByLengthAndWriteKmerAbundance(self, kmer, lengthRange, output, percentage=False):
		lengthRangeStart = int(lengthRange.split('-')[0])
		lengthRangeEnd = int(lengthRange.split('-')[1])
		lengthDict = self.separateBySequenceLength()
		os.system('rm ' + output)
		for i in range(lengthRangeStart, lengthRangeEnd):
			length = i
			headerDict = lengthDict[length]
			kmerAbundance = self.headerDict2kmerAbundance(headerDict, kmer)
			self.writeKmerAbundanceTable(kmerAbundance, output, percentage, True)

	def writeKmerAbundanceTable(self, kmerAbundanceDict, output, percentage=False, append=False):
		if append:
			writeMode = 'a'
		else:
			writeMode = 'w'
		if percentage:
			kmerAbundanceDict = self.getKmerPercentages(kmerAbundanceDict)
		separator = '\t'
		out = open(output, writeMode)
		positions = kmerAbundanceDict.keys()
		out.write('kmer')
		for position in positions:
			out.write(separator + str(position))
		out.write("\n")
		for kmer in sorted(kmerAbundanceDict[1].keys()):
			out.write(kmer)
			for position in positions:
				out.write(separator + str(kmerAbundanceDict[position][kmer]))
			out.write("\n")
		return 1

	def writeKmerAbundanceMeltedData(self, kmerAbundanceDict, output, dictionary = {}, percentage=False, append=False, headerFlag=False):
		introducedList = []
		dictionaryKeys = sorted(dictionary.keys())
		for key in sorted(dictionaryKeys):
			introducedList.append(dictionary[key])
		headerList = ['position', 'sequence', 'value'] + dictionaryKeys
		if append:
			writeMode = 'a'
		else:
			writeMode = 'w'
		if percentage:
			kmerAbundanceDict = self.getKmerPercentages(kmerAbundanceDict)
		separator = '\t'
		myDict = {}
		out = open(output, writeMode)
		if headerFlag:
			out.write(separator.join(headerList) + '\n')
		for position in kmerAbundanceDict.keys():
			for kmer in kmerAbundanceDict[position].keys():
				value = kmerAbundanceDict[position][kmer]
				ll = [str(position), kmer, str(value)] + introducedList
				out.write(separator.join(ll) + '\n')
		return 1

	def writeNucleotideAbundanceTable(self, nucleotideAbundanceDict, output, nucleotideOrder='ATGC', percentage=False):
		if percentage:
			nucleotideAbundanceDict = self.getNucleotidePercentages(nucleotideAbundanceDict)
		separator = '\t'
		nucleotides = nucleotideOrder
		out = open(output, 'w')
		positions = nucleotideAbundanceDict.keys()
		for position in positions:
			out.write(separator + str(position))
		out.write("\n")
		for nucleotide in nucleotides:
			out.write(nucleotide)
			for position in positions:
				out.write(separator + str(nucleotideAbundanceDict[position][nucleotide]))
			out.write("\n")
		return 1

	def getSequenceLengths(self):
		'''prints header and the sequence length in tab-separated format'''
		for seqObject in self.stream(4096*64*64):
			print(seqObject['h'] + '\t' + str(len(seqObject['s'])))

	def getSequenceLengthDistribution(self):
		'''prints sequence length distribution'''
		myDict = {}
		for seqObject in self.stream2(4096*64*64):
			length = seqObject.getLength()
			if length in myDict.keys():
				myDict[length] += 1
			else:
				myDict[length] = 1
		for key in sorted(myDict.keys()):
			print(str(key) + '\t' + str(myDict[key]))