import argparse
from Bio import SeqIO
import re
import bisect

def extractPep(fastaFileName, outputFileName):
	success = 0
	fail = 0
	failNoKR = 0
	failNoKRBefore = 0
	failNoKRAfter = 0
	failBetween = 0
	with open(outputFileName,'w') as newFile:
		for seq in SeqIO.parse(fastaFileName,'fasta'):
			header = seq.id
			header_sp = header.split(';')

			# parser headers
			jxnString = header_sp[1]
			jxnIndex = int(re.sub("jx_pos-","",jxnString))

			betweenString = header_sp[2]
			betweenVal = int(re.sub("between_codons-",'',betweenString))

			include5PrimeString = header_sp[3]
			include3PrimeString = header_sp[4]

			include5Prime = re.sub("includes_5'-",'',include5PrimeString)
			# assertions
			if int(include5Prime) == 0:
				include5Prime = False
			elif int(include5Prime) == 1:
				include5Prime = True
			else:
				raise Exception("not supposed to happen")

			include3Prime = re.sub("includes_3'-",'',include3PrimeString)
			# assertions
			if int(include3Prime) == 0:
				include3Prime = False
			elif int(include3Prime) == 1:
				include3Prime = True
			else:
				raise Exception("not supposed to happen")

			pep = str(seq.seq)

			kIndex = [i for i, ltr in enumerate(pep) if ltr == "K"]
			rIndex = [i for i, ltr in enumerate(pep) if ltr == "R"]

			allIndex = []
			allIndex.extend(kIndex)
			allIndex.extend(rIndex)
			allIndex.sort()

			upperIndex = bisect.bisect_right(allIndex,jxnIndex)
			lowerIndex = bisect.bisect_left(allIndex,jxnIndex)
			if len(allIndex) == 0: # no K/R
				# if K/R not found at all
				# keep if transcript has both 5' and 3' end
				if include3Prime == True and include5Prime == True:
#					print(seq.id)
#					print(seq.seq)
#					print(allIndex)
					success += 1
					newFile.write(">"+seq.id+'\n')
					newFile.write(str(seq.seq)+'\n')
				else:
					fail += 1
					failNoKR += 1
			elif jxnIndex in allIndex: # jxn is same position as K/R
				if betweenVal == 1: # jxn is between codons
					fail += 1
					failBetween += 1
				else:	
					lowerIndexVal = allIndex[allIndex.index(jxnIndex) - 1] + 1
					upperIndexVal = jxnIndex + 1
#					print(seq.id)
#					print(seq.seq)
#					print(allIndex)
#					print(lowerIndex,upperIndex)
#					print(pep[lowerIndexVal:upperIndexVal])
#					print()

					success += 1
					newFile.write(">"+seq.id+'\n')
					newFile.write(pep[lowerIndexVal:upperIndexVal]+'\n')
			elif jxnIndex > allIndex[-1]: # if last K/R is before jxn pnt
				if include3Prime: # if trascript ends at 3':
					lowerIndexVal = allIndex[-1] + 1
					upperIndexVal = len(str(seq.seq))
#					print(seq.id)
#					print(seq.seq)
#					print(allIndex)
#					print(lowerIndex,upperIndex)
#					print(pep[lowerIndexVal:upperIndexVal])
#					print()

					success += 1

					newFile.write(">"+seq.id+'\n')
					newFile.write(pep[lowerIndexVal:upperIndexVal]+'\n')
				else:
					fail += 1
					failNoKRAfter += 1
			elif jxnIndex < allIndex[0]: # if first K/R is after jxn pnt
				if include5Prime: # if transcript start at 5'
					lowerIndexVal = 0
					upperIndexVal = allIndex[0] + 1
#					print(seq.id)
#					print(seq.seq)
#					print(allIndex)
#					print(lowerIndex,upperIndex)
#					print(pep[lowerIndexVal:upperIndexVal])
#					print()

					success += 1

					newFile.write(">"+seq.id+'\n')
					newFile.write(pep[lowerIndexVal:upperIndexVal]+'\n')
				else:
					fail += 1
					failNoKRBefore += 1
			elif upperIndex == lowerIndex:
				# note that +1 are needed to get rid of prefix K/R
				# and to add suffix K/R
				lowerIndexVal = allIndex[lowerIndex-1]+1
				upperIndexVal = allIndex[upperIndex]+1
#				print(seq.id)
#				print(seq.seq)
#				print(allIndex)
#				print(lowerIndex,upperIndex)
#				print(pep[lowerIndexVal:upperIndexVal])
#				print()

				success += 1

				newFile.write(">"+seq.id+'\n')
				newFile.write(pep[lowerIndexVal:upperIndexVal]+'\n')
			else:
				raise Exception("not supposed to happen")

	print('success: ' + str(success))
	print("fail: " + str(fail))
	print("fail (No KR): " + str(failNoKR))
	print("fail (No KR after jxn point): " + str(failNoKRAfter))
	print("fail (No KR before jxn point): " + str(failNoKRBefore))
		

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='extract tryptic peptide from OHSU collab fasta')
	parser.add_argument("fastaFile",help='input fasta file')
	parser.add_argument("outputFile",help='output fasta file')
	args = parser.parse_args()
	extractPep(args.fastaFile, args.outputFile)
