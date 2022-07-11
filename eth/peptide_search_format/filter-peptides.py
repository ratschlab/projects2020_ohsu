import argparse
from Bio import SeqIO
from pathlib import Path

def filterFasta(fastaFile, outputFile):
	keep = 0
	remove = 0
	pepDic = {}
	
	with open(outputFile,'w') as newFile:
		for record in SeqIO.parse(fastaFile,'fasta'):
			if len(str(record.seq)) >=6 and len(str(record.seq)) <= 50:
				tmpStr = str(record.description)
				tmpStr_sp = tmpStr.split(';')
			
				newFile.write(">"+tmpStr_sp[0]+'\n')
				newFile.write(str(record.seq)+'\n')

				pepDic[str(record.seq)] = 1
				keep += 1
			else:
				remove += 1

	print("keep: " + str(keep))
	print('remove: ' + str(remove))
	print("# unique pep: " + str(len(pepDic)))


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Take results from extract-peptides and filter out peptides that do not meet size requirements for database search. Also cleans up description of each peptide')
	parser.add_argument("fastaFile",help='fasta file')
	parser.add_argument("outputFile",help='output fasta file')
	args = parser.parse_args()
	filterFasta(args.fastaFile, args.outputFile)
