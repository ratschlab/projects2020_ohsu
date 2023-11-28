import argparse
from Bio import SeqIO
from matplotlib import pyplot as plt
import seaborn as sns

def createCDF(fastaFileName):
	pepLenList = []
	seqDic = {}
	cnt = 0
	with open("peptide-extracted-filter-unique.fasta",'w') as newFile:
		for seq in SeqIO.parse(fastaFileName,'fasta'):
			pep = str(seq.seq)
			pepLen = len(pep)

			if pep not in seqDic:
				pepLenList.append(pepLen)
				seqDic[pep] = 1
				cnt += 1

				newFile.write(">"+seq.id+'\n')
				newFile.write(str(seq.seq)+'\n')
	print(str(cnt) + " unique sequence")
		

	fig,ax = plt.subplots()
	sns.ecdfplot(pepLenList,ax=ax)
	plt.xlim(0,50)
	plt.xlabel("length of peptide")
	plt.tight_layout()
	plt.savefig("pep-len-ecdf.png")
	plt.clf()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='create eCDF of peptide lengths')
	parser.add_argument("fasta",help='fasta file')
	args = parser.parse_args()
	createCDF(args.fasta)
