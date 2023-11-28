import argparse
from pathlib import Path
from Bio import SeqIO
from matplotlib import pyplot as plt
import seaborn as sns


def dist(fastaFolderName, mapFolder, typeFileName):
	typeDic = {} # key is file, value is cancer type
	with open(typeFileName,'r') as file1:
		header = file1.readline().strip()
		for line1 in file1:
			line1_sp = line1.split('\t')
			typeDic[line1_sp[0]] = line1_sp[1]

	for f in Path(fastaFolderName).glob("**/*/peptide-extracted-filter.fasta"):
		# set of peptides from transcripts
		pepSet = set()
		for seq in SeqIO.parse(open(f,'r'),'fasta'):
			val = str(seq.id)
			val_sp = val.split("-")
			pepSet.add(val_sp[1])

		print(f.parent.stem)
		mapFile = mapFolder + f.parent.stem + "_experiments_per_peptide.tsv"

		# key is experiment ID; value is set of peptides in experiment
		expDic = {}
		maxExp = 0
		minExp = 10000000
		with open(mapFile,'r') as file1:
			header = file1.readline().strip()
			for line1 in file1:
				line1_sp = line1.strip().split('\t')

				exp_id = line1_sp[2]
				pep_id = line1_sp[0]

				exp_id_sp = exp_id.split(";")
				for exp in exp_id_sp:
					if exp not in expDic:
						expDic[exp] = set(pep_id)
					else:
						tmpSet = expDic[exp]
						tmpSet.add(pep_id)
						expDic[exp] = tmpSet

				exp_id_sp = [int(x) for x in exp_id_sp]
				if max(exp_id_sp) > maxExp:
					maxExp = max(exp_id_sp)
				if min(exp_id_sp) < minExp:
					minExp = min(exp_id_sp)

		# find intersection between exp and peptides
		hist_list = []
		zeroCnt = 0
		for i in range(minExp,maxExp+1):
			tmpSet = expDic[str(i)]

			intersect = tmpSet.intersection(pepSet)
			hist_list.append(len(intersect))

			if len(intersect) == 0:
				zeroCnt += 1

		fig,ax = plt.subplots()
		sns.histplot(hist_list)
		plt.xlabel("# of peptides per experiment")
		plt.title(str(f.parent.stem) + " 0 cnt: " + str(zeroCnt))
		plt.tight_layout()
		plt.savefig(f.parent.stem+".png")
		plt.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='histogram of the number of peptides per experiment')
	parser.add_argument("fastaFolder",help='folder caontaining fasta files of extracted peptides')
	parser.add_argument('mapFolder',help='experiments_per_peptide folder')
	parser.add_argument("typeCancer",help='file for cancer type')
	args = parser.parse_args()
	dist(args.fastaFolder,args.mapFolder,args.typeCancer)
