import argparse
from pathlib import Path
import re


def subsetFile(tideSearchFileName,mapFolder,mapRuns,index,front):
	c_t_dic = {} # map from CPTAC file name to TCGA file name
	with open(mapRuns,'r') as file1:
		header = file1.readline()
		for line1 in file1:
			line1_sp = line1.split('\t')

			t_file = line1_sp[0]
			c_file = re.sub(".mzML.gz",'',line1_sp[1])

			c_t_dic[c_file] = t_file

	curSearchFile = Path(tideSearchFileName).parent.stem
	mapFileName = mapFolder + front + c_t_dic[curSearchFile] + "_experiments_per_peptide.tsv"

	exp_id_index = 2
	pep_id_index= 0
	pepIdDic = {} # key is peptide ID present current exp given by index
	with open(mapFileName,'r') as file1:
		header = file1.readline()
		for line1 in file1:
			line1_sp = line1.strip().split('\t')

			exp_id = line1_sp[exp_id_index]
			pep_id = line1_sp[pep_id_index]

			exp_id_sp = exp_id.split(';')

			if str(index) in exp_id_sp:
				# map file just as index
				# search file format is pepId-Number
				pepIdDic["pepID-"+pep_id] = 1

	with open(tideSearchFileName,'r') as tideSearchFile, open("tide-search-filter.txt",'w') as newFile:
		header = tideSearchFile.readline().strip()
		header_sp = header.split('\t')
		newFile.write(header+'\n')

		if header_sp[14] == 'protein id':
			prot_index = 14
		else:
			raise Exception('add to if statment')

		for line1 in tideSearchFile:
			line1_sp = line1.strip().split('\t')

			curProt = line1_sp[prot_index]

			if ',' in curProt:
				curProt_sp = curProt.split(',')
				curProt_sp = [re.sub("\(1\)",'',x) for x in curProt_sp]
				curProt_sp = [re.sub("decoy_",'',x) for x in curProt_sp]
				for word in curProt_sp:
					if word in pepIdDic:
						newFile.write(line1)
						break
			else:
				curProt = re.sub("\(1\)",'',curProt)
				curProt = re.sub("decoy_",'',curProt)
				if curProt in pepIdDic:
					newFile.write(line1)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='give tide-search, experiment mapping file, and ID. Filter tide-search output to ouput based off relevant experiment id')
	parser.add_argument('tideSearch',help='tide-search.txt')
	parser.add_argument('mapFolder',help='folder containing all pepide experiment mapping file')
	parser.add_argument('mapRuns',help="map between CPTAC and TCGA file names")
	parser.add_argument('index',help='experiment index',type=int)
	parser.add_argument("front",help='text in front of each fasta file')
	args = parser.parse_args()
	subsetFile(args.tideSearch,args.mapFolder,args.mapRuns,args.index,args.front)
