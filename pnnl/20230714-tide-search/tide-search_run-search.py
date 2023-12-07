from pathlib import Path
import sys
import argparse
import re
import subprocess

def runSearch(runFileName,mapFileName, dataFolderName, tideIndexFolder):
	t_c_dic = {} # TCGA to CPTAC dic
	c_t_dic = {} # CPTAC to TGCA dic
	with open(mapFileName,'r') as file1:
		header = file1.readline().strip()
		for line1 in file1:
			line1_sp = line1.strip().split('\t')

			tcga_sample = line1_sp[0]
			# the re.sub is needed for ovarian mzML runs. The breast samples are mgf files
			cptac_sample = re.sub(".mzML.gz",'',line1_sp[1])

			t_c_dic[tcga_sample] = cptac_sample
			c_t_dic[cptac_sample] = tcga_sample

	with open(runFileName,'r') as file1:
		for line1 in file1:
			line1 = line1.strip()

			runDir = str(Path(line1).stem)
			runDir = re.sub(".mzML",'',runDir)

			if runDir in c_t_dic:
				bashCommand = "~/projects/ohsu-eth-collab/bin/crux-toolkit/Release/src/crux tide-search --overwrite T --precursor-window 40 --concat T --top-match 1000000000 --num-threads 1 --output-dir crux-output/" + runDir + " " + dataFolderName + "/" + line1 + " " + tideIndexFolder + "/" + c_t_dic[runDir] + "/tide-indicies/final/"
				subprocess.run(bashCommand,shell=True)
			else:
				print('crap')
				print(line1)
		

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='given mapping file and run file. Run searches')
	parser.add_argument('runFile',help='file containing runs')
	parser.add_argument('mapFile',help='map between CPTAC and TCGA')
	parser.add_argument("dataFolder",help='folder containing mzML files')
	parser.add_argument("tideIndex",help='folder that contains all tide indicies')
	args = parser.parse_args()
	runSearch(args.runFile,args.mapFile,args.dataFolder, args.tideIndex)
