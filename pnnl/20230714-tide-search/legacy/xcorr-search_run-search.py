from pathlib import Path
import sys
import argparse
import re
import subprocess

def runSearch(runFileName,inst):
	with open(runFileName,'r') as file1:
		for line1 in file1:
			line1 = line1.strip()

			runDir = str(Path(line1).stem)
			bashCommand = "/home/bmx8177/projects/ohsu-eth-collab/bin/crux-toolkit/Release/src/crux tide-search --overwrite T --precursor-window 40 --concat T --top-match 1 --num-threads 5 --output-dir crux-output/" + runDir + " ~/projects/ohsu-eth-collab/data/mgf/" + line1 + " ../../20230409-tide-index/" + inst + "/tide-indicies/final/"
			subprocess.run(bashCommand,shell=True)
		

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Run database search against OHSU or ETH tide-index.')
	parser.add_argument('runFile',help='file containing LC-MS/MS runs')
	parser.add_argument('inst',help='Instution that is being analyzed. Either ohsu or eth')
	args = parser.parse_args()
	runSearch(args.runFile,args.inst)
