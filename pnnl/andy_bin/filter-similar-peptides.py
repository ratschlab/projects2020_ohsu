import argparse
import re


def filterSimilarPeptides(inputFile,bound):
	newFileName = inputFile[:-4] # removes '.txt
	dic={}
	with open(inputFile,'r') as file1, \
         open(newFileName+"_filter.fa",'w') as fastaFile, \
         open(newFileName+"_filter.txt",'w') as newFile:

		# no header line 
		line1=file1.readline().strip()
		cnt = 0
		nLine = 0
		while line1:
			line1_sp = line1.split('\t')
		
			relPep = line1_sp[0]
			neighPep = line1_sp[1]
			simScore = float(line1_sp[2])
			relPepMass = line1_sp[3]
			neighPepMass = line1_sp[4]
			massDiff = line1_sp[5]
			massDiffPpm = line1_sp[6]

			# write to file if
			# score > bound
			# ppm error < mz_thresh
			# and neighbor peptide not seen before
			if simScore > bound and neighPep not in dic:
				newFile.write(line1+'\n')

				fastaFile.write(">"+line1_sp[0]+"_"+line1_sp[1]+"_"+str(simScore)+"_"+ \
                                relPepMass+"_"+neighPepMass+"_"+massDiff+"_"+massDiffPpm+'\n')
				fastaFile.write(neighPep+'\n')

				cnt= cnt + 1
				dic[neighPep] = 1

			nLine = nLine + 1
			line1=file1.readline().strip()


	print("Printed "+str(cnt)+" lines out of "+str(nLine))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Takes the output from pepsim.py. The output of from pepsim.py is a pairwise list of peptides. For each line (ie pair of peptide), there is sequence, proportion of MS2 peaks in common, and m/z difference. This script takes that output and filters so that it only contains similar peptides above a similarity threshold (bound parameter). In addition it filters out any peptides whose masses are not within a ppm tolerance (thresh paramter). The output of this script is a fasta file that contains the neighbor peptides. Assume that the first peptide in line each is the relevant peptide and the second peptide is the neighbor peptide. Note that I assume isotope error of 0')
	parser.add_argument('inputFile',help='Output from pepsim.py')
	parser.add_argument('bound',help='Similarity threshold (lower bound)',type=float)
	args = parser.parse_args()
	filterSimilarPeptides(args.inputFile,args.bound)
