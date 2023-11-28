import argparse

def rerank(searchFileName):
	scanDic = {}
	with open(searchFileName,'r') as file1, open('tide-search-filter-rerank.txt','w') as newFile:
		header = file1.readline().strip()
		newFile.write(header+'\n')

		header_sp = header.split('\t')

		assert (header_sp[1] == "scan"),"incorrect input file"
		assert (header_sp[9] == "xcorr rank"),"incorrect input file"

		for line1 in file1:
			line1_sp = line1.split('\t')

			curScan = line1_sp[1]
			curRank = line1_sp[9]

			if curScan not in scanDic:
				line1_sp[9] = "1"
				newLine = '\t'.join(line1_sp)
				newFile.write(newLine)

				scanDic[curScan] = 1



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Given a tide-search file. Rerank PSMs so that first instance of scan is kept (which we assume is best scoring) and changed to rank=1')
	parser.add_argument("searchFile",help='tide-search file')
	args = parser.parse_args()
	rerank(args.searchFile)
