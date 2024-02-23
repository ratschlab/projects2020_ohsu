import crema
import argparse

def crema_py_API(input_files:list, outdir:str):

    psms = crema.read_tide(input_files, pairing_file_name=None)
    results = psms.assign_confidence(pep_fdr_type="psm-peptide", score_column='xcorr score', eval_fdr=0.05, threshold='q-value')
    results.to_txt(output_dir=outdir,  file_root=None, sep="\t", decoys=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes results from tide search and splits them between experimental conditions')
    parser.add_argument("--input-files", nargs='+', help='list of tide search files')
    parser.add_argument("--outdir", help='directory to store the results')
    args = parser.parse_args()
    print(args)
    crema_py_API(args.input_files, args.outdir)
