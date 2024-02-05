import crema
import argparse

def crema_py_API(input_files:list, eval_fdr:float, threshold:float, outdir:str):

    psms = crema.read_tide(input_files, pairing_file_name=None)
    results =  psms.assign_confidence(score_column=None, pep_fdr_type="peptide-only", eval_fdr=eval_fdr, threshold=threshold)
    results.to_txt(output_dir=outdir,  file_root=None, sep="\t", decoys=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes results from tide search and splits them between experimental conditions')
    parser.add_argument("--input-files", nargs='+', help='list of tide search files')
    parser.add_argument("--eval-fdr", type=float, help='fdr level used to get the best column')
    parser.add_argument("--threshold", type=float, help='fdr threshold for results')
    parser.add_argument("--outdir", help='directory to store the results')
    args = parser.parse_args()
    print(args)
    crema_py_API(args.input_files, args.eval_fdr, args.threshold, args.outdir)
