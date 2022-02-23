
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#indir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2_linked_part/SegmExpr_mtx
#outdir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2_linked_part_snappy/SegmExpr_mtx

indir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2_linked_part/JuncExpr_mtx
outdir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/output/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_1fc5828_pya.0.17.1_ref_mode2_linked_part_snappy/JuncExpr_mtx

bsub -J matrix_conversion -W 24:00 -n 12 -R "rusage[mem=19000]" "${SCRIPT_DIR}/convert_dir.py ${indir} ${outdir} 12 16000"

