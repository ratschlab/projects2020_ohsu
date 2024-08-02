# We run code from PWD directory + DATA SET
PWD='/cluster/home/prelotla/github/projects2020_ohsu/eth/MY_Master_thesis_rerun_LP/TASKS'

# This directory use for save files like (data, plots, tables and etc.)
SAVE_DIR = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/202301_myurchikova_MT_rerun_LP'

## PATH to table(s) 
# Non-filtered BRCA ETH
ETH_PATH_BRCA = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/filtering_intermediate/complete_cancer_candidates_order_r_complete.tsv.gz'
# Non-filtered BRCA OHSU
OHSU_PATH_BRCA = '/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/archive/OHSU_June2023_filter-debug_complete-annotated-shortlist.tsv.gz' #LP changed path with archive
# List of genes
BATCH_PATH_BRCA = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/batch_to_gene.txt'
# Directory from filtered BRCA ETH
FILTERING_PATH_BRCA = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102'
OHSU_BRCA_NEW=True
# Filtered BRCA OHSU (not all filtering parameters)
TAR_OHSU_BRCA= '/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_June2023_filter-debug_all_output.tar.gz'
# Filtered BRCA OHSU (all filtering parameters)
TAR_OSHU_BRCA_NEW='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/June28_renamed_kmerfiles_OHSU.tar.gz'
# One of filtering expirement
FILTERING_ID= 'filters_19May_order_5ge_wAnnot_GPstar'

# Create directory structure
# SAVE_DIR
# | DIR_CSV
# | - DIR_BRCA
#   | - NAME_TABLES
#     | - NAME_OHSU_BRCA
#     | - NAME_ETH_BRCA
#     | - NAME_NON_FILTERING_BRCA
#     | - NAME_FILTERING_BRCA
#     | - NAME_ETH_TASK_BRCA
#     | - NAME_OHSU_TASK_BRCA
#     | - NAME_FINAL_BRCA
#     | - NAME_FINAL_SORTED_FILTER_BRCA
#     | - NAME_PERCENT_BRCA
#   | - NAME_SAMPLES
#     | - *some sample*
#       | - PLOT_SORT_BY
#       | - PRETTY_PLOT_SORT_BY
#       | - plots & filters
# | - DIR_OVARIAN

# Definition constants

PATTERN='[\d+]{1,7}|Any'
FILTER_PATTERN='[\d+]{1,7}|Any'
SAMPLE_PATTERN='SampleLim[\d+]{1,7}\.'
COHORTLIM_PATTERN='Cohort[Ll]im[\d+]{1,7}\.'
ACROSS_PATTERN='Across[\d+]{1,7}|AcrossAny'        
        
# Choose necessary data from tables and concatenate filtered & non-filtered cleared tables in one format.
RESTRICTS_BRCA=['TCGAC8A12P01A11RA11507',
                 'TCGAAOA0JM01A21RA05607',
                 'TCGABHA18V01A11RA12D07',
                 'TCGAA2A0D201A21RA03407',
                 'TCGAA2A0SX01A12RA08407']

RESTRICTS_OVARIAN=['TCGA25131901A01R156513',
            'TCGA25131301A01R156513',
            'TCGA61200801A02R156813',
            'TCGA24143101A01R156613',
            'TCGA24229801A01R156913',]

DIR_CSV='DATA'
DIR_BRCA='BRCA'
DIR_OVARIAN='OVARIAN'
NAME_TABLES='TABLES'
NAME_SAMPLES='SAMPLES'

NAME_OHSU_BRCA='ohsu_BRCA.csv'
NAME_ETH_BRCA='eth_BRCA.csv'
NAME_NON_FILTERING_BRCA='out_df_non_filtering_BRCA.csv'
NAME_FILTERING_BRCA='out_df_filtering_BRCA.csv'
NAME_ETH_TASK_BRCA='ETH_task_BRCA.csv'
NAME_OHSU_TASK_BRCA='OHSU_task_BRCA.csv'
NAME_FINAL_BRCA='final_df.csv'
NAME_FINAL_SORTED_FILTER_BRCA='final_df_sorted_with_filter.csv'
NAME_PERCENT_BRCA='percent_table_BRCA.csv'

NAME_OHSU_OVARIAN='ohsu_OVARIAN.csv'
NAME_ETH_OVARIAN='eth_OVARIAN.csv'
NAME_NON_FILTERING_OVARIAN='out_df_non_filtering_OVARIAN.csv'
NAME_FILTERING_OVARIAN='out_df_filtering_OVARIAN.csv'
NAME_ETH_TASK_OVARIAN='ETH_task_OVARIAN.csv'
NAME_OHSU_TASK_OVARIAN='OHSU_task_OVARIAN.csv'
NAME_FINAL_OVARIAN='final_df_OVARIAN.csv'
NAME_FINAL_SORTED_FILTER_OVARIAN='final_df_sorted_with_filter_OVARIAN.csv'
NAME_PERCENT_OVARIAN='percent_table_OVARIAN.csv'

NAME_PLOT_PRETTY ='Pretty Ploting'
NAME_PLOT_ABSOLUT_PRETTY='Pretty Ploting_matplotlib_absolut'
NAME_PLOT_ABSOLUT_PRETTY_M1000='Pretty Ploting_matplotlib_absolut_l1000'
NAME_PLOT_ABSOLUT_PRETTY_L1000='Pretty Ploting_matplotlib_absolut_m1000'
NAME_PLOT_ABSOLUT_PRETTY_NF='Pretty Ploting_matplotlib_absolut_nf'
NAME_NORMAL_SB_NUMBERS='normalized_stacked_barplot_with_number'
NAME_PROPORTIONAL_PLOT='proportional_stacked_and_unstacked__precent_barplot'
NAME_PROPORTIONS_JP_GP='proportions_of_JP_GP'
NAME_PLOT_PERCENT_PRETTY='Pretty Ploting_matplotlib_percent'
NAME_PLOT_LOI_C1='Lost of Intersection C1'
NAME_PLOT_LOI_C2='Lost of Intersection C2'
NAME_PLOT_GTEX='GTEX'
NAME_PLOT_GP='Story of coordinate GP'
NAME_PLOT_JP='Story of coordinate JP'

LOGTHM = 'logarithm'

PDF = '.pdf'
PNG = '.png'
# 'intersection' or 'ohsu'
STORY_OF_FILTER_COOR='intersection'

# 'unstacked' or 'stacked'
PLOT_TYPE='unstacked'

# 'x-axis sorted by size of intersection and priority' or 'x-axis sorted by total' or 'x-axis sorted by size of GP-specific set and priority' or 'x-axis sorted by size of JP-specific set and priority'
PLOT_SORT_BY='x-axis sorted by size of GP-specific set and priority'

PLOT_GP_SORT_BY='x-axis sorted by size of GP-specific set' #x-axis sorted by size of intersection or x-axis sorted by size of GP-specific set

PLOT_JP_SORT_BY='x-axis sorted by size of JP-specific set' #x-axis sorted by size of intersection or x-axis sorted by size of JP-specific set
# Priority task
# PLOT_SORT_BY_PRIORITY=0

# 'percent' or None (Hight Priority in "Pretty Plotting percent")
PRETTY_PLOT_SORT_BY='percent'


LOG = False
LANG='ENG'
OHSU_COLOR = '#1d8ea9'
ETH_COLOR = '#f27700'
OHSU_ETH_COLOR = '#a3b49b'

COLORS_COOR = {
            'size_eth\ohsu':ETH_COLOR,
            'size_intersection':OHSU_ETH_COLOR,
            'size_ohsu\eth':OHSU_COLOR,
            'size_eth\ohsu_coor': ETH_COLOR,
            'size_intersection_coor':OHSU_ETH_COLOR ,
            'size_ohsu\eth_coor':OHSU_COLOR
    }


COLORS={
            'coordinates from GP\JP':COLORS_COOR['size_ohsu\eth_coor'] ,
            'coordinates from intersection':COLORS_COOR['size_intersection_coor'] ,
            'coordinates from JP\GP':COLORS_COOR['size_eth\ohsu_coor'] ,
            'kmers from GP\JP': COLORS_COOR['size_ohsu\eth'],
            'kmers from intersection':COLORS_COOR['size_intersection'],
            'kmers from JP\GP':COLORS_COOR['size_eth\ohsu']}

TEXT_SIZE = 65
TT='\n---------------------------------------------\n'

ETH_COLMNS= {
             'TCGAC8A12P01A11RA11507all':'TCGAC8A12P01A11RA11507',
             'TCGAAOA0JM01A21RA05607all':'TCGAAOA0JM01A21RA05607',
             'TCGABHA18V01A11RA12D07all':'TCGABHA18V01A11RA12D07',
             'TCGAA2A0D201A21RA03407all':'TCGAA2A0D201A21RA03407',
             'TCGAA2A0SX01A12RA08407all':'TCGAA2A0SX01A12RA08407',
            }