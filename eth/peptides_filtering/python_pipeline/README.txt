Steps 
Count number of samples assing a threshold
1.a) wrap_send_cohort_threshold.sh 
1.b) cohort_threshold.py 
1.c) (with helpers_filter.py)
# Example in send_background_threshold.sh
# Example in send_foreground_threshold_b.sh

Join the background on the foreground table
2.a) wrap_intermediate_cohorts.sh 
2.b) intermediate_cohorts.py 
# Example in send_inter_cohorts.sh

Filter depending on filter criteria defined
3.a) wrap_filtering_from_intermed.sh
3.b) filtering_from_intermediate.py (formerly
p20230216_d_filtering_table_thresholds-updated.ipynb)
3.c) With helpers_ffi.py
#Example in send_filtering_from_intermediate.py

###
Perform a similar analysis with the star junctions. Quantify their
expression, and then remove them based on the filtering criteria

Join the background on the foreground table for the STAR junctions
4.a) wrap_STAR_intermediate_cohorts.sh
4.b) star_intermediate_cohorts.py 
4.c) (with helpers_STAR.py)
# Example in send_star_inter_cohorts.sh

###
Make the intersection of the STAR and the GP filtering (on peptides)
5) p20230524_star_GP_intersection.ipynb
Renaming outputs with codes
6) p20230712_rename_outputs.ipynb

###
Developement code for filtering 
p20230216_d_filtering_table_thresholds-updated.ipynb

Development code to select some relevant genes
p20230323_extract_genes_lift_over.ipynb
