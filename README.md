# WFmicrosats
microsatellite study of winter flounder from the Saltonstall-Kennedy grant 2016-2017

**Code organization**
****

#### diversity_stats folder ####

**Microsats_20.R** does some QAQC in the form of linkage disequilibrium, cleaning up missing data, etc.  
This is for the version where the bay comparisons are based on all individuals.  
performs all the diversity stat tests, like Shannon Diversity index, etc.  
Also does isolation by distance, and an FST test. This not only contains between bay statistics but also between contingents and cohorts in Mattituck and Shinnecock. It outputs figures and info for tables to the diversity_output_files and diversity_figs folders under the subfolders "bay". 

**Microsats_BYC.R** Is a copy of Microsats_20.R, except that only 2016 YOY are used for the between bay comparisons. Also the bays where we have more than one cohort are rareified.  It outputs figures and info for tables to the diversity_output_files and diversity_figs folders under the subfolders "YOY16_bay".

**FST.R**
Computes Weir & cockerham pairwise FST for bays comparisons and cohort/contingent comparisons. The cohort/contingent comparisons are saved to the diversity_output_files folder under the name of the bay they are from. The between bay comparisons where all samples are used are saved in the bay subfolder. The results of between bay comparisons where rarified 2016 YOY are used are saved in the YOY16_bay subfolder. 


