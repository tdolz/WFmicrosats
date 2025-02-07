# WFmicrosats
microsatellite study of winter flounder from the Saltonstall-Kennedy grant 2016-2017. Code used in development of "Patterns of persistence: Genetic and behavioral population complexity of winter flounder amid population declines" https://onlinelibrary.wiley.com/doi/10.1111/jfb.15890 This code is not very organized but please feel free to reach out if you have any questions! Tara.Dolan@mass.gov

**Code organization**
****

#### diversity_stats folder ####

**Microsats_20.R**  
This is for the version where the bay comparisons are based on all individuals.  
performs all the diversity stat tests, like Shannon Diversity index, etc.  
Also does isolation by distance, and an FST test. This not only contains between bay statistics but also between contingents and cohorts in Mattituck and Shinnecock. It outputs figures and info for tables to the diversity_output_files and diversity_figs folders under the subfolders "bay". 

**datacleanup_comparison.R**   
Does the QAQC in the form of linkage disequilibrium, cleaning up missing data, deciding what loci to remove etc. We try doing this 1) using all individuals, 2) using 2016 YOY only and 3) using 2016 YOY and rareifiying the data. 
**Microsats_YOY16.R** 
Is a copy of Microsats_20.R, except that only 2016 YOY are used for the between bay comparisons. Also the bays where we have more than one cohort are NOT rareified.  It outputs figures and info for tables to the diversity_output_files and diversity_figs folders under the subfolders "YOY16_bay".

**Microsats_YOY16RARE.R** Is a copy of Microsats_YOY16.R, except that the bays where we have more than one cohort ARE rareified.  It outputs figures and info for tables to the diversity_output_files and diversity_figs folders under the subfolders "YOY16_rare".

**FST.R**
Computes Weir & cockerham pairwise FST for bays comparisons and cohort/contingent comparisons. The cohort/contingent comparisons are saved to the diversity_output_files folder under the name of the bay they are from. The between bay comparisons where all samples are used are saved in the bay subfolder. The results of between bay comparisons where rarified 2016 YOY are used are saved in the YOY16_bay subfolder. 


