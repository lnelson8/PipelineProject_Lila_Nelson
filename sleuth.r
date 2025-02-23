library(sleuth)
library(dplyr)


stab = read.table("Pipeline_Project/sleuth_input",header=TRUE)

so = sleuth_prep(stab)

so = sleuth_fit(so, ~condition, 'full')

so = sleuth_fit(so, ~1, 'reduced')

so = sleuth_lrt(so, 'reduced', 'full')

#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 

#print top 10 transcripts
head(sleuth_significant, n=10)

#write FDR < 0.05 transcripts to file
write.table(sleuth_significant, file="HCMV_results.txt",quote = FALSE,row.names = FALSE)


