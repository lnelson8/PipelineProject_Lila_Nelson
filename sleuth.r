library(sleuth)
library(dplyr)

# reading in the table that is a part of the repo that has the accession numbers and correlating conditions and paths for the kallisto output
stab = read.table("Pipeline_Project/sleuth_input",header=TRUE)

# initializing sleuth with the table 
so = sleuth_prep(stab)

# this command fits a model for the conditions
so = sleuth_fit(so, ~condition, 'full')

# this command fits the reduced model
so = sleuth_fit(so, ~1, 'reduced')

# this command performs the liklihood ratio tests for the full and reduced fit models
so = sleuth_lrt(so, 'reduced', 'full')

#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 

#print top 10 transcripts
head(sleuth_significant, n=10)

#write FDR < 0.05 transcripts to file
write.table(sleuth_significant, file="HCMV_results.txt",quote = FALSE,row.names = FALSE)


