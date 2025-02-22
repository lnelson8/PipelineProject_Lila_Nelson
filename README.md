# PipelineProject_Lila_Nelson
 
This pipeline is designed to take in FASTQ files, specifically Human herpesvirus 5 (HCMV) transcriptomes, from different patients and compare them. First TPM (transcripts per million) are calculated using kallisto, then the output from kallisto can be given to the tool sleuth and used to determine differentially expressed genes. Bowtie2 will be used to create an index for the HCMV genome and map the transcriptomes to the index. The reads that mapped to the HCMV genome can be assembled using SPAdes producing FASTA files that can be used as a blastn query to determine the most similar strains. This will output a log file (PipelineProject.log) that contains the important information obtained from each step. 

### Dependencies
**the following needs to be downloaded to use this pipeline**

[<ins>kallisto</ins>](https://pachterlab.github.io/kallisto/download)

[<ins>sleuth</ins>](https://pachterlab.github.io/sleuth/download)

[<ins>Bowtie2</ins>](https://github.com/BenLangmead/bowtie2)

[<ins>SPAdes</ins>](https://github.com/ablab/spades?tab=readme-ov-file)

[<ins>BLAST+</ins>](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
