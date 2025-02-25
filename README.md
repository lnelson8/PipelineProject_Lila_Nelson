# PipelineProject_Lila_Nelson
 
This pipeline is designed to take in FASTQ files, specifically Human herpesvirus 5 (HCMV) transcriptomes, from different patients and compare them. First TPM (transcripts per million) are calculated using kallisto, then the output from kallisto can be given to the tool sleuth and used to determine differentially expressed genes. Bowtie2 will be used to create an index for the HCMV genome and map the transcriptomes to the index. The reads that mapped to the HCMV genome can be assembled using SPAdes producing FASTA files that can be used as a blastn query to determine the most similar strains. This will output a log file (PipelineProject.log) that contains the important information obtained from each step. 

### Dependencies

**the following needs to be downloaded to use this pipeline:**

- [Python](https://www.python.org/downloads/source/) version 3.10.12+
- [Biopython](https://biopython.org/wiki/Download) version 1.83+
- [R](https://cran.r-project.org/mirrors.html) version 4.4.2+
- [<ins>kallisto</ins>](https://pachterlab.github.io/kallisto/download)version 0.51.1+ 
- [<ins>sleuth</ins>](https://pachterlab.github.io/sleuth/download)
- [<ins>Bowtie2</ins>](https://github.com/BenLangmead/bowtie2) version 2.2.4+
- [<ins>SPAdes</ins>](https://github.com/ablab/spades?tab=readme-ov-file) 
- [<ins>BLAST+</ins>](https://www.ncbi.nlm.nih.gov/books/NBK569861/)
- [<ins>samtools+</ins>](https://github.com/samtools/samtools)

## Testing/Running the wrapper.py script
the below command will run the wrapper with the test data and output a file titles Pipeline_Project.log with the output 

<code>git clone https://github.com/lnelson8/PipelineProject_Lila_Nelson.git</code>

<code>cd PipelineProject_Lila_Nelson</code>

<code>python wrapper.py -i test_data -o Pipeline_Project.log </code>

To run with your own data make sure the paired end fastq files are all in the same directory and change the input (-i) to reflect the path to that directory
