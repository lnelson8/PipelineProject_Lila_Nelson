import os
import glob
from Bio import SeqIO
import statistics
import sys
import argparse

#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="ADD TITLE OF SCRIPT HERE (shows on help)")
    parser.add_argument("-i", "--input",
    help="input file",
    required=True)
    parser.add_argument("-o", "--output",
    help="output file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

# my kallisto function takes in the input path, the output file, the accession number, and metadata
def kallisto(input, output, accession, metadata):
    # the transcriptome variable is the command for ncbi datasets that I am using to get the transcriptome for HCMV including the cds (the transcriptome)
    transcriptome = f'datasets download genome accession {accession} --include cds'
    os.system(transcriptome) # os.system calls the function and then after I unzip the dataset
    os.system('unzip ncbi_dataset.zip')
    # the below command is for creating the kallisto index from the cds in the dataset
    index_command = f'kallisto index -i HCMV_index.idx ncbi_dataset/data/{accession}/cds_from_genomic.fna'
    os.system(index_command)
    # now I am opening my output file to write to (this will write over any existing file with the name)
    with open(output, 'w') as out:
        # counting the number of cds by using SeqIO parse to make a list of the number of sequences and counting the length
        cds = list(SeqIO.parse('ncbi_dataset/data/'+accession+'/cds_from_genomic.fna', 'fasta'))
        out.write(f'The HCMV genome (NC_006273.2) has {len(cds)} CDs\n')
    # using glob to make a list of paths to my forward reads  
    reads = glob.glob(input+'/*1.fastq')
    os.system('mkdir results') # making a directory to store the results 
    for read in reads:
        name = os.path.basename(read) # this command gets just the name of the file when glob.glob was getting the entire path this is very useful for naming new files
        # the following command calculates the tpms using kallisto 
        kallisto_command = 'time kallisto quant -i HCMV_index.idx -o results/'+name[:-8]+' -b 10 -t 2 '+read+' '+read[:-7]+'2.fastq'
        os.system(kallisto_command)
    # now I am reopening the outfile as append so the previous data stays 
    with open(output, 'a') as out:
        tpm = glob.glob('results/*/abundance.tsv') # making a list of the abundance file paths from the kallisto results
        out.write('sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n') # making a header that is tab-delimited 
        for i in tpm: # for each file I open the file and make a list of each of the values 
            with open(i, 'r') as f:
                tpm_list = [] # this empty list will house the tpm values I get from each line in the results file
                abund = f.read().splitlines() # this splits the file by line
                for line in abund[1:]: # looping through starting after the header 
                    values = line.split('\t') # splitting each line into a list (I only want the last value)
                    tpm_list.append(float(values[-1])) # appending the tpm value to the list I made
                sample = i.split('/')[1] # the sample is accession number which I retreive from the path name
                condition = metadata[sample][1] # I use the sample variable to call the values that correspond in the metadata dictionary
                # below I am calculating the stats from my list of tpm values 
                mini = min(tpm_list) 
                med = statistics.median(tpm_list)
                mean = statistics.mean(tpm_list)
                maxi = max(tpm_list)
                # writing out the results for each file to the output
                out.write(sample+'\t'+condition+'\t'+str(mini)+'\t'+str(med)+'\t'+str(mean)+'\t'+str(maxi)+'\n')

# next is my sleuth command which only needs the path for the output file to write out the answers 
def sleuth(output):
    os.system('Rscript Pipeline_Project/sleuth.r') # the only command is calling the r.file that is uploaded on the github repo 
    # it reads in a file that is also uploaded on the repo 
    # this function then writes the results to the output file 
    with open(output, 'a') as out:
        out.write('target_id\ttest_stat\tpval\tqval\n')
        with open('HCMV_results.txt', 'r') as f:
            values = f.read().splitlines() # creating a list of the values int he sleuth output 
            header = values[0].split() # taking the header and making it into a list so i can find the indeces of each value I want
            # below I am getting the index position for the values I want to write out to the file
            target_id = header.index('target_id')
            test_stat = header.index('test_stat')
            pval = header.index('pval')
            qval = header.index('qval')
            # for every line in the file (except the header) I am making a list and getting the right numbers
            for q in values[1:]:
                v = q.split(' ')
                out.write(f'{v[target_id]}\t{v[test_stat]}\t{v[pval]}\t{v[qval]}\n')

# my bowtie2 function takes in the input and output paths as well as the accession and metadata 
def bt2(input, output, accession, metadata):
    # making a directory both the bt2 output and one for the index to keep the results neat
    os.system('mkdir bowtie2')
    os.system('mkdir bt2_index')
    # creating the index for the HCMV transcriptome 
    bowtie_index = 'bowtie2-build ncbi_dataset/data/'+accession+'/cds_from_genomic.fna bt2_index/HCMV_index'
    os.system(bowtie_index)
    # again I am looping through the pairs in the data folder
    reads = glob.glob(input+'/*1.fastq')
    for read in reads:
        name = os.path.basename(read)
        # for each pair I am running bowtie2 getting fastq files with the aligned reads as well as a sam file
        bowtie_command = 'bowtie2 -x bt2_index/HCMV_index -1 '+read+' -2 '+read[:-7]+'2.fastq -S bowtie2/'+name[:-8]+'HCMV_map.sam --al-conc-gz bowtie2/'+name[:-8]+'_mapped_%.fq.gz'
        os.system(bowtie_command)
        # to get the number of reads that aligned I used samtools view
        # this first command gets the total number of read pairs
        os.system('samtools view -c bowtie2/'+name[:-8]+'HCMV_map.sam -o before')
        # this second one gets the number of read pairs that did align 
        os.system('samtools view -c -F 4 bowtie2/'+name[:-8]+'HCMV_map.sam -o after')
        # getting those numbers from the files they wrote out to
        with open('before', 'r') as f:
            before = f.read().strip()
        with open('after', 'r') as f:
            after = f.read().strip()
        # writing out the before and after read values to my output file
        with open(output, 'a') as out:
            out.write(' '.join(metadata[name[:-8]])+' had '+before+' read pairs before Bowtie2 filtering and '+after+' after.\n')
        
# spades function needs the input path and output path
def spades(input, output):
    # using both 2pi and 6dpi to get a donor 1 assembly and a donor 2 assembly 
    spades_command1 = 'spades.py -k 77 -t 2 --only-assembler --pe-1 1 '+input+'/SRR5660030_1.fastq --pe-2 1 '+input+'/SRR5660030_2.fastq --pe-1 2 '+input+'/SRR5660033_1.fastq --pe-2 2 '+input+'/SRR5660033_2.fastq -o Donor1_assembly/'
    spades_command2 = 'spades.py -k 77 -t 2 --only-assembler --pe-1 1 '+input+'/SRR5660044_1.fastq --pe-2 1 '+input+'/SRR5660044_2.fastq --pe-1 2 '+input+'/SRR5660045_1.fastq --pe-2 2 '+input+'/SRR5660045_2.fastq -o Donor2_assembly/'
    os.system(spades_command1)
    os.system(spades_command2)
    # writing my commands out to my output file
    with open(output, 'a') as out:
        out.write(spades_command1+'\n')
        out.write(spades_command2+'\n')

# blast function only needs the path of where to output the results
def blast(output):
    # creating a local database for the taxon betaherpesvirinae 
    os.system('datasets download virus genome taxon Betaherpesvirinae --refseq --include genome')
    os.system('unzip ncbi_dataset.zip')
    # blast db command takes in the ncbi datasets output and makes the local db called HCMV_blast_db
    blast_db_command = 'makeblastdb -in ncbi_dataset/data/genomic.fna -out HCMV_blast_db -title HCMV_blast_db -dbtype nucl'
    os.system(blast_db_command)
    # for each assembly 
    for i in glob.glob('Donor*_assembly'):
        # I am making a list of the contigs in the fasta so i can get the longest one
        contigs = list(SeqIO.parse(i+'/K77/final_contigs.fasta', 'fasta'))
        # creating a new file that has only the first contig which is the longest one
        with open(i+'/longest_contig', 'w') as f:
            f.write(str(contigs[0].seq))
        # blasting that longest contig against the db I made
        blast_command = 'blastn -query '+i+'/longest_contig -db HCMV_blast_db -out '+i[:6]+'blast -outfmt \"6 sacc pident length qstart qend sstart send bitscore evalue stitle\"'
        os.system(blast_command)
        # writing out the header of the values that were retuned (specified in the blast command)
        with open(output, 'a') as out:
            out.write(i+':\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')
            # opening the blast output files and making alist of the matches
            with open(i[:6]+'blast', 'r') as f:
                matches = f.read().splitlines()
            # getting only the top 10 matches and writing them out to the output file
            if len(matches) > 10:
                matches = matches[0:10]
            for match in matches:
                out.write(match+'\n')

# I created functions that take in the infile (which is just a path to a directory that has all of the data in it) and any other necessary like the outfile name 
# or the accession (which I made a variable incase anyone wanted to alter this to work for another organism) and metadata which pairs up the accession number with their metadata
accession = 'GCA_000845245.1'
metadata = {'SRR5660030': ['Donor 1', '2dpi'], 'SRR5660033':['Donor 1', '6dpi'], 'SRR5660044':['Donor 2', '2dpi'], 'SRR5660045':['Donor 2', '6dpi']}

# the wrapper runs by calling each function and each step can be done separately by commenting out a function that isn't needed/wanted
kallisto(infile, outfile, accession, metadata)
sleuth(outfile)
bt2(infile, outfile, accession, metadata)
spades(infile, outfile)
blast(outfile)
