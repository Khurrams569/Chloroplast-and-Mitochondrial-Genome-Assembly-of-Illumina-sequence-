# Chloroplast-and-Mitochondrial-Genome-Assembly-of-Illumina-sequence

This method used for Illumina sequence nuclear dataset. Method used to map nuclear genome with reference genome and extract the chloroplast and mitochondrial sequence further used for assembly of chloroplast and mitochondria

* Commands and Procedure Local Server  and HCC

1.	Trimming
Command in HCC: 
>module load trimmomatic/0.36

>trimmomatic PE L4_566_GDD_1.fq L4_566_GDD_2.fq /common/yinlab/kshahzad2/L4_566_GDD_forward_paired.fq.gz /common/yinlab/kshahzad2/L4_566_GDD_forward_unpaired.fq.gz /common/yinlab/kshahzad2/L4_566_GDD_reverse_paired.fq.gz /common/yinlab/kshahzad2/L4_566_GDD_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Command in Glu: Software already installed, just follow the commands
>java -jar /home/khurrams569/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 V300031144_L4_TEArrrRAACB-572_1.fq.gz V300031144_L4_TEArrrRAACB-572_2.fq.gz L4_572_forward_paired.fq.gz L4_572_forward_unpaired.fq.gz L4_572_reverse_paired.fq.gz L4_572_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

4 new files are generated.  
*L4_566_GDD_forward_paired.fq.gz 
*L4_566_GDD_forward_unpaired.fq.gz 
*L4_566_GDD_reverse_paired.fq.gz 
*L4_566_GDD_reverse_unpaired.fq.gz

Use paired sequence files for next step. 

2.	Mapped sequence with reference genome
* Tools used Bowtie; See the commands below in HCC;
>module load bowtie/2.5

>bowtie2-build –-help

>bowtie2-build -f input_reference.fasta bowtie

>bowtie –-help
 
>bowtie2 --very-fast-local -x /common/yinlab/kshahzad2/bowtie -1 /common/yinlab/kshahzad2/L4_566_GDD_forward_paired.fq.gz -2 /common/yinlab/kshahzad2/L4_566_GDD_reverse_paired.fq.gz -S /common/yinlab/kshahzad2/L4_566_GDD_mapped.sam

*To convert .sam into .bam file use this command in samtools in HCC:

*Module load samtools

>samtools view -S -b /input_Path/filename.sam > /output_path/filename.bam

>convert .sam file into .fastq file by following command:

>samtools fastq -F 4 L3_568_mapped.sam > L3_568_mapped.fastq


Tools used Bowtie; See the commands below in Glu and Gly;
	
TO add channels

$ conda config --show channels

$ conda config --add channels bioconda

to create environments 

$ conda env list

$ conda create -n bioinformatics

$ conda activate bioinformatics

After add environment you need to install software’s. Every time when you need to use software first you have to activate the environment by above commands then use software. 

For install any software use the below command.
$ conda install -c bioconda bowtie

$ bowtie2-build -f GDD_L4_572.fasta bowtie


$ bowtie2 --very-fast-local -x /mnt/array1/khurrams/FY3/bowtie -1 /mnt/array1/khurrams/FY3/V300035165_L2_TEArrrRAABB-512_1.fq.gz -2 /mnt/array1/khurrams/FY3/V300035165_L2_TEArrrRAABB-512_2.fq.gz -S /mnt/array1/khurrams/FY3/L2_512_mito_mapped.sam

3.	Chloroplast genome assembly

$ conda create -n assembly

$ conda create -n assembly nanoplot flye bandage bwa samtools pilon

$ conda activate assembly

Convert .sam file into .fastq file by following command:

$ samtools fastq -F 4 L3_568_mapped.sam > L3_568_mapped.fastq

To install NOVOplasty in Gly Glu and HCC, use this command: 

$ git clone https://github.com/ndierckx/NOVOPlasty.git

Make a config.txt file by following link: https://github.com/ndierckx/NOVOPlasty/blob/master/config.txt

Use this command to run the NOVOPlasty in all servers: 

$ perl /home/khurrams569/NOVOPlasty/NOVOPlasty4.3.5.pl -c /mnt/array1/khurrams/FZ2/config.txt

Use this command to run the NOVOPlasty in HCC:

$ perl /common/yinlab/kshahzad2/NOVOPlasty/NOVOPlasty4.3.5.pl -c /common/yinlab/kshahzad2/assembly/config.txt

Check the sequence in Blast online: 

https://blast.ncbi.nlm.nih.gov/Blast.cgi

Use Nucleotide BLAST option. Check the box in front of “Align two or more sequence”.
Upload assembled file in “Enter Query Sequence”
Upload reference file in “Enter Subject Sequence”

Click on BLAST
The results will show 99.75% similarity sequence. 

Check the sequence in mapped file: 

$ less L4_568_mapped.sam

Extract the mapped sequence IDs from mapped sequence file L4_568_mapped.sam from following command

$ less L4_568_mapped.sam | awk '{if ($3 != "*") print $1}' > L4_568_mapped.ID &

Learn how to sort and delete the duplicate IDs from the mapped.ID file. Then extract IDs from your forward and reverse files from the following command.  

$ sort L4_568_mapped.ID | uniq -i > L4_568_map_truid

L4_568_map_truid is the file which contain all ID numbers without duplicate names. 

Use this command to check the file without duplicate names
$ cat L4_568_map_truid

To add “/1” at the end of each line
$ less L3_524_map_truid | awk '{print $0"/2"}' > L3_524_map_reverse_id




Use this command to remove all the 

Use seqtk for extract mapped sequence from original sequence files.

$ seqtk subseq L1_542_reverse_paired.fq.gz L1_542_map_reverse_id > L1_542_mapped_mito_reverse.fastq &

Spades command line in glu and gly

$ source activate spades_env

$ spades.py --careful --trusted-contigs /mnt/array1/khurrams/DZD1/mito_ref_3.fasta -o SPADES_OUT -1 L1_542_mapped_mito_forward.fastq -2 L1_542_mapped_mito_reverse.fastq

how to check the number of base pairs in Spades assembly

$ grep ‘>’ SPADES_OUT/contigs.fasta |wc -l


To blast the sequence in HCC:
$ blastn -query input_reads.fasta -db input_reads_db -out blastn_output.alignments
![image](https://github.com/Khurrams569/khurrams/assets/165841830/f78c50d6-0cad-4e22-9119-92344eb78d18)
