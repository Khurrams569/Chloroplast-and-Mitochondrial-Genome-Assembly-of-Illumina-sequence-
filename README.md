# Chloroplast-and-Mitochondrial-Genome-Assembly-of-Illumina-sequence

This method used for Illumina sequence nuclear dataset. Method used to map nuclear genome with reference genome and extract the chloroplast and mitochondrial raw sequences further used for assembly of chloroplast and mitochondria genome. 

* Commands and Procedure Local Server  and HCC

1.	Trimming
Command in HCC: 
>module load trimmomatic/0.36

>trimmomatic PE L4_566_GDD_1.fq L4_566_GDD_2.fq /common/yinlab/kshahzad2/L4_566_GDD_forward_paired.fq.gz /common/yinlab/kshahzad2/L4_566_GDD_forward_unpaired.fq.gz /common/yinlab/kshahzad2/L4_566_GDD_reverse_paired.fq.gz /common/yinlab/kshahzad2/L4_566_GDD_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Command in Glu: Software already installed, just follow the commands

>java -jar /home/khurrams569/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 V300031144_L4_TEArrrRAACB-572_1.fq.gz V300031144_L4_TEArrrRAACB-572_2.fq.gz L4_572_forward_paired.fq.gz L4_572_forward_unpaired.fq.gz L4_572_reverse_paired.fq.gz L4_572_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

4 new files are generated.  
L4_566_GDD_forward_paired.fq.gz L4_566_GDD_forward_unpaired.fq.gz L4_566_GDD_reverse_paired.fq.gz L4_566_GDD_reverse_unpaired.fq.gz

Use paired sequence files for next step. Obtain the reference chloroplast genomes (Accession: OL450428.1.  GI:2190895788) and mitochondrial genome (Accession: OM809792.1. GI: 2241228696) from NCBI. This is from a public database like NCBI or from a previous study. Make sure it's in FASTA forma

2.	Mapped sequence with reference genome
Tools used Bowtie; See the commands below in HCC;
>module load bowtie/2.5
>bowtie2-build –-help
>bowtie2-build -f input_reference.fasta bowtie
>bowtie –-help
 
>bowtie2 --very-fast-local -x /common/yinlab/kshahzad2/bowtie -1 /common/yinlab/kshahzad2/L4_566_GDD_forward_paired.fq.gz -2 /common/yinlab/kshahzad2/L4_566_GDD_reverse_paired.fq.gz -S /common/yinlab/kshahzad2/L4_566_GDD_mapped.sam

Tools used Bowtie; See the commands below in GLu and gly;
>Conda install -c bioconda bowtie
>bowtie2-build -f ref_sequence.fasta bowtie

>bowtie2 --very-fast-local -x /home/khurrams569/bowtie -1 /home/khurrams569/L3_568_forward_paired.fq.gz -2 /home/khurrams569/L3_568_reverse_paired.fq.gz -S /home/khurrams569/L3_568_mapped.sam

3.	Check and extract the mapped sequence from original sequence files: 
>less L4_568_mapped.sam

Extract the mapped sequence IDs from mapped sequence file L4_568_mapped.sam from following command:

>less L4_568_mapped.sam | awk '{if ($3 != "*") print $1}' > L4_568_mapped.ID &

To sort and delete the duplicate IDs from the mapped.ID file. Then extract IDs from your forward and reverse files from the following command.  

>sort L4_568_mapped.ID | uniq -i > L4_568_map_truid

L4_568_map_truid is the file which contain all ID numbers without duplicate names. 

To add “/1” at the end of each line
>less L3_524_map_truid | awk '{print $0"/2"}' > L3_524_map_reverse_id

Use seqtk for extract mapped sequence from original sequence files.
>seqtk subseq L1_542_reverse_paired.fq.gz L1_542_map_reverse_id > L1_542_mapped_mito_reverse.fastq &

4.	Assembly of filtered Chloroplast and Mitochondrial reads by using NOVOPlasty.
•	NOVOPlasty is a de novo assembler specialized for assembling organelle genomes (such as chloroplasts) from whole-genome data, which can also be used for mitochondrial and bacterial genomes. It utilizes a seed-and-extend algorithm, and requires a related reference sequence in fasta format and a seed sequence, which can be a part of the chloroplast genome.
•	Here is a step-by-step guide to using NOVOPlasty through an SSH terminal to assemble chloroplast sequences:
•	Install NOVOPlasty: If NOVOPlasty is not installed on the server, you can clone it from the GitHub repository:
git clone https://github.com/ndierckx/NOVOPlasty.git
•	Prepare Configuration File: NOVOPlasty requires a configuration file (e.g., config.txt) with details about the project. Here is an example of how the file format would look:
•	
•	Project:
•	-----------------------
•	Project name         	 = ChloroplastAssembly
•	Type                 	 	= chloro
•	GenBank file        	  = reference_genome.gb
•	
•	Dataset 1:
•	-----------------------
•	Read Length       	    	= 150
•	Insert size         	    	= 300
•	Platform           	    	= illumina
•	Single/Paired      	    	= PE
•	FastQ1 forward reads      = path_to_forward_reads.fastq
•	FastQ1 reverse reads       = path_to_reverse_reads.fastq
•	Optional:
•	Seed Input                       = seed_sequence.fasta
•	Chimeric reads                = path_to_chimeric_reads.fastq
•	
To install NOVOplasty in Gly Glu and HCC, use this command: 

>git clone https://github.com/ndierckx/NOVOPlasty.git
Make a config.txt file by following link: https://github.com/ndierckx/NOVOPlasty/blob/master/config.txt
Use this command to run the NOVOPlasty in all servers: 
>perl /home/khurrams569/NOVOPlasty/NOVOPlasty4.3.5.pl -c /home/khurrams569/config.txt

Check the sequence in Blast online: 
https://blast.ncbi.nlm.nih.gov/Blast.cgi
Use Nucleotide BLAST option. Check the box in front of “Align two or more sequence”.
Upload assembled file in “Enter Query Sequence”
Upload reference file in “Enter Subject Sequence”
Click on BLAST
The results will show 99.75% similarity sequence. 

5.	Spades for Mitochondria genome assembly:    
Spades command line in glu and gly:
>spades.py --careful --trusted-contigs /mnt/raid5-2/khurrams/mito_L4_524/mito_ref_3.fasta -o SPADES_OUT -1 L3_524_mapped_forward.fastq -2 L3_524_mapped_reverse.fastq

To check the number of base pairs in Spades assembly:
>grep ‘>’ SPADES_OUT/contigs.fasta |wc -l

6.	How to download dataset from NCBI with HCC:
>module load SRAtoolkit/2.11
>vdb-config -i
You will see a screen where you operate the buttons by pressing the letter highlighted in red, or by pressing the tab-key until the wanted button is reached and then pressing the space- or the enter-key.
I.	You want to enable the "Remote Access" option on the Main screen.
II.	If you would like the toolkit to default to using the smaller SRA Lite format with simplified quality scores, set the "Prefer SRA Lite files with simplified base quality scores" option on the Main screen.
III.	Proceed to the "Cache" tab where you will want to enable "local file-caching" and you want to set the "Location of user-repository".
IV.	The repository directory needs to be set to an empty folder. This is the folder where prefetch will deposit the files.
V.	Go to your cloud provider tab and accept to "report cloud instance identity".
>prefetch SRR5416919
>srapath SRR000001
>fastq-dump SRR000001
>ls -l SRR5416919.fastq
>wc -l SRR5416919.fastq
>cache-mgr --report
>cache-mgr --clear

7.	CANU assembly for PACbio data in HCC: 
First step is correction of fastq file.
>canu -correct -p camellia -d assembly genomeSize=1.1m -pacbio /common/yinlab/kshahzad2/camellia_sinensis_var_assamica/camellia_Reads.fasta.gz
Then, trim the output of the correction:
>canu -trim -p camellia -d assembly1   genomeSize=1.1m -corrected -pacbio /common/yinlab/kshahzad2/camellia_sinensis_var_assamica/camellia.correctedReads.fasta.gz
And finally, assemble the output of trimming, twice, with different stringency on which overlaps to use:
>canu -p camellia -d assembly-0.039 genomeSize=1.1m correctedErrorRate=0.039 -trimmed -corrected -pacbio /common/yinlab/kshahzad2/camellia_sinensis_var_assamica/cammelia.trimmedReads.fasta.gz
>canu -p camellia -d assembly-0.075 genomeSize=1.1m correctedErrorRate=0.075 -trimmed -corrected -pacbio /common/yinlab/kshahzad2/camellia_sinensis_var_assamica/cammelia.trimmedReads.fasta.gz

8.	Blast in Glu and HCC: 
>makeblastdb -in mito_db.fasta -out camellia -dbtype nucl -title plant -parse_seqids
>blastn -db camellia -query scaffolds.fasta -task blastn -dust no > output.txt

9.	GetOrganelle for Mitochondrial assembly in HCC:
>get_organelle_config.py --add embplant_pt,embplant_mt
>get_organelle_from_reads.py -s seed.fasta -1 forward.fq -2 reverse.fq -a chloroplast.fasta --genes MT_coding_sequnce.fasta -o mitochondria_output -R 20 -k 21,45,65,85,105 -P 1000000 -F embplant_mt
GetOrganelle in Glu commands:
>get_organelle_from_reads.py -1 R1.clean.fastq -2 R2.clean.fastq -s mitochondrion.fasta -a chloroplast.fasta --genes MT_coding_sequnce.fasta -o mitochondria_Sample1 -R 30 -t 20 -k 21,35,45,55,75,95,105,115,127 -F plant_mt

10.	Velvet for Mitochondrial assembly in HCC:
>module load velvet/1.2
>velveth output_directory/ 43 -fastq -longPaired -separate input_reads_pair_1.fastq input_reads_pair_2.fastq
>velvetg output_directory/ -min_contig_lgth 200

For mitochondrial Annotation: (http://www.1kmpg.cn/mgavas/)
![image](https://github.com/Khurrams569/Chloroplast-and-Mitochondrial-Genome-Assembly-of-Illumina-sequence-/assets/165841830/a739b958-ec3c-49d1-9f4d-8645e580f4b8)


*To convert .sam into .bam file use this command in samtools in HCC:
>module load samtools

>samtools view -S -b /input_Path/filename.sam > /output_path/filename.bam
>convert .sam file into .fastq file by following command:
>samtools fastq -F 4 L3_568_mapped.sam > L3_568_mapped.fastq

For install any software use the below command.
Tools used Bowtie; See the commands below in Glu and Gly;
	
TO add channels
>conda config --show channels
>conda config --add channels bioconda
to create environments 
>conda env list
>conda create -n bioinformatics
>conda activate bioinformatics

After add environment you need to install software’s. Every time when you need to use software first you have to activate the environment by above commands then use software. 
>conda install -c bioconda bowtie
>conda create -n assembly
>conda create -n assembly nanoplot flye bandage bwa samtools pilon
>conda activate assembly
Convert .sam file into .fastq file by following command:
>samtools fastq -F 4 L3_568_mapped.sam > L3_568_mapped.fastq

To install NOVOplasty in Gly Glu and HCC, use this command: 
>git clone https://github.com/ndierckx/NOVOPlasty.git
Make a config.txt file by following link: https://github.com/ndierckx/NOVOPlasty/blob/master/config.txt

Check the sequence in Blast online: 
https://blast.ncbi.nlm.nih.gov/Blast.cgi
Use Nucleotide BLAST option. Check the box in front of “Align two or more sequence”.
Upload assembled file in “Enter Query Sequence”
Upload reference file in “Enter Subject Sequence”
Click on BLAST
The results will show 99.75% similarity sequence. 

To blast the sequence in HCC:
>blastn -query input_reads.fasta -db input_reads_db -out blastn_output.alignments

Check the sequence in mapped file: 
>less L4_568_mapped.sam

Use seqtk for extract mapped sequence from original sequence files.
>seqtk subseq L1_542_reverse_paired.fq.gz L1_542_map_reverse_id > L1_542_mapped_mito_reverse.fastq &




![image](https://github.com/Khurrams569/khurrams/assets/165841830/f78c50d6-0cad-4e22-9119-92344eb78d18)
