#!/usr/bin/bash
#Code written by B217754-2022 and adapted from NGG Tutorial 3 to take four ONT libraries and produce a QC report,
# assembly and assembly assessment, in a mirroring of Filipović, I., Rašić, G., Hereward, J. et al. 'A high-quality 
#de novo genome assembly based on nanopore sequencing of a wild-caught coconut rhinoceros beetle (Oryctes rhinoceros) 
#BMC Genomics 23, 426 (2022). https://doi.org/10.1186/s12864-022-08628-z'

#Make and use working directory
mkdir ICA1 && cd ./ICA1

#install conda env
/localdisk/home/ubuntu-software/Anaconda3/bin/conda-env create -n NGG3 --file /localdisk/data/NGG/conda_envs/NGG3.yml
#activate conda env
source /localdisk/home/ubuntu-software/Anaconda3/bin/activate NGG3
#Adjust seaborn to enable proper plotting in NanoPlot (errors arose otherwise during NanoPlot)
conda install seaborn=0.10.1

#Get the files from SRA
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15436048/SRR15436048
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15436049/SRR15436049
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15436050/SRR15436050
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15436051/SRR15436051

#extract fastq from the SRA files, gzip it (takes a long time)
fastq-dump --gzip SRR15436048 SRR15436049 SRR15436050 SRR15436051

#Combine the SRA libraries into a single file (takes a long time)
cat SRR15436048.fastq.gz SRR15436049.fastq.gz SRR15436050.fastq.gz SRR15436051.fastq.gz > merged.fastq.gz

#Run Nanoplot for QC of reads, filtering out Phred scsores less than 8, transform data logarithmically for 
#better analysis
NanoPlot --verbose --fastq ./merged.fastq.gz -o Nanoplot_Phred_cutoff_8 --loglength --N50 -t175 #--minqual 8
#Run Readbean to assemble the reads (resource intensive)
wtdbg2 -x ont -i ./merged.fastq.gz -o ./Orhinoceros_wtdg2 -g377m -t170 
#Generated the contigs
wtpoa-cns -t170 -i ./Orhinoceros_wtdg2.ctg.lay.gz -fo ./Orhinoceros_wtdg2.ctg.fa
#Determined the quality of the assembly with BBMap's command stats.sh (gives number of contigs and N50 among others)
stats.sh ./Orhinoceros_wtdg2.ctg.fa > Assembly_assessment.txt