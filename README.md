# ReSeq-paper
Snakemake file to reproduce figures and tables in the ReSeq paper

This included the pipeline used to compare ReSeq, pIRS, NEAT and ART on 8 datasets, a bwa bowtie2 comparison with simulated data using ReSeq and various small things.

Only scripts not included in any software, screenshots and csv files that were created with modified ReSeq code are provided. The datasets, references and software has to be downloaded from their original source.

The files not included, but necessary to reproduce the figures in the paper are:
```
input/ecoli/reference/GCF_000005845.2_ASM584v2_genomic.fa
input/ecoli/SRR3191692/SRR3191692_1.fastq.gz
input/ecoli/SRR3191692/SRR3191692_2.fastq.gz
input/ecoli/SRR490124/SRR490124_2.fastq.gz
input/ecoli/SRR490124/SRR490124_1.fastq.gz
input/bcereus/reference/GCF_000007825.1_ASM782v1_genomic.fa
input/mouse/reference/GCF_000001635.26_GRCm38.p6_genomic.fa
input/mouse/ERR3085830/ERR3085830_1.fastq.gz
input/mouse/ERR3085830/ERR3085830_2.fastq.gz
input/athaliana/reference/ATgenomeTAIR9.171.fa
input/athaliana/ERR2017816/ERR2017816_2.fastq.gz
input/athaliana/ERR2017816/ERR2017816_1.fastq.gz
input/baseSpace/Nextera-Repeat-600pM-2x151-DI/* (The complete directory downloaded from Illumina BaseSpace)
input/human/ERR1955542/ERR1955542_2.fastq.gz
input/human/ERR1955542/ERR1955542_1.fastq.gz
input/human/reference/GRCh38_latest_genomic.fa
input/rsphaeroides/reference/GCF_000012905.2_ASM1290v2_genomic.fa
bin/pilon-1.21.jar
bin/neat-genReads.py (genReads.py in the NEAT git repository)
bin/neat/* (utilities/* in the NEAT git repository)
bin/quast/* (quast program folder)
```
