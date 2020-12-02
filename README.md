# ReSeq-paper
Snakemake file to reproduce figures and tables in the ReSeq paper

This included the pipeline used to compare ReSeq, pIRS, NEAT and ART on 8 datasets, a bwa bowtie2 comparison with simulated data using ReSeq and various small things.

Only scripts not included in any software, screenshots and csv files that were created with modified ReSeq code are provided. The datasets, references and software has to be downloaded from their original source.

The files not included, but necessary to reproduce the figures in the paper are:
```
bin/BEAR/adapterDB.fna
bin/BEAR/drisee.py
bin/BEAR/error_models.py
bin/BEAR/error_quality.pl
bin/BEAR/error_quality.py
bin/BEAR/generate_reads.py
bin/BEAR/parametric_abundance.pl
bin/BEAR/qiime-uclust.py
bin/BEAR/run_find_steiner.pl
bin/BEAR/seq_length_stats.py
bin/BEAR/trim_reads.pl
bin/neat-genReads.py (genReads.py in the NEAT git repository)
bin/neat/* (utilities/* in the NEAT git repository)
bin/pilon-1.21.jar
bin/quast/* (quast program folder)
input/athaliana/ERR2017816/ERR2017816_1.fastq.gz
input/athaliana/ERR2017816/ERR2017816_2.fastq.gz
input/athaliana/PRJNA562949/combined_1.fq.gz
input/athaliana/PRJNA562949/combined_2.fq.gz
input/athaliana/reference/ATgenomeTAIR9.171.fa
input/baseSpace/Nextera-Repeat-600pM-2x151-DI/* (The complete directory downloaded from Illumina BaseSpace)
input/bcereus/reference/GCF_000007825.1_ASM782v1_genomic.fa
input/ecoli/DRR058060/DRR058060_1.fastq.gz
input/ecoli/DRR058060/DRR058060_2.fastq.gz
input/ecoli/reference/GCF_000005845.2_ASM584v2_genomic.fai
input/ecoli/SRR3191692/SRR3191692_1.fastq.gz
input/ecoli/SRR3191692/SRR3191692_2.fastq.gz
input/ecoli/SRR490124/SRR490124_2.fastq.gz
input/ecoli/SRR490124/SRR490124_1.fastq.gz
input/human/ERR1955542/ERR1955542_2.fastq.gz
input/human/ERR1955542/ERR1955542_1.fastq.gz
input/human/PRJEB33197/combined_1.fq.gz (COLO829BL control samples)
input/human/PRJEB33197/combined_2.fq.gz (COLO829BL control samples)
input/human/reference/GRCh38_latest_genomic.fa
input/mouse/ERR3085830/ERR3085830_1.fastq.gz
input/mouse/ERR3085830/ERR3085830_2.fastq.gz
input/mouse/reference/GCF_000001635.26_GRCm38.p6_genomic.fa
input/rsphaeroides/reference/GCF_000012905.2_ASM1290v2_genomic.fa
```
