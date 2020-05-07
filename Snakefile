shell.prefix("sleep 5; ")

REFERENCES = {}
REFERENCES['ecoli'] = {'prefix':"GCF_000005845.2_ASM584v2_genomic", 'diploid':False, 'use_jellyfish':1}
REFERENCES['bcereus'] = {'prefix':"GCF_000007825.1_ASM782v1_genomic", 'diploid':False, 'use_jellyfish':1}
REFERENCES['rsphaeroides'] = {'prefix':"GCF_000012905.2_ASM1290v2_genomic", 'diploid':False, 'use_jellyfish':1}
REFERENCES['athaliana'] = {'prefix':"ATgenomeTAIR9.171", 'diploid':True, 'use_jellyfish':1}
REFERENCES['mouse'] = {'prefix':"GCF_000001635.26_GRCm38.p6_genomic", 'diploid':True, 'use_jellyfish':0}
REFERENCES['human'] = {'prefix':"GRCh38_latest_genomic", 'diploid':True, 'use_jellyfish':0}

REFERENCES['phix'] = {'prefix':"phix", 'diploid':False}

#REFERENCES['bovine'] = {'prefix':"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic", 'diploid':True}

# Getting coverage:
#Peak:# samtools stats -c 1,20000,1 storage/ecoli/SRR490124/real/mapping-bowtie2-s.bam | grep ^COV | cut -f 3- | awk 'BEGIN{max=0;len=0;}($2 > max && $1>20){max=$2; len=$1}END{print len}'
#Mean:# samtools stats -c 1,50000,1 storage/ecoli/SRR490124/real/mapping-bowtie2-s.bam | grep ^COV | cut -f 3- | awk 'BEGIN{sum=0;count=0;}{sum+=$1*$2; count+=$2}END{print sum/count}'
#Median:# samtools stats -c 1,50000,1 storage/ecoli/SRR490124/real/mapping-bowtie2-s.bam | grep ^COV | cut -f 3- | tac | awk '{count += $2; print $0}END{print int(count/2)}' | tac | awk '{if(NR==1){halfcount = $1}else{count += $2; if(count >= halfcount){print $1}}}' | head -n1

# DATASET       PEAK    MEAN        MEDIAN
# SRR490124     429     461.179     449
# SRR3191692    987     1031.23     1011
#   assembly    21      298.788     13
# S5L001        493     508.108     508
# S1L001        1491    1474.96     1528
# S9L001        3016    2871.54     2901
# ERR2017816    31      39.147      31
# ERR3085830    44      43.3827     43
# ERR1955542    41      39.6098     40

#Getting indel rates:
# samtools view -f 67 -q 10 storage/ecoli/SRR490124/real/mapping-bowtie2-s.bam | awk '{gsub("M","M ", $6);gsub("I","I ", $6);gsub("D","D ", $6);gsub("S","S ", $6);print $6}' | awk '{for(i=1; i<=NF; i++){count[substr($i,length($i),1)] += substr($i,1,length($i)-1)}}END{print "M=" count["M"], "I=" count["I"], "D=" count["D"], "S=" count["S"], "IR1=" count["I"]/(count["M"]+count["D"]), "DR1=" count["D"]/(count["M"]+count["D"])}'
# samtools view -f 131 -q 10 storage/ecoli/SRR490124/real/mapping-bowtie2-s.bam | awk '{gsub("M","M ", $6);gsub("I","I ", $6);gsub("D","D ", $6);gsub("S","S ", $6);print $6}' | awk '{for(i=1; i<=NF; i++){count[substr($i,length($i),1)] += substr($i,1,length($i)-1)}}END{print "M=" count["M"], "I=" count["I"], "D=" count["D"], "S=" count["S"], "IR2=" count["I"]/(count["M"]+count["D"]), "DR2=" count["D"]/(count["M"]+count["D"])}'

SAMPLES = {}
SAMPLES['SRR490124'] = {'fq1':"SRR490124_1.fastq.gz", 'fq2':"SRR490124_2.fastq.gz", 'species':"ecoli", 'readlength':100, 'insertlength':164, 'insertsd':30, 'coverage':449, 'insertion_rate1':7.6e-05, 'insertion_rate2':8.0e-05, 'deletion_rate1':1.1e-05, 'deletion_rate2':1.7e-05, 'max_kmer_count_plotted':500}
SAMPLES['SRR3191692'] = {'fq1':"SRR3191692_1.fastq.gz", 'fq2':"SRR3191692_2.fastq.gz", 'species':"ecoli", 'readlength':126, 'insertlength':408, 'insertsd':83, 'coverage':1011, 'insertion_rate1':2.0e-05, 'insertion_rate2':3.6e-05, 'deletion_rate1':2.8e-05, 'deletion_rate2':2.9e-05, 'max_kmer_count_plotted':1250, 'insertlength_asm':406, 'insertsd_asm':83, 'coverage_asm':13, 'insertion_rate1_asm':3.1e-04, 'insertion_rate2_asm':3.2e-04, 'deletion_rate1_asm':2.6e-05, 'deletion_rate2_asm':2.7e-05}
SAMPLES['S5L001'] = {'fq1':"Ecoli1_L001_S5_L001_R1_001.fastq.gz", 'fq2':"Ecoli1_L001_S5_L001_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L1", 'species':"ecoli", 'readlength':151, 'insertlength':152, 'insertsd':180, 'coverage':508, 'insertion_rate1':2.9e-05, 'insertion_rate2':2.8e-05, 'deletion_rate1':8.2e-06, 'deletion_rate2':8.6e-06, 'max_kmer_count_plotted':750}
SAMPLES['S1L001'] = {'fq1':"Bcereus1_L001_S1_L001_R1_001.fastq.gz", 'fq2':"Bcereus1_L001_S1_L001_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L1", 'species':"bcereus", 'readlength':151, 'insertlength':152, 'insertsd':180, 'coverage':1528, 'insertion_rate1':1.5e-04, 'insertion_rate2':1.4e-04, 'deletion_rate1':1.8e-04, 'deletion_rate2':1.7e-04, 'max_kmer_count_plotted':3000}
SAMPLES['S9L001'] = {'fq1':"Rsphaeroides1_L001_S9_L001_R1_001.fastq.gz", 'fq2':"Rsphaeroides1_L001_S9_L001_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L1", 'species':"rsphaeroides", 'readlength':151, 'insertlength':152, 'insertsd':180, 'coverage':2901, 'insertion_rate1':5.8e-04, 'insertion_rate2':5.3e-04, 'deletion_rate1':5.6e-05, 'deletion_rate2':5.7e-05, 'max_kmer_count_plotted':5000}
SAMPLES['ERR2017816'] = {'fq1':"ERR2017816_1.fastq.gz", 'fq2':"ERR2017816_2.fastq.gz", 'species':"athaliana", 'readlength':150, 'insertlength':319, 'insertsd':59, 'coverage':31, 'insertion_rate1':8.3e-05, 'insertion_rate2':9.0e-05, 'deletion_rate1':4.3e-05, 'deletion_rate2':5.2e-05, 'max_kmer_count_plotted':150}
SAMPLES['ERR3085830'] = {'fq1':"ERR3085830_1.fastq.gz", 'fq2':"ERR3085830_2.fastq.gz", 'species':"mouse", 'readlength':151, 'insertlength':406, 'insertsd':93, 'coverage':43, 'insertion_rate1':3.5e-04, 'insertion_rate2':3.4e-04, 'deletion_rate1':2.9e-04, 'deletion_rate2':2.8e-04, 'max_kmer_count_plotted':150}
SAMPLES['ERR1955542'] = {'fq1':"ERR1955542_1.fastq.gz", 'fq2':"ERR1955542_2.fastq.gz", 'species':"human", 'readlength':150, 'insertlength':475, 'insertsd':103, 'coverage':40, 'insertion_rate1':1.5e-04, 'insertion_rate2':1.5e-04, 'deletion_rate1':1.3e-04, 'deletion_rate2':1.2e-04, 'max_kmer_count_plotted':100}

#SAMPLES['fgcz-loh2'] = {'fq1':"20170918.A-loh2_R1.fastq.gz", 'fq2':"20170918.A-loh2_R2.fastq.gz", 'species':"athaliana", 'readlength':151, 'insertlength':673, 'insertsd':163, 'coverage':103, 'insertion_rate1':2.2e-04, 'insertion_rate2':2.5e-04, 'deletion_rate1':2.1e-04, 'deletion_rate2':2.0e-04}

#SAMPLES['SRR3191728'] = {'fq1':"SRR3191728_1.fastq.gz", 'fq2':"SRR3191728_2.fastq.gz", 'species':"ecoli", 'readlength':126, 'insertlength':396, 'insertsd':87, 'coverage':439, 'insertion_rate1':2.4e-05, 'insertion_rate2':1.2e-04, 'deletion_rate1':2.8e-05, 'deletion_rate2':3.5e-05}
#SAMPLES['fgcz-novaseqtest1'] = {'fq1':"Phix_R1.fastq.gz", 'fq2':"Phix_R2.fastq.gz", 'species':"phix", 'readlength':151, 'insertlength':361, 'insertsd':29, 'coverage':56071}
#SAMPLES['fgcz-HiSeq4000-spikein'] = {'fq1':"Undetermined_S0_L008_R1_001.fastq.gz", 'fq2':"Undetermined_S0_L008_R2_001.fastq.gz", 'species':"phix", 'readlength':151, 'insertlength':362, 'insertsd':27}
#SAMPLES['fgcz-obv31'] = {'fq1':"20180405.A-OBV_31_R1.fastq.gz", 'fq2':"20180405.A-OBV_31_R2.fastq.gz", 'species':"bovine"}

LANEREPS = {}
LANEREPS['S5L001'] = ['S17L002','S29L003']
SAMPLES['S17L002'] = {'fq1':"Ecoli1_L002_S17_L002_R1_001.fastq.gz", 'fq2':"Ecoli1_L002_S17_L002_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L2", 'species':"ecoli"}
SAMPLES['S29L003'] = {'fq1':"Ecoli1_L003_S29_L003_R1_001.fastq.gz", 'fq2':"Ecoli1_L003_S29_L003_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L3", 'species':"ecoli"}

LIBREPS = {}
LIBREPS['S5L001'] = ['S6L001','S7L001']
SAMPLES['S6L001'] = {'fq1':"Ecoli2_L001_S6_L001_R1_001.fastq.gz", 'fq2':"Ecoli2_L001_S6_L001_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L1", 'species':"ecoli"}
SAMPLES['S7L001'] = {'fq1':"Ecoli3_L001_S7_L001_R1_001.fastq.gz", 'fq2':"Ecoli3_L001_S7_L001_R2_001.fastq.gz", 'basespace':"Nextera-Repeat-600pM-2x151-DI/Nextera_L1", 'species':"ecoli"}

SIMULATORS = ['ReSeq','ART','pIRS','NEAT']
# BEAR: Second read is only reversed, not complemented
# SInC: No errors simulated

MAINSAMPLES = ['SRR490124', 'SRR3191692', 'S5L001', 'S1L001', 'S9L001', 'ERR2017816', 'ERR3085830', 'ERR1955542']

def simulatorOutputFile(template_segment, wildcards):
    sim = wildcards.simulator.split('_')[0]
    if 'ART' == sim:
        return "art{}.fq.gz".format(template_segment)
    if 'BEAR' == sim:
        return "bear-R{}.fq.gz".format(template_segment)
    if 'NEAT' == sim:
        return "neat_read{}.fq.gz".format(template_segment)
    if 'pIRS' == sim:
        return "pirs_{}.fq.gz".format(template_segment)
    if 'ReSeq' == sim:
        return "reseq-R{}.fq.gz".format(template_segment)
    #if 'SInC' == sim:
    #    return "sinc-R{}.fq.gz".format(template_segment)
    print("Simulator '{}' is unknown".format(wildcards.simulator))

jellyfish_k = 51
max_threads = 64
max_mem_mb = 419840

rule all:
    input:
        "storage/paper/tables/table_covcorrelation_spearman.txt",
        "storage/paper/figures/figure_comp_resources.pdf",
        "storage/paper/figures/figure_coverage.pdf",
        "storage/paper/figures/figure_covfit.pdf",
        "storage/paper/figures/figure_kmer.pdf",
        "storage/paper/figures/figure_kmer_cross.pdf",
        "storage/paper/figures/figure_mapper_comp.pdf",
        "storage/paper/figures/figure_preqc.pdf",
        "storage/paper/figures/figure_syserror.pdf",
        "storage/paper/additional_figures/figure_covcorrelation.pdf",
        "storage/paper/additional_figures/figure_dispersion_pars_overview.pdf"
        "storage/paper/additional_figures/figure_error_rate.pdf"
        "storage/paper/additional_figures/figure_gcbias_consistency.pdf",
        "storage/paper/additional_figures/figure_gc_logit_vs_log.pdf",
        "storage/paper/additional_figures/figure_kmer_cross_mouse_human_seqbias.pdf",
        "storage/paper/additional_figures/figure_kmer_SRR3191692_assembly_cov.pdf",
        "storage/paper/additional_figures/figure_nobias.pdf",
        "storage/paper/additional_figures/figure_preqc_SRR3191692_assembly_cov.pdf",
        "storage/paper/additional_figures/figure_quality_values.pdf"
        "storage/paper/additional_figures/figure_start_end_surrounding.pdf",
        expand("storage/summary/reseq-plots/{sample}.pdf", sample=MAINSAMPLES),
        ["storage/{0}/{1}/real/mapping-bowtie2-s-insert-length.pdf".format(SAMPLES[sample]['species'], sample) for sample in MAINSAMPLES]

###-Create figures and tables for paper-------------------------------------------------------------------------------------------------------------------------------------------------

rule figure_surbias:
    input:
        "input/paper/csv/S5L001_sum_maxlike_fit.csv",
        "input/paper/csv/S1L001_sum_maxlike_fit.csv",
        "input/paper/csv/S9L001_sum_maxlike_fit.csv"
    output:
        "storage/paper/subfigures/figure_surbias.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_surbias.R {params.loaddir} {output}
        """

rule figure_gcbias:
    input:
        "input/paper/csv/SRR490124_sum_maxlike_fit.csv"
    output:
        "storage/paper/subfigures/figure_gcbias.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_gcbias.R {params.loaddir} {output}
        """
        
rule figure_gcbias_consistency:
    input:
        "input/paper/csv/S5L001_sum_maxlike_fit.csv",
        "input/paper/csv/S1L001_sum_maxlike_fit.csv",
        "input/paper/csv/S9L001_sum_maxlike_fit.csv"
    output:
        "storage/paper/additional_figures/figure_gcbias_consistency.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_gcbias_consistency.R {params.loaddir} {output}
        """

rule figure_nobias:
    input:
        "storage/ecoli/S5L001/ReSeq/nobias/stats.pdf"
    output:
        "storage/paper/additional_figures/figure_nobias.pdf"
    params:
        workdir="work/paper/figures/figure_nobias"
    shell:
        """
        mkdir -p {params.workdir}
        pdftk {input} cat 14 output {params.workdir}/surbias_full.pdf
        pdfcrop {params.workdir}/surbias_full.pdf {params.workdir}/surbias.pdf 1>/dev/null
        pdftk {input} cat 11 output {params.workdir}/gcbias_full.pdf
        pdfcrop {params.workdir}/gcbias_full.pdf {params.workdir}/gcbias.pdf 1>/dev/null
        Rscript bin/figure_nobias.R {params.workdir}
        pdfcrop {params.workdir}/figure_nobias.pdf {output} 1>/dev/null
        """

rule figure_samplemean:
    input:
        "input/paper/csv/SRR490124_sum_dispersion_fit.csv",
        "input/paper/csv/SRR490124_mult_dispersion_fit.csv",
        "input/paper/csv/S5L001_sum_dispersion_fit.csv",
        "input/paper/csv/S5L001_mult_dispersion_fit.csv"
    output:
        "storage/paper/subfigures/figure_samplemean.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_samplemean.R {params.loaddir} {output}
        """
        
rule figure_dispersion:
    input:
        "input/paper/csv/SRR490124_sum_dispersion_fit.csv",
        "input/paper/csv/S5L001_sum_dispersion_fit.csv"
    output:
        "storage/paper/subfigures/figure_dispersion_a.pdf",
        "storage/paper/subfigures/figure_dispersion_b.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_dispersion.R {params.loaddir} {output}
        """
        
rule figure_covfit:
    input:
        surbias="storage/paper/subfigures/figure_surbias.pdf",
        gcbias="storage/paper/subfigures/figure_gcbias.pdf",
        samplemean="storage/paper/subfigures/figure_samplemean.pdf",
        dispa="storage/paper/subfigures/figure_dispersion_a.pdf",
        dispb="storage/paper/subfigures/figure_dispersion_b.pdf"
    output:
        "storage/paper/figures/figure_covfit.pdf"
    params:
        workdir="work/paper/figures/figure_covfit"
    shell:
        """
        mkdir -p {params.workdir}
        pdfcrop --margins '0 0 1 0' {input.surbias} {params.workdir}/surbias.pdf 1>/dev/null
        pdfcrop {input.gcbias} {params.workdir}/gcbias.pdf 1>/dev/null
        pdfcrop {input.samplemean} {params.workdir}/samplemean.pdf 1>/dev/null
        pdfcrop {input.dispa} {params.workdir}/dispa.pdf 1>/dev/null
        pdfcrop {input.dispb} {params.workdir}/dispb.pdf 1>/dev/null
        Rscript bin/figure_covfit.R {params.workdir}
        pdfcrop {params.workdir}/figure_covfit.pdf {output} 1>/dev/null
        """

rule figure_coverage:
    input:
        SRR490124="storage/summary/reseq-plots/SRR490124.pdf",
        S5L001="storage/summary/reseq-plots/S5L001.pdf"
    output:
        "storage/paper/figures/figure_coverage.pdf"
    params:
        workdir="work/paper/figures/figure_coverage"
    shell:
        """
        mkdir -p {params.workdir}
        pdftk {input.SRR490124} cat 1 output {params.workdir}/SRR490124-cov-full.pdf
        pdfcrop {params.workdir}/SRR490124-cov-full.pdf {params.workdir}/SRR490124-cov.pdf 1>/dev/null
        pdftk {input.SRR490124} cat 5 output {params.workdir}/SRR490124-dup-full.pdf
        pdfcrop {params.workdir}/SRR490124-dup-full.pdf {params.workdir}/SRR490124-dup.pdf 1>/dev/null
        pdftk {input.S5L001} cat 1 output {params.workdir}/S5L001-cov-full.pdf
        pdfcrop {params.workdir}/S5L001-cov-full.pdf {params.workdir}/S5L001-cov.pdf 1>/dev/null
        pdftk {input.S5L001} cat 5 output {params.workdir}/S5L001-dup-full.pdf
        pdfcrop {params.workdir}/S5L001-dup-full.pdf {params.workdir}/S5L001-dup.pdf 1>/dev/null
        Rscript bin/figure_coverage.R {params.workdir}
        pdfcrop {params.workdir}/figure_coverage.pdf {output} 1>/dev/null
        """
        
rule table_covcorrelation:
    input:
        lane="storage/ecoli/S5L001/real/correlation/{method}-coverage-lanereps.txt",
        library="storage/ecoli/S5L001/real/correlation/{method}-coverage-libreps.txt",
        simreal=expand("storage/ecoli/S5L001/{simulator}/correlation/{{method}}-coverage-real.txt",simulator=SIMULATORS),
        simself=expand("storage/ecoli/S5L001/{simulator}/correlation/{{method}}-coverage-simreps.txt",simulator=SIMULATORS)
    output:
        "storage/paper/tables/table_covcorrelation_{method}.txt"
    params:
        simulator=SIMULATORS
    shell:
        """
        echo "\\begin{{tabular}}{{l|cccc}}" > {output}
        echo "& Lane & Library & Real & Simulation\\\\\\\\" >> {output}
        echo "\\hline" >> {output}
        cat {input.lane} | awk 'BEGIN{{sum=0;count=0}}{{sum+=$3;count+=1}}END{{printf "Real & $%.2f$", sum/count}}' >> {output}
        cat {input.library} | awk 'BEGIN{{sum=0;count=0}}{{sum+=$3;count+=1}}END{{printf " & $%.2f$ & - & -\\\\\\\\", sum/count; print ""}}' >> {output}
        for SIM in {params.simulator}; do
            cat storage/ecoli/S5L001/$SIM/correlation/{wildcards.method}-coverage-real.txt | awk -v sim=$SIM 'BEGIN{{sum=0;count=0}}{{sum+=$3;count+=1}}END{{printf "%s & - & - & $%.2f$", sim, sum/count}}' >> {output}
            cat storage/ecoli/S5L001/$SIM/correlation/{wildcards.method}-coverage-simreps.txt | awk 'BEGIN{{sum=0;count=0}}{{sum+=$3;count+=1}}END{{printf " & $%.2f$\\\\\\\\", sum/count; print ""}}' >> {output}
        done
        echo "\\end{{tabular}}" >> {output}
        """

rule figure_covcorrelation_summary:
    input:
        lane="storage/ecoli/S5L001/real/correlation/pearson-coverage-lanereps.pdf",
        lib="storage/ecoli/S5L001/real/correlation/pearson-coverage-libreps.pdf",
        simreal=expand("storage/ecoli/S5L001/{simulator}/correlation/pearson-coverage-real.pdf",simulator=SIMULATORS),
        simself=expand("storage/ecoli/S5L001/{simulator}/correlation/pearson-coverage-simreps.pdf",simulator=SIMULATORS)
    output:
        "storage/paper/additional_figures/figure_covcorrelation.pdf"
    params:
        workdir="work/paper/figures/figure_covcorrelation",
        simulator=SIMULATORS
    shell:
        """
        mkdir -p {params.workdir}
        pdfcrop {input.lane} {params.workdir}/cor1.pdf 1>/dev/null
        pdfcrop {input.lib} {params.workdir}/cor2.pdf 1>/dev/null
        i=3
        for SIM in {params.simulator}; do
            pdfcrop storage/ecoli/S5L001/$SIM/correlation/pearson-coverage-real.pdf {params.workdir}/cor${{i}}.pdf 1>/dev/null
            i=$((i + 1))
            pdfcrop storage/ecoli/S5L001/$SIM/correlation/pearson-coverage-simreps.pdf {params.workdir}/cor${{i}}.pdf 1>/dev/null
            i=$((i + 1))
        done
        Rscript bin/figure_covcorrelation_summary.R {params.workdir}
        pdfcrop {params.workdir}/figure_covcorrelation.pdf {output} 1>/dev/null
        """
        
rule figure_syserror:
    input:
        SRR490124="storage/summary/reseq-plots/SRR490124.pdf",
        S5L001="storage/summary/reseq-plots/S5L001.pdf",
        SRR490124screen="input/paper/screenshots/SRR490124-NC_000913.3_pilon-1677350-1677450.png",
        S5L001screen="input/paper/screenshots/S5L001-NC_000913.3_pilon-3768350-3768450.png",
        repscreen="input/paper/screenshots/S5L001-S17L002-S29L003-S6L001-S7L001-NC_000913.3_pilon-3768350-3768450-mod.png"
    output:
        "storage/paper/figures/figure_syserror.pdf"
    params:
        workdir="work/paper/figures/figure_syserror"
    shell:
        """
        mkdir -p {params.workdir}
        pdftk {input.SRR490124} cat 42 output {params.workdir}/SRR490124-errors-full.pdf
        pdfcrop {params.workdir}/SRR490124-errors-full.pdf {params.workdir}/SRR490124-errors.pdf 1>/dev/null
        pdftk {input.S5L001} cat 42 output {params.workdir}/S5L001-errors-full.pdf
        pdfcrop {params.workdir}/S5L001-errors-full.pdf {params.workdir}/S5L001-errors.pdf 1>/dev/null
        ln -srf {input.SRR490124screen} {params.workdir}/SRR490124-screen.png
        ln -srf {input.S5L001screen} {params.workdir}/S5L001-screen.png
        ln -srf {input.repscreen} {params.workdir}/rep-screen.png
        Rscript bin/figure_syserror.R {params.workdir}
        pdfcrop {params.workdir}/figure_syserror.pdf {output} 1>/dev/null
        """
        
rule figure_kmer:
    input:
        expand("storage/summary/kmer/{sample}.pdf",sample=["SRR490124", "SRR3191692", "SRR3191692_assembly", "S5L001", "S1L001", "S9L001","ERR2017816","ERR3085830","ERR1955542"])
    output:
        "storage/paper/figures/figure_kmer.pdf"
    params:
        workdir="work/paper/figures/figure_kmer"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_kmer.pdf {params.workdir}/figure[1-9].pdf
        pdfcrop {params.workdir}/figure_kmer.pdf {output} 1>/dev/null
        """
        
rule figure_kmer_cross:
    input:
        expand("storage/summary/kmer_cross/{sample}.pdf",sample=["ecoli_S5L001_rsphaeroides_S9L001", "bcereus_S1L001_ecoli_S5L001", "rsphaeroides_S9L001_ecoli_S5L001", "athaliana_ERR2017816_mouse_ERR3085830", "mouse_ERR3085830_human_ERR1955542", "human_ERR1955542_mouse_ERR3085830"])
    output:
        "storage/paper/figures/figure_kmer_cross.pdf"
    params:
        workdir="work/paper/figures/figure_kmer_cross"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_kmer_cross.pdf {params.workdir}/figure[1-9].pdf
        pdfcrop {params.workdir}/figure_kmer_cross.pdf {output} 1>/dev/null
        """

rule figure_kmer_cov:
    input:
        "storage/summary/kmer/SRR3191692_assembly_cov.pdf"
    output:
        "storage/paper/additional_figures/figure_kmer_SRR3191692_assembly_cov.pdf"
    shell:
        "pdfcrop {input} {output} 1>/dev/null"

rule figure_kmer_cross_seqbias:
    input:
        "storage/summary/kmer_cross/mouse_ERR3085830_human_ERR1955542_seqbias.pdf"
    output:
        "storage/paper/additional_figures/figure_kmer_cross_mouse_human_seqbias.pdf"
    shell:
        "pdfcrop {input} {output} 1>/dev/null"

rule figure_preqc:
    input:
        expand("storage/summary/preqc/{sample}.pdf",sample=["SRR490124", "SRR3191692", "SRR3191692_assembly", "S5L001", "S1L001", "S9L001","ERR2017816","ERR3085830","ERR1955542"])
    output:
        "storage/paper/figures/figure_preqc.pdf"
    params:
        workdir="work/paper/figures/figure_preqc"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdftk $fi cat 4 output {params.workdir}/figure${{i}}-full.pdf
            pdfcrop {params.workdir}/figure${{i}}-full.pdf {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_preqc.pdf {params.workdir}/figure[1-9].pdf
        pdfcrop {params.workdir}/figure_preqc.pdf {output} 1>/dev/null
        """

rule figure_preqc_cov:
    input:
        "storage/summary/preqc/SRR3191692_assembly_cov.pdf"
    output:
        "storage/paper/additional_figures/figure_preqc_SRR3191692_assembly_cov.pdf"
    params:
        workdir="work/paper/figures/figure_preqc_cov"
    shell:
        """
        mkdir -p {params.workdir}
        pdftk {input} cat 4 output {params.workdir}/preqc_SRR3191692_assembly_cov.pdf
        pdfcrop {params.workdir}/preqc_SRR3191692_assembly_cov.pdf {output} 1>/dev/null
        """
        
rule figure_comp_res:
    input:
        "storage/summary/comp_resources/cpu_training.pdf",
        "storage/summary/comp_resources/elapsed_training.pdf",
        "storage/summary/comp_resources/memory_training.pdf",
        "storage/summary/comp_resources/cpu_simulation.pdf",
        "storage/summary/comp_resources/elapsed_simulation.pdf",
        "storage/summary/comp_resources/memory_simulation.pdf"
    output:
        "storage/paper/figures/figure_comp_resources.pdf"
    params:
        workdir="work/paper/figures/figure_comp_resources"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_comp_resources.pdf {params.workdir}/figure[1-6].pdf
        pdfcrop {params.workdir}/figure_comp_resources.pdf {output} 1>/dev/null
        """

rule figure_map_stats:
    input:
        "work/paper/figures/figure_map_stats/map_stats.csv"
    output:
        "storage/summary/mapper_comp/{sample}_stats.pdf"
    shell:
        "Rscript bin/figure_map_stats.R {input} {wildcards.sample} {output}"

rule figure_map_stats_prepare:
    input:
        real_map=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_real_mapping-{mapper}-s.bam.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bt2_map=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_eval_mapping-{mapper}-s.bam.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bwa_map=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_bwa_eval_mapping-{mapper}-s.bam.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        real_cigar=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_real_mapping-{mapper}-s.bam_cigar.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bt2_cigar=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_eval_mapping-{mapper}-s.bam_cigar.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bwa_cigar=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_bwa_eval_mapping-{mapper}-s.bam_cigar.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"])
    output:
        "work/paper/figures/figure_map_stats/map_stats.csv"
    shell:
        """
        echo "Variable, Simulator, SRR490124_bowtie2, SRR490124_bwa, SRR3191692_bowtie2, SRR3191692_bwa, S5L001_bowtie2, S5L001_bwa" > {output}
        printf "Unmapped pairs, real" >> {output}
        for f in {input.real_map}; do printf ", %d" $(awk '(1==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Unmapped pairs, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_map}; do printf ", %d" $(awk '(1==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Unmapped pairs, ReSeq-bwa" >> {output}
        for f in {input.bwa_map}; do printf ", %d" $(awk '(1==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Single reads, real" >> {output}
        for f in {input.real_map}; do printf ", %d" $(awk '(2==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Single reads, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_map}; do printf ", %d" $(awk '(2==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Single reads, ReSeq-bwa" >> {output}
        for f in {input.bwa_map}; do printf ", %d" $(awk '(2==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Mapped pairs, real" >> {output}
        for f in {input.real_map}; do printf ", %d" $(awk '(3==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Mapped pairs, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_map}; do printf ", %d" $(awk '(3==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Mapped pairs, ReSeq-bwa" >> {output}
        for f in {input.bwa_map}; do printf ", %d" $(awk '(3==NR){{print $1}}' $f) >> {output}; done
        echo >> {output}
        
        printf "Matches, real" >> {output}
        for f in {input.real_cigar}; do printf ", %d" $(awk '{{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Matches, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf ", %d" $(awk '{{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Matches, ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf ", %d" $(awk '{{print $1}}' $f) >> {output}; done
        echo >> {output}
        printf "Soft trimmed, real" >> {output}
        for f in {input.real_cigar}; do printf ", %d" $(awk '{{print $2}}' $f) >> {output}; done
        echo >> {output}
        printf "Soft trimmed, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf ", %d" $(awk '{{print $2}}' $f) >> {output}; done
        echo >> {output}
        printf "Soft trimmed, ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf ", %d" $(awk '{{print $2}}' $f) >> {output}; done
        echo >> {output}
        printf "Insertions, real" >> {output}
        for f in {input.real_cigar}; do printf ", %d" $(awk '{{print $3}}' $f) >> {output}; done
        echo >> {output}
        printf "Insertions, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf ", %d" $(awk '{{print $3}}' $f) >> {output}; done
        echo >> {output}
        printf "Insertions, ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf ", %d" $(awk '{{print $3}}' $f) >> {output}; done
        echo >> {output}
        printf "Deletions, real" >> {output}
        for f in {input.real_cigar}; do printf ", %d" $(awk '{{print $4}}' $f) >> {output}; done
        echo >> {output}
        printf "Deletions, ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf ", %d" $(awk '{{print $4}}' $f) >> {output}; done
        echo >> {output}
        printf "Deletions, ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf ", %d" $(awk '{{print $4}}' $f) >> {output}; done
        echo >> {output}
        """ 

rule table_map_stats:
    input:
        real_map=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_real_mapping-{mapper}-s.bam.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bt2_map=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_eval_mapping-{mapper}-s.bam.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bwa_map=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_bwa_eval_mapping-{mapper}-s.bam.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        real_cigar=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_real_mapping-{mapper}-s.bam_cigar.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bt2_cigar=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_eval_mapping-{mapper}-s.bam_cigar.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"]),
        bwa_cigar=expand(expand("work/paper/tables/table_map_stats/storage_ecoli_{{sample}}_ReSeq_bwa_eval_mapping-{mapper}-s.bam_cigar.txt", mapper=["bowtie2","bwa"]), sample=["SRR490124","SRR3191692","S5L001"])
    output:
        "storage/paper/tables/table_map_stats.txt"
    shell:
        """
        echo "\\begin{{tabular}}{{ll|ccccccccc}}" > {output}
        echo "  & & SRR490124 & & SRR3191692 & & S5L001 & \\\\\\\\" >> {output}
        echo "  & & bowtie2 & bwa & bowtie2 & bwa & bowtie2 & bwa \\\\\\\\" >> {output}
        echo "  \hline" >> {output}
        
        printf "  Unmapped pairs & real" >> {output}
        for f in {input.real_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(1==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(1==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bwa" >> {output}
        for f in {input.bwa_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(1==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "  Single reads & real" >> {output}
        for f in {input.real_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(2==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(2==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bwa" >> {output}
        for f in {input.bwa_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(2==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "  Mapped pairs & real" >> {output}
        for f in {input.real_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(3==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(3==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bwa" >> {output}
        for f in {input.bwa_map}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{sum+=$1;if(3==NR){{val=$1}}}}END{{print val, val*100/sum}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        
        echo "  \hline" >> {output}
        printf "  Matches & real" >> {output}
        for f in {input.real_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $1, $1*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $1, $1*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $1, $1*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "  Soft trimmed & real" >> {output}
        for f in {input.real_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $2, $2*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $2, $2*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $2, $2*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "  Insertion & real" >> {output}
        for f in {input.real_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $3, $3*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $3, $3*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $3, $3*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "  Deletion & real" >> {output}
        for f in {input.real_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $4, $4*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & ReSeq-bowtie2" >> {output}
        for f in {input.bt2_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $4, $4*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}
        printf "   & reseq-bwa" >> {output}
        for f in {input.bwa_cigar}; do printf "%.2e}} (%.0f\\\\%%)" $(awk '{{print $4, $4*100/($1+$2+$3+$4)}}' $f) | awk '{{gsub("e[+]0"," \\\\cdot 10^{{", $0); printf " & $%s$", $0}}' >> {output}; done
        echo " \\\\\\\\" >> {output}

        echo "\\end{{tabular}}" >> {output}
        """
        
rule table_map_stats_prepare:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/mapping-{3}-s.bam".format(wildcards.species, wildcards.sample, "/".join(wildcards.simulator.rsplit("_", 1)), wildcards.mapper)
    output:
        "work/paper/tables/table_map_stats/storage_{species,[0-9A-Za-z]*}_{sample,[0-9A-Za-z]*}_{simulator}_mapping-{mapper}-s.bam.txt"
    shell:
        """
        samtools view -c -f 76 -F 2304 {input} > {output}
        samtools view -c -f 8 -F 2308 {input} >> {output}
        samtools view -c -f 64 -F 2316 {input} >> {output}
        """
        
rule table_map_stats_prepare2:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/mapping-{3}-s.bam".format(wildcards.species, wildcards.sample, "/".join(wildcards.simulator.rsplit("_", 1)), wildcards.mapper)
    output:
        "work/paper/tables/table_map_stats/storage_{species,[0-9A-Za-z]*}_{sample,[0-9A-Za-z]*}_{simulator}_mapping-{mapper}-s.bam_cigar.txt"
    shell:
        """
        samtools view -F 2308 {input} | awk '{{gsub("M","M ", $6);gsub("I","I ", $6);gsub("D","D ", $6);gsub("S","S ", $6);print $6}}' | awk 'BEGIN{{count["M"]=0;count["S"]=0;count["I"]=0;count["D"]=0}}{{for(i=1; i<=NF; i++){{count[substr($i,length($i),1)] += substr($i,1,length($i)-1)}}}}END{{print count["M"], count["S"], count["I"], count["D"]}}' > {output}
        """
        
rule table_map_stats_prepare3_bowtie2:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/mapping-bowtie2-s.bam".format(wildcards.species, wildcards.sample, "/".join(wildcards.simulator.rsplit("_", 1)))
    output:
        "work/paper/tables/table_map_stats/storage_{species,[0-9A-Za-z]*}_{sample,[0-9A-Za-z]*}_{simulator}_mapping-bowtie2-s.bam_errorrate.txt"
    shell:
        """
        samtools view -F 2308 {input} | awk '{{if("NM:i:" == substr($17, 1, 5)){{errors += substr($17, 6, length($17)); bases += length($10)}}else if("NM:i:" == substr($18, 1, 5)){{errors += substr($18, 6, length($18)); bases += length($10)}}else{{print "failure"}}}}END{{print errors, bases}}' > {output}
        """
        
rule table_map_stats_prepare3_bwa:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/mapping-bwa-s.bam".format(wildcards.species, wildcards.sample, "/".join(wildcards.simulator.rsplit("_", 1)))
    output:
        "work/paper/tables/table_map_stats/storage_{species,[0-9A-Za-z]*}_{sample,[0-9A-Za-z]*}_{simulator}_mapping-bwa-s.bam_errorrate.txt"
    shell:
        """
        samtools view -F 2308 {input} | awk '{{if("NM:i:" == substr($12, 1, 5)){{errors += substr($12, 6, length($12)); bases += length($10)}}else{{print "failure"}}}}END{{print errors, bases}}' > {output}
        """

rule figure_mapper_comp_summary:
    input:
        expand("storage/summary/mapper_comp/{sample}_stats.pdf",sample=["SRR490124", "SRR3191692", "S5L001"]),
        expand("storage/summary/mapper_comp/{sample}_mapq.pdf",sample=["SRR490124", "SRR3191692", "S5L001"]),
        expand("storage/summary/mapper_comp/{sample}_correctness.pdf",sample=["SRR490124", "SRR3191692", "S5L001"])
    output:
        "storage/paper/figures/figure_mapper_comp.pdf"
    params:
        workdir="work/paper/figures/figure_mapper_comp"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_mapper_comp.pdf {params.workdir}/figure[1-9].pdf
        pdfcrop {params.workdir}/figure_mapper_comp.pdf {output} 1>/dev/null
        """

rule figure_mapper_comp:
    input:
        "work/paper/figures/figure_mapper_comp/storage_ecoli_{sample}_ReSeq_eval_mapping-bowtie2-s_mapping_correctness.csv",
        "work/paper/figures/figure_mapper_comp/storage_ecoli_{sample}_ReSeq_eval_mapping-bwa-s_mapping_correctness.csv",
        "work/paper/figures/figure_mapper_comp/storage_ecoli_{sample}_ReSeq_bwa_eval_mapping-bowtie2-s_mapping_correctness.csv",
        "work/paper/figures/figure_mapper_comp/storage_ecoli_{sample}_ReSeq_bwa_eval_mapping-bwa-s_mapping_correctness.csv",
        "work/paper/figures/figure_mapper_comp/storage_ecoli_{sample}_real_mapping-bowtie2-s_mapping_qualities.csv",
        "work/paper/figures/figure_mapper_comp/storage_ecoli_{sample}_real_mapping-bwa-s_mapping_qualities.csv"
    output:
        "storage/summary/mapper_comp/{sample}_correctness.pdf",
        "storage/summary/mapper_comp/{sample}_mapq.pdf"
    params:
        workdir="work/paper/figures/figure_mapper_comp",
        prefix="storage/summary/mapper_comp/{sample}"
    shell:
        "Rscript bin/figure_mapper_comp.R {params.workdir} {wildcards.sample} {params.prefix}"

rule figure_mapper_comp_prepare2:
    input:
        "storage/{species}/{sample}/real/mapping-{mapper}-s.bam"
    output:
        "work/paper/figures/figure_mapper_comp/storage_{species}_{sample}_real_mapping-{mapper}-s_mapping_qualities.csv"
    shell:
        """
        echo "mapq, unmapped, count" > {output}
        samtools view -F 2304 {input} | awk '{{print $5, int($2%8/4)}}' | sort -k1,1rn -k2,2n | uniq -c | awk '{{print $2 ", " $3 ", " $1}}' >> {output}
        """

rule figure_mapper_comp_prepare:
    input:
        "storage/{species}/{sample}/ReSeq{mapper2}/eval/mapping-{mapper}-s.bam"
    output:
        "work/paper/figures/figure_mapper_comp/storage_{species}_{sample}_ReSeq{mapper2,|_bwa}_eval_mapping-{mapper}-s_mapping_correctness.csv"
    shell:
        """
        echo "mapq, correctness, unmapped, negative, count" > {output}
        samtools view -F 2304 {input} | awk '{{
            split($1, id, ":");
            if(1==int($2%8/4)){{
                if("Adapter" == id[3]){{
                    print $5, 3, 1, 1
                }}
                else{{
                    print $5, 0, 1, 0
                }}
            }}
            else{{
                if(id[3]!=$3){{
                    if("Adapter" == id[3]){{
                        print $5, 0, 0, 1
                    }}
                    else{{
                        print $5, 0, 0, 0
                    }}
                }}
                else{{
                    soft=0;
                    len=0;
                    s=1;
                    for(i=1; i<=length($6); i++){{
                        if(substr($6, i, 1) ~ /^[A-Z]+$/){{
                            if(substr($6, i, 1) ~ /^[MD]+$/){{
                                len+=substr($6,s,i-s)
                            }};
                            if(substr($6, i, 1) == "S"){{
                                soft=1
                            }};
                            s=i+1
                        }}
                    }};
                    rev = int($2%32/16);
                    if(rev){{
                        start=$4+len-1
                    }}
                    else{{
                        start=$4
                    }};
                    second = int($2%256/128);
                    cor_or=0;
                    if(id[2]==id[4]){{
                        cor_or=1
                    }}
                    else{{
                        if(id[2]<id[4]){{
                            if(second == rev){{
                                cor_or=1
                            }}
                        }}
                        else{{
                            if(second != rev){{
                                cor_or=1
                            }}
                        }}
                    }};
                    if(id[2]<=id[4]){{
                        from=id[2];
                        to=id[4]
                    }}
                    else{{
                        from=id[4];
                        to=id[2]
                    }};
                    overlap=0;
                    if((cor_or && rev) || (!cor_or && !rev)){{
                        if($4+len-1 >= to-length($10) && $4 <= to){{
                            overlap=1
                        }}
                    }}
                    else{{
                        if($4+len-1 >= from && $4 <= from+length($10)){{
                            overlap=1
                        }}
                    }};
                    if(!cor_or){{
                        if(overlap){{
                            print $5, 1, 0, 0
                        }}
                        else{{
                            print $5, 0, 0, 0
                        }}
                    }}
                    else{{
                        if(!overlap){{
                            print $5, 0, 0, 0
                        }}
                        else{{
                            if((rev && start != to) || (!rev && start != from)){{
                                print $5, 1, 0, 0
                            }}
                            else{{
                                if((rev && ($4 < from || ($4 > from && soft))) || (!rev && (start+len-1 > to || (start+len-1 < to && soft)))){{
                                    print $5, 2, 0, 0
                                }}
                                else{{
                                    print $5, 3, 0, 0
                                }}
                            }}
                        }}
                    }}
                }}
            }}
        }}' | sort -k1,1nr -k2,2nr -k3,3n -k4,4n | uniq -c | awk '{{print $2 ", " $3 ", " $4 ", " $5 ", " $1}}' >> {output}
        """

rule figure_start_end_surrounding:
    input:
        "storage/ecoli/SRR490124/real/surrounding.pdf",
        "storage/ecoli/S5L001/real/surrounding.pdf"
    output:
        "storage/paper/additional_figures/figure_start_end_surrounding.pdf"
    params:
        workdir="work/paper/figures/figure_start_end_surrounding"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp4.R {params.workdir}/figure_start_end_surrounding.pdf {params.workdir}/figure[1-4].pdf
        pdfcrop {params.workdir}/figure_start_end_surrounding.pdf {output} 1>/dev/null
        """
        
rule figure_dispersion_pars:
    input:
        "input/paper/csv/SRR490124_sum_maxlike_fit.csv",
        "input/paper/csv/SRR3191692_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_sum_maxlike_fit.csv",
        "input/paper/csv/S1L001_sum_maxlike_fit.csv",
        "input/paper/csv/S9L001_sum_maxlike_fit.csv",
        "input/paper/csv/ERR2017816_sum_maxlike_fit.csv",
        "input/paper/csv/ERR3085830_sum_maxlike_fit.csv",
        "input/paper/csv/ERR1955542_sum_maxlike_fit.csv"
    output:
        "storage/paper/subfigures/figure_dispersion_pars.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_dispersion_pars.R {params.loaddir} {output}
        """

rule figure_dispersion_start_pars:
    input:
        "input/paper/csv/S5L001_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_b0.5_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_b0.75_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_b2.0_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_b10.0_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a10.0_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a0.5_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a0.2_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a0.1_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a0.5_b0.5_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a0.1_b0.5_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_a0.314_b1.24_sum_maxlike_fit.csv"
    output:
        "storage/paper/subfigures/figure_dispersion_start_pars1.pdf",
        "storage/paper/subfigures/figure_dispersion_start_pars2.pdf",
        "storage/paper/subfigures/figure_dispersion_start_pars3.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_dispersion_start_pars.R {params.loaddir} {output}
        """

rule figure_dispersion_pars_overview:
    input:
        "storage/paper/subfigures/figure_dispersion_pars.pdf",
        "storage/paper/subfigures/figure_dispersion_start_pars1.pdf",
        "storage/paper/subfigures/figure_dispersion_start_pars2.pdf",
        "storage/paper/subfigures/figure_dispersion_start_pars3.pdf"
    output:
        "storage/paper/additional_figures/figure_dispersion_pars_overview.pdf"
    params:
        workdir="work/paper/figures/figure_dispersion_pars_overview"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp4.R {params.workdir}/figure_dispersion_pars_overview.pdf {params.workdir}/figure[1-4].pdf
        pdfcrop {params.workdir}/figure_dispersion_pars_overview.pdf {output} 1>/dev/null
        """

rule figure_gc_logit_vs_log:
    input:
        "input/paper/csv/SRR490124_logit_maxlike_fit.csv",
        "input/paper/csv/SRR490124_sum_maxlike_fit.csv",
        "input/paper/csv/SRR3191692_logit_maxlike_fit.csv",
        "input/paper/csv/SRR3191692_sum_maxlike_fit.csv",
        "input/paper/csv/S5L001_logit_maxlike_fit.csv",
        "input/paper/csv/S5L001_sum_maxlike_fit.csv",
        "input/paper/csv/S1L001_logit_maxlike_fit.csv",
        "input/paper/csv/S1L001_sum_maxlike_fit.csv",
        "input/paper/csv/S9L001_logit_maxlike_fit.csv",
        "input/paper/csv/S9L001_sum_maxlike_fit.csv",
        "input/paper/csv/ERR2017816_logit_maxlike_fit.csv",
        "input/paper/csv/ERR2017816_sum_maxlike_fit.csv",
        "input/paper/csv/ERR3085830_logit_maxlike_fit.csv",
        "input/paper/csv/ERR3085830_sum_maxlike_fit.csv",
        "input/paper/csv/ERR1955542_logit_maxlike_fit.csv",
        "input/paper/csv/ERR1955542_sum_maxlike_fit.csv"
    output:
        "storage/paper/additional_figures/figure_gc_logit_vs_log.pdf"
    params:
        loaddir="input/paper/csv/"
    shell:
        """
        Rscript bin/figure_gc_logit_vs_log.R {params.loaddir} {output}
        """

rule figure_error_rate:
    input:
        expand("storage/summary/preqc/{sample}.pdf",sample=["SRR490124", "SRR3191692", "SRR3191692_assembly", "S5L001", "S1L001", "S9L001","ERR2017816","ERR3085830","ERR1955542"])
    output:
        "storage/paper/additional_figures/figure_error_rate.pdf"
    params:
        workdir="work/paper/figures/figure_error_rate"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdftk $fi cat 5 output {params.workdir}/figure${{i}}-full.pdf
            pdfcrop {params.workdir}/figure${{i}}-full.pdf {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_error_rate.pdf {params.workdir}/figure[1-9].pdf
        pdfcrop {params.workdir}/figure_error_rate.pdf {output} 1>/dev/null
        """
        
rule figure_quality_values:
    input:
        expand("storage/summary/preqc/{sample}.pdf",sample=["SRR490124", "SRR3191692", "SRR3191692_assembly", "S5L001", "S1L001", "S9L001","ERR2017816","ERR3085830","ERR1955542"])
    output:
        "storage/paper/additional_figures/figure_quality_values.pdf"
    params:
        workdir="work/paper/figures/figure_quality_values"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdftk $fi cat 10 output {params.workdir}/figure${{i}}-full.pdf
            pdfcrop {params.workdir}/figure${{i}}-full.pdf {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        Rscript bin/figure_comp9.R {params.workdir}/figure_quality_values.pdf {params.workdir}/figure[1-9].pdf
        pdfcrop {params.workdir}/figure_quality_values.pdf {output} 1>/dev/null
        """

###-Create summary plots---------------------------------------------------------------------------------------------------------------------------------------------------------------
rule comp_res_summary:
    input:
        "input/paper/csv/comp_resources.csv"
    output:
        "storage/summary/comp_resources/cpu_training.pdf",
        "storage/summary/comp_resources/cpu_simulation.pdf",
        "storage/summary/comp_resources/elapsed_training.pdf",
        "storage/summary/comp_resources/elapsed_simulation.pdf",
        "storage/summary/comp_resources/memory_training.pdf",
        "storage/summary/comp_resources/memory_simulation.pdf"
    params:
        outdir="storage/summary/comp_resources/"
    shell:
        "Rscript bin/figure_comp_resources.R {input} {params.outdir}"
        
rule cross_kmer_plot_summary_seqbias:
    input:
       "work/summary/kmer/{sample}/real.histo",
       "work/summary/kmer_cross_seqbias/{species}_{sample}_{crossspecies}_{crosssample}/ReSeq.histo",
       lambda wildcards: expand("work/summary/kmer_cross/{0}_{1}_{2}_{3}/{{simulator}}.histo".format(wildcards.species, wildcards.sample, wildcards.crossspecies, wildcards.crosssample), simulator=[sim for sim in SIMULATORS if sim != "ReSeq"] )
    output:
        "storage/summary/kmer_cross/{species}_{sample}_{crossspecies}_{crosssample}_seqbias.pdf"
    params:
        maxcount=lambda wildcards: SAMPLES[wildcards.sample]['max_kmer_count_plotted'],
        k=jellyfish_k
    shell:
        "./bin/plotKmerSpec.py -k {params.k} -x {params.maxcount} -o {output} {input}"

rule cross_kmer_plot_summary_seqbias_prep:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}_{3}_{4}/eval/kmers_seqbias/{2}-k{5}.histo".format(wildcards.species,wildcards.sample,wildcards.simulator,wildcards.crossspecies,wildcards.crosssample,jellyfish_k)
    output:
        "work/summary/kmer_cross_seqbias/{species}_{sample}_{crossspecies}_{crosssample}/{simulator}.histo",
    shell:
        "ln -srf {input} {output}"

rule cross_kmer_plot_summary:
    input:
       "work/summary/kmer/{sample}/real.histo",
       lambda wildcards: expand("work/summary/kmer_cross/{0}_{1}_{2}_{3}/{{simulator}}.histo".format(wildcards.species, wildcards.sample, wildcards.crossspecies, wildcards.crosssample), simulator=SIMULATORS )
    output:
        "storage/summary/kmer_cross/{species}_{sample}_{crossspecies}_{crosssample}.pdf"
    params:
        maxcount=lambda wildcards: SAMPLES[wildcards.sample]['max_kmer_count_plotted'],
        k=jellyfish_k
    shell:
        "./bin/plotKmerSpec.py -k {params.k} -x {params.maxcount} -o {output} {input}"

rule cross_kmer_plot_summary_prep:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}_{3}_{4}/eval/kmers/{2}-k{5}.histo".format(wildcards.species,wildcards.sample,wildcards.simulator,wildcards.crossspecies,wildcards.crosssample,jellyfish_k)
    output:
        "work/summary/kmer_cross/{species}_{sample}_{crossspecies}_{crosssample}/{simulator}.histo"
    shell:
        "ln -srf {input} {output}"

rule kmer_plot_summary:
    input:
       "work/summary/kmer/{sample}/real.histo",
       lambda wildcards: expand("work/summary/kmer/{0}{1}{2}/{{simulator}}.histo".format(wildcards.sample, wildcards.asm, wildcards.cov), simulator=SIMULATORS)
    output:
        "storage/summary/kmer/{sample,[0-9A-Za-z]*}{asm,|_assembly}{cov,|_cov}.pdf"
    params:
        maxcount=lambda wildcards: SAMPLES[wildcards.sample]['max_kmer_count_plotted'],
        k=jellyfish_k
    shell:
        "./bin/plotKmerSpec.py -k {params.k} -x {params.maxcount} -o {output} {input}"

rule kmer_plot_summary_prep:
    input:
        lambda wildcards: "storage/{0}/{1}/real/kmers/{1}-k{2}.histo".format(SAMPLES[wildcards.sample]['species'],wildcards.sample,jellyfish_k)
    output:
        "work/summary/kmer/{sample}/real.histo"
    shell:
        "ln -srf {input} {output}"

rule kmer_plot_summary_prep2:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/eval/kmers{3}/{2}-k{4}.histo".format(SAMPLES[wildcards.sample]['species'],wildcards.sample+wildcards.asm,wildcards.simulator,wildcards.cov,jellyfish_k)
    output:
        "work/summary/kmer/{sample,[0-9A-Za-z]*}{asm,|_assembly}{cov,|_cov}/{simulator}.histo"
    shell:
        "ln -srf {input} {output}"

rule sga_preqc_plot_summary:
    input:
        "work/summary/preqc/{sample}/real.preqc",
        lambda wildcards: expand("storage/{0}/{1}/{{simulator}}/eval/preqc{2}/{{simulator}}.preqc".format(SAMPLES[wildcards.sample]['species'], wildcards.sample+wildcards.asm, wildcards.cov), simulator=SIMULATORS)
    output:
        "storage/summary/preqc/{sample,[0-9A-Za-z]*}{asm,|_assembly}{cov,|_cov}.pdf"
    log:
        "logs/sga_preqc_plot_summary/{sample}{asm}{cov}.log"
    shell:
        "./bin/sga-preqc-report.py -p -o $(dirname {output})/$(basename {output} .pdf) {input} 1>{log} 2>&1"

rule sga_preqc_plot_summary_prep:
    input:
        lambda wildcards: "storage/{0}/{1}/real/preqc/{1}.preqc".format(SAMPLES[wildcards.sample]['species'], wildcards.sample)
    output:
        "work/summary/preqc/{sample}/real.preqc"
    shell:
        "ln -srf {input} {output}"

rule plot_eval_with_reseq_summary:
    input:
        "work/summary/reseq-plots/{sample}/real.reseq",
        expand("work/summary/reseq-plots/{{sample}}/{simulator}.reseq", simulator=SIMULATORS)
    output:
        "storage/summary/reseq-plots/{sample}.pdf"
    log:
        "logs/summary/{sample}.log"
    shell:
        """
        plotDataStats.py -o {output} {input} >{log} 2>&1
        """
        
rule plot_eval_with_reseq_summary_prep:
    input:
        lambda wildcards: "storage/{0}/{1}/ReSeq/original.reseq".format(SAMPLES[wildcards.sample]['species'],wildcards.sample)
    output:
        "work/summary/reseq-plots/{sample}/real.reseq"
    shell:
        "ln -srf {input} {output}"

rule plot_eval_with_reseq_summary_prep2:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/eval/mapping-bowtie2-s.bam.reseq".format(SAMPLES[wildcards.sample]['species'],wildcards.sample,wildcards.simulator)
    output:
        "work/summary/reseq-plots/{sample}/{simulator}.reseq"
    shell:
        "ln -srf {input} {output}"

###-Correlation analysis---------------------------------------------------------------------------------------------------------------------------------------------------------------
rule coverage_correlation_sim_real_plot:
    input:
        real="work/{species}/{sample}/real/correlation/{sample}_q0.depth",
        sim="work/{species}/{sample}/{simulation}/replicates/sim1_q0.depth"
    output:
        "storage/{species}/{sample}/{simulation}/correlation/pearson-coverage-real.pdf"
    params:
        name1="Coverage original",
        name2="Coverage {simulation}",
        legendpos=lambda wildcards: "lr" if wildcards.simulation == "ART" or wildcards.simulation == "NEAT" else "ul"
    shell:
        "Rscript bin/figure_covcorrelation.R {input.real} {input.sim} '{params.name1}' '{params.name2}' {output} {params.legendpos}"

rule coverage_correlation_sim_real_spearman:
    input:
        real="work/{species}/{sample}/real/correlation/{sample}_q0.depth",
        sim=lambda wildcards: expand("work/{0}/{1}/{2}/replicates/sim{{simnum}}_q0.depth".format(wildcards.species, wildcards.sample, wildcards.simulation), simnum=["1","2","3"])
    output:
        "storage/{species}/{sample}/{simulation}/correlation/spearman-coverage-real.txt"
    shell:
        """
        Rscript bin/coverage_spearman.R {input.real} {input.sim} {output}
        """

rule coverage_correlation_sim_real:
    input:
        real="work/{species}/{sample}/real/correlation/{sample}_q0.depth",
        sim=lambda wildcards: expand("work/{0}/{1}/{2}/replicates/sim{{simnum}}_q0.depth".format(wildcards.species, wildcards.sample, wildcards.simulation), simnum=["1","2","3"])
    output:
        "storage/{species}/{sample}/{simulation}/correlation/pearson-coverage-real.txt"
    params:
        tmp="work/{species}/{sample}/{simulation}/correlation/tmp_depth_diff_coverage_real.txt"
    shell:
        """
        mkdir -p $(dirname "{params.tmp}")
        rm -f {output}
        for S1 in {input.sim}; do
            join -j 4 -e "-" -o 1.1 1.2 2.1 2.2 1.3 2.3 -a 1 -a 2 <(cat $S1 | awk '{{print $1, $2, $3, $1 "_" $2}}' | sort -k4) <(cat {input.real} | awk '{{print $1, $2, $3, $1 "_" $2}}' | sort -k4) | awk '($1 == $3){{print $1, $2, $5, $6}}($1 == "-"){{print $3, $4, 0, $6}}($3 == "-"){{print $1, $2, $5, 0}}' | sort -k1,1 -k2,2n > {params.tmp}
            cat {params.tmp} | awk 'BEGIN{{sum1=0;sum2=0;count=0}}{{sum1 += $3; sum2 += $4; count++}}END{{print sum1/count, sum2/count}}' | awk -v f1="$S1" -v f2="{input.real}" 'BEGIN{{sum1=0; sum2=0; sum_m=0}}{{if(1==NR){{mean1=$1; mean2=$2}}else{{sum1 += ($3-mean1)^2; sum2 += ($4-mean2)^2; sum_m += ($3-mean1)*($4-mean2)}}}}END{{print f1, f2, sum_m/sqrt(sum1)/sqrt(sum2)}}' - {params.tmp} >> {output}
        done
        rm -f {params.tmp}
        """
        
rule coverage_correlation_sim_plot:
    input:
        lambda wildcards: expand("work/{0}/{1}/{2}/replicates/sim{{simnum}}_q0.depth".format(wildcards.species, wildcards.sample, wildcards.simulation), simnum=["1","2"])
    output:
        "storage/{species}/{sample}/{simulation}/correlation/pearson-coverage-simreps.pdf"
    params:
        name1="Coverage {simulation} sim 1",
        name2="Coverage {simulation} sim 2"
    shell:
        "Rscript bin/figure_covcorrelation.R {input} '{params.name1}' '{params.name2}' {output} ul"

rule coverage_correlation_sim_spearman:
    input:
        lambda wildcards: expand("work/{0}/{1}/{2}/replicates/sim{{simnum}}_q0.depth".format(wildcards.species, wildcards.sample, wildcards.simulation), simnum=["1","2","3"])
    output:
        "storage/{species}/{sample}/{simulation}/correlation/spearman-coverage-simreps.txt"
    shell:
        """
        Rscript bin/coverage_spearman.R {input} {output}
        """

rule coverage_correlation_sim:
    input:
        lambda wildcards: expand("work/{0}/{1}/{2}/replicates/sim{{simnum}}_q0.depth".format(wildcards.species, wildcards.sample, wildcards.simulation), simnum=["1","2","3"])
    output:
        "storage/{species}/{sample}/{simulation}/correlation/pearson-coverage-simreps.txt"
    params:
        tmp="work/{species}/{sample}/{simulation}/correlation/tmp_depth_diff_coverage_simreps.txt"
    shell:
        """
        mkdir -p $(dirname "{params.tmp}")
        rm -f {output}
        for S1 in {input}; do
            for S2 in {input}; do
                if [[ "$S1" < "$S2" ]]; then
                    join -j 4 -e "-" -o 1.1 1.2 2.1 2.2 1.3 2.3 -a 1 -a 2 <(cat $S1 | awk '{{print $1, $2, $3, $1 "_" $2}}' | sort -k4) <(cat $S2 | awk '{{print $1, $2, $3, $1 "_" $2}}' | sort -k4) | awk '($1 == $3){{print $1, $2, $5, $6}}($1 == "-"){{print $3, $4, 0, $6}}($3 == "-"){{print $1, $2, $5, 0}}' | sort -k1,1 -k2,2n > {params.tmp}
                    cat {params.tmp} | awk 'BEGIN{{sum1=0;sum2=0;count=0}}{{sum1 += $3; sum2 += $4; count++}}END{{print sum1/count, sum2/count}}' | awk -v f1="$S1" -v f2="$S2" 'BEGIN{{sum1=0; sum2=0; sum_m=0}}{{if(1==NR){{mean1=$1; mean2=$2}}else{{sum1 += ($3-mean1)^2; sum2 += ($4-mean2)^2; sum_m += ($3-mean1)*($4-mean2)}}}}END{{print f1, f2, sum_m/sqrt(sum1)/sqrt(sum2)}}' - {params.tmp} >> {output}
                fi
            done
        done
        rm -f {params.tmp}
        """

rule coverage_correlation_plot:
    input:
        lambda wildcards: "work/{0}/{1}/real/correlation/{1}_q0.depth".format(wildcards.species, wildcards.sample),
        lambda wildcards: "work/{0}/{1}/real/correlation/{2}_q0.depth".format(wildcards.species, wildcards.sample, LANEREPS[wildcards.sample][0] if 'lane' == wildcards.type else LIBREPS[wildcards.sample][0])
    output:
        "storage/{species}/{sample}/real/correlation/pearson-coverage-{type,lane|lib}reps.pdf"
    params:
        name1="Coverage original",
        name2=lambda wildcards: "Coverage {0} replicate".format("lane" if "lane" == wildcards.type else "library")
    shell:
        "Rscript bin/figure_covcorrelation.R {input} '{params.name1}' '{params.name2}' {output} ul"

rule coverage_correlation_spearman:
    input:
        lambda wildcards: "work/{0}/{1}/real/correlation/{1}_q0.depth".format(wildcards.species, wildcards.sample),
        lambda wildcards: expand("work/{0}/{1}/real/correlation/{{sample2}}_q0.depth".format(wildcards.species, wildcards.sample), sample2=LANEREPS[wildcards.sample] if 'lane' == wildcards.type else LIBREPS[wildcards.sample] )
    output:
        "storage/{species}/{sample}/real/correlation/spearman-coverage-{type,lane|lib}reps.txt"
    shell:
        """
        Rscript bin/coverage_spearman.R {input} {output}
        """

rule coverage_correlation:
    input:
        lambda wildcards: "work/{0}/{1}/real/correlation/{1}_q0.depth".format(wildcards.species, wildcards.sample),
        lambda wildcards: expand("work/{0}/{1}/real/correlation/{{sample2}}_q0.depth".format(wildcards.species, wildcards.sample), sample2=LANEREPS[wildcards.sample] if 'lane' == wildcards.type else LIBREPS[wildcards.sample] )
    output:
        "storage/{species}/{sample}/real/correlation/pearson-coverage-{type,lane|lib}reps.txt"
    params:
        tmp="work/{species}/{sample}/real/correlation/tmp_depth_diff_coverage_{type}.txt"
    shell:
        """
        rm -f {output}
        for S1 in {input}; do
            for S2 in {input}; do
                if [[ "$S1" < "$S2" ]]; then
                    join -j 4 -e "-" -o 1.1 1.2 2.1 2.2 1.3 2.3 -a 1 -a 2 <(cat $S1 | awk '{{print $1, $2, $3, $1 "_" $2}}' | sort -k4) <(cat $S2 | awk '{{print $1, $2, $3, $1 "_" $2}}' | sort -k4) | awk '($1 == $3){{print $1, $2, $5, $6}}($1 == "-"){{print $3, $4, 0, $6}}($3 == "-"){{print $1, $2, $5, 0}}' | sort -k1,1 -k2,2n > {params.tmp}
                    cat {params.tmp} | awk 'BEGIN{{sum1=0;sum2=0;count=0}}{{sum1 += $3; sum2 += $4; count++}}END{{print sum1/count, sum2/count}}' | awk -v f1="$S1" -v f2="$S2" 'BEGIN{{sum1=0; sum2=0; sum_m=0}}{{if(1==NR){{mean1=$1; mean2=$2}}else{{sum1 += ($3-mean1)^2; sum2 += ($4-mean2)^2; sum_m += ($3-mean1)*($4-mean2)}}}}END{{print f1, f2, sum_m/sqrt(sum1)/sqrt(sum2)}}' - {params.tmp} >> {output}
                fi
            done
        done
        rm -f {params.tmp}
        """

rule samtools_depth_sim:
    input:
        lambda wildcards: "storage/{0}/{1}/{2}/eval/mapping-bowtie2-s.bam".format(wildcards.species, wildcards.sample, wildcards.simulator) if "1" == wildcards.simnum else "work/{0}/{1}/{2}/replicates/sim{3}.bam".format(wildcards.species, wildcards.sample, wildcards.simulator,wildcards.simnum)
    output:
        "work/{species}/{sample}/{simulator}/replicates/sim{simnum}_q{qual}.depth"
    shell:
        "samtools depth -a -d 0 -Q {wildcards.qual} {input} > {output}"

rule samtools_depth:
    input:
        lambda wildcards: "storage/{0}/{1}/real/mapping-bowtie2-s.bam".format(wildcards.species, wildcards.sample) if wildcards.sample == wildcards.sample2 else "storage/{0}/{1}/real/replicates/{2}.bam".format(wildcards.species, wildcards.sample, wildcards.sample2)
    output:
        "work/{species}/{sample}/real/correlation/{sample2}_q{qual}.depth"
    shell:
        "samtools depth -a -d 0 -Q {wildcards.qual} {input} > {output}"

rule bowtie2_map_sim_replicates:        
    input:
        "work/{species}/reference/{sample}/index_bowtie2/pilon-corrected.1.bt2",
        fq1=lambda wildcards: "storage/{0}/{1}/{2}/{3}{4}".format(wildcards.species,wildcards.sample,wildcards.simulator,("" if "1" == wildcards.simnum else "replicates/sim{0}/".format(wildcards.simnum)) if "" == wildcards.type else "syserror-replicates/sim{0}/".format(wildcards.simnum),simulatorOutputFile(1,wildcards)),
        fq2=lambda wildcards: "storage/{0}/{1}/{2}/{3}{4}".format(wildcards.species,wildcards.sample,wildcards.simulator,("" if "1" == wildcards.simnum else "replicates/sim{0}/".format(wildcards.simnum)) if "" == wildcards.type else "syserror-replicates/sim{0}/".format(wildcards.simnum),simulatorOutputFile(2,wildcards))
    output:
        "work/{species}/{sample}/{simulator}/{type,|syserror-}replicates/sim{simnum}.bam"
    params:
        refind="work/{species}/reference/{sample}/index_bowtie2/pilon-corrected",
        tmpdir="work/{species}/{sample}/{simulator}/{type}replicates/sim{simnum}"
    log:
        "logs/bowtie2_map_sim_{type}replicates/{sample}/{simulator}-{simnum}.log"
    threads: max_threads
    shell:
        "mkdir -p {params.tmpdir}; "
        "bowtie2 -p {threads} -X 2000 -x {params.refind} -1 {input.fq1} -2 {input.fq2}"
        " 2>{log} | samtools sort -m 10G -@ 4 -T {params.tmpdir}/_tmp -o {output} -"

rule bowtie2_map_replicates:
    input:
        "work/{species}/reference/{sample}/index_bowtie2/pilon-corrected.1.bt2",
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample2]['basespace'],SAMPLES[wildcards.sample2]['fq1']) if 'basespace' in SAMPLES[wildcards.sample2] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample2,SAMPLES[wildcards.sample2]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample2]['basespace'],SAMPLES[wildcards.sample2]['fq2']) if 'basespace' in SAMPLES[wildcards.sample2] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample2,SAMPLES[wildcards.sample2]['fq2'])
    output:
        "storage/{species}/{sample}/real/replicates/{sample2}.bam"
    params:
        refind="work/{species}/reference/{sample}/index_bowtie2/pilon-corrected",
        tmp="work/{species}/{sample}/real/tmp_{sample2}_"
    log:
        "logs/bowtie2_map_replicates/{sample}-{sample2}.log"
    threads: max_threads
    shell:
        "bowtie2 -p {threads} -X 2000 -x {params.refind} -1 {input.fq1} -2 {input.fq2}"
        " 2>{log} | samtools sort -m 10G -@ 4 -T {params.tmp} -o {output} -"

###-Evaluate simulators----------------------------------------------------------------------------------------------------------------------------------------------------------------
rule bwa_map:
    input:
        lambda wildcards: "work/{0}/reference/{1}/index_bwa/{2}.bwt".format(wildcards.species,wildcards.sample, "pilon-corrected" if "" == wildcards.asm else "sga-contigs"),
        fq1=lambda wildcards: "storage/{0}/{1}/{2}/{4}{3}".format(wildcards.species,wildcards.sample+wildcards.asm,wildcards.simulator,simulatorOutputFile(1,wildcards),"nobias/" if wildcards.folder == "nobias" else ""),
        fq2=lambda wildcards: "storage/{0}/{1}/{2}/{4}{3}".format(wildcards.species,wildcards.sample+wildcards.asm,wildcards.simulator,simulatorOutputFile(2,wildcards),"nobias/" if wildcards.folder == "nobias" else "")
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/{simulator}/{folder}/mapping-bwa-s.bam"
    params:
        refind=lambda wildcards: "work/{0}/reference/{1}/index_bwa/{2}".format(wildcards.species,wildcards.sample, "pilon-corrected" if "" == wildcards.asm else "sga-contigs"),
        tmpdir="work/{species}/{sample}{asm}/{simulator}/eval"
    log:
        "logs/bwa_map/{species}/{sample}{asm}/{folder}_{simulator}.log"
    threads: max_threads
    shell:
        "mkdir -p {params.tmpdir}; "
        "bwa mem -t {threads} {params.refind} {input.fq1} {input.fq2}"
        " 2>{log} | samtools sort -m 10G -@ 4 -T {params.tmpdir}/tmp_bwa_ -o {output} -"

rule plot_eval_with_reseq_nobias:
    input:
        "storage/{species}/{sample}/ReSeq/original.reseq",
        "storage/{species}/{sample}/{simulator}/nobias/nobias.reseq"
    output:
        "storage/{species}/{sample}/{simulator}/nobias/stats.pdf"
    log:
        "logs/eval_with_reseq/{species}/{sample}/plot_nobias_{simulator}.log"
    shell:
        "plotDataStats.py -o {output} {input} >{log} 2>&1"
        
rule plot_eval_with_reseq:
    input:
        "storage/{species}/{sample}/ReSeq/original.reseq",
        "storage/{species}/{sample}/{simulator}/eval/mapping-bowtie2-s.bam.reseq"
    output:
        "storage/{species}/{sample}/{simulator}/eval/stats.pdf"
    log:
        "logs/eval_with_reseq/{species}/{sample}/plot_{simulator}.log"
    shell:
        "plotDataStats.py -o {output} {input} >{log} 2>&1"

rule eval_with_reseq_nobias:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        bam="storage/{species}/{sample}/{simulator}/nobias/mapping-bowtie2-s.bam"
    output:
        "storage/{species}/{sample}/{simulator}/nobias/nobias.reseq"
    log:
        "logs/eval_with_reseq/{species}/{sample}/nobias_{simulator}.log"
    threads: max_threads
    shell:
        "reseq illuminaPE -j {threads} -r {input.ref} -b {input.bam} -S {output} --statsOnly >{log} 2>&1"

rule eval_with_reseq:
    input:
        ref=lambda wildcards: "storage/{0}/reference/{1}/{2}".format(wildcards.species, wildcards.sample, "sga-contigs.fa" if wildcards.asm != "" else "pilon-corrected.fa"),
        bam="storage/{species}/{sample}{asm}/{simulator}/eval/mapping-bowtie2-s.bam"
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/{simulator}/eval/mapping-bowtie2-s.bam.reseq"
    params:
        adapter=lambda wildcards: "" if wildcards.simulator == "ReSeq" else "--adapterFile TruSeq_single"
    log:
        "logs/eval_with_reseq/{species}/{sample}{asm}/{simulator}.log"
    threads: max_threads
    shell:
        "reseq illuminaPE -j {threads} -r {input.ref} -b {input.bam} -S {output} {params.adapter} --statsOnly >{log} 2>&1"

rule bowtie2_map:
    input:
        lambda wildcards: "work/{0}/reference/{1}/index_bowtie2/{2}.1.bt2".format(wildcards.species,wildcards.sample, "pilon-corrected" if "" == wildcards.asm else "sga-contigs"),
        fq1=lambda wildcards: "storage/{0}/{1}/{2}/{4}{3}".format(wildcards.species,wildcards.sample+wildcards.asm,wildcards.simulator,simulatorOutputFile(1,wildcards),"nobias/" if wildcards.folder == "nobias" else ""),
        fq2=lambda wildcards: "storage/{0}/{1}/{2}/{4}{3}".format(wildcards.species,wildcards.sample+wildcards.asm,wildcards.simulator,simulatorOutputFile(2,wildcards),"nobias/" if wildcards.folder == "nobias" else "")
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/{simulator}/{folder}/mapping-bowtie2-s.bam"
    params:
        refind=lambda wildcards: "work/{0}/reference/{1}/index_bowtie2/{2}".format(wildcards.species,wildcards.sample, "pilon-corrected" if "" == wildcards.asm else "sga-contigs"),
        tmpdir="work/{species}/{sample}{asm}/{simulator}/eval"
    log:
        "logs/bowtie2_map/{species}/{sample}{asm}/{folder}_{simulator}.log"
    threads: max_threads
    shell:
        "mkdir -p {params.tmpdir}; "
        "bowtie2 -p {threads} -X 2000 -x {params.refind} -1 {input.fq1} -2 {input.fq2} 2>{log}"
        " | samtools sort -m 10G -@ 4 -T {params.tmpdir}/tmp_ -o {output} -"

rule sga_preqc_plot_sim:
    input:
        "storage/{species}/{sample}/{simulator}/eval/preqc/{simulator}.preqc"
    output:
        "storage/{species}/{sample}/{simulator}/eval/preqc/{simulator}.pdf"
    params:
        output="storage/{species}/{sample}/{simulator}/eval/preqc/{simulator}"
    log:
        "logs/sga_preqc_plot/{simulator}-{sample}.log"
    shell:
        "./bin/sga-preqc-report.py -p -o {params.output} {input} 1>{log} 2>&1"

rule sga_preqc_sim:
    input:
        fq1=lambda wildcards: "storage/{0}/{1}/{2}/{3}".format(wildcards.species,wildcards.sample,wildcards.simulator,simulatorOutputFile(1,wildcards) if "" == wildcards.cov else simulatorOutputFile(1,wildcards).replace(wildcards.simulator.lower(), wildcards.simulator.lower()+wildcards.cov)),
        fq2=lambda wildcards: "storage/{0}/{1}/{2}/{3}".format(wildcards.species,wildcards.sample,wildcards.simulator,simulatorOutputFile(2,wildcards) if "" == wildcards.cov else simulatorOutputFile(2,wildcards).replace(wildcards.simulator.lower(), wildcards.simulator.lower()+wildcards.cov))
    output:
        "storage/{species}/{sample}/{simulator}/eval/preqc{cov,|_cov}/{simulator}.preqc"
    params:
        output="storage/{species}/{sample}/{simulator}/eval/preqc{cov}/{simulator}"
    log:
        "logs/sga_preqc/{simulator}-{sample}{cov}.log"
    threads: max_threads
    shell:
        "bin/preqc.sh -j {threads} -l {log} -o {params.output} {input.fq1} {input.fq2}"

rule plot_kmer:
    input:
        real="storage/{species}/{sample}/real/kmers/{sample}-k{kmer_length}.histo",
        sim="storage/{species}/{sample}{asm}/{simulator}/eval/kmers/{simulator}-k{kmer_length}.histo"
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/{simulator}/eval/kmers/{simulator}-k{kmer_length}.pdf"
    params:
        maxcount=lambda wildcards: SAMPLES[wildcards.sample]['max_kmer_count_plotted']
    shell:
        "./bin/plotKmerSpec.py -k {wildcards.kmer_length} -x {params.maxcount} -o {output} {input.real} {input.sim}"

rule kmer_cross_sim:
    input:
        fq1=lambda wildcards: "storage/{0}/{1}/{2}_{3}_{4}/{5}".format(wildcards.species,wildcards.sample,wildcards.simulator,wildcards.crossspecies,wildcards.crosssample,simulatorOutputFile(1,wildcards) if "" == wildcards.tag else simulatorOutputFile(1,wildcards).replace(wildcards.simulator.lower(), wildcards.simulator.lower()+wildcards.tag)),
        fq2=lambda wildcards: "storage/{0}/{1}/{2}_{3}_{4}/{5}".format(wildcards.species,wildcards.sample,wildcards.simulator,wildcards.crossspecies,wildcards.crosssample,simulatorOutputFile(2,wildcards) if "" == wildcards.tag else simulatorOutputFile(2,wildcards).replace(wildcards.simulator.lower(), wildcards.simulator.lower()+wildcards.tag))
    output:
        "storage/{species}/{sample}/{simulator}_{crossspecies}_{crosssample}/eval/kmers{tag,|_seqbias}/{simulator}-k{kmer_length}.histo"
    params:
        output="storage/{species}/{sample}/{simulator}_{crossspecies}_{crosssample}/eval/kmers{tag}/{simulator}-k{kmer_length}",
        tmp_dir="work/{species}/{sample}/{simulator}_{crossspecies}_{crosssample}/eval/kmers{tag}/",
        jellyfish=lambda wildcards: REFERENCES[wildcards.species]['use_jellyfish']
    log:
        "logs/kmc/{species}/{sample}/{simulator}_{crossspecies}_{crosssample}{tag}.log"
    threads: max_threads
    resources:
        mem_gb=256
    shell:
        """
        if [ 1 -eq {params.jellyfish} ]
        then
            jellyfish count -m {wildcards.kmer_length} -s 12G -t {threads} -o {params.output}.jf -C -F 2 <(zcat {input.fq1}) <(zcat {input.fq2})
            jellyfish histo -o {output} {params.output}.jf
            rm -f {params.output}.jf
        else
            mkdir -p {params.tmp_dir}
            printf "{input.fq1}\n{input.fq2}\n" > {params.tmp_dir}/input_list.txt
            kmc -k{wildcards.kmer_length} -m{resources.mem_gb} -ci1 -cs65000 -t{threads} @{params.tmp_dir}/input_list.txt {params.output} {params.tmp_dir} 1>{log} 2>&1
            kmc_tools transform {params.output} histogram {output} 1>>{log} 2>&1
            rm -f {params.output}.kmc_pre {params.output}.kmc_suf
        fi
        """

rule kmer_sim:
    input:
        fq1=lambda wildcards: "storage/{0}/{1}/{2}/{3}".format(wildcards.species,wildcards.sample+wildcards.asm,wildcards.simulator,simulatorOutputFile(1,wildcards) if "" == wildcards.cov else simulatorOutputFile(1,wildcards).replace(wildcards.simulator.lower(), wildcards.simulator.lower()+wildcards.cov)),
        fq2=lambda wildcards: "storage/{0}/{1}/{2}/{3}".format(wildcards.species,wildcards.sample+wildcards.asm,wildcards.simulator,simulatorOutputFile(2,wildcards) if "" == wildcards.cov else simulatorOutputFile(2,wildcards).replace(wildcards.simulator.lower(), wildcards.simulator.lower()+wildcards.cov))
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/{simulator}/eval/kmers{cov,|_cov}/{simulator}-k{kmer_length}.histo"
    params:
        output="storage/{species}/{sample}{asm}/{simulator}/eval/kmers{cov}/{simulator}-k{kmer_length}",
        tmp_dir="work/{species}/{sample}{asm}/{simulator}/eval/kmers{cov}/",
        jellyfish=lambda wildcards: REFERENCES[wildcards.species]['use_jellyfish']
    log:
        "logs/kmc/{species}/{sample}{asm}/{simulator}{cov}.log"
    threads: max_threads
    resources:
        mem_gb=256
    shell:
        """
        if [ 1 -eq {params.jellyfish} ]
        then
            jellyfish count -m {wildcards.kmer_length} -s 12G -t {threads} -o {params.output}.jf -C -F 2 <(zcat {input.fq1}) <(zcat {input.fq2})
            jellyfish histo -o {output} {params.output}.jf
            rm -f {params.output}.jf
        else
            mkdir -p {params.tmp_dir}
            printf "{input.fq1}\n{input.fq2}\n" > {params.tmp_dir}/input_list.txt
            kmc -k{wildcards.kmer_length} -m{resources.mem_gb} -ci1 -cs65000 -t{threads} @{params.tmp_dir}/input_list.txt {params.output} {params.tmp_dir} 1>{log} 2>&1
            kmc_tools transform {params.output} histogram {output} 1>>{log} 2>&1
            rm -f {params.output}.kmc_pre {params.output}.kmc_suf
        fi
        """

###-Information on real data to compare simulators to----------------------------------------------------------------------------------------------------------------------------------

rule sga_preqc_plot:
    input:
        "storage/{species}/{sample}/real/preqc/{sample}.preqc"
    output:
        "storage/{species}/{sample}/real/preqc/{sample}.pdf"
    params:
        output="storage/{species}/{sample}/real/preqc/{sample}"
    log:
        "logs/sga_preqc_plot/real-{sample}.log"
    shell:
        "./bin/sga-preqc-report.py -p -o {params.output} {input} 1>{log} 2>&1"

rule sga_preqc:
    input:
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        protected("storage/{species}/{sample}/real/preqc/{sample}.preqc")
    params:
        output="storage/{species}/{sample}/real/preqc/{sample}"
    log:
        "logs/sga_preqc/real-{sample}.log"
    threads: max_threads
    shell:
        "bin/preqc.sh -j {threads} -l {log} -o {params.output} {input.fq1} {input.fq2}"

rule kmer:
    input:
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        "storage/{species}/{sample}/real/kmers/{sample}-k{kmer_length}.histo"
    params:
        output="storage/{species}/{sample}/real/kmers/{sample}-k{kmer_length}",
        tmp_dir="work/{species}/{sample}/real/kmers/",
        jellyfish=lambda wildcards: REFERENCES[wildcards.species]['use_jellyfish']
    log:
        "logs/kmc/{species}/{sample}/real.log"
    threads: max_threads
    resources:
        mem_gb=256
    shell:
        """
        if [ 1 -eq {params.jellyfish} ]
        then
            jellyfish count -m {wildcards.kmer_length} -s 12G -t {threads} -o {params.output}.jf -C -F 2 <(zcat {input.fq1}) <(zcat {input.fq2})
            jellyfish histo -o {output} {params.output}.jf
            rm -f {params.output}.jf
        else
            mkdir -p {params.tmp_dir}
            printf "{input.fq1}\n{input.fq2}\n" > {params.tmp_dir}/input_list.txt
            kmc -k{wildcards.kmer_length} -m{resources.mem_gb} -ci1 -cs65000 -t{threads} @{params.tmp_dir}/input_list.txt {params.output} {params.tmp_dir} 1>{log} 2>&1
            kmc_tools transform {params.output} histogram {output} 1>>{log} 2>&1
            rm -f {params.output}.kmc_pre {params.output}.kmc_suf
        fi
        """

rule surrounding_pdf:
    input:
        expand("storage/{{species}}/{{sample}}/real/surrounding_{strand}_{segment}.csv", strand=["forward","reverse"], segment=["first","second"])
    output:
        "storage/{species}/{sample}/real/surrounding.pdf"
    params:
        workdir="storage/{species}/{sample}/real/"
    shell:
        "Rscript bin/figure_surrounding.R {params.workdir}"

def SurroundingCsvFlags(strand, segment):
    if "forward" == strand:
        if "first" == segment:
            return "-f 99 -F 3868"
        elif "second" == segment:
            return "-f 163 -F 3868"
        else:
            raise KeyError("Unknown segment specified: {0}".format(segment))
    elif "reverse" == strand:
        if "first" == segment:
            return "-f 83 -F 3884"
        elif "second" == segment:
            return "-f 147 -F 3884"
        else:
            raise KeyError("Unknown segment specified: {0}".format(segment))
    else:
        raise KeyError("Unknown strand specified: {0}".format(strand))

rule surrounding_csv_forward:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        bam="storage/{species}/{sample}/real/mapping-bowtie2-s.bam"
    output:
        "storage/{species}/{sample}/real/surrounding_forward_{segment}.csv"
    params:
        flags=lambda wildcards: SurroundingCsvFlags("forward", wildcards.segment)
    shell:
        "echo 'base,pos,count' > {output}; samtools faidx -r <(samtools view -q 10 {params.flags} {input.bam} | awk '($4>50){{print $3 \":\" $4-5 \"-\" $4+9}}') {input.ref} | awk -F '' '(NR%2==0){{ for(p=1; p <= NF; p++){{ stor[$p][p] += 1 }} }}END{{ for(base in stor){{ for(p in stor[base]){{ print base \",\" p-6 \",\" stor[base][p]}}}} }}' >> {output}"
        

rule surrounding_csv_reverse:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        bam="storage/{species}/{sample}/real/mapping-bowtie2-s.bam"
    output:
        "storage/{species}/{sample}/real/surrounding_reverse_{segment}.csv"
    params:
        flags=lambda wildcards: SurroundingCsvFlags("reverse", wildcards.segment)
    shell:
        "echo 'base,pos,count' > {output}; samtools faidx --reverse-complement -r <(samtools view -q 10 {params.flags} {input.bam} | awk '($8>50){{end=$8-$9-1; print $3 \":\" end-9 \"-\" end+50}}') {input.ref} 2>/dev/null | awk -F '' '((NR%2==0) && (60 == length($0))){{ for(p=46; p <= NF; p++){{ stor[$p][p-45] += 1 }} }}END{{ for(base in stor){{ for(p in stor[base]){{ print base \",\" p-6 \",\" stor[base][p]}}}} }}' >> {output}"

###-ReSeq------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

rule gzip:
    input:
        "{file}"
    output:
        "{file}.gz"
    shell:
        "gzip -f {input}"

rule reseq_simulate_nobias:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{species}/{sample}/ReSeq{mapper}/original.reseq",
        ipf="storage/{species}/{sample}/ReSeq{mapper}/original.reseq.ipf",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else []
    output:
        fq1="storage/{species}/{sample}/ReSeq{mapper,|_bwa}/nobias/reseq-R1.fq",
        fq2="storage/{species}/{sample}/ReSeq{mapper,|_bwa}/nobias/reseq-R2.fq"
    params:
        diploid=lambda wildcards: "-V " if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/reseq{mapper}/{species}/{sample}_reseq_simulate_nobias.log"
    threads: max_threads
    shell:
        "reseq illuminaPE -j {threads} -R {input.ref} --noBias -s {input.stats} -p {input.ipf} --ipfIterations 0 {params.diploid}{input.vcf} -1 {output.fq1} -2 {output.fq2} >{log} 2>&1"

rule reseq_cross_simulate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{crossspecies}/{crosssample}/ReSeq/original.reseq",
        ipf="storage/{crossspecies}/{crosssample}/ReSeq/original.reseq.ipf",
        seqbias=lambda wildcards: [] if "" == wildcards.seqbias else "work/{0}/{1}/ReSeq/original_seqbias.tsv".format(wildcards.species, wildcards.sample)
    output:
        fq1="storage/{species}/{sample}/ReSeq_{crossspecies}_{crosssample}/reseq{seqbias,|_seqbias}-R1.fq",
        fq2="storage/{species}/{sample}/ReSeq_{crossspecies}_{crosssample}/reseq{seqbias,|_seqbias}-R2.fq",
        timelog="storage/{species}/{sample}/ReSeq_{crossspecies}_{crosssample}/benchmark/reseq_simulate{seqbias,|_seqbias}.log"
    params:
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        seqbias=lambda wildcards: "no" if "" == wildcards.seqbias else "file --refBiasFile work/{0}/{1}/ReSeq/original_seqbias.tsv".format(wildcards.species, wildcards.sample)
    log:
        "logs/reseq/{species}/{sample}_reseq_{crossspecies}_{crosssample}_simulate{seqbias}.log"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} reseq illuminaPE -j {threads} -R {input.ref} --refBias {params.seqbias} -c {params.coverage} -s {input.stats} -p {input.ipf} --ipfIterations 0 -1 {output.fq1} -2 {output.fq2} >{log} 2>&1"

rule reseq_extract_seqbias:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{species}/{sample}/ReSeq/original.reseq"
    output:
        "work/{species}/{sample}/ReSeq/original_seqbias.tsv"
    shell:
        "reseq getRefSeqBias -o {output} -r {input.ref} -s {input.stats}"
    
rule reseq_simulate_syserror_replicate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{species}/{sample}/ReSeq/original.reseq",
        ipf="storage/{species}/{sample}/ReSeq/original.reseq.ipf",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        syserr="storage/{species}/{sample}/ReSeq/syserror-replicates/syserror.fq"
    output:
        fq1="storage/{species}/{sample}/ReSeq/syserror-replicates/sim{simnum,(?!1).*}/reseq-R1.fq",
        fq2="storage/{species}/{sample}/ReSeq/syserror-replicates/sim{simnum,(?!1).*}/reseq-R2.fq"
    params:
        diploid=lambda wildcards: "-V " if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/reseq/{species}/{sample}_reseq_simulate_syserror_sim{simnum}.log"
    threads: max_threads
    shell:
        "reseq illuminaPE -j {threads} -R {input.ref} --refBias keep -s {input.stats} -p {input.ipf} --ipfIterations 0 {params.diploid}{input.vcf} --readSysError {input.syserr} -1 {output.fq1} -2 {output.fq2} >{log} 2>&1"

rule reseq_simulate_syserror:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{species}/{sample}/ReSeq/original.reseq",
        ipf="storage/{species}/{sample}/ReSeq/original.reseq.ipf",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else []
    output:
        syserr="storage/{species}/{sample}/ReSeq/syserror-replicates/syserror.fq",
        fq1="storage/{species}/{sample}/ReSeq/syserror-replicates/sim1/reseq-R1.fq",
        fq2="storage/{species}/{sample}/ReSeq/syserror-replicates/sim1/reseq-R2.fq"
    params:
        diploid=lambda wildcards: "-V " if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/reseq/{species}/{sample}_reseq_simulate_syserror_sim1.log"
    threads: max_threads
    shell:
        "reseq illuminaPE -j {threads} -R {input.ref} --refBias keep -s {input.stats} -p {input.ipf} --ipfIterations 0 {params.diploid}{input.vcf} --writeSysError {output.syserr} -1 {output.fq1} -2 {output.fq2} >{log} 2>&1"

rule reseq_simulate_replicate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{species}/{sample}/ReSeq/original.reseq",
        ipf="storage/{species}/{sample}/ReSeq/original.reseq.ipf",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else []
    output:
        fq1="storage/{species}/{sample}/ReSeq/replicates/sim{simnum}/reseq-R1.fq",
        fq2="storage/{species}/{sample}/ReSeq/replicates/sim{simnum}/reseq-R2.fq"
    params:
        diploid=lambda wildcards: "-V " if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/reseq/{species}/{sample}_reseq_simulate_rep{simnum}.log"
    threads: max_threads
    shell:
        "reseq illuminaPE -j {threads} -R {input.ref} --refBias keep -s {input.stats} -p {input.ipf} --ipfIterations 0 {params.diploid}{input.vcf} -1 {output.fq1} -2 {output.fq2} >{log} 2>&1"

rule reseq_simulate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        stats="storage/{species}/{sample}{asm}/ReSeq{mapper}/original.reseq",
        ipf="storage/{species}/{sample}{asm}/ReSeq{mapper}/original.reseq.ipf",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else []
    output:
        fq1="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ReSeq{mapper,|_bwa}/reseq{cov,|_cov}-R1.fq",
        fq2="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ReSeq{mapper,|_bwa}/reseq{cov,|_cov}-R2.fq",
        timelog="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ReSeq{mapper,|_bwa}/benchmark/reseq_simulate{cov,|_cov}.log"
    params:
        diploid=lambda wildcards: "-V " if REFERENCES[wildcards.species]['diploid'] else "",
        coverage=lambda wildcards: "" if "" == wildcards.asm or "" == wildcards.cov else "-c {}".format(SAMPLES[wildcards.sample]['coverage'])
    log:
        "logs/reseq{mapper}/{species}/{sample}{asm}{cov}_reseq_simulate.log"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} reseq illuminaPE -j {threads} -R {input.ref} --refBias keep -s {input.stats} -p {input.ipf} {params.coverage} --ipfIterations 0 {params.diploid}{input.vcf} -1 {output.fq1} -2 {output.fq2} >{log} 2>&1"

rule reseq_iterative_proportional_fitting:
    input:
        stats="storage/{species}/{sample}/ReSeq{mapper}/original.reseq"
    output:
        ipf="storage/{species}/{sample}/ReSeq{mapper,|_bwa}/original.reseq.ipf",
        timelog="storage/{species}/{sample}/ReSeq{mapper,|_bwa}/benchmark/reseq_iterative_proportional_fitting.log"
    log:
        "logs/reseq{mapper}/{species}/{sample}_reseq_iterative_proportional_fitting.log"
    threads: max_threads
    shell:
        """
        env time -v -o {output.timelog} reseq illuminaPE -j {threads}  -s {input.stats} --stopAfterEstimation >{log} 2>&1
        """

rule reseq_create_profile_bwa:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        bam="storage/{species}/{sample}/real/mapping-bwa-s.bam",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else []
    output:
        stats="storage/{species}/{sample}/ReSeq_bwa/original.reseq",
        timelog="storage/{species}/{sample}/ReSeq_bwa/benchmark/reseq_create_profile.log"
    params:
        diploid=lambda wildcards: "-v " if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/reseq/{species}/{sample}_reseq_create_profile_bwa.log"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} reseq illuminaPE -j {threads} -r {input.ref} -b {input.bam} -S {output.stats} --statsOnly --noTiles {params.diploid}{input.vcf} >{log} 2>&1"

rule reseq_create_profile:
    input:
        ref=lambda wildcards: "storage/{0}/reference/{1}/{2}".format(wildcards.species, wildcards.sample, "sga-contigs.fa" if wildcards.asm != "" else "pilon-corrected.fa"),
        bam="storage/{species}/{sample}{asm}/real/mapping-bowtie2-s.bam",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/{2}".format(wildcards.species, wildcards.sample, "sga-contigs.vcf" if wildcards.asm != "" else "pilon-corrected.vcf") if REFERENCES[wildcards.species]['diploid'] else []
    output:
        stats="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ReSeq/original.reseq",
        timelog="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ReSeq/benchmark/reseq_create_profile.log"
    params:
        diploid=lambda wildcards: "-v " if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/reseq/{species}/{sample}{asm}_reseq_create_profile.log"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} reseq illuminaPE -j {threads} -r {input.ref} -b {input.bam} -S {output.stats} --statsOnly --noTiles {params.diploid}{input.vcf} >{log} 2>&1"

###-pIRS-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rule pirs_cross_simulate:
    input:
        ref=lambda wildcards: "work/{0}/reference/{1}/pilon-corrected-h1.fa".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else "storage/{0}/reference/{1}/pilon-corrected.fa".format(wildcards.species, wildcards.sample),
        ref2=lambda wildcards: "work/{0}/reference/{1}/pilon-corrected-h2.fa".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        basecall="work/{crossspecies}/{crosssample}/pIRS/profiles/baseCalling.count.matrix",
        gcbias=lambda wildcards: "work/{0}/{1}/pIRS/profiles/gc_bias_{2}.dat".format(wildcards.crossspecies, wildcards.crosssample, SAMPLES[wildcards.crosssample]['insertlength']),
        indel="work/{crossspecies}/{crosssample}/pIRS/profiles/stats.InDel.matrix"
    output:
        "storage/{species}/{sample}/pIRS_{crossspecies}_{crosssample}/pirs_1.fq",
        "storage/{species}/{sample}/pIRS_{crossspecies}_{crosssample}/pirs_2.fq",
        timelog="storage/{species}/{sample}/pIRS_{crossspecies}_{crosssample}/pirs-simulate-time.log"
    params:
        outdir="storage/{species}/{sample}/pIRS_{crossspecies}_{crosssample}",
        readlength=lambda wildcards: SAMPLES[wildcards.crosssample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.crosssample]['insertlength'],
        insertsd=lambda wildcards: SAMPLES[wildcards.crosssample]['insertsd'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        diploid=lambda wildcards: "-d" if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/pIRS/{species}/{sample}_pirs_{crossspecies}_{crosssample}_simulate.log"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} pirs simulate -t {threads} -B {input.basecall} -G {input.gcbias} -I {input.indel} -l {params.readlength} -x {params.coverage} -m {params.insertlength} -v {params.insertsd} -o {params.outdir} -s pirs {params.diploid} {input.ref} {input.ref2} >{log} 2>&1"

rule pirs_simulate_replicate:
    input:
        ref=lambda wildcards: "work/{0}/reference/{1}/pilon-corrected-h1.fa".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else "storage/{0}/reference/{1}/pilon-corrected.fa".format(wildcards.species, wildcards.sample),
        ref2=lambda wildcards: "work/{0}/reference/{1}/pilon-corrected-h2.fa".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        basecall="work/{species}/{sample}/pIRS/profiles/baseCalling.count.matrix",
        gcbias=lambda wildcards: "work/{0}/{1}/pIRS/profiles/gc_bias_{2}.dat".format(wildcards.species, wildcards.sample, SAMPLES[wildcards.sample]['insertlength']),
        indel="work/{species}/{sample}/pIRS/profiles/stats.InDel.matrix"
    output:
        "storage/{species}/{sample}/pIRS/replicates/sim{simnum}/pirs_1.fq",
        "storage/{species}/{sample}/pIRS/replicates/sim{simnum}/pirs_2.fq"
    params:
        outdir="storage/{species}/{sample}/pIRS/replicates/sim{simnum}",
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength'],
        insertsd=lambda wildcards: SAMPLES[wildcards.sample]['insertsd'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        diploid=lambda wildcards: "-d" if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/pIRS/{species}/{sample}_pirs_simulate_rep{simnum}.log"
    threads: max_threads
    shell:
        "pirs simulate -t {threads} -B {input.basecall} -G {input.gcbias} -I {input.indel} -l {params.readlength} -x {params.coverage} -m {params.insertlength} -v {params.insertsd} -o {params.outdir} -s pirs {params.diploid} {input.ref} {input.ref2} >{log} 2>&1"

rule pirs_simulate:
    input:
        ref=lambda wildcards: "work/{0}/reference/{1}/pilon-corrected-h1.fa".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else "storage/{0}/reference/{1}/pilon-corrected.fa".format(wildcards.species, wildcards.sample),
        ref2=lambda wildcards: "work/{0}/reference/{1}/pilon-corrected-h2.fa".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        basecall="work/{species}/{sample}{asm}/pIRS/profiles/baseCalling.count.matrix",
        gcbias=lambda wildcards: "work/{0}/{1}/pIRS/profiles/gc_bias_{2}.dat".format(wildcards.species, wildcards.sample+wildcards.asm, SAMPLES[wildcards.sample]['insertlength'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertlength_asm']),
        indel="work/{species}/{sample}{asm}/pIRS/profiles/stats.InDel.matrix"
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/pirs{cov,|_cov}_1.fq",
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/pirs{cov,|_cov}_2.fq",
        timelog="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/pirs-simulate{cov,|_cov}-time.log"
    params:
        outdir="storage/{species}/{sample}{asm}/pIRS",
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertlength_asm'],
        insertsd=lambda wildcards: SAMPLES[wildcards.sample]['insertsd'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertsd_asm'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'] if "" == wildcards.asm or "_cov" == wildcards.cov else SAMPLES[wildcards.sample]['coverage_asm'],
        diploid=lambda wildcards: "-d" if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/pIRS/{species}/{sample}{asm}{cov}_pirs_simulate.log"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} pirs simulate -t {threads} -B {input.basecall} -G {input.gcbias} -I {input.indel} -l {params.readlength} -x {params.coverage} -m {params.insertlength} -v {params.insertsd} -o {params.outdir} -s pirs{wildcards.cov} {params.diploid} {input.ref} {input.ref2} >{log} 2>&1"

rule pirs_profile_gc_bias:
    input:
        ref=lambda wildcards: "storage/{0}/reference/{1}/{2}".format(wildcards.species, wildcards.sample, "sga-contigs.fa" if wildcards.asm != "" else "pilon-corrected.fa"),
        bam="storage/{species}/{sample}{asm}/real/mapping-bowtie2-s.bam"
    output:
        "work/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/gc_bias_{insertlength}.dat",
        timelog_soap="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/soap-coverage-{insertlength}-time.log",
        timelog_gcbias="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/gc-coverage-bias-{insertlength}-time.log"
    params:
        outdir="work/{species}/{sample}{asm}/pIRS/profiles",
        tmp="work/{species}/{sample}{asm}/real/mapping-bowtie2-s.sam"
    log:
        soap="logs/pIRS/{species}/{sample}{asm}_soap_coverage.log",
        gcbias="logs/pIRS/{species}/{sample}{asm}_gc_bias_profile.log"
    threads: max_threads
    shell:
        """
        mkdir -p $(dirname "{params.tmp}")
        samtools view -h {input.bam} > {params.tmp}
        env time -v -o {output.timelog_soap} soap.coverage -p {threads} -sam -cvg -onlyuniq -i {params.tmp} -refsingle {input.ref} -o {params.outdir}/soap_coverage.dresult -depthsingle {params.outdir}/soap_coverage.depth >{log.soap} 2>&1
        env time -v -o {output.timelog_gcbias} gc_coverage_bias -r {input.ref} -o {params.outdir}/gc_bias -w {wildcards.insertlength} {params.outdir}/soap_coverage.depth >{log.gcbias} 2>&1
        rm -f {params.tmp}
        """

rule pirs_profile_indel:
    input:
        "storage/{species}/{sample}/real/mapping-bowtie2-s.bam"
    output:
        "work/{species}/{sample}/pIRS/profiles/stats.InDel.dat",
        "work/{species}/{sample}/pIRS/profiles/stats.InDel.matrix",
        timelog="storage/{species}/{sample}/pIRS/profiles/indelstat-time.log"
    params:
        outdir="work/{species}/{sample}/pIRS/profiles"
    log:
        "logs/pIRS/{species}/{sample}_indelstat_profile.log"
    shell:
        "env time -v -o {output.timelog} indelstat_sam_bam {input} {params.outdir}/stats >{log} 2>&1"

rule pirs_profile_base_calling:
    input:
        ref=lambda wildcards: "storage/{0}/reference/{1}/{2}".format(wildcards.species, wildcards.sample, "sga-contigs.fa" if wildcards.asm != "" else "pilon-corrected.fa"),
        bam="storage/{species}/{sample}{asm}/real/mapping-bowtie2-s.bam"
    output:
        "work/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/baseCalling.count.matrix",
        "work/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/baseCalling.ratio.matrix",
        "work/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/baseCalling.stat",
        "work/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/baseCalling.info",
        timelog="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/pIRS/profiles/baseCalling-time.log"
    params:
        outdir="work/{species}/{sample}{asm}/pIRS/profiles"
    log:
        "logs/pIRS/{species}/{sample}{asm}_baseCalling_profile.log"
    shell:
        "env time -v -o {output.timelog} baseCalling_Matrix_calculator -m 42 -r {input.ref} -i {input.bam} -o {params.outdir}/baseCalling >{log} 2>&1"

###-ART--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rule art_cross_simulate:
    input:
        q1="work/{crossspecies}/{crosssample}/ART/profiles/qualityR1.txt",
        q2="work/{crossspecies}/{crosssample}/ART/profiles/qualityR2.txt",
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa"
    output:
        "storage/{species}/{sample}/ART_{crossspecies}_{crosssample}/art1.fq",
        "storage/{species}/{sample}/ART_{crossspecies}_{crosssample}/art2.fq",
        timelog="storage/{species}/{sample}/ART_{crossspecies}_{crosssample}/art-simulate-time.log"
    params:
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        ins1=lambda wildcards: SAMPLES[wildcards.crosssample]['insertion_rate1'],
        ins2=lambda wildcards: SAMPLES[wildcards.crosssample]['insertion_rate2'],
        del1=lambda wildcards: SAMPLES[wildcards.crosssample]['deletion_rate1'],
        del2=lambda wildcards: SAMPLES[wildcards.crosssample]['deletion_rate2'],
        readlength=lambda wildcards: SAMPLES[wildcards.crosssample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.crosssample]['insertlength'],
        outprefix="storage/{species}/{sample}/ART_{crossspecies}_{crosssample}/art",
        insertsd=lambda wildcards: SAMPLES[wildcards.crosssample]['insertsd'],
    log:
        "logs/ART/{species}/{sample}_art_{crossspecies}_{crosssample}_simulate.log"
    shell:
        "env time -v -o {output.timelog} art_illumina -1 {input.q1} -2 {input.q2} -f {params.coverage} -i {input.ref} -ir {params.ins1} -ir2 {params.ins2} -dr {params.del1} -dr2 {params.del2} "
        "-l {params.readlength} -m {params.insertlength} -o {params.outprefix} -p -s {params.insertsd} -sp >{log} 2>&1"

rule art_simulate_replicate:
    input:
        q1="work/{species}/{sample}/ART/profiles/qualityR1.txt",
        q2="work/{species}/{sample}/ART/profiles/qualityR2.txt",
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
    output:
        "storage/{species}/{sample}/ART/replicates/sim{simnum}/art1.fq",
        "storage/{species}/{sample}/ART/replicates/sim{simnum}/art2.fq"
    params:
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        ins1=lambda wildcards: SAMPLES[wildcards.sample]['insertion_rate1'],
        ins2=lambda wildcards: SAMPLES[wildcards.sample]['insertion_rate2'],
        del1=lambda wildcards: SAMPLES[wildcards.sample]['deletion_rate1'],
        del2=lambda wildcards: SAMPLES[wildcards.sample]['deletion_rate2'],
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength'],
        outprefix="storage/{species}/{sample}/ART/replicates/sim{simnum}/art",
        insertsd=lambda wildcards: SAMPLES[wildcards.sample]['insertsd']
    log:
        "logs/ART/{species}/{sample}_art_simulate_rep{simnum}.log"
    shell:
        "art_illumina -1 {input.q1} -2 {input.q2} -f {params.coverage} -i {input.ref} -ir {params.ins1} -ir2 {params.ins2} -dr {params.del1} -dr2 {params.del2} "
        "-l {params.readlength} -m {params.insertlength} -o {params.outprefix} -p -s {params.insertsd} -sp >{log} 2>&1"

rule art_simulate:
    input:
        q1="work/{species}/{sample}/ART/profiles/qualityR1.txt",
        q2="work/{species}/{sample}/ART/profiles/qualityR2.txt",
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa"
    output:
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ART/art{cov,|_cov}1.fq",
        "storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ART/art{cov,|_cov}2.fq",
        timelog="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/ART/art-simulate{cov,|_cov}-time.log"
    params:
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'] if "" == wildcards.asm or "_cov" == wildcards.cov else SAMPLES[wildcards.sample]['coverage_asm'],
        ins1=lambda wildcards: SAMPLES[wildcards.sample]['insertion_rate1'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertion_rate1_asm'],
        ins2=lambda wildcards: SAMPLES[wildcards.sample]['insertion_rate2'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertion_rate2_asm'],
        del1=lambda wildcards: SAMPLES[wildcards.sample]['deletion_rate1'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['deletion_rate1_asm'],
        del2=lambda wildcards: SAMPLES[wildcards.sample]['deletion_rate2'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['deletion_rate2_asm'],
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertlength_asm'],
        outprefix="storage/{species}/{sample}{asm}/ART/art{cov}",
        insertsd=lambda wildcards: SAMPLES[wildcards.sample]['insertsd'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertsd_asm']
    log:
        "logs/ART/{species}/{sample}{asm}{cov}_art_simulate.log"
    shell:
        "env time -v -o {output.timelog} art_illumina -1 {input.q1} -2 {input.q2} -f {params.coverage} -i {input.ref} -ir {params.ins1} -ir2 {params.ins2} -dr {params.del1} -dr2 {params.del2} "
        "-l {params.readlength} -m {params.insertlength} -o {params.outprefix} -p -s {params.insertsd} -sp >{log} 2>&1"

rule art_profile:
    input:
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
        #f1="input/{species}/{sample}/{sample}_1.fastq.gz",
        #f2="input/{species}/{sample}/{sample}_2.fastq.gz"
    output:
        "work/{species}/{sample}/ART/profiles/qualityR1.txt",
        "work/{species}/{sample}/ART/profiles/qualityR2.txt",
        timelog="storage/{species}/{sample}/ART/profiles/quality.log"
    log:
        "logs/ART/{species}/{sample}_quality_profile.log"
    params:
        profile="work/{species}/{sample}/ART/profiles/quality",
        workdir="work/{species}/{sample}/ART/profiles/"
    threads: max_threads
    shell:
        """
        cp {input.fq1} {params.workdir}/seq_1.fq.gz
        cp {input.fq2} {params.workdir}/seq_2.fq.gz
        env time -v -o {output.timelog} art_profiler_illumina {params.profile} {params.workdir} fq.gz {threads} >{log} 2>&1
        rm -f {params.workdir}/seq_1.fq.gz {params.workdir}/seq_2.fq.gz
        """

###-BEAR-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rule bear_simulate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2']),
        abundance="work/{species}/{sample}/BEAR/profiles/abundance.txt",
        per1="work/{species}/{sample}/BEAR/profiles/dri_read1.txt.per",
        per2="work/{species}/{sample}/BEAR/profiles/dri_read2.txt.per",
        errqual1="work/{species}/{sample}/BEAR/profiles/reads1.fastq.err.qual",
        errqual2="work/{species}/{sample}/BEAR/profiles/reads2.fastq.err.qual",
        errmat1="work/{species}/{sample}/BEAR/profiles/reads1.fastq.err.matr",
        errmat2="work/{species}/{sample}/BEAR/profiles/reads2.fastq.err.matr"
    output:
        fq1="storage/{species}/{sample}/BEAR/bear-R1.fq",
        fq2="storage/{species}/{sample}/BEAR/bear-R2.fq",
        uniformtl="storage/{species}/{sample}/BEAR/benchmark/uniform-time.log",
        trim1tl="storage/{species}/{sample}/BEAR/benchmark/trimming1-time.log",
        trim2tl="storage/{species}/{sample}/BEAR/benchmark/trimming2-time.log"
    log:
        r1="logs/BEAR/{species}/{sample}_trim1.log",
        r2="logs/BEAR/{species}/{sample}_trim2.log"
    params:
        workdir="work/{species}/{sample}/BEAR/",
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength'],
        insertsd=lambda wildcards: SAMPLES[wildcards.sample]['insertsd']
    shell:
        """
        mkdir -p {params.workdir}
        cd bin/BEAR/
        gzip -cd ../../{input.fq1} > ../../{params.workdir}/reads1.fastq
        gzip -cd ../../{input.fq2} > ../../{params.workdir}/reads2.fastq
        env time -v -o ../../{output.uniformtl} ./generate_reads.py -d -r ../../{input.ref} -a ../../{input.abundance} -o ../../{params.workdir}/uniform-reads -t $(($(wc -l ../../{params.workdir}/reads1.fastq | awk '{{print $1}}')/4)) -l {params.readlength} -i {params.insertlength} -s {params.insertsd}
        env time -v -o ../../{output.trim1tl} ./trim_reads.pl -i ../../{params.workdir}/reads1.fastq -f ../../{params.workdir}/uniform-reads.1.fasta -o ../../{params.workdir}/tmp_bear1.fq -r ../../{input.per1} -q ../../{input.errqual1} -m ../../{input.errmat1} 1>../../{log.r1} 2>&1
        env time -v -o ../../{output.trim2tl} ./trim_reads.pl -i ../../{params.workdir}/reads2.fastq -f ../../{params.workdir}/uniform-reads.2.fasta -o ../../{params.workdir}/tmp_bear2.fq -r ../../{input.per2} -q ../../{input.errqual2} -m ../../{input.errmat2} 1>../../{log.r2} 2>&1
        cat ../../{params.workdir}/tmp_bear1.fq | awk '{{if(1==NR%4){{print "@bear_" (NR-1)/4+1, substr($0, 2, length($0)-1)}}else{{print}}}}' > ../../{output.fq1}
        cat ../../{params.workdir}/tmp_bear2.fq | awk '{{if(1==NR%4){{print "@bear_" (NR-1)/4+1, substr($0, 2, length($0)-1)}}else{{print}}}}' > ../../{output.fq2}
        rm -f ../../{params.workdir}/reads{{1,2}}.fastq ../../{params.workdir}/tmp_bear{{1,2}}.fq
        cd -
        """

rule bear_profile:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        abundance="work/{species}/{sample}/BEAR/profiles/abundance.txt",
        abundancetl="storage/{species}/{sample}/BEAR/benchmark/abundance-time.log",
        per1="work/{species}/{sample}/BEAR/profiles/dri_read1.txt.per",
        drisee1tl="storage/{species}/{sample}/BEAR/benchmark/drisee1-time.log",
        per2="work/{species}/{sample}/BEAR/profiles/dri_read2.txt.per",
        drisee2tl="storage/{species}/{sample}/BEAR/benchmark/drisee2-time.log",
        errqual1="work/{species}/{sample}/BEAR/profiles/reads1.fastq.err.qual",
        errqual2="work/{species}/{sample}/BEAR/profiles/reads2.fastq.err.qual",
        errmat1="work/{species}/{sample}/BEAR/profiles/reads1.fastq.err.matr",
        errmat2="work/{species}/{sample}/BEAR/profiles/reads2.fastq.err.matr",
        err1tl="storage/{species}/{sample}/BEAR/benchmark/error-profile1-time.log",
        err2tl="storage/{species}/{sample}/BEAR/benchmark/error-profile2-time.log",
    log:
        drisee1="logs/BEAR/{species}/{sample}_drisee1.log",
        drisee2="logs/BEAR/{species}/{sample}_drisee2.log",
        err1="logs/BEAR/{species}/{sample}_errorqual1.log",
        err2="logs/BEAR/{species}/{sample}_errorqual2.log"
    params:
        workdir="work/{species}/{sample}/BEAR/profiles"
    threads: max_threads
    shell:
        """
        cd bin/BEAR/
        env time -v -o ../../{output.abundancetl} ./parametric_abundance.pl ../../{input.ref} low > ../../{output.abundance}
        gzip -cd ../../{input.fq1} > ../../{params.workdir}/reads1.fastq
        gzip -cd ../../{input.fq2} > ../../{params.workdir}/reads2.fastq
        env time -v -o ../../{output.drisee1tl} python2 drisee.py -p {threads} -t fastq -l ../../{log.drisee1} --percent ../../{params.workdir}/reads1.fastq ../../{params.workdir}/dri_read1.txt
        env time -v -o ../../{output.drisee2tl} python2 drisee.py -p {threads} -t fastq -l ../../{log.drisee2} --percent ../../{params.workdir}/reads2.fastq ../../{params.workdir}/dri_read2.txt
        env time -v -o ../../{output.err1tl} ./error_quality.pl ../../{params.workdir}/reads1.fastq ../../{params.workdir}/reads1.fastq.uc 1>../../{log.err1} 2>&1
        env time -v -o ../../{output.err2tl} ./error_quality.pl ../../{params.workdir}/reads2.fastq ../../{params.workdir}/reads2.fastq.uc 1>../../{log.err2} 2>&1
        rm -f ../../{params.workdir}/reads{{1,2}}.fastq
        cd -
        """

###-SInC-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rule sinc_simulate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        profile1=lambda wildcards: "work/{0}/{1}/SInC/profiles/{2}_bp_read_1_profile.txt".format(wildcards.species, wildcards.sample, SAMPLES[wildcards.sample]['readlength']),
        profile2=lambda wildcards: "work/{0}/{1}/SInC/profiles/{2}_bp_read_2_profile.txt".format(wildcards.species, wildcards.sample, SAMPLES[wildcards.sample]['readlength'])
    output:
        fq1="storage/{species}/{sample}/SInC/sinc-R1.fq",
        fq2="storage/{species}/{sample}/SInC/sinc-R2.fq",
        timelog="storage/{species}/{sample}/SInC/benchmark/sinc-simulation-time.log"
    log:
        "logs/SInC/{species}/{sample}_sinc_simulate.log"
    params:
        workdir="work/{species}/{sample}/SInC",
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength']-2*SAMPLES[wildcards.sample]['readlength'], # Inner insert length
        insertsd=lambda wildcards: SAMPLES[wildcards.sample]['insertsd'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage']
    threads: max_threads
    shell:
        """
        cat {input.ref} | awk '{{print $1}}' > {params.workdir}/sinc_refcp.fa
        env time -v -o {output.timelog} SInC_readGen -D "{params.insertlength}" -S {params.insertsd} -C {params.coverage} -T {threads} -R {params.readlength} {params.workdir}/sinc_refcp.fa {input.profile1} {input.profile2} 1>{log} 2>&1
        cat {params.workdir}/sinc_refcp.fa_1_{params.insertlength}_{params.insertsd}_{params.coverage}.0_{params.readlength}.fq | awk '{{if(1==NR%4){{print "@sinc_" (NR-1)/4+1, substr($0, 2, length($0)-1)}}else{{print}}}}' > {output.fq1}
        cat {params.workdir}/sinc_refcp.fa_2_{params.insertlength}_{params.insertsd}_{params.coverage}.0_{params.readlength}.fq | awk '{{if(1==NR%4){{print "@sinc_" (NR-1)/4+1, substr($0, 2, length($0)-1)}}else{{print}}}}' > {output.fq2}
        rm -f {params.workdir}/*sinc_refcp*
        """

rule sinc_profile:
    input:
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        profile1="work/{species}/{sample}/SInC/profiles/{readlength}_bp_read_1_profile.txt",
        prof1tl="storage/{species}/{sample}/SInC/profiles/benchmark/{readlength}_profile1-time.log",
        profile2="work/{species}/{sample}/SInC/profiles/{readlength}_bp_read_2_profile.txt",
        prof2tl="storage/{species}/{sample}/SInC/profiles/benchmark/{readlength}_profile2-time.log"
    params:
        workdir="work/{species}/{sample}/SInC/profiles/"
    shell:
        """
        mkdir -p {params.workdir}
        gzip -cd {input.fq1} > {params.workdir}/reads1.fq
        echo {params.workdir}/reads1.fq > {params.workdir}/list1.txt
        env time -v -o {output.prof1tl} genProfile -R 1 -l {wildcards.readlength} {params.workdir}/list1.txt
        mv {wildcards.readlength}_bp_read_1_profile.txt {output.profile1}
        gzip -cd {input.fq2} > {params.workdir}/reads2.fq
        echo {params.workdir}/reads2.fq > {params.workdir}/list2.txt
        env time -v -o {output.prof2tl} genProfile -R 2 -l {wildcards.readlength} {params.workdir}/list2.txt
        mv {wildcards.readlength}_bp_read_2_profile.txt {output.profile2}
        rm -rf {params.workdir}/reads1.fq {params.workdir}/reads2.fq
        """

###-NEAT-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rule neat_cross_simulate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        cp="work/{crossspecies}/{crosssample}/NEAT/profiles/gc_bias.p",
        fp="work/{crossspecies}/{crosssample}/NEAT/profiles/fraglen.p",
        ep="work/{crossspecies}/{crosssample}/NEAT/profiles/errors.p"
    output:
        fq1="storage/{species}/{sample}/NEAT_{crossspecies}_{crosssample}/neat_read1.fq",
        fq2="storage/{species}/{sample}/NEAT_{crossspecies}_{crosssample}/neat_read2.fq",
        timelog="storage/{species}/{sample}/NEAT_{crossspecies}_{crosssample}/benchmark/simulation-time.log"
    params:
        prefix="storage/{species}/{sample}/NEAT_{crossspecies}_{crosssample}/neat",
        readlength=lambda wildcards: SAMPLES[wildcards.crosssample]['readlength'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        diploid=lambda wildcards: "-p 2 -v storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else "-p 1"
    log:
        "logs/NEAT/{species}/{sample}_simulate_{crossspecies}_{crosssample}.log"
    shell:
        "env time -v -o {output.timelog} python2 bin/neat-genReads.py -r {input.ref} -R {params.readlength} -o {params.prefix} -c {params.coverage} -e {input.ep} -M 0 "
        "--pe-model {input.fp} --gc-model {input.cp} {params.diploid} 1>{log} 2>&1"

rule neat_simulate_replicate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        cp="work/{species}/{sample}/NEAT/profiles/gc_bias.p",
        fp="work/{species}/{sample}/NEAT/profiles/fraglen.p",
        ep="work/{species}/{sample}/NEAT/profiles/errors.p"
    output:
        "storage/{species}/{sample}/NEAT/replicates/sim{simnum}/neat_read1.fq",
        "storage/{species}/{sample}/NEAT/replicates/sim{simnum}/neat_read2.fq"
    params:
        prefix="storage/{species}/{sample}/NEAT/replicates/sim{simnum}/neat",
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'],
        diploid=lambda wildcards: "-p 2 -v storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else "-p 1"
    log:
        "logs/NEAT/{species}/{sample}_simulate_rep{simnum}.log"
    shell:
        "python2 bin/neat-genReads.py -r {input.ref} -R {params.readlength} -o {params.prefix} -c {params.coverage} -e {input.ep} -M 0 "
        "--pe-model {input.fp} --gc-model {input.cp} {params.diploid} 1>{log} 2>&1"

rule neat_simulate:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        vcf=lambda wildcards: "storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else [],
        cp="work/{species}/{sample}{asm}/NEAT/profiles/gc_bias.p",
        fp="work/{species}/{sample}{asm}/NEAT/profiles/fraglen.p",
        ep="work/{species}/{sample}/NEAT/profiles/errors.p"
    output:
        fq1="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/NEAT/neat{cov,|_cov}_read1.fq",
        fq2="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/NEAT/neat{cov,|_cov}_read2.fq",
        timelog="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/NEAT/benchmark/simulation{cov,|_cov}-time.log"
    params:
        prefix="storage/{species}/{sample}{asm}/NEAT/neat{cov}",
        readlength=lambda wildcards: SAMPLES[wildcards.sample]['readlength'],
        coverage=lambda wildcards: SAMPLES[wildcards.sample]['coverage'] if "" == wildcards.asm or "_cov" == wildcards.cov else SAMPLES[wildcards.sample]['coverage_asm'],
        diploid=lambda wildcards: "-p 2 -v storage/{0}/reference/{1}/pilon-corrected.vcf".format(wildcards.species, wildcards.sample) if REFERENCES[wildcards.species]['diploid'] else "-p 1"
    log:
        "logs/NEAT/{species}/{sample}{asm}{cov}_simulate.log"
    shell:
        "env time -v -o {output.timelog} python2 bin/neat-genReads.py -r {input.ref} -R {params.readlength} -o {params.prefix} -c {params.coverage} -e {input.ep} -M 0 "
        "--pe-model {input.fp} --gc-model {input.cp} {params.diploid} 1>{log} 2>&1"

rule neat_coverage_profile:
    input:
        ref=lambda wildcards: "storage/{0}/reference/{1}/{2}".format(wildcards.species, wildcards.sample, "sga-contigs.fa" if wildcards.asm != "" else "pilon-corrected.fa"),
        bam="storage/{species}/{sample}{asm}/real/mapping-bowtie2-s.bam"
    output:
        bias="work/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/NEAT/profiles/gc_bias.p",
        timelog_bedtools="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/NEAT/benchmark/bedtools-time.log",
        timelog_computegc="storage/{species}/{sample,[0-9A-Za-z]*}{asm,|_assembly}/NEAT/benchmark/computegc-time.log"
    params:
        workdir="work/{species}/{sample}{asm}/NEAT/profiles",
        insertlength=lambda wildcards: SAMPLES[wildcards.sample]['insertlength'] if "" == wildcards.asm else SAMPLES[wildcards.sample]['insertlength_asm']
    log:
        "logs/NEAT/{species}/{sample}{asm}_computegc.log"
    shell:
        """
        env time -v -o {output.timelog_bedtools} bedtools genomecov -d -ibam {input.bam} -g {input.ref} 1>{params.workdir}/genome_coverage.txt
        env time -v -o {output.timelog_computegc} python2 bin/neat/computeGC.py -r <(cat {input.ref} | awk '{{if(">" == substr($1,1,1)){{print $1}}else{{print $0}}}}') -i {params.workdir}/genome_coverage.txt -w {params.insertlength} -o {output.bias} 1>{log} 2>&1
        rm -f {params.workdir}/genome_coverage.txt
        """
        
rule neat_fraglen_profile:
    input:
        "storage/{species}/{sample}/real/mapping-bowtie2-s.bam"
    output:
        bias="work/{species}/{sample}/NEAT/profiles/fraglen.p",
        timelog="storage/{species}/{sample}/NEAT/benchmark/computeFraglen-time.log",
        timelog2="storage/{species}/{sample}/NEAT/benchmark/samtools-time.log"
    params:
        workdir="work/{species}/{sample}/NEAT/profiles"
    log:
        "logs/NEAT/{species}/{sample}_computeFraglen.log"  
    shell:
        """
        export org_path=$(pwd)
        cd {params.workdir}
        env time -v -o ${{org_path}}/{output.timelog2} samtools view ${{org_path}}/{input} | env time -v -o ${{org_path}}/{output.timelog} python2 ${{org_path}}/bin/neat/computeFraglen.py 1>${{org_path}}/{log} 2>&1
        cd -
        """
        
rule neat_error_profile:
    input:
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        bias="work/{species}/{sample}/NEAT/profiles/errors.p",
        timelog="storage/{species}/{sample}/NEAT/benchmark/genSeqErrorModel-time.log"
    log:
        "logs/NEAT/{species}/{sample}_genSeqErrorModel.log"  
    shell:
        "env time -v -o {output.timelog} python2 bin/neat/genSeqErrorModel.py -i {input.fq1} -i2 {input.fq2} -o {output.bias} --plot 1>{log} 2>&1"

###-Prepare input files for simulators-------------------------------------------------------------------------------------------------------------------------------------------------

rule apply_vcf:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        vcf="work/{species}/reference/{sample}/pilon-corrected.vcf.gz"
    output:
        "work/{species}/reference/{sample}/pilon-corrected-h{haplotype}.fa",
    log:
        "logs/apply_diploid_vcf/{species}/{sample}-h{haplotype}.log"
    shell:
        "bcftools consensus -f {input.ref} -H {wildcards.haplotype} -o {output} {input.vcf} 1>{log} 2>&1"
        
rule prepare_vcf_for_bcftools:
    input:
        "storage/{species}/reference/{sample}/pilon-corrected.vcf"
    output:
        "work/{species}/reference/{sample}/pilon-corrected.vcf.gz"
    log:
        "logs/prepare_vcf_for_bcftools/{species}/{sample}.log"
    shell:
        """
        bcftools view -Oz -o {output} {input}
        bcftools index {output}
        """

rule freebayes:
    input:
        ref="storage/{species}/reference/{sample}/pilon-corrected.fa",
        bam="storage/{species}/{sample}/real/mapping-bowtie2-s.bam",
    output:
        "storage/{species}/reference/{sample}/pilon-corrected.vcf"
    log:
        "logs/freebayes/{species}/{sample}.log"
    shell:
        "freebayes -f {input.ref} {input.bam} 2>{log} | awk '(\"#\" == substr($1, 1, 1) || (10 <= $6 && \"0/0\" != substr($10,1,3) && 0 == gsub(\"N\",\"N\",$4) + gsub(\"N\",\"N\",$5)) )' > {output}" 

rule getInsertLength:
    input:
        bam="storage/{species}/{sample}/real/mapping-bowtie2-s.bam",
        bai="storage/{species}/{sample}/real/mapping-bowtie2-s.bam.bai"
    output:
        pdf="storage/{species}/{sample}/real/mapping-bowtie2-s-insert-length.pdf",
        txt="storage/{species}/{sample}/real/mapping-bowtie2-s-insert-length.txt"
    shell:
        "./bin/plotInsertSize.py -o {output.pdf} {input.bam} 1>{output.txt}"

rule sga_assembly:
    input:
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        version="storage/{species}/reference/{sample}/sga.version",
        assembly="storage/{species}/reference/{sample}/sga-contigs.fa"
    params:
        tmp_dir="work/{species}/reference/{sample}/sga",
        out_prefix="storage/{species}/reference/{sample}/sga"
    log:
        preprocess="logs/sga_assembly/{species}/{sample}-preprocess.log",
        index="logs/sga_assembly/{species}/{sample}-index.log",
        correct="logs/sga_assembly/{species}/{sample}-correct.log",
        filter="logs/sga_assembly/{species}/{sample}-filter.log",
        overlap="logs/sga_assembly/{species}/{sample}-overlap.log",
        assemble="logs/sga_assembly/{species}/{sample}-assemble.log"
    threads: max_threads
    shell:
        """
        mkdir -p {params.tmp_dir}
        sga --version > {output.version}
        sga preprocess -p 1 -o {params.tmp_dir}/reads.fq {input.fq1} {input.fq2} 1>{log.preprocess} 2>&1
        sga index -t {threads} -a ropebwt -p {params.tmp_dir}/reads {params.tmp_dir}/reads.fq 1>{log.index} 2>&1
        sga correct -t {threads} -o {params.tmp_dir}/reads.ec.fq -p {params.tmp_dir}/reads {params.tmp_dir}/reads.fq 1>{log.correct} 2>&1
        sga index -t {threads} -a ropebwt -p {params.tmp_dir}/reads.ec {params.tmp_dir}/reads.ec.fq 1>{log.index} 2>&1
        sga filter -t {threads} -o {params.tmp_dir}/reads.ec.filter.pass.fa {params.tmp_dir}/reads.ec.fq 1>{log.filter} 2>&1
        sga overlap -t {threads} {params.tmp_dir}/reads.ec.filter.pass.fa 1>{log.overlap} 2>&1
        mv reads.ec.filter.pass.asqg.gz {params.tmp_dir}/reads.ec.filter.pass.asqg.gz
        sga assemble -o {params.out_prefix} {params.tmp_dir}/reads.ec.filter.pass.asqg.gz 1>{log.assemble} 2>&1
        """

rule pilon:
    input:
        ref=lambda wildcards: "input/{0}/reference/{1}.fa".format(wildcards.species, REFERENCES[wildcards.species]['prefix']),
        bam="work/{species}/{sample}/real/mapping-bowtie2-s-uncorrected.bam",
        bai="work/{species}/{sample}/real/mapping-bowtie2-s-uncorrected.bam.bai"
    output:
        protected("storage/{species}/reference/{sample}/pilon-corrected.fa"),
        protected("storage/{species}/reference/{sample}/pilon-corrected.changes")
    params:
        outdir="storage/{species}/reference/{sample}",
        diploid=lambda wildcards: "--diploid" if REFERENCES[wildcards.species]['diploid'] else ""
    log:
        "logs/pilon/{species}/{sample}.log"
    threads: max_threads
    resources:
        mem_mb=max_mem_mb
    shell:
        """
        java -Xmx{resources.mem_mb}M -jar bin/pilon-1.21.jar --genome {input.ref} --frags {input.bam} --output pilon-corrected --outdir {params.outdir} {params.diploid} --changes --threads {threads} 1>{log} 2>&1
        mv {params.outdir}/pilon-corrected.fasta {params.outdir}/pilon-corrected.fa 
        """

rule samtools_index:
    input:
        "{bam}"
    output:
        "{bam}.bai"
    shell:
        "samtools index {wildcards.bam}"

rule bwa_map_mod_names:
    input:
        lambda wildcards: "work/{0}/reference/{1}index_bwa/{2}.bwt".format(wildcards.species, "" if "work" == wildcards.folder else "{}/".format(wildcards.sample), REFERENCES[wildcards.species]['prefix'] if "work" == wildcards.folder else "pilon-corrected"),
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        "{folder}/{species}/{sample}/real/mapping-bwa-s{infotag,.*}.bam"
    params:
        refind=lambda wildcards: "work/{0}/reference/{1}index_bwa/{2}".format(wildcards.species, "" if "work" == wildcards.folder else "{}/".format(wildcards.sample), REFERENCES[wildcards.species]['prefix'] if "work" == wildcards.folder else "pilon-corrected"),
        tmp="work/{species}/{sample}/real/tmp_"
    log:
        "logs/bwa_map_mod_names/{sample}{infotag}.log"
    threads: max_threads
    shell:
        "bwa mem -t {threads} {params.refind}"
        " <(reseq-prepare-names.py {input.fq1} {input.fq2})"
        " <(reseq-prepare-names.py {input.fq2} {input.fq1})"
        " 2>{log} | samtools sort -m 10G -@ 4 -T {params.tmp} -o {output} -"

rule bwa_index:
    input:
        lambda wildcards: "{0}/{1}/reference/{2}{3}.fa".format("storage" if len(wildcards.sampledir) else "input",wildcards.species,wildcards.sampledir,wildcards.refname)
    output:
        "work/{species}/reference/{sampledir,.*}index_bwa/{refname}.bwt"
    params:
        refind="work/{species}/reference/{sampledir}index_bwa/{refname}"
    log:
        "logs/bwa_index/{species}/{sampledir}{refname}.log"
    shell:
        "bwa index -p {params.refind} {input} >{log} 2>&1"

rule bowtie2_map_mod_names:
    input:
        lambda wildcards: "work/{0}/reference/{1}index_bowtie2/{2}.1.bt2".format(wildcards.species, "" if "work" == wildcards.folder else "{}/".format(wildcards.sample), REFERENCES[wildcards.species]['prefix'] if "work" == wildcards.folder else ("sga-contigs" if "_assembly" == wildcards.sample_tag else "pilon-corrected")),
        fq1=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq1']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq1']),
        fq2=lambda wildcards: "work/baseSpace/{0}/{1}".format(SAMPLES[wildcards.sample]['basespace'],SAMPLES[wildcards.sample]['fq2']) if 'basespace' in SAMPLES[wildcards.sample] else "input/{0}/{1}/{2}".format(wildcards.species,wildcards.sample,SAMPLES[wildcards.sample]['fq2'])
    output:
        "{folder}/{species}/{sample,[0-9A-Za-z]*}{sample_tag,|_assembly}/real/mapping-bowtie2-s{infotag,.*}.bam"
    params:
        refind=lambda wildcards: "work/{0}/reference/{1}index_bowtie2/{2}".format(wildcards.species, "" if "work" == wildcards.folder else "{}/".format(wildcards.sample), REFERENCES[wildcards.species]['prefix'] if "work" == wildcards.folder else ("sga-contigs" if "_assembly" == wildcards.sample_tag else "pilon-corrected")),
        tmp="work/{species}/{sample}/real/tmp_"
    log:
        "logs/bowtie2_map_mod_names/{species}/{folder}_{sample}{sample_tag}{infotag}.log"
    threads: max_threads
    shell:
        "bowtie2 -p {threads} -X 2000 -x {params.refind}"
        " -1 <(reseq-prepare-names.py {input.fq1} {input.fq2})"
        " -2 <(reseq-prepare-names.py {input.fq2} {input.fq1})"
        " 2>{log} | samtools sort -m 10G -@ 4 -T {params.tmp} -o {output} -"

rule bowtie2_index:
    input:
        lambda wildcards: "{0}/{1}/reference/{2}{3}.fa".format("storage" if len(wildcards.sampledir) else "input",wildcards.species,wildcards.sampledir,wildcards.refname)
    output:
        "work/{species}/reference/{sampledir,.*}index_bowtie2/{refname}.1.bt2"
    params:
        refind="work/{species}/reference/{sampledir}index_bowtie2/{refname}"
    log:
        "logs/bowtie2_index/{species}/{sampledir}{refname}.log"
    shell:
        "bowtie2-build {input} {params.refind} >{log} 2>&1"

rule bcl2fastq:
    input:
        "input/baseSpace/{project}/SampleSheet.csv"
    output:
        "work/baseSpace/{project}/DeliverySummary.csv"
    params:
        indir="input/baseSpace/{project}/",
        outdir="work/baseSpace/{project}/",
        sample_sheet="work/baseSpace/{project}_SampleSheet.csv"
    log:
        "logs/bcl2fastq/{project}.log"
    threads: max_threads
    shell:
        """
        bcl2fastq --version 1>{params.outdir}/bcl2fastq.version 2>$1
        cat {input} | awk '{{ if(""==$2){{ if(NR!=16 && NR!=17){{print $0}}}}else{{print $1 "_" $2}} }}' | awk 'BEGIN{{ FS=","; OFS = "," }}{{if(""!=$10 && "Lane" != $1){{print $1, $2 "_L00" $1, $3 "_L00" $1, $4, $5, $6, $7, $8, $9, $10, $11}}else{{print $0}}}}' > {params.sample_sheet}
        bcl2fastq -p {threads} --ignore-missing-positions -R {params.indir} -o {params.outdir} --sample-sheet {params.sample_sheet} 1>{log} 2>&1
        """
