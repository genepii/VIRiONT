#!usr/bin/en python3
import os
import re
import glob
import string
import pandas as pd

configfile : "config/params.yaml"

datapath=config['data_loc']
if (datapath[-1] != "/"):
	datapath=datapath+"/"
resultpath=config['result_loc']
if (resultpath[-1] != "/"):
	resultpath=resultpath+"/"
refpath=config['ref_loc']
trim_min=config['min_length']
trim_max=config['max_length']
trim_head=config['head_trim']
trim_tail=config['tail_trim']
quality_read=config['min_qual_ONT']
MI_cutoff=config['MI_cutoff']
mincov_cons=config['min_cov']
variant_frequency=config['Vfreq']
mpileup_depth=config['max_depth']
mpileup_basequal=config['basequality']
mutation_research=config['HBV_mut']
mutation_table_path=config['path_table']
if (mutation_table_path[-1] != "/"):
	mutation_table_path=mutation_table_path+"/"
min_freq=config['freq_min']
window=config['window_pos']
#slurm needs
pipeline_loc=config['pipeline_loc']

#CANU
max_threads=int(int(config['thread_number'])/2)
max_mem=int(int(config['mem_cost'])/2)
do_correction=config['correction']
coverage_correction=int(config['cov_correction'])
error_rate=config['error_rate']

#Amplicon removal
do_amplicon_removal=config['do_amplicon_removal']
bedfile_path=config['bedfilepath']

#Read MI results from VIRiONT_MI1.py
data_multiinf = resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv"
try:
    multiinf_table = pd.read_csv(resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv",sep="\t",names=["barcode", "reference","readcount","ratio"])
except:
    sys.exit("The '04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv' file produced during the blast step is missing. Exiting.")
sample_list=list(multiinf_table['barcode'])
reference_list=list(multiinf_table['reference'])
assoc_sample_ref_number=len(sample_list)
BARCODE = list(dict.fromkeys(sample_list))

if (len(BARCODE)==0):
    sys.exit("The '04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv' file is empty. Any data available for analysis. Exiting.")


#get database name
filename=os.path.basename(refpath)
list_split=filename.split(".")
database_name=list_split[0]

######## INTERMEDIATE FILES ########
#produce read list
readlist=[]
for i in range(0,assoc_sample_ref_number):
	readlist.append(resultpath +"READLIST/"+sample_list[i]+"/"+reference_list[i]+"_readlist.txt")
#produce filtered fastq
fastq_filtered=[]
for i in range(0,assoc_sample_ref_number):
	fastq_filtered.append(resultpath +"05_REFILTERED_FASTQ/"+sample_list[i]+"/"+reference_list[i]+"_filtered.fastq.gz")
######## CORRECTION ########
#produce corrected fastq if enabled
corrected_fastq=[]
for i in range(0,assoc_sample_ref_number):
	corrected_fastq.append(resultpath +"05_REFILTERED_FASTQ/"+sample_list[i]+"/"+reference_list[i]+"_filtered_corrected.fastq.gz")
######## PRECONSENSUS FILES ########
#produce bam
bamfile_precons=[]
for i in range(0,assoc_sample_ref_number):
	bamfile_precons.append(resultpath +"06_PRECONSENSUS/BAM/"+sample_list[i]+"/"+reference_list[i]+"_sorted.bam")
#produce vcf files
vcffile_precons=[]
for i in range(0,assoc_sample_ref_number):
	vcffile_precons.append(resultpath +"06_PRECONSENSUS/VCF/"+sample_list[i]+"/"+reference_list[i]+"_sammpileup.vcf")
#produce consensus files
consfile=[]
for i in range(0,assoc_sample_ref_number):
	consfile.append(resultpath +"06_PRECONSENSUS/SEQUENCES/"+sample_list[i]+"/"+reference_list[i]+"_cons.fasta")
######## FINAL CONSENSUS FILES ########
#produce bam aligned on preconsensus
bamfile=[]
for i in range(0,assoc_sample_ref_number):
	bamfile.append(resultpath +"07_BAM/"+sample_list[i]+"/"+reference_list[i]+"_sorted.bam")
#produce vcf files on preconsensus
vcffile=[]
for i in range(0,assoc_sample_ref_number):
	vcffile.append(resultpath +"08_VCF/"+sample_list[i]+"/"+reference_list[i]+"_sammpileup.vcf")
#produce final consensus files
consfile_final=[]
for i in range(0,assoc_sample_ref_number):
	consfile_final.append(resultpath +"09_CONSENSUS/"+sample_list[i]+"/"+reference_list[i]+"_cons.fasta")
######## MUTATION HBV ########
#copy vcf result from reference alignment for mutation screening
vcf_mut=[]
for i in range(0,assoc_sample_ref_number):
	vcf_mut.append(resultpath +"13_MUTATION_SCREENING/"+sample_list[i]+"/"+reference_list[i]+"_sammpileup.vcf_variants.txt")
#produce table with mutation screening if enabled
table_mut=[]
if (mutation_research==True):
    for i in range(0,assoc_sample_ref_number):
	    table_mut.append(resultpath +"13_MUTATION_SCREENING/"+sample_list[i]+"/all_results/"+sample_list[i]+"_"+reference_list[i]+"_PreCore.csv")
######## COVERAGE ANALYSIS ######## 
#produce cov files
covfile=[]
for i in range(0,assoc_sample_ref_number):
	covfile.append(resultpath +"12_COVERAGE/"+sample_list[i]+"/"+reference_list[i]+".cov")
######## QC METRIC FILES ########
#produce filterseq files
filtseqfilecount=[]
for i in range(0,assoc_sample_ref_number):
	filtseqfilecount.append(resultpath +"10_QC_ANALYSIS/"+sample_list[i]+"/"+reference_list[i]+"_filterseqcount.csv")


rule pipeline_output:
    input:
        ######## MUTATION HBV ########
        #vcf_mut,
        table_mut,
        ######## INTERMEDIATE FILES ########
        #readlist,
        #filtseqfile,
        #fastq_filtered,
        ######## CORRECTION ########
        #corrected_fastq,
        ######## PRECONSENSUS FILES ########
        #bamfile_precons, 
        #vcffile_precons,
        #consfile_precons,
        ######## FINAL CONSENSUS FILES ########  
        #bamfile,
        #vcffile,
        #consfile,    
        cons_seq = resultpath+"09_CONSENSUS/all_cons.fasta",
        ######## COVERAGE ANALYSIS ########  
        #covfile,
        cov_plot = resultpath+"12_COVERAGE/cov_plot.pdf",
        ######## QC PHYLOGENETIC TREE FILES ########  
        #allseq = resultpath+"11_PHYLOGENETIC_TREE/allseq.fasta",
        #align_seq = resultpath+"11_PHYLOGENETIC_TREE/align_seq.fasta",
        #NWK_tree = resultpath+"11_PHYLOGENETIC_TREE/IQtree_analysis.treefile",
        tree_pdf = resultpath+"11_PHYLOGENETIC_TREE/RADIAL_tree.pdf",
        ######## QC METRIC FILES ########  
        full_summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_table.csv",
        report = resultpath+"param_file.txt",

######## RULES ########
rule write_param_used:
    message:
        "write an analysis report with all used parameters."
    output:
        report= resultpath+"param_file.txt"
    run:
        textfile=open(resultpath+"param_file.txt","w")
        textfile.write("######################\n")
        textfile.write("#### PARAMS USED #####\n")
        textfile.write("######################\n")
        textfile.write("data repository:"+datapath+"\n")
        textfile.write("result repository:"+resultpath+"\n")
        textfile.write("database used:"+refpath+"\n")
        textfile.write("read minlength:"+str(trim_min)+"\n")
        textfile.write("read maxlength:"+str(trim_max)+"\n")
        textfile.write("quality filtering:"+str(quality_read)+"\n")
        textfile.write("read 5' trimming length:"+str(trim_head)+"\n")
        textfile.write("read 3' trimming length:"+str(trim_tail)+"\n")
        textfile.write("min coverage for consensus generation:"+str(mincov_cons)+"\n")
        textfile.write("variant calling depth:"+str(mpileup_depth)+"\n")
        textfile.write("variant calling base quality threeshold:"+str(mpileup_basequal)+"\n")
        textfile.write("multi-infection cutoff:"+str(MI_cutoff)+"\n")
        textfile.write("variant frequency:"+str(variant_frequency)+"\n")
        textfile.write("mutation research:"+str(mutation_research)+"\n")
        textfile.write("read correction enabled:"+str(do_correction)+"\n")
        textfile.write("maximum of threads usable for CANU:"+str(max_threads)+"\n")
        textfile.write("maximum of memory usable for CANU:"+str(max_mem)+"\n")
        textfile.write("coverage correction:"+str(coverage_correction)+"\n")
        textfile.write("corrected error rate for CANU:"+str(error_rate)+"\n")
        textfile.close()

rule getfastqlist:
    message:
        "Get read names from {wildcards.barcode} matching with {wildcards.reference} using R."
    input:
        blastn_result = resultpath+"04_BLASTN_ANALYSIS/{barcode}_blastnR.tsv"
    output:
        fastqlist = temp(resultpath +"READLIST/{barcode}/{reference}_readlist.txt" )
    #conda:
    #    "env/Renv.yaml"  
    params:
        VIRiONT_path = pipeline_loc
    shell:
        """
        Rscript {params.VIRiONT_path}script/get_read_list.R {input.blastn_result} {wildcards.reference} {output.fastqlist}
        """

rule extract_matching_read:
	message:
		"Extract reads from {wildcards.barcode} matching with the '{wildcards.reference}' reference using seqkit."
	input:
		read_list = rules.getfastqlist.output.fastqlist,
		nonhuman_fastq = resultpath + '03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq.gz'      
	output:
		merged_filtered = temp(resultpath +"05_REFILTERED_FASTQ/{barcode}/{reference}_filtered.fastq"),
		merged_filtered_compressed = resultpath +"05_REFILTERED_FASTQ/{barcode}/{reference}_filtered.fastq.gz" 
	#conda:
	#	"env/seqkit.yaml"          
	shell:
		"""
		seqkit grep --pattern-file {input.read_list} {input.nonhuman_fastq} > {output.merged_filtered}
        gzip -c {output.merged_filtered} > {output.merged_filtered_compressed}
		"""   

rule read_correction:
    message:
        "correct the {wildcards.barcode}-{wildcards.reference} reads for the ONT error rate using CANU."
    input:
        merged_filtered = rules.extract_matching_read.output.merged_filtered_compressed ,
        split_ref_path = resultpath+"00_SUPDATA/REFSEQ/" ,
    output:
        corrected_filtered_uncompressed = temp(resultpath +"05_REFILTERED_FASTQ/{barcode}/{reference}_filtered_corrected.fastq"), 
        corrected_filtered = resultpath +"05_REFILTERED_FASTQ/{barcode}/{reference}_filtered_corrected.fastq.gz" 
    #conda:
    #    "env/canu.yaml"
    threads: 
        max_threads
    resources: 
        mem_mb= max_mem
    params:
        correction_rep = resultpath+"00_SUPDATA/read_correction/{barcode}/{reference}/"
    shell:
        """
        #compute reference size
        size_ref=`tail -n 1 {input.split_ref_path}/{wildcards.reference}.fasta | wc -c`
        #launch CANU in correction mode only
        canu -correct \
            -nanopore {input.merged_filtered} \
            -p {wildcards.reference} \
            -d {params.correction_rep} \
            genomeSize=${{size_ref}} \
            maxMemory={max_mem} \
            maxThreads={threads} \
            correctedErrorRate={error_rate} \
            maxInputCoverage={coverage_correction} 
        #process final modifications
        seqtk seq -F "?" {params.correction_rep}{wildcards.reference}.correctedReads.fasta.gz > {output.corrected_filtered_uncompressed}
        gzip -c {output.corrected_filtered_uncompressed} > {output.corrected_filtered}
        """

rule filtered_fastq_alignemnt:
    message:
        "Align the {wildcards.barcode}-{wildcards.reference} reads on '{wildcards.reference}' reference using minimap2."
    input:   
        merged_filtered = rules.extract_matching_read.output.merged_filtered_compressed if (do_correction==False) else rules.read_correction.output.corrected_filtered ,
        split_ref_path = resultpath+"00_SUPDATA/REFSEQ/" ,
    output:
        spliced_bam = resultpath+"06_PRECONSENSUS/BAM/{barcode}/{reference}_sorted.bam" 
    #conda:
    #    "env/minimap2.yaml"              
    shell:
        """
        minimap2 -ax splice {input.split_ref_path}/{wildcards.reference}.fasta {input.merged_filtered} | samtools sort > {output.spliced_bam}
        samtools index {output.spliced_bam}  
        """
rule remove_amplicon:
    message:
        "Remove amplicon bam if bed supplied."
    input:
        spliced_bam = rules.filtered_fastq_alignemnt.output.spliced_bam,
        bedfile = bedfile_path
    output:
        ampliconclip_bam = resultpath+"06_PRECONSENSUS/BAM/{barcode}/{reference}_ampliconclip_sorted.bam" 
    #conda:
    #    "env/samtools.yaml"  
    shell:
        """
        samtools ampliconclip --both-ends -b {input.bedfile} {input.spliced_bam} | samtools sort > {output.ampliconclip_bam}
        samtools index {output.ampliconclip_bam}
        """    
rule compute_coverage:
    message:
        "Compute read coverage from {wildcards.barcode} for the '{wildcards.reference}' reference using bedtools."
    input:
        sorted_bam = rules.filtered_fastq_alignemnt.output.spliced_bam if (do_amplicon_removal==False) else rules.remove_amplicon.output.ampliconclip_bam
    output:
        coverage = resultpath+"12_COVERAGE/{barcode}/{reference}.cov" 
    #conda:
    #    "env/bedtools.yaml"          
    shell:
        """
        bedtools genomecov -ibam {input} -d -split | sed 's/$/\t{wildcards.barcode}\tMINION/' > {output.coverage} 
        """

rule plot_coverage:
    message:
        "Plot coverage summary for all barcodes using R."
    input:
        cov = covfile
    output:
        cov_sum = resultpath+"12_COVERAGE/cov_sum.cov" ,
        cov_plot = resultpath+"12_COVERAGE/cov_plot.pdf"
    #conda:
    #    "env/Renv.yaml"
    params:
        VIRiONT_path = pipeline_loc 
    shell:
        """
        cat {input.cov} > {output.cov_sum}
        Rscript {params.VIRiONT_path}script/plot_cov_MI.R {output.cov_sum} {output.cov_plot}
        """

rule variant_calling:
    message:
        "Count nucleotidic base repartition from {wildcards.barcode} aligned reads for each position of the '{wildcards.reference}' reference  using samtools."
    input:
        split_ref_path = resultpath+"00_SUPDATA/REFSEQ/"  ,
        sorted_bam = rules.filtered_fastq_alignemnt.output.spliced_bam if (do_amplicon_removal==False) else rules.remove_amplicon.output.ampliconclip_bam ,
    output:
        vcf = resultpath+"06_PRECONSENSUS/VCF/{barcode}/{reference}_sammpileup.vcf"
    #conda:
    #    "env/samtools.yaml"  
    shell:
        """
        samtools mpileup -d {mpileup_depth} -Q {mpileup_basequal} -f {input.split_ref_path}/{wildcards.reference}.fasta {input.sorted_bam}  > {output}
        """   

rule generate_consensus:
    message:
        "Generate pre-consensus sequence for the {wildcards.barcode} based on the '{wildcards.reference}' reference."
    input:
        vcf = rules.variant_calling.output.vcf
    output:
        fasta_cons_temp = temp(resultpath+"06_PRECONSENSUS/SEQUENCES/{barcode}/{reference}_cons_temp.fasta") ,
        fasta_cons = resultpath+"06_PRECONSENSUS/SEQUENCES/{barcode}/{reference}_cons.fasta",
        vcf = resultpath+"06_PRECONSENSUS/VCF/{barcode}/{reference}_sammpileup.vcf_variants.txt"
    params:
        VIRiONT_path = pipeline_loc
    shell:
        """
        perl {params.VIRiONT_path}script/pathogen_varcaller_MINION.PL {input} 0.5 {output.fasta_cons_temp} {mincov_cons}
        sed  's/>.*/>{wildcards.barcode}_PRECONS/' {output.fasta_cons_temp} > {output.fasta_cons}
        """

rule precons_alignemnt:
    message:
        "Align the '{wildcards.reference}' selected reads from the {wildcards.barcode} on the pre-consensus sequence."
    input:
        fasta_cons = rules.generate_consensus.output.fasta_cons,
        merged_filtered = rules.extract_matching_read.output.merged_filtered_compressed if (do_correction==False) else rules.read_correction.output.corrected_filtered ,
    output:
        consbam = resultpath +"07_BAM/{barcode}/{reference}_sorted.bam"
    #conda:
    #    "env/minimap2.yaml"              
    shell:
        """
        minimap2 -ax splice {input.fasta_cons} {input.merged_filtered} | samtools sort > {output.consbam}
        samtools index {output.consbam}  
        """

rule variant_calling_precons:
    message:
        "Nucleotidic base count for the {wildcards.barcode} aligned on the pre-consensus sequence."
    input:
        fasta_cons = rules.generate_consensus.output.fasta_cons,
        consbam = rules.precons_alignemnt.output.consbam ,
    output:
        vcf = resultpath+"08_VCF/{barcode}/{reference}_sammpileup.vcf"
    #conda:
    #    "env/samtools.yaml"  
    shell:
        """
        samtools mpileup -d {mpileup_depth} -Q {mpileup_basequal} -f {input.fasta_cons} {input.consbam}  > {output.vcf}
        """   

rule generate_finalconsensus:
    message:
        "Generate the final consensus sequence for the {wildcards.barcode}."
    input:
        vcf = rules.variant_calling_precons.output.vcf
    output:
        fasta_cons_temp = temp(resultpath+"09_CONSENSUS/{barcode}/{reference}_cons_temp.fasta") ,
        fasta_cons = resultpath+"09_CONSENSUS/{barcode}/{reference}_cons.fasta"
    params:
        VIRiONT_path = pipeline_loc
    shell:
        """
        perl {params.VIRiONT_path}script/pathogen_varcaller_MINION.PL {input} {variant_frequency} {output.fasta_cons_temp} {mincov_cons}
        sed  's/>.*/>{wildcards.barcode}_{wildcards.reference}_ONT{variant_frequency}/' {output.fasta_cons_temp} > {output.fasta_cons}
        """


rule read_metric:
    message:
        "Extract read sequence from  raw / dehosted / filtered {wildcards.barcode} reads for metric computation."
    input:
        raw_fastq = resultpath+"01_MERGED/{barcode}_merged.fastq.gz"  ,
        trimmed_fastq = resultpath+"03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq.gz" ,
        dehosted_fastq = resultpath + '02_DEHOSTING/{barcode}_meta.fastq.gz' ,
    output:
        raw_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_rawseq.txt"),
        trimmed_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_trimmseq.txt"),
        dehosted_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}_dehostseq.txt"),
        raw_count = temp(resultpath+"10_QC_ANALYSIS/{barcode}_rawcount.csv"),
        trimm_count = temp(resultpath+"10_QC_ANALYSIS/{barcode}_trimmcount.csv"),
        dehost_count = temp(resultpath+"10_QC_ANALYSIS/{barcode}_dehostcount.csv"),
    run:
        shell("zcat {input.raw_fastq} | awk '(NR%4==2)'  > {output.raw_read}")
        shell("zcat {input.trimmed_fastq} | awk '(NR%4==2)'  > {output.trimmed_read}")
        shell("zcat {input.dehosted_fastq} | awk '(NR%4==2)'  > {output.dehosted_read}")
        ######
        rawlengthfile=open(output.raw_count,'w')
        rawRfile=open(output.raw_read,'r')
        for line in rawRfile:
            rawlengthfile.write(str(len(line.rstrip()))+";"+wildcards.barcode+";01_RAW_FASTQ;NONE\n")
        rawRfile.close()
        rawlengthfile.close()
        ######
        trimlengthfile=open(output.trimm_count,'w')
        trimRfile=open(output.trimmed_read,'r')
        for line in trimRfile:
            trimlengthfile.write(str(len(line.rstrip()))+";"+wildcards.barcode+";03_TRIMM_FASTQ;NONE\n")
        trimRfile.close()
        trimlengthfile.close()
        ######
        dehostlengthfile=open(output.dehost_count,'w')
        dehostRfile=open(output.dehosted_read,'r')
        for line in dehostRfile:
            dehostlengthfile.write(str(len(line.rstrip()))+";"+wildcards.barcode+";02_DEHOST_FASTQ;NONE\n")
        dehostRfile.close()
        dehostlengthfile.close()

rule read_metric_MI:
    message:
        "Extract read sequence from '{wildcards.reference}' matched reads in the {wildcards.barcode} for metric computation."
    input:
        bestmatched_fastq = rules.extract_matching_read.output.merged_filtered_compressed if (do_correction==False) else rules.read_correction.output.corrected_filtered ,
    output:
        bestmatched_read = temp(resultpath+"10_QC_ANALYSIS/{barcode}/{reference}_filterseq.txt"),
        bestmatched_count = temp(resultpath+"10_QC_ANALYSIS/{barcode}/{reference}_filterseqcount.csv"),
    run:
        shell("zcat {input.bestmatched_fastq} | awk '(NR%4==2)'  > {output.bestmatched_read}")
        ######
        BMlengthfile=open(output.bestmatched_count,'w')
        BMRfile=open(output.bestmatched_read,'r')
        for line in BMRfile:
            BMlengthfile.write(str(len(line.rstrip()))+";"+wildcards.barcode+";05_REFILTERED_FASTQ;"+wildcards.reference+"\n")
        BMRfile.close()
        BMlengthfile.close()

rule concat_datametric:
    input:
        rawcount = expand(rules.read_metric.output.raw_count,barcode=BARCODE),
        trimcount = expand(rules.read_metric.output.trimm_count,barcode=BARCODE),
        dehostcount = expand(rules.read_metric.output.dehost_count,barcode=BARCODE),
        refilteredcount = expand(filtseqfilecount),
    output:
        allraw = temp(resultpath+"10_QC_ANALYSIS/all_rawcount.csv"),
        alldehost = temp(resultpath+"10_QC_ANALYSIS/all_dehostcount.csv"),
        alltrim = temp(resultpath+"10_QC_ANALYSIS/all_trimcount.csv"),
        allrefilter = temp(resultpath+"10_QC_ANALYSIS/all_refcount.csv"),
    run:
        shell("cat {input.rawcount} > {output.allraw}")
        shell("cat {input.dehostcount} > {output.alldehost}")
        shell("cat {input.trimcount} > {output.alltrim}")
        shell("cat {input.refilteredcount} > {output.allrefilter}")

rule compute_metric:
    message:
        "Produce final metric table."
    input:
        allraw = rules.concat_datametric.output.allraw,
        alldehost = rules.concat_datametric.output.alldehost,
        alltrim = rules.concat_datametric.output.alltrim,
        allrefilter = rules.concat_datametric.output.allrefilter,
    output:
        full_summ_table = resultpath+"10_QC_ANALYSIS/METRIC_summary_table.csv",
    #conda:
    #    "env/Renv.yaml"
    params:
        VIRiONT_path = pipeline_loc
    shell:
        """
        Rscript {params.VIRiONT_path}script/compute_metrics.R {input.allraw} {input.alldehost} \
            {input.alltrim} {input.allrefilter} {output.full_summ_table}
        """

rule filterIncompleteSeq:
    message:
        "Filter consensus sequences containing N bases above 10%."
    input:
        allcons = expand(consfile_final),
    output:
        cons_seq = resultpath+"09_CONSENSUS/all_cons.fasta",
        filteredcons = resultpath+"11_PHYLOGENETIC_TREE/filteredcons.fasta",
    run:
        shell("cat {input.allcons} > {output.cons_seq}")
        fastas={}
        fasta_file=open(output.cons_seq,"r")
        for line in fasta_file:
            if (line[0]=='>'):
                header=line
                fastas[header]=''
            else:
                fastas[header]+=line
        fasta_file.close()
        #
        filtered_fasta=open(output.filteredcons,'w')
        for header,sequence in fastas.items():
            sequence_ok=sequence.translate({ord(c): None for c in string.whitespace}).upper()
            seqlen=len(sequence_ok)
            Ncount=sequence_ok.count("N")
            CompPerc=(seqlen-Ncount)/seqlen*100
            if (CompPerc >= 90):
                filtered_fasta.write(header.rstrip("\n")+"\n")
                filtered_fasta.write(sequence_ok+"\n")
        filtered_fasta.close()         

rule prepareSEQ:
    message:
        "Concatenate all consensus sequences with all reference sequences."
    input:
        filteredcons = rules.filterIncompleteSeq.output.filteredcons,
        reference_file = refpath
    output:
        allseq = temp(resultpath+"11_PHYLOGENETIC_TREE/allseq.fasta"),
    shell:
        """
        cat {input.reference_file} {input.filteredcons} > {output.allseq}
        """

rule sequenceAlign:
    message:
        "Multiple Alignment on all sequences using Muscle."
    input:
        allseq = rules.prepareSEQ.output.allseq
    output:
        align_seq = temp(resultpath+"11_PHYLOGENETIC_TREE/align_seq.fasta")
    #conda:
    #    "env/muscle.yaml"
    shell:
        "muscle -align {input} -output {output}"

rule buildTree:
    message:
        "Build Newick tree from alignment using iqtree."
    input:
        align_seq = rules.sequenceAlign.output.align_seq,
    output:
        iqtree = expand(resultpath+"11_PHYLOGENETIC_TREE/IQtree_analysis"+".{ext}", ext=["contree","bionj","ckp.gz","iqtree","log","mldist","splits.nex","treefile"])
    #conda:
    #    "env/iqtree.yaml"
    shell:
        "iqtree -s {input.align_seq} -m K80 -B 1000 -T AUTO --prefix {resultpath}11_PHYLOGENETIC_TREE/IQtree_analysis "

rule plotTree:
    message:
        "Plot radial tree using ETE 3."
    input:
        iqtree = rules.buildTree.output.iqtree ,
        NWK_data = resultpath+"11_PHYLOGENETIC_TREE/IQtree_analysis.treefile"
    output:
        tree_pdf = resultpath+"11_PHYLOGENETIC_TREE/RADIAL_tree.pdf"
    #conda:
    #    "env/ETE3.yaml"
    params:
        VIRiONT_path = pipeline_loc
    shell:
        "python3 {params.VIRiONT_path}script/makeTREE.py {input.NWK_data} {output.tree_pdf} "

rule copy_vcf_result:
    message:
        "Copying vcf results into the mutation folder."
    input:
        vcf = rules.generate_consensus.output.vcf ,
    output:
        vcf = resultpath+"13_MUTATION_SCREENING/{barcode}/{reference}_sammpileup.vcf_variants.txt"
    shell:
        "cp {input.vcf} {output.vcf}"

rule search_HBV_mutation:
    message:
        "Search mutation when a whole genome HBV reference is provided."
    input:
        vcf = rules.copy_vcf_result.output.vcf,
        #table_mut = mutation_table_path + "mutation_{reference}.csv"
    params:
        table_mut = mutation_table_path,
        VIRiONT_path = pipeline_loc
    output:
        #raw
        result_PC = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_PreCore.csv",
        result_DS = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_DomaineS.csv",
        result_RT = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_DomaineRT.csv",
        result_BCP = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_BCP.csv",
        result_DPS1 = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_DomainePreS1.csv",
        result_DPS2 = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_DomainePreS2.csv",
        result_DHBx = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_DomaineHBx.csv",
        result_C = resultpath+"13_MUTATION_SCREENING/{barcode}/all_results/{barcode}_{reference}_Core.csv",
        #filtered
        result_PC_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_PreCore.csv",
        result_DS_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_DomaineS.csv",
        result_RT_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_DomaineRT.csv",
        result_BCP_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_BCP.csv",
        result_DPS1_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_DomainePreS1.csv",
        result_DPS2_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_DomainePreS2.csv",
        result_DHBx_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_DomaineHBx.csv",
        result_C_F = resultpath+"13_MUTATION_SCREENING/{barcode}/filtered/{barcode}_{reference}_Core.csv",
    #conda:
    #    "env/Renv.yaml"
    shell:
        """
        Rscript {params.VIRiONT_path}script/search_mutation.R {input.vcf} {params.table_mut} {min_freq} {window} \
            {output.result_PC} {output.result_BCP} {output.result_DS} {output.result_RT} \
            {output.result_DPS1} {output.result_DPS2} {output.result_DHBx} {output.result_C} \
            {output.result_PC_F} {output.result_BCP_F} {output.result_DS_F} {output.result_RT_F} \
            {output.result_DPS1_F} {output.result_DPS2_F} {output.result_DHBx_F} {output.result_C_F} 
        """