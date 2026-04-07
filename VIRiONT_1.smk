#!usr/bin/en python3
import os
import glob
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
bitscore=config['bitscore_min']
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

#Check if input data are present / get all barcodes in a list after demultiplexing 
barcode_list = glob.glob(datapath+"barcode*")
if (len(barcode_list) < 1):
	sys.exit("No barcode repositories found in the indicated path. Please check the 'data_loc' parameter and be sure that barcode* directories are here. Exiting.")
else:
	BARCODE=[]
	IGNORED_BARCODE=[]
	BAD_BARCODE=[]
	list_invalid=[]
	for BC in barcode_list:
		barcode=str(os.path.basename(BC))
		#barcode rep empty?
		if len(os.listdir(BC) ) == 0:
			IGNORED_BARCODE.append(barcode)
		else:
			valid_BC=True
			for file in os.listdir(BC):
				test_file=BC+file
				if test_file.lower().endswith(('.fastq','.gz')) == False:
					valid_BC=False
					list_invalid.append(test_file)
			if valid_BC == True:
				BARCODE.append(barcode)
			else:
				BAD_BARCODE.append(barcode)

#Check reference fasta file sanity / get database name
def read_fasta(file):   
    fastas={}
    fasta_file=open(file,"r")
    for line in fasta_file:
        if (line[0]=='>'):
            header=line
            fastas[header]=''
        else:
            fastas[header]+=line
    fasta_file.close()
    return fastas
try:
	ref_list=read_fasta(refpath) 
except:
	sys.exit("No fasta reference file found in the indicated path. Please check the 'ref_loc' parameter. Exiting.")
if (len(ref_list) < 1):
	sys.exit("No fasta sequence can be detected in the indicated 'ref_loc' file. Please check your format file. Exiting.")
else:
	filename=os.path.basename(refpath)
	list_split=filename.split(".")
	database_name=list_split[0]

#Check parameter sanity
#trim_min
try:
	trim_min=int(trim_min)
except:
	sys.exit("A numeric value is expected for the 'min_length' parameter. Please check this parameter. Exiting.")
#trim_max
try:
	trim_max=int(trim_max)
except:
	sys.exit("A numeric value is expected for the 'max_length' parameter. Please check this parameters. Exiting.")
#trim_head
try:
	trim_head=int(trim_head)
except:
	sys.exit("A numeric value is expected for the 'head_trim' parameter. Please check this parameters. Exiting.")
#trim_tail
try:
	trim_tail=int(trim_tail)
except:
	sys.exit("A numeric value is expected for the 'tail_trim' parameter. Please check this parameters. Exiting.")
#quality_read
try:
	quality_read=int(quality_read)
except:
	sys.exit("A numeric value is expected for the 'min_qual_ONT' parameter. Please check this parameters. Exiting.")
#MI_cutoff
try:
	MI_cutoff=float(MI_cutoff)
except:
	sys.exit("A numeric value is expected for the 'MI_cutoff' parameter. Please check this parameters. Exiting.")
if (MI_cutoff==0 or MI_cutoff>100):
	sys.exit("The Multi Infection cutoff should be set over 0 and not exceed 100. Exiting.")

#mpileup_depth
try:
	mpileup_depth=int(mpileup_depth)
except:
	sys.exit("A numeric value is expected for the 'max_depth' parameter. Please check this parameters. Exiting.")
#mpileup_basequal
try:
	mpileup_basequal=int(mpileup_basequal)
except:
	sys.exit("A numeric value is expected for the 'basequality' parameter. Please check this parameters. Exiting.")
#mincov_cons
try:
	mincov_cons=int(mincov_cons)
except:
	sys.exit("A numeric value is expected for the 'min_cov' parameter. Please check this parameters. Exiting.")
#variant_frequency
try:
	variant_frequency=float(variant_frequency)
except:
	sys.exit("A numeric value is expected for the 'Vfreq' parameter. Please check this parameters. Exiting.")
if ((variant_frequency==0) or (variant_frequency>=1)):
	sys.exit("The Variant frequency should be set between 0 and 1. Please check the 'Vfreq' parameters. Exiting.")

#Check if mutation table are present 
if (mutation_research=="TRUE"):
	mut_table_list = glob.glob(mutation_table_path+"mutation*.csv")
	if (len(mut_table_list) < 1):
		sys.exit("Mutation research enabled and no mutation table found in the indicated path. Please check the 'path_table' parameter. Exiting.")
	try:
		min_freq=float(min_freq)
	except:
		sys.exit("A numeric value is expected for the 'freq_min' parameter. Please check this parameters. Exiting.")
	if (min_freq==0 or min_freq>100):
		sys.exit("The 'freq_min' cutoff should be set over 0 and not exceed 100. Exiting.")
	try:
		window=int(window)
	except:
		sys.exit("A numeric value is expected for the 'window' parameter. Please check this parameters. Exiting.")
	

#NanoFilt Command construction with input parameters
if (trim_min>0):
	lengmin="--minlength "+str(trim_min)
else:
	lengmin=""
if (trim_max>0):
	lengmax="--maxlength "+str(trim_max)
else:
	lengmax=""
if (trim_head>0):
	head="--headcrop "+str(trim_head)
else:
	head=""
if (trim_tail>0):
	tail="--tailcrop "+str(trim_tail)
else:
	tail=""
if (quality_read>0):
	qfiltering="--quality "+str(quality_read)
else:
	qfiltering=""


#final output
rule pipeline_ending:
	input:
		######## INTERMEDIATE FILES ########
		#database = expand(resultpath+"00_SUPDATA/DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"]),
		#merged_file = expand(resultpath+'01_MERGED/{barcode}_merged.fastq',barcode=BARCODE),
		#human_bam = expand(resultpath+"02_DEHOSTING/{barcode}_human.bam",barcode=BARCODE),
		#meta_fastq = expand(resultpath + '02_DEHOSTING/{barcode}_meta.fastq',barcode=BARCODE),
		#trimmed_file = expand(resultpath+'03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq',barcode=BARCODE),
		#converted_fastq = expand(resultpath+"FASTA/{barcode}.fasta", barcode=BARCODE),
		#ref_rep=resultpath+"00_SUPDATA/REFSEQ/",
		#R_data = expand(resultpath+"04_BLASTN_ANALYSIS/{barcode}_fmt.txt" ,barcode=BARCODE),
		#merged_data = resultpath+"04_BLASTN_ANALYSIS/ALL_refcount.tsv"
		######## FINAL OUTPUTS ########		
		#blastn_result=expand(resultpath+"04_BLASTN_ANALYSIS/{barcode}_blastnR.tsv",barcode=BARCODE),
		summ_multiinf = resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv",
		fastq_content = resultpath+"fastq_content.txt"

rule fastq_analysis:
	message:
		"Describe barcode repositories content."
	input:
	output:
		fastq_content = resultpath+"fastq_content.txt"
	run:
		textfile=open(resultpath+"fastq_content.txt","w")
		textfile.write("##########################\n")
		textfile.write("##### FASTQ ANALYSIS #####\n")
		textfile.write("##########################\n")
		textfile.write("barcode repository containing fastq/gz files and used for analysis:\n")
		for bc in BARCODE:
			textfile.write(bc+"\n")
		textfile.write("##########################\n")
		textfile.write("barcode repository containing other files than fastq/gz and ignored for analysis:\n")
		for bc in IGNORED_BARCODE:
			textfile.write(bc+"\n")
		textfile.write("##########################\n")
		textfile.write("list of problematic files:\n")
		for file in list_invalid:
			textfile.write(file+"\n")
		textfile.write("##########################\n")
		textfile.write("empty barcode repositories and ignored for analysis:\n")
		for bc in IGNORED_BARCODE:
			textfile.write(bc+"\n")
		textfile.close()

rule merging_fastq:
	message:
		"Merging the fastq/fastq.gz files from the {wildcards.barcode}/ folder."
	input: 
		barcode_rep = datapath+"{barcode}/",
	output: 
		merged_fastq = temp(resultpath+"01_MERGED/{barcode}_merged.fastq"),
		merged_fastq_compressed = resultpath+"01_MERGED/{barcode}_merged.fastq.gz"  
	shell: 
		"""
		gz_file_num=`find {input.barcode_rep} -name "*.gz" | wc -l`
		if [ $gz_file_num -gt 0 ]
		then
        	zcat {input.barcode_rep}/* > {output.merged_fastq}
		else
        	cat {input.barcode_rep}/* > {output.merged_fastq}
		fi
		gzip -c {output.merged_fastq} > {output.merged_fastq_compressed}
		"""

#rule get_hg19:
#	message:
#		"Download if necessary the hg19 reference genome. Stored in the VIRiONT/ref/ folder. Executed once per VIRiONT installation."
#	output:
#		hg19_ref =  "ref/hg19.fa"
#	shell:
#		"""
#		wget -P ref/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#		gunzip ref/hg19.fa.gz       
#		"""

#rule index_hg19:
#	message:
#		"Indexing the hg19 reference genome for a quicker dehosting. Executed once per VIRiONT installation."
#	input:
#		hg19 = rules.get_hg19.output.hg19_ref
#	output:
#		hg19_index = "ref/hg19.mmi"
	#conda:
	#	"env/minimap2.yaml"
#	shell:
#		"minimap2 -d {output} {input}"

rule hg19_dehosting:
	message:
		"Aligning reads from {wildcards.barcode} fastq on human genome for identifying host reads using minimap2."
	input:
		merged_fastq = rules.merging_fastq.output.merged_fastq_compressed ,
		ref_file=  pipeline_loc +"ref/hg19.mmi"
	output:
		human_bam = temp(resultpath+"02_DEHOSTING/{barcode}_human.bam")
	#conda:
	#	"env/minimap2.yaml" 
	threads: 4
	shell:
		"""
		minimap2 -t {threads} -ax splice {input.ref_file} {input.merged_fastq}  | samtools view -b > {output.human_bam}
		"""    

rule nonhuman_read_extract:
	message:
		"Extract unaligned reads from {wildcards.barcode} using samtools."
	input:
		human_bam = rules.hg19_dehosting.output.human_bam
	output:
		nonhuman_bam = temp(resultpath+"02_DEHOSTING/{barcode}_meta.bam")
	#conda:
	#	"env/samtools.yaml" 
	shell: 
		"samtools view -b -f 4 {input.human_bam} > {output.nonhuman_bam} "   

rule converting_bam_fastq:
	message:
		"Convert the {wildcards.barcode} bam into fastq containing only viral reads using bedtools."
	input:
		nonhuman_bam = rules.nonhuman_read_extract.output.nonhuman_bam 
	output:
		nonhuman_fastq = temp(resultpath + '02_DEHOSTING/{barcode}_meta.fastq'),
		nonhuman_fastq_compressed = resultpath + '02_DEHOSTING/{barcode}_meta.fastq.gz',
	#conda:
	#	"env/bedtools.yaml"
	shell:
		"""
		bedtools bamtofastq  -i {input.nonhuman_bam} -fq {output.nonhuman_fastq} 
		gzip -c {output.nonhuman_fastq} > {output.nonhuman_fastq_compressed} 
		"""

rule trimming_fastq:
	message:
		"Filtering and trimming viral reads from {wildcards.barcode} using chopper with input parameters."
	input:
		meta_fastq = rules.converting_bam_fastq.output.nonhuman_fastq_compressed
	output:
		trimmed_fastq = temp(resultpath+"03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq"),
		trimmed_fastq_compressed = resultpath+"03_FILTERED_TRIMMED/{barcode}_filtered_trimmed.fastq.gz"

	#conda:
	#	"env/chopper.yaml"
	shell:
		"""
		gunzip -c {input.meta_fastq} | chopper {lengmin} {lengmax} {head} {tail} {qfiltering}  > {output.trimmed_fastq} 
		gzip -c {output.trimmed_fastq} > {output.trimmed_fastq_compressed} 
		"""

rule converting_fastq_fasta:
	message:
		"Converting the {wildcards.barcode} reads into fasta sequences for blastn research using seqtk."
	input:
		nonhuman_fastq = rules.trimming_fastq.output.trimmed_fastq
	output:
		converted_fastq = temp(resultpath+"FASTA/{barcode}.fasta" ),
		converted_fastq_temp = temp(resultpath+"FASTA/{barcode}_temp.fasta" )
	#conda:
	#	"env/seqtk.yaml"        
	shell:
		"""
		seqtk seq -A {input.nonhuman_fastq} > {output.converted_fastq}
		sed '/^>/d' {output.converted_fastq} > {output.converted_fastq_temp}
		"""    

rule split_reference:
	message:
		"Spliting the input reference for isolate each genotype sequence if multiple reference are given."
	input:
		ref_file= refpath
	output:
		ref_rep=directory(resultpath+"00_SUPDATA/REFSEQ/")
	params:
		VIRiONT_path = pipeline_loc
	shell:
		"""
		mkdir -p {output} 
		{params.VIRiONT_path}script/split_reference.py {input} {output} 
		"""

rule make_db:
	message:
		"Build blast database from input reference {refpath} using makeblastdb."
	input:
		ref_fasta_file = refpath 
	output:
		database = expand(resultpath+"00_SUPDATA/DB/"+database_name+".{ext}", ext=["nhr", "nin", "nsq"])
	params:
		database_path = resultpath+"00_SUPDATA/DB/"+database_name
	#conda:
	#	"env/blast.yaml"          
	shell:
		"""
		makeblastdb -in {input.ref_fasta_file} -out {params.database_path} -input_type fasta -dbtype nucl
		"""

rule blastn_ref:
	message:
		"Blasting {wildcards.barcode} reads on the custom database."
	input: 
		fasta_file = rules.converting_fastq_fasta.output.converted_fastq,
		database = rules.make_db.output.database ,
	output:
		R_data = resultpath+"04_BLASTN_ANALYSIS/{barcode}_fmt.txt" 
	threads: 4
	params:
		database_path = resultpath+"00_SUPDATA/DB/"+database_name    
	#conda:
	#	"env/blast.yaml"               
	shell:
		"""
		blastn -db {params.database_path} -query {input.fasta_file} \
			-outfmt "6 qseqid sseqid bitscore slen qlen length pident" \
			-out {output.R_data} -num_threads {threads}
		"""   

rule count_refmatching:
	message:
		"Count matching reads for each references in the {wildcards.barcode}."
	input:
		R_data = rules.blastn_ref.output.R_data ,
		ref_table = rules.split_reference.output.ref_rep,
	output:
		ref_count = temp(resultpath+"04_BLASTN_ANALYSIS/{barcode}_refcount.tsv"),
		blastn_result = resultpath+"04_BLASTN_ANALYSIS/{barcode}_blastnR.tsv"
	params:
		VIRiONT_path = pipeline_loc
	#conda:
	#	"env/Renv.yaml"     
	shell:
		"""
		if [ -s {input.R_data} ] 
		then
			Rscript {params.VIRiONT_path}script/count_ref.R {input.R_data} \
				{input.ref_table}/R_table_analysis.csv \
				{wildcards.barcode} \
				{output.ref_count} \
				{output.blastn_result} \
				{bitscore}
		else
			touch {output.ref_count}
			echo "Empty: Any blast results." > {output.blastn_result}
		fi
		"""  

rule MI_analysis:
	message:
		"Summarize matched references for each barcode into the '04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv' file, mandatory for the following."
	input:
		count_ref_data = expand(rules.count_refmatching.output.ref_count,barcode=BARCODE)
	output:
		merged_data = temp(resultpath+"04_BLASTN_ANALYSIS/ALL_refcount.tsv"),
		plot_pdf = resultpath+"04_BLASTN_ANALYSIS/read_repartition.pdf",
		summ_multiinf = resultpath+"04_BLASTN_ANALYSIS/SUMMARY_Multi_Infection.tsv"
	#conda:
	#	"env/Renv.yaml" 
	params:
		VIRiONT_path = pipeline_loc  
	shell:
		"""
		cat {input.count_ref_data} > {output.merged_data}
		Rscript {params.VIRiONT_path}script/MI_analysis.R {output.merged_data} \
			{MI_cutoff} \
			{output.plot_pdf} \
			{output.summ_multiinf}
		"""