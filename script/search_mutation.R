library(stringr)
library(Biostrings)
library(tidyr)

#setwd("C:/Users/regueex/Desktop/FIGURE_PAPIER_CARO")
argv <- commandArgs(TRUE)

#Parametres
vcf_input<-argv[1]
tablemut_input<-argv[2]
filter_threshold<-as.numeric(as.character(argv[3]))
window_pos<-as.numeric(as.character(argv[4]))

file_output_PC<-argv[5]
file_output_BCP<-argv[6]
file_output_DS<-argv[7]
file_output_RT<-argv[8]
file_output_DPS1<-argv[9]
file_output_DPS2<-argv[10]
file_output_DHBx<-argv[11]
file_output_C<-argv[12]

file_output_PC_F<-argv[13]
file_output_BCP_F<-argv[14]
file_output_DS_F<-argv[15]
file_output_RT_F<-argv[16]
file_output_DPS1_F<-argv[17]
file_output_DPS2_F<-argv[18]
file_output_DHBx_F<-argv[19]
file_output_C_F<-argv[20]

#Functions
getAA<-function(codon){
  #codon<-"AC"
  if (str_detect(codon,"N")==T) {
    AA<-"ININT"
  } else if (nchar(codon)==3) {
    #coder un truc pour quand les codons n'existent pas?
    list_codon<-names(GENETIC_CODE)
    if(codon %in% list_codon) AA<-GENETIC_CODE[[codon]] else AA<-"ININT"
  } else AA<-"ININT"
  return(AA)
}
searchMUT_CODON<-function(vcf,mut_table,mutation){
  #mutation<-"Domaine RT"
  mut_table<-subset(table_mut,Region==mutation)
  mut_table<-mut_table[,c("Position_EcoR1","NUC1_P3","NUC2_P3","NUC3_P3")]
  mut_table<-pivot_longer(data=mut_table,
                          col=c("NUC1_P3","NUC2_P3","NUC3_P3"),
                          names_to="POS_TYPE",
                          values_to = "REF_POS")
  mut_table$NUM_CODON<-mut_table$Position_EcoR1
  mut_table<-mut_table[,c("REF_POS","POS_TYPE","NUM_CODON")]
  #move window if no primer in reference
  mut_table$REF_POS<-mut_table$REF_POS+window_pos
  table_vcf_mut<-merge(vcf,mut_table,by="REF_POS")
  if(nrow(table_vcf_mut)==0) {
      return(table_vcf_mut)
  } else {
  CODON_REF<-aggregate(table_vcf_mut$REF,list(table_vcf_mut$base_status,table_vcf_mut$NUM_CODON),paste,collapse="")

  colnames(CODON_REF)<-c("base_status","NUM_CODON","CODON_REF")
  table_vcf_mut<-merge(table_vcf_mut,CODON_REF,by=c("NUM_CODON","base_status"))
  
  table_vcf_mut$CODON_ALT<-table_vcf_mut$CODON_REF
  table_vcf_mut$CODON_ALT<-ifelse(table_vcf_mut$POS_TYPE=="NUC1_P3",paste0(table_vcf_mut$base,str_sub(table_vcf_mut$CODON_REF,start=2,end=3)),table_vcf_mut$CODON_ALT)
  table_vcf_mut$CODON_ALT<-ifelse(table_vcf_mut$POS_TYPE=="NUC2_P3",paste0(str_sub(table_vcf_mut$CODON_REF,start=1,end=1),table_vcf_mut$base,str_sub(table_vcf_mut$CODON_REF,start=3,end=3)),table_vcf_mut$CODON_ALT)
  table_vcf_mut$CODON_ALT<-ifelse(table_vcf_mut$POS_TYPE=="NUC3_P3",paste0(str_sub(table_vcf_mut$CODON_REF,start=1,end=2),table_vcf_mut$base),table_vcf_mut$CODON_ALT)
  
  table_vcf_mut$REF_AA<-sapply(table_vcf_mut$CODON_REF,getAA)
  table_vcf_mut$ALT_AA<-sapply(table_vcf_mut$CODON_ALT,getAA)
  table_vcf_mut$warning<-ifelse(table_vcf_mut$REF_AA==table_vcf_mut$ALT_AA,"OK","alerte!")
  
  table_vcf_mut$Mutation_name<-paste0(table_vcf_mut$REF_AA,as.character(table_vcf_mut$NUM_CODON),table_vcf_mut$ALT_AA)
  
  table_vcf_mut<-table_vcf_mut[,c("REFERENCE","REF_POS","REF","total_count",
                                  "base_status","base","count","NUM_CODON","POS_TYPE","REF_AA","ALT_AA","Mutation_name","freq","warning")]
  table_vcf_mut<-table_vcf_mut[with(table_vcf_mut, order(REF_POS,POS_TYPE)),]
  return(table_vcf_mut)
}
}
searchMUT_NT<-function(vcf,mut_table,mutation){
  mut_table<-subset(table_mut,Region==mutation)
  mut_table<-mut_table[,c("Position_EcoR1","NUC1_P3")]
  #move window if no primer in reference
  mut_table$NUC1_P3<-mut_table$NUC1_P3+window_pos
  table_vcf_mut<-merge(vcf,mut_table,by.x="REF_POS",by.y="NUC1_P3")

  table_vcf_mut$Mutation_name<-paste0(table_vcf_mut$REF,as.character(table_vcf_mut$Position_EcoR1),table_vcf_mut$base)
  table_vcf_mut<-table_vcf_mut[,c("REFERENCE","REF_POS","Position_EcoR1","REF","total_count",
                                  "base_status","base","count","Mutation_name","freq")]
  table_vcf_mut<-table_vcf_mut[with(table_vcf_mut, order(REF_POS,base_status)),]
  return(table_vcf_mut)
}

#Read VCF
table_vcf<-read.delim(vcf_input,sep=" ",header=F)
table_vcf<-table_vcf[,-c(3,10,11)]
colnames(table_vcf)<-c("REFERENCE","REF_POS","COVERAGE","REF","majo","mino_1st","mino_2nd","mino_3d")
table_vcf<-subset(table_vcf,REFERENCE!="EOF@EOF")

table_vcf$REF<-str_sub(table_vcf$REF,start = 1,end = 1)

table_vcf<-pivot_longer(data=table_vcf,
                        col=c("majo","mino_1st","mino_2nd","mino_3d"),
                        names_to="base_status",
                        values_to = "base_count")

table_vcf$base<-str_sub(table_vcf$base_count,start = 1,end = 1)
table_vcf$count<-as.numeric(str_sub(table_vcf$base_count,start = 3))
total_count<-aggregate(table_vcf$count,by=list(table_vcf$REF_POS),sum)
colnames(total_count)<-c("REF_POS","total_count")

table_vcf<-merge(table_vcf,total_count,by="REF_POS")
table_vcf$freq<-round(table_vcf$count/table_vcf$total_count*100,digits=2)

#get table to use
reference<-unique(as.character(table_vcf$REFERENCE))
reference<-str_replace(reference,"P3_GT","")
reference<-str_replace(reference,"GT","")
if(grepl("A",reference)) table<-"mutation_GTA.csv"
if(grepl("B",reference)) table<-"mutation_GTB.csv"
if(grepl("C",reference)) table<-"mutation_GTC.csv"
if(grepl("D",reference)) table<-"mutation_GTD.csv"
if(grepl("E",reference)) table<-"mutation_GTE.csv"
if(grepl("F",reference)) table<-"mutation_GTF.csv"
if(grepl("G",reference)) table<-"mutation_GTG.csv"
if(grepl("H",reference)) table<-"mutation_GTH.csv"
if(grepl("I",reference)) table<-"mutation_GTI.csv"

if(grepl("GTG|GTH",table)==TRUE){
  write.table("Mutation table not available for this genotype",file_output_PC)
  write.table("Mutation table not available for this genotype",file_output_BCP)
  write.table("Mutation table not available for this genotype",file_output_DS)
  write.table("Mutation table not available for this genotype",file_output_RT)
  write.table("Mutation table not available for this genotype",file_output_DPS1)
  write.table("Mutation table not available for this genotype",file_output_DPS2)
  write.table("Mutation table not available for this genotype",file_output_DHBx)
  write.table("Mutation table not available for this genotype",file_output_C)
  write.table("Mutation table not available for this genotype",file_output_PC_F)
  write.table("Mutation table not available for this genotype",file_output_BCP_F)
  write.table("Mutation table not available for this genotype",file_output_DS_F)
  write.table("Mutation table not available for this genotype",file_output_RT_F)
  write.table("Mutation table not available for this genotype",file_output_DPS1_F)
  write.table("Mutation table not available for this genotype",file_output_DPS2_F)
  write.table("Mutation table not available for this genotype",file_output_DHBx_F)
  write.table("Mutation table not available for this genotype",file_output_C_F)
  
  quit(save = "no")
}

#read mutation table
table_mut<-read.csv2(paste0(tablemut_input,table))
table_mut$Region<-as.character(table_mut$Region)


############# Nucleotid mut research #############
table_vcf_PC<-searchMUT_NT(table_vcf,table_mut,"PreCore")
table_vcf_BCP<-searchMUT_NT(table_vcf,table_mut,"BCP")

############# Protein mut research #############
table_vcf_RT<-searchMUT_CODON(table_vcf,table_mut,"Domaine RT")
table_vcf_C<-searchMUT_CODON(table_vcf,table_mut,"Core")
table_vcf_DS<-searchMUT_CODON(table_vcf,table_mut,"Domaine S")
table_vcf_DPS1<-searchMUT_CODON(table_vcf,table_mut,"Domaine PreS1")
table_vcf_DPS2<-searchMUT_CODON(table_vcf,table_mut,"Domaine PreS2")
table_vcf_DHBx<-searchMUT_CODON(table_vcf,table_mut,"Domaine HBx")

#filtered tables
############# Nucleotid mut research #############
table_vcf_PC_filter<-subset(table_vcf_PC,freq>=filter_threshold)
table_vcf_BCP_filter<-subset(table_vcf_BCP,freq>=filter_threshold)

############# Protein mut research #############
table_vcf_RT_filter<-subset(table_vcf_RT,freq>=filter_threshold)
table_vcf_C_filter<-subset(table_vcf_C,freq>=filter_threshold)
table_vcf_DS_filter<-subset(table_vcf_DS,freq>=filter_threshold)
table_vcf_DPS1_filter<-subset(table_vcf_DPS1,freq>=filter_threshold)
table_vcf_DPS2_filter<-subset(table_vcf_DPS2,freq>=filter_threshold)
table_vcf_DHBx_filter<-subset(table_vcf_DHBx,freq>=filter_threshold)

############# Write tables #############
write.csv2(table_vcf_BCP,file_output_BCP,row.names = F,quote = F)
write.csv2(table_vcf_PC,file_output_PC,row.names = F,quote = F)
write.csv2(table_vcf_RT,file_output_RT,row.names = F,quote = F)
write.csv2(table_vcf_DS,file_output_DS,row.names = F,quote = F)
write.csv2(table_vcf_DPS1,file_output_DPS1,row.names = F,quote = F)
write.csv2(table_vcf_DPS2,file_output_DPS2,row.names = F,quote = F)
write.csv2(table_vcf_DHBx,file_output_DHBx,row.names = F,quote = F)
write.csv2(table_vcf_C,file_output_C,row.names = F,quote = F)

############# Write tables #############
write.csv2(table_vcf_BCP_filter,file_output_BCP_F,row.names = F,quote = F)
write.csv2(table_vcf_PC_filter,file_output_PC_F,row.names = F,quote = F)
write.csv2(table_vcf_RT_filter,file_output_RT_F,row.names = F,quote = F)
write.csv2(table_vcf_DS_filter,file_output_DS_F,row.names = F,quote = F)
write.csv2(table_vcf_DPS1_filter,file_output_DPS1_F,row.names = F,quote = F)
write.csv2(table_vcf_DPS2_filter,file_output_DPS2_F,row.names = F,quote = F)
write.csv2(table_vcf_DHBx_filter,file_output_DHBx_F,row.names = F,quote = F)
write.csv2(table_vcf_C_filter,file_output_C_F,row.names = F,quote = F)
