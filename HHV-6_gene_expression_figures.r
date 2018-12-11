###trying to fix this shit by rewriting it all#####
#####taking advantage of local blast here#####
library(Biostrings)
library(Rsamtools)
library(GenomicAlignments)
library(edgeR)
library(rtracklayer)
library(xlsx)
library(ggplot2)
library(gtools)
library(plotrix)
library(foreach)
library(doParallel)
library(seqinr)
library(scales)
library(RColorBrewer)
library(svMisc)

setwd('/Users/gerbix/Documents/vikas/Thesis_from_other_mac/hhv6b_rna/')

data <- read.xlsx("gtex_icihhv6_hhv6b_positives.xlsx", sheetIndex = 2)
datatrimmed <- data[, -c(5:9),(12:13)]
data<-data[complete.cases(data$Run), ]
data<-data[data$LibrarySelection=='cDNA',]
#reads z29 gff3, trims out repeats
gff3<-readGFF('hhv6breference.gff3',version=3,columns=NULL, tags=NULL, filter=NULL, nrows=-1, raw_data=FALSE)
for (j in 1:3) {
  for( i in 1:nrow(gff3)){ 
    if(( gff3$rpt_type[i]=='TANDEM')==TRUE && is.na(gff3$rpt_type[i])==FALSE){ 
      print('ok')
      gff3<-gff3[-i,]
      
    }
  }
}

for (j in 1:3) {
  for( i in 1:nrow(gff3)){ 
    
    if(( gff3$rpt_type[i]=='TERMINAL')==TRUE && is.na(gff3$rpt_type[i])==FALSE){ 
      print('ok')
      gff3<-gff3[-i,]
      
    }
  }
}

#deletes U86
gff3<-gff3[!is.na(gff3$gene),]
gff3<-gff3[!(gff3$gene=='U86'),]

#turns sam files into bam files. Merges all bam files into one big file
system('for i in *.sam; do samtools view -Sb $i > $i.bam; done ')

srr<-c()
read_id<-c()
read_seq<-c()
body_site<-c()
read_pos<-c()
sample_identifier<-c()

filelist<-list.files(pattern = "*.bam")

###takes forever but let it run against localblast later
##more reproducible
##read_ids have srr number in them so you can sort later
### run this whole shit against localblast 
for( i in 1:length(filelist)){ 
  progress(i)
  #print(filelist[i])
  tempbam<-scanBam(filelist[i])
  for(j in 1:length(tempbam[[1]]$qname)){
    #progress(j)
    ####setting thresholds for where to pull reads from. This avoids a bunch of repeat reads####
    if( tempbam[[1]]$pos[j] > 10000 && tempbam[[1]]$pos[j] < 155000) { 
      #print('ok')
  srr<-append(srr,strsplit(filelist[i],split = '.sam')[[1]][1])
  read_id<-append(read_id,tempbam[[1]]$qname[j])
  read_seq<-append(read_seq,tempbam[[1]]$seq[j])
  read_pos<-append(read_pos,tempbam[[1]]$pos[j])
    }
  }
}

alldata<-data.frame(srr,read_id,read_seq,read_pos)

fastaFile<-alldata[,c(2,3)]
fastaFile$read_id<-as.character(fastaFile$read_id)

system('mkdir fastas')
setwd('fastas')
for (i in 1:nrow(fastaFile)){ 
  print(i/nrow(fastaFile)*100)
  write.fasta(fastaFile$read_seq[i], fastaFile$read_id[i], paste0(fastaFile$read_id[i],'.fasta'), open = "w", nbchar = 200, as.string = TRUE)
  
  }

system("find . -path './*fasta' -prune -type f -exec cat {} + > big-fasta.txt")
system("blastn -query /Users/gerbix/Documents/vikas/Thesis_from_other_mac/hhv6b_rna/fastas/big-fasta.txt -db /Users/gerbix/Downloads/blast_viruses/all_virus.fasta -num_threads 8 -out blast_results.txt")
system("python /Users/gerbix/Documents/vikas/Thesis_from_other_mac/hhv6b_rna/testing/blast_hits.py blast_results.txt herpesvirus 6B")
####run local blast on that file### 
###run the verison of blast_hits for the local database on the imac### 
###make dataframe for any reads that survive only### 
##pull info like gene location etc for those positive files only###

blast_hits <- read.csv('blast_hits.csv')
#####anything >1 is 1 real hit

colnames(blast_hits)[1]<-'read_id'
colnames(blast_hits)[2]<-'hits'

blast_hits<-blast_hits[blast_hits$hits>1,]
blast_hits$sample<-NA
for(i in 1:nrow(blast_hits)){
  progress(i)
blast_hits$sample[i]<-strsplit(as.character(blast_hits$read_id[i]), "[.]")[[1]][1]
}

positivebams<-unique(blast_hits$sample)

setwd('/Users/gerbix/Documents/vikas/Thesis_from_other_mac/hhv6b_rna/')

bamslist<-c()
for (i in positivebams){ 
  tempbamname<-paste0(i, '.sam.bam')
  temp<-scanBam(tempbamname)
  bamslist<-append(bamslist, temp)
  }

total_list<-blast_hits
total_list$start = NA
total_list$end = NA

trimmed_alldata<-alldata




####really really damn slow (>1week)###
# for( i in 1:nrow(blast_hits)) { 
#   print(i)
#   if( blast_hits$read_id[i] %in% alldata) { } 
#     else { 
#       trimmed_alldata[i]<-NA
#       }
#     }







# #dataframe for SRR
# datadf<-as.data.frame(data$Run)
# colnames(datadf)[1]<-'SRR'
# datadf$body_site<-data$body_site
# datadf$sampleID<-data$submitted_subject_id
# datadf<-datadf[complete.cases(datadf), ]
# 
# #dataframe for gene info
# countdf<-data.frame(gff3$gene)
# colnames(countdf)[1]<-'gene'
# countdf[2]<-gff3$start
# colnames(countdf)[2]<-'start'
# countdf[3]<-gff3$end
# colnames(countdf)[3]<-'end'
# countdf<-na.omit(countdf)
# countdf[4]<-0
# colnames(countdf)[4]<-'count'
# 
# 
# count<-0
# df_length<-nrow(countdf) * nrow(datadf)
# totaldf<-data.frame(matrix(ncol = 7, nrow = df_length))
# for (i in 1:nrow(datadf)) { 
#   for (j in 1:(nrow(countdf))){
#     count<-count+1
#     totaldf[count,1]<-as.character(datadf$SRR[i])
#     colnames(totaldf)[1]<-'SRR'
#     totaldf[count,2]<-as.character(datadf$sampleID[i])
#     colnames(totaldf)[2]<-'sample_id'
#     totaldf[count,3]<-as.character(datadf$body_site[i])
#     colnames(totaldf)[3]<-'body_site'
#     totaldf[count,4]<-as.character(countdf$gene[j])
#     colnames(totaldf)[4]<-'gene'
#     totaldf[count,5]<-as.character(countdf$start[j])
#     colnames(totaldf)[5]<-'start'
#     totaldf[count,6]<-as.character(countdf$end[j])
#     colnames(totaldf)[6]<-'end'
#     totaldf[count,7]<-as.character(countdf$count[j])
#     colnames(totaldf)[7]<-'count'
#   }
#   print(100*i/nrow(countdf))
# }
# totaldf$count<-0
# totaldf$pos<-0
# ###
# totaldf$read_id<-0
# for(i in 1:nrow(totaldf)){ 
#   tempbamname<-paste0(as.character(totaldf$SRR[i]),'.sam.bam')
#   if( file.exists(tempbamname)) { 
#     tempbam<-scanBam(tempbamname)
#     for (j in 1:length(tempbam[[1]]$pos)) {
#       if(as.integer(totaldf$start[i])<tempbam[[1]]$pos[j] & tempbam[[1]]$pos[j]<as.integer(totaldf$end[i])) {
#         totaldf$count[i]<-(totaldf$count[i] + 1)
#         ###
#         #print(tempbam[[1]]$pos[j])
#         totaldf$pos[i]<-j
#         totaldf$read_id[i]<-tempbam[[1]]$qname[j]
#         ###
#         #  print(totaldf$SRR[i])
#       }
#     } 
#     print(100*i/nrow(totaldf))
#   }
#   else  {
#     print (paste0(tempbamname, ' does not exist'))
#   }
# }
# 
# #for parallelization
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload computer
# registerDoParallel(cl)
# 
# #trims duplicates from totaldf
# tempdf<-data.frame(matrix(ncol = 9, nrow = 0))
# srrlist<-datadf$SRR
# srrlist <- srrlist[!is.na(srrlist)]
# srrlist<-as.factor(srrlist)
# dflist<-c()
# x<-c()
# 
# makedfs<-function (SRRnumber, df){ 
#   df<-as.data.frame(df)
#   if(is.na(SRRnumber)==TRUE){ }
#   else {
#     templist<-c()
#     x<-c()
#     for (j in 1:nrow(df)){ 
#       # tempdf<-data.frame(matrix(ncol = 7, nrow = 0))
#       x<-append(x,df$SRR[j] == as.character(SRRnumber))
#       #tempdf
#     } 
#     
#     # print(x)
#     #print(i)
#     
#     tempdf<-NULL
#     print('ok')
#     print(class(df))
#     tempdf<-df[x,]
#     trimmedtempdf<-NULL
#     x<-c()
#     trimmedtempdf<-aggregate(tempdf$count~tempdf$gene,data=tempdf,FUN=sum)
#     colnames(trimmedtempdf)[1]<-'gene_name'
#     colnames(trimmedtempdf)[2]<-'count'
#     
#     trimmedtempdf$SRR<-as.character(SRRnumber)
#     
#     #adds in body site
#     bodysite_temp<-tempdf$body_site[SRRnumber==trimmedtempdf$SRR[1]]
#     bodysite_temp<-bodysite_temp[1]
#     trimmedtempdf$body_site<-bodysite_temp[1]
#     #trimmedtempdf$body_site<-tempdf$body_site[srrlist==trimmedtempdf$SRR[1]]
#     
#     #adds in sample id 
#     sample_id_temp<-tempdf$sample_id[SRRnumber==trimmedtempdf$SRR[1]]
#     sample_id_temp<-sample_id_temp[1]
#     trimmedtempdf$sample_id<-sample_id_temp[1]
#     #dflist[[i]] <-trimmedtempdf
#     assign(paste0(SRRnumber),trimmedtempdf)
#     print(as.character(SRRnumber))
#     #print(paste0(i,"/",length(srrlisttest)))
#     
#     #adds in pos 
#     pos_temp<-tempdf$pos[SRRnumber==trimmedtempdf$SRR[1]]
#     pos_temp<-pos_temp[1]
#     trimmedtempdf$pos<-pos_temp[1]
#     
#     #adds in read_id
#     read_id_temp<-tempdf$read_id[SRRnumber==trimmedtempdf$SRR[1]]
#     read_id_temp<-read_id_temp[1]
#     trimmedtempdf$read_id<-read_id_temp[1]
#     
#     return(trimmedtempdf)
#   }
# }
# 
# dflist[[i]]<-foreach(i=1:length(srrlist)) %dopar% {
#   #x<-c()
#   makedfs(as.character(srrlist[i]),totaldf)
#   #print(i)
# }
# dflist[sapply(dflist, is.null)] <- NULL
# final_df_trimmed<- do.call(rbind, dflist[[1]])
# 
# #stop cluster
# stopCluster(cl)
# 
# #removes all 0 count values 
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$count>0,]
# 
# #fills in start values 
# final_df_trimmed$start<-0
# for (i in 1:nrow(gff3)){ 
#   for (j in 1:nrow(final_df_trimmed)){ 
#     if (final_df_trimmed$gene_name[j]==gff3$gene[i]) {
#       final_df_trimmed$start[j]<-gff3$start[i]
#     }
#   }
# }
# 
# #fills in end values 
# final_df_trimmed$end<-0
# for (i in 1:nrow(gff3)){ 
#   for (j in 1:nrow(final_df_trimmed)){ 
#     if (final_df_trimmed$gene_name[j]==gff3$gene[i]) {
#       final_df_trimmed$end[j]<-gff3$end[i]
#     }
#   }
# }
# #####resetting position values 
# final_df_trimmed$pos<-0
# totaldf_0_trimmed<-totaldf[totaldf$count>0,]
# for (i in 1:nrow(final_df_trimmed)){ 
#   for(j in 1:nrow(totaldf_0_trimmed)){ 
#     if(final_df_trimmed$SRR[i]==totaldf_0_trimmed$SRR[j] & final_df_trimmed$start[i]==totaldf_0_trimmed$start[j]){ 
#       final_df_trimmed$pos[i]<-totaldf_0_trimmed$pos[j]
#     }
#   }
#   print(100*i/nrow(final_df_trimmed))
# }
# 
# ######  
# 
# 
# #helps keep x axis with all genes later on
# genedf<-as.data.frame(gff3$gene)
# colnames(genedf)[1]<-'gene_name'
# genedf[2]<-0
# colnames(genedf)[2]<-'count'
# genedf<-aggregate(genedf$count~genedf$gene_name,data=genedf,FUN=sum)
# colnames(genedf)[1]<-'gene_name'
# colnames(genedf)[2]<-'count'
# genedf[3]<-0
# colnames(genedf)[3]<-'SRR'
# genedf[4]<-0
# colnames(genedf)[4]<-'body_site'
# genedf[5]<-0
# colnames(genedf)[5]<-'sample_id'
# genedf[6]<-0
# colnames(genedf)[6]<-'start'
# genedf[7]<-0
# colnames(genedf)[7]<-'end'
# genedf[8]<-0
# colnames(genedf)[8]<-'pos'
# genedf[9]<-0
# colnames(genedf)[9]<-'read_id'
# final_df_trimmed<-rbind(final_df_trimmed,genedf)
# 
# 
# #gets rid of DR1, DR3, DR6
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='DR1',]
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='DR6',]
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='DR3',]
# 
# #remove B4
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='B4',]
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='B5',]
# 
# #remove U3 
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='U3',]
# 
# #manually setting U4 srr1343859 counts to 2 because the rest is repeat region. The rest of the U4s are repeat so deleting them.
# final_df_trimmed$count[final_df_trimmed$SRR=='SRR1343859']<-2
# u4row<-final_df_trimmed[final_df_trimmed$SRR=='SRR1343859',]
# final_df_trimmed<-final_df_trimmed[final_df_trimmed$gene_name!='U4' & final_df_trimmed$SRR!='SRR1343859',]
# final_df_trimmed<-rbind(final_df_trimmed,u4row)
# 
# 
# #calculating rpkm
# final_df_trimmed$rpkm<-0
# total_reads_per_sample<-0  
# 
# for (i in 1:nrow(final_df_trimmed)) { 
#   if( final_df_trimmed$count[i]>0){
#     total_reads_per_sample<-(data$MBases[data$Run==final_df_trimmed$SRR[i]])
#     print(total_reads_per_sample)
#     transcript_length<-(final_df_trimmed$end[i]-final_df_trimmed$start[i])
#     print(transcript_length)
#     final_df_trimmed$rpkm[i]<-final_df_trimmed$count[i]/(((transcript_length)/1000)*(total_reads_per_sample*1000000))
#   } 
#   
# }
# final_df_trimmed$rpkm<-as.numeric(final_df_trimmed$rpkm)
# #######
# 
# 
# #reorder x axis
# final_df_trimmed$gene_name<-as.character(final_df_trimmed$gene_name)
# ordered_final_df<-final_df_trimmed[mixedorder(final_df_trimmed$gene_name),]
# ordered_final_df$order<-1:nrow(ordered_final_df)
# 
# #color palette (still needs work) 
# coloring = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
#              "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
#              "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
#              "#8A7C64", "#599861","#0000ff","#ff0000")
# 
# #sets all genes in x axis withot plotting 0s
# for( i in 1:nrow(ordered_final_df)){ 
#   if(ordered_final_df$count[i]==0){ 
#     ordered_final_df$count[i]<-NA
#     ordered_final_df$body_site[i]<-NA
#     ordered_final_df$sample_id[i]<-NA
#     ordered_final_df$rpkm[i]<-NA
#     #ordered_final_df$pos[i]<-NA
#   }
# }
# 
# #plotting graph
# graph_1<-ggplot(ordered_final_df,aes(x=reorder(gene_name, order), y=rpkm)) + geom_point(aes(shape=sample_id, color=body_site), size=4)
# graph_1 +scale_colour_manual(values=coloring) +
#   coord_trans(y = "log10") + 
#   annotation_logticks(sides="l", scaled = FALSE) + 
#   theme(legend.key.size = unit(.05, "cm")) + 
#   theme_bw() +
#   scale_shape_manual(values=c(15, 16, 17, 18))+ 
#   theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1.1)) 
