#Given fasta and gff3 file, takes 10 biggest genes and creates csv with accession, gene, and start/end pos

library (ape)
library(seqinr)
library(Biostrings)
library(rtracklayer)


#import gff3 file to find where genes are 
#sequence.gff3 is a reference genome 
generef<-readGFF("sequence.gff3", version=3,
                 columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                 raw_data=FALSE)
#removes cds0 from data frame 
generef<-generef[-c(1), ]

#finds sizes of each gene
genesizelist<-list()
for (i in 1:nrow(generef)) {
  x<-generef$end[i]-generef$start[i]
  genesizelist[i]<-x
}
rm(i)
rm(x)

#creates data frame for gene size with CDS
geneID<-subset.data.frame(generef,select=c("ID"))
geneID$new<-1:nrow(geneID)
genesizelistdf <- data.frame(matrix(unlist(genesizelist), nrow=125, byrow=T),stringsAsFactors=FALSE)
genesizelistdf$new<- 1:nrow(genesizelistdf)
genewithsize<-merge.data.frame(geneID, genesizelistdf,by="new")
genewithsize$new<-NULL
colnames(genewithsize)[2] <- "size"

#splits genewithsize into chunks of 10 rows 
splitgenelist<- split(genewithsize,rep(1:13,each=10))

#creates individual data frames for each name in list
for (i in 1:13)
{
  r<-as.data.frame(splitgenelist[i])
  colnames(r)[1] <- "gene"
  colnames(r)[2] <- "size"
  assign(paste0("split",i), r)
  
}
rm(r)

#finds biggest gene in each list
dflist<-c('split1','split2','split3','split4','split5','split6','split7','split8','split9','split10','split11','split12','split13')
largestgenes<- c()
for (i in 1:13) {
  z<-get(dflist[i])
  q<-order(z$size, decreasing = TRUE)
  f<-q[1]
  t<-z$gene[f]
  largestgenes<- append(largestgenes, t)  
}

#pulls gene name from notes column
notelist<-list() 

for( i in 1:13) {
  q<- generef$ID==largestgenes[i]
  z<-generef$Note[q]
  z<-as.character(z)
  if (is.na(z))
  {
    notelist[i]<-largestgenes[i]
    
  }
  else {
    notelist[i]<-z
  }
}
notelisttrimmed<-gsub(',', '', notelist)
#notelisttrimmed<-as.list(notelisttrimmed)
write.table(notelisttrimmed, file="notelist.txt", quote = FALSE, row.names = FALSE,col.names = FALSE)

#start and end of biggest genes (to be used in for loop)
genestart= generef$start[which(generef$ID==largestgenes[1])]
geneend= generef$end[which(generef$ID==largestgenes[1])]

#imports .gff3 with all 130 organisms 
allorganisms<-readGFF("allsequences.gff3", version=3,
                      columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                      raw_data=FALSE)
#trimming columns
allorganisms <- allorganisms[,c('seqid','start','end','ID','Note')]

#adds z29 strain to root tree in future steps 
z29root<- readGFF("z29.gff3",version=3,
                  columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                  raw_data=FALSE)
#trimming columns
z29root <- z29root[,c('seqid','start','end','ID','Note')] 

#combines data frames
allorganismsrooted<-rbind(z29root,allorganisms)

allorganismsfiltered<-subset(allorganismsrooted,
                             
for (i in 1:13) {
  if ( startsWith(as.character(notelist[i]), 'U')) {
    allorganismsrooted$Note==notelist[i]
    } 
  if (startsWith(as.character(notelist[i]), 'cds')) {
      allorganismsrooted$ID==notelist[i]
                                 
    } 
  }  
)

#loop to make individual data frames for each CDS 
#separates cds from U 
for (i in 1:13) {
  if (startsWith(as.character(notelist[i]),'U'))
  {
    r<-subset(allorganismsrooted,allorganismsrooted$Note==notelist[i])
    assign(paste0(notelist[i]), r)  
    
  }
  if (startsWith(as.character(notelist[i]),'cds'))
  {
    r<-subset(allorganismsrooted,allorganismsrooted$ID==notelist[i])
    assign(paste0(notelist[i]), r)  
    
  }
}
#makes DNAstrings object from fasta and includes root
hhv6sequences <- readDNAStringSet("hhv6sequences.fasta")
z29sequence <- readDNAStringSet("z29sequence.fasta")
hhv6strings<- c(z29sequence,hhv6sequences)
#makes list of data frames
dflist=(notelist)

#makes fasta files and puts them in individual folders for each cds
#will print any error files to console

startpositionslist<-list()
endpositionlist<-list()
genenames<-list()
genewidth<-list()
dataframenames<-list()
geneids<-list()

for (df in dflist){
  cwd <- getwd()        
  newdir <- paste0(df)
  dir.create(newdir)      
  setwd(newdir) 
  for(i in 1:131)
  {tryCatch( {
    z<-subseq(hhv6strings[i], start=get(df)$start[i], end=get(df)$end[i], width=NA)
    q<-names(hhv6strings[i])
    q<-sub('\\..*', '.fasta', q)
    startpositionslist[[i]]<-get(df)$start[i]
    endpositionlist[[i]]<-get(df)$end[i]
    genenames[[i]]<-get(df)$seqid[i]
    genewidth[i]<-unlist(endpositionlist[i])-unlist(startpositionslist[i])
    writeXStringSet(z, q, format='fasta')  
  }
  , error=function(e) print (q)
  
  )
    
  }
  geneids[[i]]<-df
  r<- data.frame(unlist(geneids),unlist(genenames),unlist(startpositionslist),unlist(endpositionlist),unlist(genewidth))
  z<-paste0("dataframe",df)
  assign(z, r)
  dataframenames[df]<-z
  setwd(cwd)
  
}

dataframenames<-as.character(dataframenames)
geneinfotable<-rbind(dataframeB1,dataframecds120,dataframecds24,dataframecds9,dataframeU28,dataframeU35,dataframeU45,dataframeU52,dataframeU59,dataframeU75,dataframeU8,dataframeU94)
colnames(geneinfotable)[1]<-'gene'
colnames(geneinfotable)[2]<-'accession number'
colnames(geneinfotable)[3]<-'start position'
colnames(geneinfotable)[4]<-'end position'
colnames(geneinfotable)[5]<-'width'
rm(dataframeB1,dataframecds120,dataframecds24,dataframecds9,dataframeU28,dataframeU35,dataframeU45,dataframeU52,dataframeU59,dataframeU75,dataframeU8,dataframeU94)
write.csv(geneinfotable, file="geneinfo.csv")




setwd(cwd)
