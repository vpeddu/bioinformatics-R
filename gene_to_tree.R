#Takes sample gff3 + fasta, root gff3 + fasta, gene, and makes gene tree 
library (ape)
library(seqinr)
library(Biostrings)
library(rtracklayer)

#input gene name, sequence file names 
#args[1]=gff file name
#args[2]= gene name 
#args[3]= fasta file name 
#args[4]= root gff3
#args[5]=root fasta

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

#imports gff files
  importedgffs<-readGFF(args[1], version=3,
                        columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                        raw_data=FALSE)
  importedgffcondensed <- importedgffs[,c('seqid','start','end','ID','Note')] 
  
#imports z29 root gff3 and fasta
  root<- readGFF(args[4],version=3,
                    columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                    raw_data=FALSE)
  root <- root[,c('seqid','start','end','ID','Note')] 
  z29rootfasta <- readDNAStringSet(args[5])
  
  
#attaches z29 root to other files 
  rooted<-rbind(root,importedgffcondensed)
  rootedfiltered<-subset(rooted, rooted$Note==args[2])
  rownumber<-nrow(rootedfiltered)
  sequences <- readDNAStringSet(args[3])
  
#makes new directory for fasta files 
  newdir <- paste0(args[2])
  cwd=getwd()
  dir.create(newdir)
  setwd(newdir)
  
#blank lists for geneinfo table
  startpositionslist<-list()
  endpositionlist<-list()
  genenames<-list()
  genewidth<-list()
  dataframenames<-list()
  geneids<-list()
  
  for (i in 1:rownumber)
  {
    {tryCatch( {
      z<-subseq(sequences[i], start=rootedfiltered$start[i], end=rootedfiltered$end[i], width=NA)
      q<-names(sequences[i])
      q<-sub('\\..*', '.fasta', q)
      startpositionslist[[i]]<-rootedfiltered$start[i]
      endpositionlist[[i]]<-rootedfiltered$end[i]
      genenames[[i]]<-rootedfiltered$seqid[i]
      genewidth[i]<-unlist(endpositionlist[i])-unlist(startpositionslist[i])
      writeXStringSet(z, q, format='fasta')
    }
    , error=function(e) print (q)
    )
    }
  }
  #makes geneinfo table
  geneinfotable<- data.frame(unlist(genenames),unlist(startpositionslist),unlist(endpositionlist),unlist(genewidth))
  filename<- paste(args[2] , ".csv")
  colnames(geneinfotable)[1]<-'accession number'
  colnames(geneinfotable)[2]<-'start position'
  colnames(geneinfotable)[3]<-'end position'
  colnames(geneinfotable)[4]<-'width'
  write.csv(geneinfotable, file=filename)
  
  
#system commands to concatenate, align, and run fasttree on fasta files 
  system("
         lastdir=`ls -tr &lt;parentdir&gt; | tail -1`
         cd lastdir
         cat *.fasta >${PWD##*/}combined.fasta
         mafft ${PWD##*/}combined.fasta >  ${PWD##*/}aligned.fasta
         echo ${PWD##*/}aligned.fasta
         PATH=$PATH:~/Desktop/fasttree/
         FastTree -gtr -nt < ${PWD##*/}aligned.fasta > ${PWD##*/}tree.newick
         echo ${PWD##*/}tree.newick
         cd -
         ")
  setwd(cwd)
  
  quit(save="no")
