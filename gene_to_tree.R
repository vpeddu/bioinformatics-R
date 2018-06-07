library (ape)
library(seqinr)
library(Biostrings)
library(rtracklayer)

#input gene name, sequence file names 
#requires z29 genome 
geneanalyze<- function(gfffile, genename, fastafile )
{
 #make data frames from gff3 start/end to pull data from
  #bind z29 to fastas 
  
  importedgffs<-readGFF(gfffile, version=3,
                   columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                   raw_data=FALSE)
  importedgffcondensed <- importedgffs[,c('seqid','start','end','ID','Note')] 
  
  z29root<- readGFF("z29.gff3",version=3,
                    columns=NULL, tags=NULL, filter=NULL, nrows=-1,
                    raw_data=FALSE)
  z29root <- z29root[,c('seqid','start','end','ID','Note')] 
  
  rooted<-rbind(z29root,importedgffcondensed)
  rootedfiltered<-subset(rooted, rooted$Note==genename)
  rownumber<-nrow(rootedfiltered)
  sequences <- readDNAStringSet(fastafile)
  
  #makes new directory
  newdir <- paste0(genename)
  cwd=getwd()
  dir.create(newdir)
  setwd(newdir)
  
  #for csv file
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
geneinfotable<- data.frame(unlist(genenames),unlist(startpositionslist),unlist(endpositionlist),unlist(genewidth))
filename<- paste(genename , ".csv")
colnames(geneinfotable)[1]<-'accession number'
colnames(geneinfotable)[2]<-'start position'
colnames(geneinfotable)[3]<-'end position'
colnames(geneinfotable)[4]<-'width'
write.csv(geneinfotable, file=filename)

 # foldername<-c(genename)
 # write.csv(foldername, file=foldername)
 # write.table(foldername, file=foldername)

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
}

