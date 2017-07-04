
#//////Function for visualization of the quality for aligned sequences and calculating the number of removed positional from the alignment under the given conditions of the trimming
#
# file_q - file name format fastq with initial consensus sequences from paired-end reads
# file_f - file name with aligned sequences of fasta format (sequences in a fasta format can be part of the original consensus sequences from paired-end reads)
# fastq and fsata fale must have identical sequence names
# qs - quality threshold for the trimming procedure
# dp - maximum portion of nucleotides per position with the quality smaller then qs after trimming

plot_quality <- function(file_q, file_f, qs, dp)
 {
  
  reads<-readFastq(file_q)
  pas_m<-readFasta(file_f)

  sr<-as.character(id(reads))
  sp<-as.character(id(pas_m))
  reads<-reads[which(sr %in% sp)]
  sr<-as.character(id(reads))
  sp<-as.character(id(pas_m))
  reads<-reads[match(sr,sp)]


  pas<-sread(pas_m)
  str<-as.character(pas)
  vstr<-strsplit(str, split="")
  
  p<-1
  pb<-txtProgressBar(min = 0, max=2*length(vstr)-1+length(vstr[[1]]), style = 3)
  
  q<-quality(reads)
  mq<-as(q, "matrix")
  mq_l<-list(mq[1,!is.na(mq[1,])])
  for(i in 2:nrow(mq))
   {
    mq_l<-c(mq_l, list(mq[i,!is.na(mq[i,])]))
    setTxtProgressBar(pb, p)
    p<-p+1
   }
  mq<-mq_l
 
  

  matrix_q<-matrix(, ncol=length(vstr[[1]]), nrow=length(vstr))

  for(i in 1:nrow(matrix_q))
   { 
    matrix_q[i,which(vstr[[i]]!="-" & vstr[[i]]!=".")]<-mq_l[[i]]
    setTxtProgressBar(pb, p)
    p<-p+1
   }

  q_row<-numeric(ncol(matrix_q))

  for(i in 1:ncol(matrix_q))
   {
    x<-matrix_q[,i][!is.na(matrix_q[,i])]
    q_row[i]<-length(x[x<qs])/length(x)
    setTxtProgressBar(pb, p)
    p<-p+1
   }

  plot(q_row,typ="n")
  lines(q_row)
  
  abline(h=dp , col="red", lwd=3, lty=2)

  n<-which(q_row>=dp) 
  
  cat("\n Alignment length = ")
  cat(length(vstr[[1]]))

  cat("\n Number of deleted items = ")
  cat(length(n))

  cat("\n Share of lost alignment items = ")
  cat(length(n)/length(vstr[[1]]))
  
  cat("\n")
  
 }


#//////Function for removed sequences with poor start and end quality
#
# file_q - file name format fastq with initial consensus sequences from paired-end reads
# file_f - file name with aligned sequences of fasta format (sequences in a fasta format can be part of the original consensus sequences from paired-end reads)
# file_trem - file name for saving the result of trimming
# qs - quality threshold for the trimming procedure(the minimum average quality of reading at the beginning or end of a sequence)
# begin_lt - the length of the beginning part of the sequence for which determine qs
# bend_lt - the length of the end part of the sequence for which determine qs	


begin_end_trimm <- function(file_q, file_f, file_trem, qs, begin_lt, end_lt)
 {

  reads<-readFastq(file_q)
  pas_m<-readFasta(file_f)

  sr<-as.character(id(reads))
  sp<-as.character(id(pas_m))
  reads<-reads[which(sr %in% sp)]
  sr<-as.character(id(reads))
  sp<-as.character(id(pas_m))
  reads<-reads[match(sr,sp)]


  pas<-sread(pas_m)
  str<-as.character(pas)
  vstr<-strsplit(str, split="")
  
  p<-1
  pb<-txtProgressBar(min = 0, max=2*length(vstr)-1+length(vstr[[1]]), style = 3)
  
  q<-quality(reads)
  mq<-as(q, "matrix")
  mq_l<-list(mq[1,!is.na(mq[1,])])
  for(i in 2:nrow(mq))
   {
    mq_l<-c(mq_l, list(mq[i,!is.na(mq[i,])]))
    setTxtProgressBar(pb, p)
    p<-p+1
   }
  mq<-mq_l

  matrix_q<-matrix(, ncol=length(vstr[[1]]), nrow=length(vstr))

  for(i in 1:nrow(matrix_q))
   { 
    matrix_q[i,which(vstr[[i]]!="-" & vstr[[i]]!=".")]<-mq_l[[i]]
    setTxtProgressBar(pb, p)
    p<-p+1
   }

  q_row<-rep(0, nrow(matrix_q))

  
  for(i in 1:nrow(matrix_q))
   { 
    if(begin_lt>0)
     {
      x<-matrix_q[i,c(1:begin_lt)]
      x<-x[which(!is.na(x))]
      if(length(x)>0)
       if(mean(x)<=qs) q_row[i]<-1
     } 
 
    if(end_lt>0)
     {
      y<-matrix_q[i,c((ncol(matrix_q)-end_lt):ncol(matrix_q))]
      y<-y[which(!is.na(x))]    
      if(length(y)>0)
       if(mean(y)<=qs) q_row[i]<-1
     }
    setTxtProgressBar(pb, p)
    p<-p+1
   } 
 
  n<-which(q_row==1) 
  n
  
  if(length(n)>0) pas_m<-pas_m[-n]
 
  writeFasta(pas_m, file_trem, mode="w")
  cat("\n")
 
 }

#//////Function to remove positions from the alignment under the given conditions of the trimming
#
# file_q - filename format fastq with initial consensus sequences from paired-end reads
# file_f - file name with aligned sequences of fasta format (sequences in a fasta format can be part of the original consensus sequences from paired-end reads)
# file_trem - filename to output results of trimming
# qs - quality threshold for the trimming procedure
# dp - maximum portion of nucleotides per position with the quality smaller then qs after trimming


total_trimm <- function(file_q, file_f, file_trem, pas_m, qs, dp)
 {
  
  reads<-readFastq(file_q)
  pas_m<-readFasta(file_f)

  sr<-as.character(id(reads))
  sp<-as.character(id(pas_m))
  reads<-reads[which(sr %in% sp)]
  sr<-as.character(id(reads))
  sp<-as.character(id(pas_m))
  reads<-reads[match(sr,sp)]


  pas<-sread(pas_m)
  str<-as.character(pas)
  vstr<-strsplit(str, split="")
  
  p<-1
  pb<-txtProgressBar(min = 0, max=2*length(vstr)-1+length(vstr[[1]])+length(vstr), style = 3)
  
  q<-quality(reads)
  mq<-as(q, "matrix")
  mq_l<-list(mq[1,!is.na(mq[1,])])
  for(i in 2:nrow(mq))
   {
    mq_l<-c(mq_l, list(mq[i,!is.na(mq[i,])]))
    setTxtProgressBar(pb, p)
    p<-p+1
   }
  mq<-mq_l
 
  

  matrix_q<-matrix(, ncol=length(vstr[[1]]), nrow=length(vstr))

  for(i in 1:nrow(matrix_q))
   { 
    matrix_q[i,which(vstr[[i]]!="-" & vstr[[i]]!=".")]<-mq_l[[i]]
    setTxtProgressBar(pb, p)
    p<-p+1
   }

  q_row<-numeric(ncol(matrix_q))

  for(i in 1:ncol(matrix_q))
   {
    x<-matrix_q[,i][!is.na(matrix_q[,i])]
    q_row[i]<-length(x[x<qs])/length(x)
    setTxtProgressBar(pb, p)
    p<-p+1
   }

  n<-which(q_row>=dp)

  id<-as.character(id(pas_m))
  id<-paste(">", id, sep="")

  pastrem<-character()

  for(i in 1:nrow(matrix_q))
   {
    pastrem<-c(pastrem,paste(vstr[[i]][-n],sep="", collapse=""))
    setTxtProgressBar(pb, p)
    p<-p+1
   }
  
  lines<-character(length(id)+length(pastrem))
  
  lines[seq(1, 2*length(id), 2)]<-id
  lines[seq(0, 2*length(pastrem), 2)[-1]]<-pastrem

  writeLines(lines,file_trem, sep = "\n")
  
  cat("\n")
  
 }

#Example of using functions

#connecting the necessary library

library(ShortRead)

#set of working directory
setwd("/home/my_data/")

#visualization of reading quality
plot_quality(file_q="example_data.fastq", file_f="example_data.fasta", qs=20, dp=0.1)

#removed sequences with poor start and end quality
begin_end_trimm(file_q="example_data.fastq", file_f="example_data.fasta", file_trem="intermediate_result.fasta", qs=37, begin_lt=60, end_lt=0)

#visualization of reading quality after removing the sequences
plot_quality(file_q="example_data.fastq", file_f="intermediate_result.fasta", qs=20, dp=0.1)

#removed positional from the alignment under the given conditions of the trimming
total_trimm(file_q="example_data.fastq", file_f="intermediate_result.fasta", file_trem="example_trem.fasta", qs=20, dp=0.1)
