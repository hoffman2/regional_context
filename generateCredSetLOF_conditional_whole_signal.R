library(DBI)
library(odbc)
options(gtx.dbConnection = dbConnect(odbc(), dsn = 'impaladsn'))
library(data.table)
#Location of credset function
library(dplyr)
library(parallel)
source("dependencies/abf.R")
source("dependencies/coloc.R")
library(stringr)
library(sparklyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  myInfile <- args[1]
  myOutfile <- args[2]
} else {
  stop("<inFile> <outfile>")
}


#Read in file
myVariants<-fread(myInfile)

generateCredSet_conditional <- function(searchVariants,distFromVariant,outFile){

  calcPosterior <- function(myDT){
  newDT <- myDT
  newDT <- newDT[!is.na(beta_cond)&!is.na(se_cond),]
  newDT[,traitABF:=with(newDT,abf.Wakefield(beta=beta_cond,se=se_cond,priorsd=1,log=TRUE))]
  newDT[,posteriorTrait:=norm1(newDT[["traitABF"]], log = TRUE)]
  newDT<-newDT[order(-posteriorTrait)]
  newDT[,cumSumTrait:=cumsum(x=newDT[["posteriorTrait"]])]
  newDT[,credSet:=credset(newDT[["posteriorTrait"]])]
  return(newDT)
}
  
  
  #Start by  making a copy of the table with variant info. We will create additional columns for new minimum pvalue variant info and cred set information.
  myTable <- searchVariants
  myTable[,c("SIGNAL","MINP_CHROM","MINP_POS","MINP_PVALUE","MINP_lnABF","MINP_POSTERIOR","MINP_MAF","TARGET_PVALUE","TARGET_MAF","TARGET_lnABF","TARGET_POSTERIOR","TARGET_CREDSET","TARGET_RANK_PVALUE","TARGET_RANK_ABF","CREDSET_SIZE"):=NA]
  myTable[,MINP_PVALUE:=as.double(MINP_PVALUE)][,TARGET_PVALUE:=as.double(TARGET_PVALUE)][,TARGET_MAF:=as.double(TARGET_MAF)][,MINP_CHROM:=as.character(MINP_CHROM)][,MINP_POS:=as.integer(MINP_POS)][,TARGET_POSTERIOR:=as.double(TARGET_POSTERIOR)][,TARGET_lnABF:=as.double(TARGET_lnABF)][,MINP_POSTERIOR:=as.double(MINP_POSTERIOR)][,MINP_lnABF:=as.double(MINP_lnABF)][,TARGET_CREDSET:=as.character(TARGET_CREDSET)][,TARGET_RANK_PVALUE:=as.double(TARGET_RANK_PVALUE)][,TARGET_RANK_ABF:=as.double(TARGET_RANK_ABF)][,CREDSET_SIZE:=as.integer(CREDSET_SIZE)][,MINP_MAF:=as.double(MINP_MAF)][,SIGNAL:=as.integer(SIGNAL)][,uniqueComb:=paste(chrom,pos,ref,alt,VARIANT_TYPE,analysis,signal,sep=".")]
  #Need to avoid issue of duplicated 500K analysis since we are linking by phenotype. To skirt this I am appending a unique identifier to "analysis" column when duplicated. Shouldnt be necessary with script has ben redesigned
  myTable[,uniqueComb:=make.unique(myTable[["uniqueComb"]],sep = "_")]
  numbRows=nrow(myTable)
  myCount=0
  #Create progress bar for tracking
  progressBar=txtProgressBar(min = 0, max = nrow(myTable), initial = 0)
  for(i in myTable[["uniqueComb"]]){
    myChrom=as.character(myTable[uniqueComb %in% i,][["chrom"]])
    myPos=as.integer(myTable[uniqueComb %in% i,][["pos"]])
    myRef=as.character(myTable[uniqueComb %in% i,][["ref"]])
    myAlt=as.character(myTable[uniqueComb %in% i,][["alt"]])
    myType=as.character(myTable[uniqueComb %in% i,][["VARIANT_TYPE"]])
    myPhenotype=as.character(myTable[uniqueComb %in% i,][["analysis"]])
    mySignal=as.integer(myTable[uniqueComb %in% i,][["signal"]])
    myUnique=as.character(myTable[uniqueComb %in% i,][["uniqueComb"]])
    part1="SELECT * FROM gene_23andme_l2_shared.gwas_results_cond_pub where analysis="
    #mySQL_query=part1 + information taken from each variant in annotation file
    #Originally had a freq filter here, but decided to remove so that we can make sure all variants can be looked up before cred set
    #Not implementing the bp window flag, just pulling the entire signal that the variant belongs to
    mySQL_query=paste0(part1,"'",myPhenotype,"' AND chrom=","'",myChrom,"'"," AND signal=",mySignal)
    #Save the results of the query as a data.table for downstream analysis
    myQueryResults <- as.data.table(dbGetQuery(getOption("gtx.dbConnection"), mySQL_query))[order(pval_cond),]
    #We dont have frequency information in conditional results
    #myQueryResults[,MAF:=ifelse(freq>0.5,as.double(1-freq),as.double(freq))]
    #Add in 500K variant summary stats, but need to make sure it is present
    if(nrow(myQueryResults[pos==myPos&ref==myRef&alt==myAlt,])==0){
      myTable[uniqueComb %in% i,][["TARGET_PVALUE"]]="NA"
      #myTable[uniqueComb %in% i,][["TARGET_MAF"]]="NA"
    } else{
      myTable[uniqueComb %in% i,][["TARGET_PVALUE"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["pval_cond"]][1]
      #myTable[uniqueComb %in% i,][["TARGET_MAF"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["MAF"]][1]
      myTable[uniqueComb %in% i,][["SIGNAL"]]=myQueryResults[pos == myPos,][["signal"]][1]
      }
    #Set allele frequency frequency after capturing variant in 500K so that we can perform cred set analysis
    #myQueryResults<-myQueryResults[MAF>=0.0001,]
    myQueryResults <-calcPosterior(myQueryResults)
    myQueryResults[,credSet:=as.character(credSet)]
    #Put in rank of variants based on pvalue
    myQueryResults[,variantRank:=frank(myQueryResults[["pval_cond"]])]
    #Put in rank based on ABF
    myQueryResults[,variantRank_cred:=frank(-myQueryResults[["traitABF"]])]
    #If variant is multi allelic I am taking variant with lowest pvalue hence the addition of [1] at end of each query.This shouldnt be an issue now that ref/alt added in to unique name
    myTable[uniqueComb %in% i,][["MINP_CHROM"]]=myQueryResults[which.min(pval_cond),][["chrom"]][1]
    myTable[uniqueComb %in% i,][["MINP_POS"]]=myQueryResults[which.min(pval_cond),][["pos"]][1]
    myTable[uniqueComb %in% i,][["MINP_PVALUE"]]=myQueryResults[which.min(pval_cond),][["pval_cond"]][1]
    myTable[uniqueComb %in% i,][["MINP_POSTERIOR"]]=myQueryResults[which.min(pval_cond),][["posteriorTrait"]][1]
    myTable[uniqueComb %in% i,][["MINP_lnABF"]]=myQueryResults[which.min(pval_cond),][["traitABF"]][1]
    #myTable[uniqueComb %in% i,][["MINP_MAF"]]=myQueryResults[which.min(pval_cond),][["MAF"]][1]
    myTable[uniqueComb %in% i,][["CREDSET_SIZE"]]=length(myQueryResults[credSet=="TRUE",credSet])
    #Check if LOF exists in Toby's database, if not, dont add in additional info about LOF variant from GWAS.
    if(nrow(myQueryResults[pos==myPos&ref==myRef&alt==myAlt,])==0){
      myTable[uniqueComb %in% i,][["TARGET_CREDSET"]]="TARGET_Not_Found"
        } else{
          #If LOF exists in GWAS, then run perform following query
          myTable[uniqueComb %in% i,][["TARGET_RANK_PVALUE"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["variantRank"]][1]
          myTable[uniqueComb %in% i,][["TARGET_RANK_ABF"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["variantRank_cred"]][1]
          myTable[uniqueComb %in% i,][["TARGET_POSTERIOR"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["posteriorTrait"]][1]
          myTable[uniqueComb %in% i,][["TARGET_lnABF"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["traitABF"]][1]
          myTable[uniqueComb %in% i,][["TARGET_CREDSET"]]=myQueryResults[pos==myPos&ref==myRef&alt==myAlt,][["credSet"]][1]
        }
    myCount=myCount+1
    if(myCount%%100==0) print(paste0("Finished with number:",myCount," of ",nrow(myTable)))}
  
  myTable[,TARGET_IS_MOST_SIG:=ifelse(MINP_POS==pos,"YES","NO")]
  write.table(myTable,row.names = F,quote=F,sep="\t",file=paste0(outFile,".txt"))
  return(myTable)
  print("Im Done!")}
  
#Export as text files
system.time(generateCredSet_conditional(myVariants,1e6,myOutfile))

##To run, locally
##Create time stamp
#myTime=gsub(x=Sys.time(),"-|:|[[:space:]]","_")
#system(paste0("rm -r my_temp_",myTime))
#system(paste0("mkdir my_temp_",myTime))
##Split into chunks
#lapply(1:length(myList),function(x) write.table(myList[[x]],paste0("my_temp_",myTime,"/target_variant_analysis.",names(myList)[x],".txt"),row.names=F))
##Run the script to perform regional context analysis
##Test in background
#for(i in 1:8){system(paste0("Rscript generateCredSetLOF_conditional_whole_signal.R ","my_temp_",myTime,"/target_variant_analysis.",i,".txt ","my_temp_",myTime,"/target_regional_results.set",i," &"))}  

