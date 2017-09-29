#!/usr/bin/env Rscript
##this script is designed to calculate SSMD and SSMD score
#Update 20160714 
#change the read count cutoff from 10 to a variable which can be specified via command line
#the Score fomula from SSMD to abs(SSMD) 
#update 2016/07/31, add help function 
#update 2016/08/02 fix the bug where arugments are taken as string
#update 2016/08/09, sort(decreasing) according to the last column Score, only works for the pairwise comparison.
#update 2017/01/11, give the option to use the sum( control+sample) > threshold to filter out data
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
stop("4 arguments should be provided.
Usage: RScript --vanilla script.R count_table output_prefix count_treshold [0(default)|1].n", call.=FALSE)
} else if(length(args)==1){
 #default values for the later two arguments
 args[2]="output"
 args[3]= 10
 args[4]=0	
} else if( length(args)==2){
	args[3]=10
	args[4]=0
}else if(length(args)==3){
 args[4]=0
}

library(reshape2)
library(plyr)
#It reads the table with headers, please make sure you have the correct headers: RSL.guide,guide.set,Control,Sample1,Sample2...
tab<-read.csv(args[1],header=T)
outnormTab <- paste(args[2],"_normalized.csv",sep="")
#outMESMADSamp <- paste(args[2],"_MES_MAD_Sample.csv",sep="")
#outMESMADCtrl <- paste(args[2],"_MES_MAD_CTRL.csv",sep="")
outfilename<-paste(args[2],"_SSMD_Score.csv",sep="")
taboutfilename=paste(args[2],"_SSMD_Score.tsv",sep="")
## 1st get values for each samples,  2nd, sum up each columns 3rd, avg column values
sample_name<-colnames(tab)
if(as.numeric(args[4])> 0){
tabfilt=subset(tab,(Control+tab[,eval(sample_name[4])]) >=as.numeric(args[3]))
write.csv(file=paste(args[2],"filtDataBaseonSum.csv",sep="_"),tabfilt,quote=F,row.names=F)
} else{
	tabfilt=subset(tab,Control >=as.numeric(args[3]))
write.csv(file=paste(args[2],"filtDataBaseonCtrl.csv",sep="_"),tabfilt,quote=F,row.names=F)
}
value<-tabfilt[,3:ncol(tabfilt)]
tabrow<-colSums(value)
avgTotal <- sum(tabrow)/ncol(value)
norTab<-mapply(function(x, y) x / y, value*avgTotal ,tabrow)
norTabpre<-data.frame(cbind(tabfilt[,1:2],norTab),stringsAsFactors=TRUE)
sam_num<-length(sample_name)-3
sample <- tail(sample_name,sam_num)
write.csv(file=outnormTab, norTabpre, quote=F,row.names=F)
## get subset: to remove guides less than N normalized reads
#Non-target guides
norTabFiltNont <-norTabpre[grep("NONT",norTabpre$RSL.guide,perl=T,invert=F),]
enorTabFiltNont <- log((norTabFiltNont[,c(sample)]+1) /(norTabFiltNont[,3]+1))/log(2)
enorTabFiltNontF <- data.frame(cbind(norTabFiltNont[,1:2],enorTabFiltNont),stringsAsFactors=TRUE)
colnames(enorTabFiltNontF) <- c(head(sample_name,2),sample)

#Real guides
norTabFiltGuide <-norTabpre
#norTabFilt[grep("NONT",norTabpre$RSL.guide,perl=T,invert=T),]
enorTabFiltGuide<-log((norTabFiltGuide[,c(sample)]+1) /(norTabFiltGuide[,3]+1))/log(2)
enorTabFiltGuideF <- data.frame(cbind(norTabFiltGuide[,1:2],enorTabFiltGuide),stringsAsFactors=TRUE)
colnames(enorTabFiltGuideF) <- c(head(sample_name,2),sample)

###Section 1: for each treatments, get ESij
enorTabFiltGuideFM<-aggregate(enorTabFiltGuideF[, c(sample)], list(enorTabFiltGuideF$guide.set), median)
colnames(enorTabFiltGuideFM) <- c("guide.set",paste(sample,".MES",sep=""))
#Absolute Deviation
enorTabFiltGuideMerge<-join(enorTabFiltGuideF,enorTabFiltGuideFM)
enorTabFiltGuideAD<- data.frame(enorTabFiltGuideMerge[,1:2], 1.4826* abs(enorTabFiltGuideMerge[,c(sample)]- enorTabFiltGuideMerge[,c(paste(sample,".MES",sep=""))]),stringsAsFactors=TRUE)
#MAD of Samples
enorTabFiltGuideMAD<-aggregate(enorTabFiltGuideAD[, 3:ncol(enorTabFiltGuideAD)], list(enorTabFiltGuideAD$guide.set), median)
colnames(enorTabFiltGuideMAD) <- c("guide.set",paste(sample,".MAD",sep=""))
MESMADsamp<-join(enorTabFiltGuideFM,enorTabFiltGuideMAD)
#write.csv(file=outMESMADSamp,MESMADsamp,quote=F,row.names=F)


###Section 2: calculating from the non-targeting guies MEScon and MADcont , pooling all the non-targeting guides, get median of effective size
NonMedianPre<-melt(enorTabFiltNontF, id.vars="guide.set", measure.vars=c(sample))
#effective size and median (MEScon) for each control
NonMedian<-aggregate(NonMedianPre$value,list(NonMedianPre$variable),median)
colnames(NonMedian) <- c("variable","MedianC")
NonMedianPreN<-join(NonMedianPre,NonMedian)
#(1)Absolute Deviation for pooled non-targeting ES; (2) MADcontrol
NonMedianPreN[,"AD"] <- 1.4826 * abs(NonMedianPreN$value - NonMedianPreN$MedianC)
MADCON<-aggregate(NonMedianPreN$AD,list(NonMedianPreN$variable),median)
colnames(MADCON) <- c("variable","MAD")
MESMADcon<-join(NonMedian,MADCON)
#write.csv(file=outMESMADCtrl,MESMADcon,quote=F,row.names=F)

###Section 3: Merge tables from Samples and Control and Calculating SSMD and Score

# merge MES MAD from sample with MEScon and MADcon
SamCon<-data.frame(MESMADsamp,t(MESMADcon$MedianC),t(MESMADcon$MAD),stringsAsFactors=TRUE)
colnames(SamCon)<-c("guide.set",paste(sample,".MES",sep=""),paste(sample,".MAD",sep=""),paste(sample,".MESc",sep=""),paste(sample,".MADc",sep=""))
SSMD<-(SamCon[,2:(sam_num+1)]-SamCon[,(2*sam_num+2):(3*sam_num+1)])/sqrt(SamCon[,(sam_num+2):(2*sam_num+1)] ^ 2 + SamCon[,(sam_num*3+2):ncol(SamCon)]^2)
Score<-SamCon[,2:(sam_num+1)] * abs(SamCon[,2:(sam_num+1)]-SamCon[,(2*sam_num+2):(3*sam_num+1)])/sqrt(SamCon[,(sam_num+2):(2*sam_num+1)] ^ 2 + SamCon[,(sam_num*3+2):ncol(SamCon)]^2)
SSMDScore<-data.frame(SSMD,Score)
colnames(SSMDScore) <-c(paste(sample,".SSMD",sep=""),paste(sample,".Score",sep=""))
outtable<-data.frame(SamCon,SSMDScore,stringsAsFactors=TRUE)
deouttable <- outtable[order(outtable[ncol(outtable)],decreasing=F),]
deouttable["Rank"]=seq(1:nrow(deouttable))
write.table(file=taboutfilename,deouttable,quote=F,sep="\t",row.names=F)
