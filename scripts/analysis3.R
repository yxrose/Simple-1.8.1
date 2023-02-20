setwd("./output")
###setwd("../output")
library("ggplot2")
library("ggrepel")
library("reshape2")
###library(qmrparser)
###library("utf8latex")
###JEN inputting the line information and setting filenames
args <- commandArgs(trailingOnly=TRUE)
line <- args[1]
###line<-"PH"
plotfilename <- paste(line, ".allSNPs.txt", sep="")
candidatefilename <- paste(line, ".candidates.txt", sep="")


mntn=read.delim(plotfilename, header=T,as.is=T)#In my previous version of R, but NOT on the bioross server, R adds an extra NA column
cnds=read.delim(candidatefilename, header =T, sep = "\t",as.is=T)

tf<-apply(mntn[,10:13],1,function(x) any(x=="."))
mntn<-mntn[!tf,]
mntn$mut.ref<-as.numeric(mntn$mut.ref)
mntn$mut.alt<-as.numeric(mntn$mut.alt)
mntn$wt.ref<-as.numeric(mntn$wt.ref)
mntn$wt.alt<-as.numeric(mntn$wt.alt)

mntn= mntn[!(mntn$chr=="ChrSy" | mntn$chr=="ChrUn"),]
cnds=cnds[!(cnds$chr=="ChrSy" | cnds$chr=="ChrUn"),]
mntn[,1]<-as.numeric(gsub("Chr","",mntn[,1]))
cnds[,1]<-as.numeric(gsub("Chr","",cnds[,1]))
mntn<-mntn[order(mntn[,1],mntn[,2]),]
cnds<-cnds[order(cnds[,1],cnds[,2]),]
wt=mntn$wt.ref/(mntn$wt.ref+mntn$wt.alt)
mut=mntn$mut.ref/(mntn$mut.ref+mntn$mut.alt)
ratio=wt-mut

allSNPs=data.frame(mntn, ratio)
#head(allSNPs)
write.table(allSNPs, paste0(line, ".allSNPs.txt"), sep="\t", row.names=F, quote=F)


tbl1=data.frame(At_num=mntn$At_num, gene=mntn$gene, chr=mntn$chr, pos=mntn$pos, mut.ref=mntn$mut.ref, mut.alt=mntn$mut.alt, wt.ref=mntn$wt.ref, wt.alt=mntn$wt.alt, mut.ratio=mut, wt.ratio=wt, ratio)
tbl1=tbl1[complete.cases(tbl1),]

#tbl1.cands=tbl1[(tbl1[,3] %in% cnds[,1]) & (tbl1[,4] %in% cnds[,2]),]
breaks=seq(0, max(tbl1$pos), round(max(tbl1$pos)/3, digits=-7))

#########################################################################
##separated chromosomes original data; after filtering (tbl2)-LOESS fitted
##removing genes with ratio below 0.1
#########################################################################
tbl2=tbl1[(tbl1[,11]>0.1),]
t2_s=split(tbl2, tbl2$chr)
lll=lapply(t2_s, function(x) {loess(x$ratio~x$pos, degree=2, span=0.3)})
mmm=lapply(lll, '[[', 'fitted')
fitted=Reduce(c, mmm)
tbl3=data.frame(tbl2, fitted)
tbl3.cands=tbl3[(tbl3[,3] %in% cnds[,1]) & (tbl3[,4] %in% cnds[,2]),]

#making the x-axes labels less messy
fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     l <- gsub("0e\\+00","0",l)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e\\+","e",l)
     l <- gsub("e", "%*%10^", l)
     l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
     l <- gsub("\\.0", "", l)
     #l <- gsub("\\'[\\.0]'", "", l)
     # return this as an expression
     parse(text=l)
}

##getting the loess fitted plot
##x11()
xxx=ggplot(tbl3, aes(pos, fitted)) + geom_point(aes(color=chr),size=0.3)+ facet_grid (.~ chr, scales = "free_x", space = "free_x")+geom_point(data=tbl3.cands, aes(x=pos, y=fitted), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=fitted, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks, labels=fancy_scientific)+labs(x="position", y="ratio")

####JEN changed file name
Rplot_loess1_file <- paste(line, ".Rplot.loess.1.pdf", sep="")
ggsave(Rplot_loess1_file, plot=xxx)

#########################################################################
##separated chromosomes original data; after filtering (tbl2)-LOESS fitted
##removing genes with ratio below 0.3
#########################################################################
tbl2=tbl1[(tbl1[,11]>0.3),]
t2_s=split(tbl2, tbl2$chr)
lll=lapply(t2_s, function(x) {loess(x$ratio~x$pos, degree=2, span=0.3)})
mmm=lapply(lll, '[[', 'fitted')
fitted=Reduce(c, mmm)
tbl3=data.frame(tbl2, fitted)
tbl3.cands=tbl3[(tbl3[,3] %in% cnds[,1]) & (tbl3[,4] %in% cnds[,2]),]

###making the x-axes labels less messy
fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     l <- gsub("0e\\+00","0",l)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e\\+","e",l)
     l <- gsub("e", "%*%10^", l)
     l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
     l <- gsub("\\.0", "", l)
     #l <- gsub("\\'[\\.0]'", "", l)
     # return this as an expression
     parse(text=l)
}

###getting the loess fitted plot
###x11()
yyy=ggplot(tbl3, aes(pos, fitted)) + geom_point(aes(color=chr),size=0.3)+ facet_grid (.~ chr, scales = "free_x", space = "free_x")+geom_point(data=tbl3.cands, aes(x=pos, y=fitted), shape=5)+geom_text_repel(data=tbl3.cands, aes(x=pos, y=fitted, label=gene), size=3)+theme(legend.position="none")+scale_x_continuous(breaks=breaks, labels=fancy_scientific)+labs(x="position", y="ratio")

### JEN changed file name
Rplot_loess3_file <- paste(line, ".Rplot.loess.3.pdf", sep="")
ggsave(Rplot_loess3_file, plot=yyy)

###and getting the allele frequency
tbl3.m=melt(tbl3, id.vars=c('At_num', 'gene', 'chr', 'pos', 'mut.ref', 'mut.alt', 'wt.ref', 'wt.alt', 'ratio', 'fitted'))

tbl3.cands.m=melt(tbl3.cands, id.vars=c('At_num', 'gene', 'chr', 'pos', 'mut.ref', 'mut.alt', 'wt.ref', 'wt.alt', 'ratio', 'fitted'))


###x11()
zzz=ggplot(tbl3.m, aes(pos, value)) + geom_point(aes(color=variable),size=0.3)+ facet_grid (.~ chr, scales = "free_x", space = "free_x")+geom_point(data=tbl3.cands.m, aes(x=pos, y=value), shape=5)+geom_text_repel(data=tbl3.cands.m, aes(x=pos, y=value, label=gene), size=3)+theme(legend.title = element_text(size = 0))+scale_x_continuous(breaks=breaks, labels=fancy_scientific)+labs(x="position", y="allele frequency")

###JEN changed file name
Rplot_allele_file <- paste(line, ".Rplot_allele.pdf", sep="")
ggsave(Rplot_allele_file, plot=zzz)
