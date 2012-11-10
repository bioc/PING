### R code from vignette source 'PING-PE.Rnw'

###################################################
### code chunk number 1: PING-PE.Rnw:26-27
###################################################
options(continue=" ")


###################################################
### code chunk number 2: Loading-PING
###################################################
library(PING)


###################################################
### code chunk number 3: Read-data
###################################################
yeastBam<- system.file("extdata/yeastChrI_M.bam",package="PING")


###################################################
### code chunk number 4: bam2gr
###################################################
gr<-bam2gr(bamFile=yeastBam, PE=TRUE)


###################################################
### code chunk number 5: Cluster-initialization
###################################################
library(parallel)


###################################################
### code chunk number 6: subset-GR
###################################################
grM<-gr[seqnames(gr)=="chrM"]


###################################################
### code chunk number 7: Genome-segmentation
###################################################
segPE<-segmentPING(grM, PE=TRUE)


###################################################
### code chunk number 8: PING-analysis
###################################################
ping<-PING(segPE, nCores=2)


###################################################
### code chunk number 9: Post-process-PING-result
###################################################
{sigmaB2=3600; rho2=15; alpha2=98; beta2=200000}
PS=postPING(ping, segPE, rho2=rho2, alpha2=alpha2, beta2=beta2, sigmaB2=sigmaB2)


###################################################
### code chunk number 10: makeRangedDataOutput (eval = FALSE)
###################################################
## rdBed<-makeRangedDataOutput(PS, type="bed")
## library(rtracklayer)
## export(rdBed, "nucPrediction.bed")


###################################################
### code chunk number 11: plotSummary-PE
###################################################
plotSummary(PS, ping,  grM, chr="chrM", from=1000, to=4000)


