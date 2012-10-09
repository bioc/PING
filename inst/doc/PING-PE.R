### R code from vignette source 'PING-PE.Rnw'

###################################################
### code chunk number 1: PING-PE.Rnw:24-25
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
### code chunk number 4: bam2gr-no-save
###################################################
gr<-bam2gr(bamFile=yeastBam)


###################################################
### code chunk number 5: Cluster-initialization
###################################################
library(parallel)


###################################################
### code chunk number 6: Genome-segmentation
###################################################
segPE<-segmentPING(gr,  chr="chrM", islandDepth=3, min_cut=50, max_cut=1000)


###################################################
### code chunk number 7: PING-analysis
###################################################
ping<-PING(segPE, PE=TRUE)


###################################################
### code chunk number 8: Post-process-PING-result
###################################################
{sigmaB2=3600; rho2=15; alpha2=98; beta2=200000}
PS=postPING(ping, segPE, rho2=rho2, alpha2=alpha2, beta2=beta2, sigmaB2=sigmaB2, PE=TRUE)


###################################################
### code chunk number 9: Display-result
###################################################
head(PS)


###################################################
### code chunk number 10: makeRangedDataOutput (eval = FALSE)
###################################################
## rdBed<-makeRangedDataOutput(PS, type="bed")
## library(rtracklayer)
## export(rdBed, "nucPrediction.bed")


###################################################
### code chunk number 11: plotSummary-PE
###################################################
plotSummary(PS, gr, chr="chrM", from=1000, to=4000, PE=TRUE)


