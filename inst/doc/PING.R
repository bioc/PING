### R code from vignette source 'PING.Rnw'

###################################################
### code chunk number 1: Loading-PING
###################################################
library(PING)


###################################################
### code chunk number 2: Reading-the-data
###################################################

path<- system.file("extdata",package="PING")
#Read the experiment : 
dataIP<-read.table(file.path(path,"GSM351492_R4_chr1.bed"),header=TRUE)
dataIP<-as(dataIP,"GRanges")


###################################################
### code chunk number 3: Cluster-initialization
###################################################
library(parallel)


###################################################
### code chunk number 4: Genome-segmentation
###################################################
seg<-segmentPING(dataIP,  minReads=NULL, maxLregion=1200,minLregion=80, jitter=TRUE)


###################################################
### code chunk number 5: PING-analysis
###################################################
ping<-PING(seg)


###################################################
### code chunk number 6: Post-process-PING-result
###################################################
PS=postPING(ping, seg)


###################################################
### code chunk number 7: makeRangedDataOutput
###################################################
rdBed<-makeRangedDataOutput(ping, type="bed")
rdFix<-makeRangedDataOutput(PS, type="fixed")


###################################################
### code chunk number 8: export-to-file (eval = FALSE)
###################################################
## library(rtracklayer)
## export(rdBed, file="ping.bed")
## export(rdFix, file="postPING.bed")


###################################################
### code chunk number 9: plotSummary
###################################################
plotSummary(PS, dataIP, "chr1", "gen", from=149000, to=153000)


###################################################
### code chunk number 10: custom-plot
###################################################
library(Gviz)
cTrack<-CoverageTrack(dataIP, "chr1", "gen")
rTrack<-RawReadsTrack(dataIP, "chr1", "gen")
nTrack<-NucleosomeTrack(PS, "chr1", "gen")


###################################################
### code chunk number 11: Gviz-tracks
###################################################
gTrack<-GenomeAxisTrack(add53=TRUE, add35=TRUE)
aTrack<-AnnotationTrack(range=IRanges(start=149500, end=151000), showFeatureId=TRUE, id="random annotation", col.title="orange", chr="chr1", gen="gen", name="custom")
plotTracks(trackList=c(gTrack, cTrack, aTrack, rTrack, nTrack), main="Custom plot", from=149000, to=153000 )


