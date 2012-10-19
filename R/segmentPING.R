segmentPING<-function(data, dataC=NULL, map=NULL,
		minReads=2, minReadsInRegion=3, jitter=FALSE, maxLregion=1200, minLregion=80, step=NULL, width=NULL, dataType="H",
		islandDepth=NULL, min_cut=NULL, max_cut=NULL, chr=NULL)
{
  #get reads type: SE/PE
  if(is.null(islandDepth))
  {
	  cat("Performing segmentation for single end reads\n")
	  if(dataType=="TF")
	  {
		  if(is.null(step))
		  {
			  step<-20L
		  }
		  if(is.null(width))
		  {
			  width<-150L
		  }
	  }
	  if(dataType=="H")
	  {
		  if(is.null(step))
		  {
			  step<-2L
		  }
		  if(is.null(width))
		  {
			  width<-150L
		  }
	  }
	  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=minReadsInRegion, jitter=jitter, maxLregion=maxLregion,minLregion=minLregion, step=step, width=width, package="PING")
	  
	  
  }
  else #TODO: Find a safe way to determine SE/PE.
  {
	  cat("Performing segmentation for paired-end reads\n")
	  if(!is.character(chr))
		  stop("Argument chr should be a character\n")
	  if(length(chr)>1)
		  stop("Paired-end sequencing data segmentation does not support multiple chromosomes\n") 
	  if(!is.numeric(islandDepth))
		  stop("Argument islandDepth should be an integer\n")
	  if(!is.numeric(min_cut) | !is.numeric(max_cut))
		  stop("Arguments min_cut and max_cut should be integers, provided: ", class(min_cut),", ", class(max_cut), "\n")
	  
	  #with GRanges as input
	  data<-data[seqnames(data)==chr]
	  #Save xi for later
	  xi<-mean(width(data))
	  PE.RD<-IRanges(start=start(data), end=end(data))
	  PE.RD<-PE.RD[which(width(PE.RD)>min_cut | width(PE.RD)<max_cut)]		
	  PE.RD<-IRanges(start=start(PE.RD), end=end(PE.RD))
	  candidate_RD<-candidate.region(PE.RD, islandDepth, min_cut, max_cut)
	  
	  #Formating data for candidate.region
	  #idx <- ( data$P$"pos.-" >  data$P$"pos.+")
	  #PE_data <- data$P[idx, ]
	  #PE.RD.temp <- IRanges(start=PE_data$"pos.+" , end=PE_data$"pos.-")
	  #PE.RD <- PE.RD.temp[width(PE.RD.temp)<max_cut, ]	## Exclude too wide fragment ##
	  
	  #Get the candidate regions (IRanges)
	  #candidate_RD <- candidate.region(PE.RD, islandDepth, min_cut, max_cut)

	  #Getting parameters for segmentation
	  #PEMF.RD <- IRanges(start=data$yFm$"pos.+", end=data$yFm$"pos.+"+1)
	  #PEMR.RD <- IRanges(start=data$yRm$"pos.-"-1, end=data$yRm$"pos.-")
	  PEMF.RD<-IRanges()
	  PEMR.RD<-IRanges()

	  #TODO: Add the treatment in case !is.null(dataC)
	  N   <- length(PE.RD)
	  NFm <- length(PEMF.RD)
	  NRm <- length(PEMR.RD)
	  Nc <- NcFm <- NcRm <- as.integer(0)#integer(0)
	  
	  #Perform the segmentation
	  if(is.null(map)){
		  segmentation <- segChrRead(candidate_RD, PE.RD, PEMF.RD, PEMR.RD, PEC.RD=NULL, PECMF.RD=NULL, 
				  PECMR.RD=NULL, map.Start=NULL, map.End=NULL, chr=chr)
	  }else{ 
		  segmentation <- segChrRead(candidate_RD, PE.RD, PEMF.RD, PEMR.RD, PEC.RD=NULL, PECMF.RD=NULL, PECMR.RD=NULL, 
				  map.Start=start(map[chr]), map.End=end(map[chr]), chr=chr)
	  }
	  paraSW<-list(islandDepth=islandDepth, min_cut=min_cut, max_cut=max_cut, xi=xi)
	  newSet<-segReadsListPE(segmentation, paraSW=paraSW, N=N, NFm, NRm, Nc, NcFm, NcRm)
  }
  

  return(newSet)
}



###########################################
## Pre-process bam files for segmentation
## 
#prePING<-function(bamFile, outpath="./", save=TRUE)
#{
	#paras <- ScanBamParam(what=c("qname", "rname", "strand", "pos", "mapq", "qwidth"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
	#rawData <- scanBam(paste(bamFile), param=paras)[[1]]
	#rawData$rname <- as.character(rawData$rname)
	#rawData$qname=substring(as.character(rawData$qname),15)
	#rawData$qname=gsub(pattern="/3", replacement="", rawData$qname)
	#addon <- rawData$qwidth
	#addon[rawData$strand!="-"] <- 0  # to get end of "-" reads, seq length need to be added to their positions
	#retained <-  rawData$mapq > 10  # filter out reads with bad quality scores
	#temp1 <- data.frame(qname=rawData$qname[retained], 
			#rname=rawData$rname[retained], 
			#strand=rawData$strand[retained], 
			#pos=(rawData$pos+addon)[retained])
	#temp1 <- split(temp1[,c( "qname", "strand", "pos" )],temp1$rname) # split data into chromosomes
#	
	#PEreads<-vector('list', length(temp1))
	#names(PEreads)<-names(temp1)
	#for(chrs in names(temp1))
	#{
		#print(chrs)
	#	
		#dat <- temp1[[chrs]]
		#temp3=reshape(dat, timevar="strand", idvar="qname", direction="wide")
		#idx <- is.na(temp3[,c("pos.+","pos.-")])
		#reads <- list(P=temp3[!(idx[,1]|idx[,2]),], yFm=temp3[idx[,2],], yRm=temp3[idx[,1],])
		#if(isTRUE(save))
		#{
			#fName=paste(outpath,"/",bamFile,"_",chrs,".rda",sep="")
			#cat("Saving reads for chromosome", chrs, "in the file", fName,"\n")
			#save(reads, file=fName)
		#}
		#else
		#{
			#PEreads[[chrs]]<-reads
		#}
	#}
#	
	#return(invisible(PEreads))
#}

###
# bam2gr:
#   INPUT: A bam file with paired-end sequencingdata
#          An optional chr if only selected chr are needed
#   OUTPUT: A GRanges object that can be used in segmentPING
###
bam2gr<-function(bamFile, chr=NULL, PE=FALSE)
{
	paras <- ScanBamParam(what=c("qname", "rname", "strand", "pos", "mapq", "qwidth"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
	bga<-readBamGappedAlignments(bamFile, use.names=TRUE, param=paras)
	gr<-GRanges()

	if(!is.null(chr))
	  chrs<-chr
	else
	  chrs<-unique(as.character(runValue(seqnames(bga))))#
	for(i in 1:length(chrs))
	{ 
	  cat("Chromosome ", chrs[[i]], "\n")
	  bga2<-bga[seqnames(bga)==chrs[i]]
	  bga2<-bga2[elementMetadata(bga2)$mapq>10]#filter out elt with low quality scores
	  if(isTRUE(PE))
	  {
	    #change names
	    qname<-elementMetadata(bga2)$qname
	    qname<-substring(qname,15)
	    qname<-gsub(pattern="/3", replacement="", qname)
	    elementMetadata(bga2)$qname<-qname
	    #merge pairs
	    asdf<-as(bga2, "data.frame")
	    df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide")
	    df2<-df[,c("start.+", "end.-")]
	    rownames(df2)<-df[,"qname"]
	    #Split PE and SE
	    idx <- is.na(df2[,c("start.+","end.-")])
	    reads <- list(P=df2[!(idx[,1]|idx[,2]),], yFm=df2[idx[,2],], yRm=df2[idx[,1],])
	    #Build object
	    ir<-IRanges(start=reads$P$`start.+`, end=reads$P$`end.-`)
	    gr<-c(gr,GRanges(ranges=ir, seqnames=chrs[i])) #GR with PE (no SE)
	  }
	  else
	  {
	    ir<-IRanges(start=elementMetadata(bga2)$pos, width=elementMetadata(bga2)$qwidth)
	    gr<-c(gr,GRanges(ranges=ir, strand=elementMetadata(bga2)$strand, seqnames=chrs[i]))
	  }
	}
	return(gr)
}
	
