segmentPING<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, maxLregion=1200, minLregion=80, step=NULL, width=NULL, dataType="H", 
		islandDepth=NULL, min_cut=NULL, max_cut=NULL, chr=NULL)
{
  #get reads type: SE/PE
  if(class(data)=="GRanges")
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
	  if(!is.numeric(islandDepth))
		  stop("Argument islandDepth should be an integer\n")
	  if(!is.numeric(min_cut) | !is.numeric(max_cut))
		  stop("Arguments min_cut and max_cut should be integers, provided: ", class(min_cut),", ", class(max_cut), "\n")
	  
	  #Formating data for candidate.region
	  idx <- ( data$P$"pos.-" >  data$P$"pos.+")
	  PE_data <- data$P[idx, ]
	  PE.RD.temp <- IRanges(start=PE_data$"pos.+" , end=PE_data$"pos.-")
	  PE.RD <- PE.RD.temp[width(PE.RD.temp)<max_cut, ]	## Exclude too wide fragment ##
	  
	  #Get the candidate regions (IRanges)
	  candidate_RD <- candidate.region(PE.RD, islandDepth, min_cut, max_cut)

	  #Getting parameters for segmentation
	  PEMF.RD <- IRanges(start=data$yFm$"pos.+", end=data$yFm$"pos.+"+1)
	  PEMR.RD <- IRanges(start=data$yRm$"pos.-"-1, end=data$yRm$"pos.-")
	  N   <- length(PE.RD)
	  NFm <- length(PEMF.RD)
	  NRm <- length(PEMR.RD)
	  Nc <- NcFm <- NcRm <- integer(0)
	  
	  #Perform the segmentation
	  if(is.null(map)){
		  segmentation <- segChrRead(candidate_RD, PE.RD, PEMF.RD, PEMR.RD, PEC.RD=NULL, PECMF.RD=NULL, 
				  PECMR.RD=NULL, map.Start=NULL, map.End=NULL, chr=chr)
	  }else{ 
		  segmentation <- segChrRead(candidate_RD, PE.RD, PEMF.RD, PEMR.RD, PEC.RD=NULL, PECMF.RD=NULL, PECMR.RD=NULL, 
				  map.Start=start(map[chr]), map.End=end(map[chr]), chr=chr)
	  }
	  paraSW<-list(islandDepth=islandDepth, min_cut=min_cut, max_cut=max_cut)
	  newSet<-segReadsListPE(segmentation, paraSW=paraSW, N=N, NFm, NRm, Nc, NcFm, NcRm)
  }
  

  return(newSet)
}
