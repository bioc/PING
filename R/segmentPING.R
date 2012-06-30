segmentPING<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, maxLregion=1200, minLregion=80, step=NULL, width=NULL, dataType="H")
{
  #dataType="H"

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
  cat("step=",step,"| width=",width,"| dataType=",dataType,"\n")
  ## Perform the segmentation
  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=minReadsInRegion, jitter=jitter, maxLregion=maxLregion,minLregion=minLregion, step=step, width=width, package="PING")
  return(newSet)
}
