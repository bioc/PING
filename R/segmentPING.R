segmentPING<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, maxLregion=1200,minLregion=80)
{
  dataType="H"

  if(dataType=="TF")
  {
    step<-paraSegTF$step
    width<-paraSegTF$width
  }
  if(dataType=="H")
  {
    step<-paraSegH$step
    width<-paraSegH$width
  }
  ## Perform the segmentation
  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=minReadsInRegion, jitter=jitter, maxLregion=maxLregion,minLregion=minLregion,
		  step=step, width=width, package="PING")
  return(newSet)
}
