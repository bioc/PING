paraSegH<-list(step=2L,width=150L)
paraSegTF<-list(step=20L,width=150L)


setParaSeg<-function(step=10,width=75,dataType="H")
{
  if(dataType!="TF" & dataType!="H")
  {
    stop("Object 'dataType' must be either 'TF' or 'H'")
  }
  if(!is.integer(width) & width<=0)
  {
    stop("'width' must be a positive integer")
  }
  if(!is.integer(step) & step<=0)
  {
    stop("'step' must be a positive integer")
  }
  
  if(dataType=="TF")
  {
    unlockBinding("paraSegTF", environment(PING))
    assign("paraSegTF",list(step=step,width=width),envir=environment(PING))
    assign("paraSegTF",list(step=step,width=width),envir=.GlobalEnv)
    lockBinding("paraSegTF", environment(PING))
  }
  else if(dataType=="H")
  {
    unlockBinding("paraSegH", environment(PING))
    assign("paraSegH",list(step=step,width=width),envir=environment(PING))
    assign("paraSegH",list(step=step,width=width),envir=.GlobalEnv)
    lockBinding("paraSegH", environment(PING))
  }
}
