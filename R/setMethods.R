setMethod("unique", "GenomeData",
function(x,incomparables = FALSE, ...)
{
  GenomeData(lapply(x,function(x){lapply(x,unique)}))
  })

setAs("pingList", "RangedData",
      function(from) {
            makeRangedDataOutput(from, type="bed", filter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500),score=c(1,Inf)),length=100)
      })

setAs("RangedData", "GenomeData",
      function(from) {
            readStart <- ifelse(strand(from) == "-",end(from),start(from))
            alignLocs <-
                split(data.frame(position = readStart, strand = strand(from)),
                      space(from)[drop=TRUE])
            GenomeData(lapply(alignLocs,function(df) with(df, split(position, strand))[c("-", "+")]))
      })

      setAs("data.frame", "GenomeData",            
            function(from) {
                  from<-as(from,"RangedData")
                  readStart <- ifelse(strand(from) == "-",end(from),start(from))
                  alignLocs <-
                      split(data.frame(position = readStart, strand = strand(from)),
                            space(from)[drop=TRUE])
                  GenomeData(lapply(alignLocs,function(df) with(df, split(position, strand))[c("-", "+")]))
            })

setAs("pingList", "data.frame",
function(from)
{
  ans <- data.frame(ID=rep(1:length(from),K(from)),chr=chromosome(from),w=w(from), mu=mu(from),
  					delta=delta(from), sigmaSqF=sigmaSqF(from), sigmaSqR=sigmaSqR(from),se=se(from),
					score=score(from), scoreF=scoreForward(from),scoreR=scoreReverse(from),
					minRange=minRange(from), maxRange=maxRange(from), seF=seF(from), seR=seR(from) )
  ans$chr	<- as.character(ans$chr)
  ans		<- ans[is.finite(ans$mu),]
  return(ans)
})

#setAs("segReadsListPE", "segReadsList",
#function(from)
#{
#	temp.List <- from@List
#	for(i in 1:length(temp.List)) temp.List[[i]] <- as(temp.List[[i]], "segReads")
#	ans <- segReadsList(temp.List,from@paraSW,from@N,from@Nc)
#  return(ans)
#})
#
#setAs("segReadsPE", "segReads",
#function(from)
#{
#  ans <- segReads(from@yF, from@yR, from@cF, from@cR, from@map, from@chr)
#  return(ans)
#})

## show and summary methods
setMethod("show", "segReads",
          function(object)
      {
          cat("Object of class 'segReads'","\n")
          cat("This object has the following slots: \n")
          cat("yR, yF, cF, cR, map\n")
      })

setMethod("show", "segReadsList",
          function(object)
      {
          cat("Object of class 'segReadsList'","\n")
          cat("This object has the following slots: \n")
          cat("List, paraSW, N, Nc\n")
          cat("List is a list of 'segReads' ojects, each of which has the following slots:\n")
          cat("yR, yF, cR, cF, map, chr\n")
      })

setMethod("show", "ping",
      function(object)
      {
        cat("Object of class 'ping'","\n")
        cat("This object has the following slots: \n")
        cat("estimates, score, scoreF, scoreR, Nmerged, converge, chr, range\n")     
        })

setMethod("show", "pingError",
          function(object)
          {
            cat("Object of class 'pingError'","\n")
            cat("This object has the following slot: \n")
            cat("errorCode\n")     
          })

setMethod("show", "pingList",
          function(object)
          {
            cat("Object of class 'pingList'","\n")
            cat("This object has the following slots: \n")
            cat("List, paraEM, paraPrior, minReads, N, Nc\n")
            cat("List is a list of 'ping' or pingError ojects\n")
          })


# setGeneric("score", function(x, ...) standardGeneric("score"))
setMethod("score", "ping",
          function(x)
          {
            return(x@score)
})
          
setMethod("score", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("score", "pingList",
          function(x)
          {
            ans<-.Call("getScore", x@List, PACKAGE="PING");
            return(ans)
            
          })


setGeneric("minRange", function(x, ...) standardGeneric("minRange"))
setMethod("minRange", "ping",
          function(x)
          {
            return(x@range[1])
})
          
setMethod("minRange", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("minRange", "pingList",
          function(x)
          {
            ans<-.Call("getMin", x@List, PACKAGE="PING");
            return(ans)
          })

setGeneric("maxRange", function(x, ...) standardGeneric("maxRange"))
setMethod("maxRange", "ping",
          function(x)
          {
            return(x@range[2])
})
          
setMethod("maxRange", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("maxRange", "pingList",
          function(x)
          {
            ans<-.Call("getMax", x@List, PACKAGE="PING");
            return(ans)
          })


setGeneric("scoreReverse", function(x, ...) standardGeneric("scoreReverse"))
setMethod("scoreReverse", "ping",
          function(x)
          {
            return(x@scoreR)
})
          
setMethod("scoreReverse", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("scoreReverse", "pingList",
          function(x)
          {
            ans<-.Call("getScoreR", x@List, PACKAGE="PING");
            return(ans)
            
          })
          
setGeneric("scoreForward", function(x, ...) standardGeneric("scoreForward"))
setMethod("scoreForward", "ping",
          function(x)
          {
            return(x@scoreF)
})
          
setMethod("scoreForward", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("scoreForward", "pingList",
          function(x)
          {
            ans<-.Call("getScoreF", x@List, PACKAGE="PING");
            return(ans)
            
          })

setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
setMethod("chromosome", "ping",
          function(x)
          {
            return(rep(x@chr,length(x@estimates$w)))
})
setMethod("chromosome", "pingError",
          function(x)
          {
            return(NULL)
          }
)
setMethod("chromosome", "pingList",
          function(x)
          {
            ans<-.Call("getChr", x@List, PACKAGE="PING");
            return(ans)
          }
)

setGeneric("map", function(x, ...) standardGeneric("map"))
setMethod("map", "segReads",
          function(x)
          {
            if(is.null(x@map) | (nrow(x@map)==0))
            {
              return(0);
            }
            else
            {
              n<-nrow(x@map)
              m<-min(x@yF[1],x@yR[1],x@map[1,1]);M<-max(tail(x@yF,1),tail(x@yR,1),x@map[n,2]);
              return(sum(diff(t(x@map)))/max(M-m,1));
            }
})

setMethod("map", "segReadsList",
          function(x)
          {
            ans<-.Call("getMap", x@List, PACKAGE="PING");
            return(ans)
          }
)


setGeneric("se", function(x, ...) standardGeneric("se"))
setMethod("se", "ping",
          function(x)
          {
            return(x@estimates$seMu)
})

setMethod("se", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("se", "pingList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(5), PACKAGE="PING");
              return(ans)
          }
)

setGeneric("seF", function(x, ...) standardGeneric("seF"))
setMethod("seF", "ping",
          function(x)
          {
            return(x@estimates$seMuF)
})

setMethod("seF", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("seF", "pingList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(6), PACKAGE="PING");
              return(ans)
          }
)

setGeneric("seR", function(x, ...) standardGeneric("seR"))
setMethod("seR", "ping",
          function(x)
          {
            return(x@estimates$seMuR)
})

setMethod("seR", "pingError",
          function(x)
          {
            return(NULL)
})


setMethod("seR", "pingList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(7), PACKAGE="PING");
              return(ans)
          }
)


setGeneric("sigmaSqF", function(x, ...) standardGeneric("sigmaSqF"))
setMethod("sigmaSqF", "ping",
          function(x)
          {
            return(x@estimates$sigmaSqF)
})

setMethod("sigmaSqF", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("sigmaSqF", "pingList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(3), PACKAGE="PING");
              return(ans)
          }
)

setGeneric("sigmaSqR", function(x, ...) standardGeneric("sigmaSqR"))
setMethod("sigmaSqR", "ping",
          function(x)
          {
            return(x@estimates$sigmaSqR)
})

setMethod("sigmaSqR", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("sigmaSqR", "pingList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(4), PACKAGE="PING");
              return(ans)
          }
)

setGeneric("delta", function(x, ...) standardGeneric("delta"))
setMethod("delta", "ping",
          function(x)
          {
            return(x@estimates$delta)
})

setMethod("delta", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("delta", "pingList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(2), PACKAGE="PING");
              return(ans)
          }
)


setGeneric("mu", function(x, ...) standardGeneric("mu"))
setMethod("mu", "ping",
          function(x)
          {
            return(x@estimates$mu)
})

setMethod("mu", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("mu", "pingList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(1), PACKAGE="PING");
              return(ans)
          }
)

setGeneric("w", function(x, ...) standardGeneric("w"))
setMethod("w", "ping",
          function(x)
          {
            return(x@estimates$w)
})

setMethod("w", "pingError",
          function(x)
          {
            return(NULL)
})

setMethod("w", "pingList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(0), PACKAGE="PING");
              return(ans)
          }
)

setGeneric("K", function(x, ...) standardGeneric("K"))
setMethod("K", "ping",
          function(x)
          {
            return(length(x@estimates$w))
})

setMethod("K", "pingError",
          function(x)
          {
            return(0)
})

setMethod("K", "pingList",
          function(x)
          {
              ans<-.Call("getK", x@List, PACKAGE="PING");
              return(ans)
          }
)

setGeneric("code", function(x, ...) standardGeneric("code"))
setMethod("code", "ping",
          function(x)
          {
            return("")
})

setMethod("code", "pingError",
          function(x)
          {
            return(x@errorCode)
})

setMethod("code", "pingList",
          function(x)
          {
              temp<-lapply(x@List,"code")
              return(unlist(temp))
          }
)

setMethod("length", "pingList",
          function(x)
          {
            return(length(x@List))
})

setMethod("length", "segReadsList",
          function(x)
          {
            return(length(x@List))
})

#setMethod("length", "segReadsListPE",
#          function(x)
#          {
#            return(length(x@List))
#})

# setGeneric("density", function(x, ...) standardGeneric("density"))
setMethod("density", "ping",
          function(x,strand="+",step=10,sum=FALSE,filter=NULL,scale=TRUE)
          {
            
            # Check that all filters are passed
            missingNames<-!c("delta","sigmaSqF","sigmaSqR","se","score")%in%names(filter)
            filter[c("delta","sigmaSqF","sigmaSqR","se","score")[missingNames]]<-list(c(0,Inf))
            if(strand=="+")
            {
              strand<-1
            }
            else if(strand=="-")
            {
              strand<--1
            }
            else if(strand=="*")
            {
              strand<-0
            }
            else
            {
              stop("Strand must be either '+', '-' or '*'")
            }            
            strand<-as.double(paste(strand,"1",sep=""))
            ans<-.Call("getDensity", x, strand, step, filter, sum, scale, PACKAGE="PING")
            return(ans);
          }
)

setMethod("density", "pingList",
          function(x,strand="+",step=10,sum=FALSE,filter=NULL,scale=TRUE)
          {
            # Check that all filters are passed
            missingNames<-!c("delta","sigmaSqF","sigmaSqR","se","score")%in%names(filter)
            filter[c("delta","sigmaSqF","sigmaSqR","se","score")[missingNames]]<-list(c(0,Inf))

            if(strand=="+")
            {
              strand<-1
            }
            else if(strand=="-")
            {
              strand<--1
            }
            else if(strand=="*")
            {
              strand<-0
            }
            else
            {
              stop("Strand must be either '+', '-' or '*'")
            }
            ans<-.Call("getDensityList", x, strand, step, filter, sum, scale, PACKAGE="PING")
            return(ans);
          }
)

setMethod("density", "pingError",
          function(x,strand=NULL,step=NULL,sum=NULL,filter=NULL)
          {
            return(NULL)
          }
)
setMethod("summary", "segReadsList",
          function(object)
      {
          cat("** Experiment information ** \n")
          cat("Chromosomes interogated: ")
          cat(unique(unlist(lapply(object@List,function(obj){obj@chr}))),"\n")
          cat("Number of reads")
          cat(" in IP: ",object@N," and in control: ",object@Nc,"\n")
          cat("** Segmentation parameters ** \n")
          cat("The following settings were used:\n")          
          cat("  Sliding window half width: ", object@paraSW$width,"\n")
          cat("  Step size: ", object@paraSW$step,"\n")          
          cat("  Minimum number of reads: ", object@paraSW$minReads,"\n")
          cat("** Segmentation summary ** \n")                    
          cat("Number of segmented regions:",length(object@List),"\n")
          cat("Summary on the number of Forward/Reverse reads per region:\n")
          cat("  Forward:\n") 
          cat("  ")
          tempF<-lapply(object@List,function(obj){length(obj@yF)})
          print(summary(as.integer(unlist(tempF))))
          cat("  Reverse:\n") 
          cat("  ")
          tempR<-lapply(object@List,function(obj){length(obj@yR)})
          print(summary(as.integer(unlist(tempR))))
          cat("Summary on the number of control Forward/Reverse reads per region:\n")
          cat("  Forward:\n") 
          cat("  ")
          tempF<-lapply(object@List,function(obj){length(obj@cF)})
          print(summary(as.integer(unlist(tempF))))
          cat("  Reverse:\n") 
          cat("  ")
          tempR<-lapply(object@List,function(obj){length(obj@cR)})
          print(summary(as.integer(unlist(tempR))))                    
          tempMap<-map(object)
          cat("** Mappability summary **\n")
          cat("Non mappable intervals cover an average ", mean(unlist(tempMap)),"% of all regions \n")                  
      })


setMethod("summary", "segReads",
      function(object)
      {
        m<-min(object@yF[1],object@yR[1])
        M<-max(tail(object@yF,1),tail(object@yR,1))        
        cat("** Region summary ** \n")                    
        cat("Summary on Forward reads:\n",summary(object@yF,digits=100),"\n")
        cat("Summary on Reverse reads:\n",summary(object@yR,digits=100),"\n")
        cat("Summary on control Forward reads:\n",summary(object@cF,digits=100),"\n")
        cat("Summary on control Reverse reads:\n",summary(object@cR,digits=100),"\n")        
        cat("Non mappable intervals cover ", sum(diff(t(object@map)))/(M-m),"% of the region \n")
      })


setMethod("[","segReadsList",
		function(x,i, j,..., drop=FALSE)
		{
			if(missing(i))
			{
				return(x)
			}
			if(!missing(j))
			{
			  stop("incorrect number of dimensions")
			}
      else
      {
        segReadsList(x@List[i],x@paraSW,x@N,x@Nc)
      }
		})

setMethod("[[","segReadsList",
    function(x, i, j, ..., exact = TRUE)
    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@List[[i]]
})


#setMethod("[","segReadsListPE",
#		function(x,i, j,..., drop=FALSE)
#		{
#			if(missing(i))
#			{
#				return(x)
#			}
#			if(!missing(j))
#			{
#			  stop("incorrect number of dimensions")
#			}
#      else
#      {
#        segReadsListPE(x@List[i],x@paraSW,x@N,x@NFm,x@NRm,x@Nc,x@NcFm,x@NcRm)
#      }
#		})
#
#setMethod("[[","segReadsListPE",
#    function(x, i, j, ..., exact = TRUE)
#    {
#      if(length(i) != 1)
#      {
#        stop("subscript out of bounds (index must have length 1)")
#      }
#      if(missing(i))
#      {
#        return(x)
#      }
#      if(!missing(j))
#      {
#        stop("incorrect number of dimensions")
#      }
#      x@List[[i]]
#})


setMethod("[","pingList",
		function(x,i, j,..., drop=FALSE)
		{
			if(missing(i))
			{
				return(x)
			}
			if(!missing(j))
			{
			  stop("incorrect number of dimensions")
			}
      else
      {
        newPingList(x@List[i], x@paraEM, x@paraPrior, x@minReads, x@N, x@Nc)        
      }
		})

setMethod("[[","pingList",
    function(x, i, j, ..., exact = TRUE)
    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@List[[i]]
})


setMethod("summary", "pingList",
          function(object)
          {
            cat("** Experiment information ** \n")
            cat("Chromosomes interogated: ")
            cat(unique(chromosome(object)),"\n")
            cat("Number of reads:")
            cat("In IP: ",object@N," in control: ",object@Nc,"\n")
            cat("** Prior parameters ** \n")
            cat("The following settings were used:\n")          
            cat("  Hyper parameters for the fragment length distribution:\n")
            cat("  xi, rho, alpha, beta: ", object@paraPrior$xi,",", object@paraPrior$rho, ",", object@paraPrior$alpha, ",", object@paraPrior$beta,"\n")          
            cat("** Score summary ** \n")                    
            print(summary(score(object)))
            cat("** Fragment length distribution summary ** \n")                    
            print(summary(delta(object)))
            cat("** Summary on the number of binding events per candidate region** \n")
            summary(K(object))
      })


setMethod("summary", "ping",
          function(object)
          {
            cat("** Score ** \n")                    
            cat(score(object),"\n")
            cat("** Fragment length estimate ** \n")                    
            cat(delta(object),"\n")
            cat("** Number of binding events in the candidate region** \n")
            cat(K(object),"\n")
          })



setMethod("plot", signature("ping", "segReads"),
          function(x, y, addKernel=FALSE, addNucleosome=FALSE, addSe=TRUE, main=NULL, rescale=1, ...)
      {
         #Set outer and figure margins to reduce gap between plots
        if(addNucleosome)
        {
          nG<-3
          par(oma=c(2.5,5,5,5),mar=c(0,5,0,0),cex.lab=2)
          layout(matrix(1:nG,ncol=1), heights = c(.4,.5,.1))
        }
        else
        {
          nG<-2
          par(oma=c(2.5,5,5,5),mar=c(0,5,0,0),cex.lab=2)
          layout(matrix(1:nG,ncol=1), heights = c(.5,.5))
        }
        
        step<-5/rescale
        .densityMix<-function(x,para)
        {
          v<-4
          w<-para$w
          mu<-para$mu
          sigmaSq<-para$sigmaSq
          sigma<-sqrt(sigmaSq)
          xNorm<-outer(-mu,x,"+")/sigma 
          return(colSums(w*dt(xNorm,df=v)/sigma))
        }

        yF<-y@yF
        yR<-y@yR
        cF<-y@cF
        cR<-y@cR        
        map<-y@map
        m<-min(yF[1],yR[1])-100/rescale
        M<-max(tail(yF,1),tail(yR,1))+100/rescale

        paraR<-list(w=x@estimates$w, mu=x@estimates$mu+x@estimates$delta/2, sigmaSq=x@estimates$sigmaSqR)
        paraF<-list(w=x@estimates$w, mu=x@estimates$mu-x@estimates$delta/2, sigmaSq=x@estimates$sigmaSqF)

        dR<-.densityMix(seq(m,M,step),paraR)
        dF<-.densityMix(seq(m,M,step),paraF)
        maxRange<-max(c(dF,dR))
        plot(seq(m,M,step),dF,xlim=c(m,M),ylim=c(0,maxRange),lty=2,type="l",xlab="",ylab="density",xaxt='n',axes=FALSE)
        title(main=main,outer=TRUE,cex.main=2)
        axis(2)
        axis(1)
        
        
        lines(seq(m,M,step),dR,lty=2,col=2)

        # if(length(map)>0)
        # {
        #   nMap<-nrow(map)
        #   for(i in 1:nMap)
        #   {
        #     segments(map[i,1], 0, map[i,2], 0,lwd=3,col=3)
        #   }
        # }
        
        # Add kernel density estimate
        if((addKernel==TRUE) & (length(yF)>1 & length(yR)>1))
        {
          dkF<-density(yF)
          dkR<-density(yR)
          lines(dkF,lty=3)
          lines(dkR,col=2,lty=3)
        }

        #Add single components and se's
        K<-length(x@estimates$w)
        for(k in 1:K)
        {
          paraR<-list(w=x@estimates$w[k], mu=x@estimates$mu[k]+x@estimates$delta[k]/2, sigmaSq=x@estimates$sigmaSqR[k])
          paraF<-list(w=x@estimates$w[k], mu=x@estimates$mu[k]-x@estimates$delta[k]/2, sigmaSq=x@estimates$sigmaSqF[k])

          dsR<-.densityMix(seq(m,M,step),paraR)
          dsF<-.densityMix(seq(m,M,step),paraF)

          lines(seq(m,M,step),dsF,lty=1)
          lines(seq(m,M,step),dsR,col=2,lty=1)
        }

        stripchart(yF[1],pch=">",method="stack",cex=2,at=.55,add=FALSE,axes=FALSE,xlim=c(m,M),ylim=c(0,1))        
        if(length(map)>0)
        {
          nMap<-nrow(map)
          symbols((map[,1]+map[,2])/2,rep(.35,nMap),rectangle=cbind(map[,2]-map[,1],rep(.6,nMap)), inches=FALSE, bg=grey(.6), fg=0, add=TRUE,xlim=c(m,M),ylim=c(0,1))
        }

        stripchart(yF,pch=">",method="stack",cex=2,at=.55,axes=FALSE,xlim=c(m,M),ylim=c(0,1),add=TRUE)
        mtext("IP",cex=1.2,side=2,las=2,at=.45)
        stripchart(yR,pch="<",method="stack",cex=2,at=.35,col=2,add=TRUE,offset=-1/3)
        
        abline(h=.45,lty=1,col="grey")
        if(addSe)
        {
          points(x@estimates$mu,rep(.45,K),pch="+",cex=2)
          if (any(x@estimates$seMu!=0))
          {
            points(x@estimates$mu-2*x@estimates$seMu,rep(.45,K),pch="[",cex=1)
            points(x@estimates$mu+2*x@estimates$seMu,rep(.45,K),pch="]",cex=1)
            segments(x@estimates$mu-2*x@estimates$seMu,rep(.45,K),x@estimates$mu+2*x@estimates$seMu,rep(.45,K),lwd=1,lty=rep(1,K))
          }
        }

        if(length(cF)>0)
        {
          stripchart(cF,pch=">",method="stack",at=0.15,cex=2,add=TRUE,xlim=c(m,M),ylab="Cont.",axes=FALSE)
          abline(h=.1,lty=1,col="grey")
        }
        if(length(cR)>0)
        {
          stripchart(cR,pch="<",method="stack",at=0.05,cex=2,col=2,add=TRUE,offset=-1/3)
        }
        if((length(cR)==0) & (length(cF)==0))
        {
        }
        mtext("Cont.",cex=1.2,side=2,las=2,at=.1)
        
        if(addNucleosome)
        {
          plot(c(m,M),c(0,1),axes=FALSE,col=0,ylim=c(0,1),xlim=c(m,M),ylab="")
          seq<-c(seq(-3,0,.1),seq(3,0,-.1))
          sapply(seq,function(shift,x,m,M,K){symbols(x@estimates$mu+shift*se(x),rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg=0, fg=grey(abs(shift)*se(x)/(3*(se(x)))), add=TRUE,xlim=c(m,M),ylim=c(0,1),lwd=2)},x,m,M,K)
          symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg="white", fg=grey(abs(0)), add=TRUE,xlim=c(m,M),ylim=c(0,1))
          
          # if(addSe)
          # {
          #   symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg=grey(.5*pmin(se(x)/50,1)), fg=0, add=TRUE,xlim=c(m,M),ylim=c(0,1))
          # }
          # else
          # {
          #   
          #   symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, fg=0, bg=1, add=TRUE,xlim=c(m,M),ylim=c(0,1))
          # }
          mtext("Nucl.",cex=1.2,side=2,las=2,at=.5)
        }
})

setMethod("plot", signature("pingError", "segReads"),
          function(x, y, addKernel=FALSE, main=NULL, rescale=1,  ...)
      {
        par(oma=c(2.5,5,5,5),mar=c(0,5,0,0))
        layout(matrix(1:2,ncol=1), heights = c(.2,.1))
        
        yF<-y@yF
        yR<-y@yR
        cF<-y@cF
        cR<-y@cR
        map<-y@map
        m<-min(yF[1],yR[1])-100/rescale
        M<-max(tail(yF,1),tail(yR,1))+100/rescale

        stripchart(yF,pch=">",method="stack",cex=2,at=.5,add=FALSE,axes=FALSE,xlim=c(m,M),ylab="Cont | Inp.",ylim=c(0,1))
        stripchart(yR,pch="<",method="stack",cex=2,at=.5,col=2,add=TRUE,offset=-1/3)
        title(main=main,outer=TRUE,cex.main=2)
        
        abline(h=.35,lty=3)

        # Add kernel density estimate
        if((addKernel==TRUE) & (length(yF)>1 & length(yR)>1))
        {
          dkF<-density(yF,bw=75/rescale)
          dkR<-density(yR,bw=75/rescale)
          plot(dkF,lty=3)
          lines(dkR,col=2,lty=3)
        }

        if(length(cF)>0)
        {
          stripchart(cF,pch=">",method="stack",at=0.2,cex=2,add=TRUE,xlim=c(m,M),ylab="Cont.",axes=FALSE)
        }
        if(length(cR)>0)
        {
          stripchart(cR,pch="<",method="stack",at=0.2,cex=2,col=2,add=TRUE)
        }
})

setMethod("plot", signature("pingList", "segReadsList"),
          function(x, y, regionIndex=NULL, addKernel=FALSE, addNucleosome=FALSE, addSe=TRUE,main=NULL, rescale=1, ...)
{
  setMain<-is.null(main)
  if(is.null(regionIndex))
  {
    regionIndex<-1:length(x@List)
  }
  for(i in regionIndex)
  {    
    if(setMain)
    {
      main<-paste(as.character(i)," (",y@List[[i]]@chr,")",sep="")
    }
    if(class(x@List[[i]])!="pingError")
    {
      plot(x@List[[i]],y@List[[i]],addKernel=addKernel, addNucleosome=addNucleosome, addSe=addSe,main=main, rescale=rescale, ...)
    }
    else
    {
      plot(x@List[[i]],y@List[[i]],addKernel=addKernel, main=paste(as.character(i)," (",y@List[[i]]@chr,")",sep=""), rescale=rescale, ...)
      warning("Object of class pingError, no PING density displayed")
    }
  }
})

setMethod("plot", signature("pingList", "pingList"),
          function(x, y, filter=NULL, h=.1, ...)
{
  FDR<-pingFDR(x,y,filter=filter)
  plot(FDR[,2],FDR[,1],xlab="score",ylab="FDR",panel.first=grid(nx=50),...)
  # points(FDR[,2],FDR[,3]/max(FDR[,3]),xaxt="n",yaxt="n",lty=3,col=3,pch=2)
  # axis(4,at=seq(0,1,.05),labels=max(FDR[,3])*seq(0,1,.05))
  FDRex<-FDR[FDR[,1]>0,]
  notDup<-rev(!duplicated(rev(FDRex[,1])))
  lines(FDRex[notDup,2],FDRex[notDup,1],col=2,lty=2,lwd=1.5)
  abline(h=h,lw=1.5,col="grey")
})

# plot function for two ping results in data.frame format
setMethod("plot", signature("data.frame", "data.frame"),
          function(x, y, h=.1, logscale=F, ...)
{
  FDR<-pingFDR2(x,y)
  if(logscale) {FDR$score=log(FDR$score); xlab="log(score)"} else xlab="score"
  plot(FDR[,"score"],FDR[,"FDR"],xlab=xlab,ylab="FDR",panel.first=grid(nx=50), 
  		 ylim=range(tail(head(FDR$FDR,-1),-1)), ...)
  FDRex<-FDR[FDR[,"FDR"] > 0,]
  notDup<-rev(!duplicated(rev(FDRex[,"FDR"])))
  lines(FDRex[notDup,"score"],FDRex[notDup,"FDR"],col=2,lty=2,lwd=1.5)
  lines(FDR[,"score"],FDR[,"FDR"],lty=1,lwd=1.5)
  abline(h=h,lw=1.5,col="grey")
})
