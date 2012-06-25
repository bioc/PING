###############################
## Convert GenomeData object into RangedData objects, default shift each reads into the center of nuc
###############################
## input:
# gd: the genomeData object to be converted
# entend: the length of nucleosome
# shift: shift distance
# lengths: length of each reads in RD object
## output: a rangedData object, each read are shifted to the center of nucleosome
###############################
genomeData2RD <- function(gd, extend=146, shift=73, lengths=1)
{
	spaces <- names(gd)
	dfs <- NULL
	for(sp in spaces)
	{
		tmp <- gd[[sp]]
		pos <- c(tmp$"+", tmp$"-" - extend)
		pos <- pos + shift
		strands <- rep(c("+","-"),sapply(tmp,length))
		dfs <- rbind(dfs, data.frame(start=pos, end = pos +lengths, strand=strands, space=sp))
	}
	return(as(dfs,"RangedData"))
}

###############################
## count number of reads in a genomeData object, return a matrix row for chr and column for strand
###############################
## input:
# dateset: the genomeData object to be counted
## output:
# result: counts of reads in the genomeData object stratified by "chr" and "strand"
###############################
Reads.counts <- function(dataset)
{
	chrnames=names(dataset)
	result = matrix(0,length(chrnames),2)
	rownames(result)=chrnames
	for(idx.chr in chrnames)
	{
		result[idx.chr,] = sapply(dataset[[idx.chr]],length)
 	}
 	result=data.frame(result)
	colnames(result)=names(dataset[[1]])
	return(result)
}

################
# Plot coverage curve of a RangedData object
################
## inputs:
# RD: the RangedData object to be ploted
# xlim: xlim of the plots, if missing calculated from range of RD
# ... : other options for "plot" function
################
plotRD <- function(RD,xlim=NULL,shift=0, add=F, scales=1, fbasis=NULL, plotit=T, ...)
{
	if (length(unique(RD$space)) > 1) 	{print("data in multiple space!"); return(NULL)}

	RD <- RD[as.character(unique(RD$space))]
	if(is.null(RD$score)) 
	{
		covRD <- coverage(RD)
	}else
	{
		scoreRD <- score(RD)
		scale2 <- 10^6 / max(scoreRD) 
		scales <- scales * scale2
		scoreRD <- scoreRD * scale2
		covRD <- coverage(RD,weight=list(scoreRD)) 
	}
	#if(length(names(RD))>1) covRD <- covRD[[1]]
	covRD <- as(covRD[[1]],"vector")
	mat <- cbind(seq_along(covRD) - 0.5, covRD)
	d <- diff(as.numeric(covRD)) != 0 
	mat <- rbind(cbind(mat[d, 1] + 1, mat[d, 2]), mat) 
	mat <- mat[order(mat[, 1]), ] 
	mat[,2] <- mat[,2] / scales 
	mat[,1] <- mat[,1] - shift
	if(is.null(xlim)) 
	{
		xlim <- range(mat[,1])
	}else 
	{
		mat <- mat[mat[,1]>xlim[1] & mat[,1]<xlim[2], ]
	}
	if(!plotit) return(range(mat[,2])) # if do not plot the figure, return the range of yaxis
	
	if(!is.null(fbasis)) plottype="p" else plottype="l"
	if(add)
	{
		lines(mat, type=plottype, ...)
	}else
	{
		plot(mat, type=plottype, xlim=xlim, ...)
	}
	if(!is.null(fbasis))
	{
		fit.sm <- fda:::smooth.basis(mat[,1], mat[,2], fbasis)
		lines(fit.sm$fd,...)
	}
}
#plotRD(RangedData(c(bb[[1]],bb[[2]])),xlim=c(-1000,1000),shift=bb[[3]])
#plotRD(RangedData(c(aa[[1]],aa[[2]])),xlim=c(-1000,1000),shift=aa[[3]],add=T,col=2)

RD2covmat <- function(RD, xlim=NULL, shift=0,  scales=1)
{
	if (length(unique(RD$space)) > 1) 	{print("data in multiple space!"); return(NULL)}

	RD <- RD[as.character(unique(RD$space))]
	if(is.null(RD$score)) 
	{
		covRD <- coverage(RD)
	}else
	{
		scoreRD <- score(RD)
		scale2 <- 10^6 / max(scoreRD) 
		scales <- scales * scale2
		scoreRD <- scoreRD * scale2
		covRD <- coverage(RD,weight=list(scoreRD)) 
	}
	#if(length(names(RD))>1) covRD <- covRD[[1]]
	covRD <- as(covRD[[1]],"vector")
	mat <- cbind(seq_along(covRD) - 0.5, covRD)
	d <- diff(covRD) != 0 
	mat <- rbind(cbind(mat[d, 1] + 1, mat[d, 2]), mat) 
	mat <- mat[order(mat[, 1]), ] 
	mat[,2] <- mat[,2] / scales 
	mat[,1] <- mat[,1] - shift
	if(is.null(xlim)) 
	{
		xlim <- range(mat[,1])
	}else 
	{
		mat <- mat[mat[,1]>xlim[1] & mat[,1]<xlim[2], ]
	}
	return(mat)
}

 
################
# Take subsets of "GD" whose distance to centers give in "subs" is less than 
#	a threshold given in "distance", the subsets of reads is shifted to the centers if "shift=T"
################
## inputs:
# GD: genomeData object to be subseted
# subs: a dataframe of centers used to subset GD, it must has column "chr" and "score" and "mu" or "centers" or "position"
# distance: max distance between Reads givin in GD and centers given in subs
# extend: extend to 3' and 5' direction
# shift: indicator of whether or not shift subseted GD to the centers
################
subsetGD <- function(GD, subs, distance=1000, extend=c(0,146), shift=T)
{
	if(is.null(subs$centers) & !is.null(subs$mu)) subs$centers <- subs$mu
	spaces <- names(GD)
	posP <- posM <- vector("list", length(spaces))
	names(posP) <- names(posM) <- spaces
	extra.shift <- distance + max(extend)
	for(sp in spaces)
	{
		tmpCts <- subset(subs, subs$chr == sp)
		tmpGD <- GD[[sp]]
		if(nrow(tmpCts)==0) next
		resP <- resM <- vector("list", nrow(tmpCts))
		for(i in 1:nrow(tmpCts))
		{
			tmpVecP <- tmpGD$"+"
			tmpVecP <- tmpVecP[(tmpVecP + extend[2]) > (tmpCts$centers[i] - distance)]
			tmpVecP <- tmpVecP[(tmpVecP - extend[1]) < (tmpCts$centers[i] + distance)]
			resP[[i]] <- tmpVecP + (extra.shift - tmpCts$centers[i]) * shift 
			tmpVecM <- tmpGD$"-"
			tmpVecM <- tmpVecM[(tmpVecM + extend[1]) > (tmpCts$centers[i] - distance)]
			tmpVecM <- tmpVecM[(tmpVecM - extend[2]) < (tmpCts$centers[i] + distance)]
			resM[[i]] <- tmpVecM + (extra.shift - tmpCts$centers[i]) * shift
		}
		posP[[sp]] <- unlist(resP)
		posM[[sp]] <- unlist(resM)
	}
	posP <- unlist(posP)
	irP <- IRanges(start = posP - extend[1], end = posP + extend[2])
	posM <- unlist(posM)
	irM <- IRanges(start = posM - extend[2], end = posM + extend[1])
	return(list("+"=irP, "-"=irM, extra.shift=extra.shift))
}
 
 
subsetGD2 <- function(GD, subs, distance=1000, extend=c(0,146), shift=T)
{
	extra.shift <- distance + max(extend)
	if(is.null(subs$centers) & !is.null(subs$mu)) subs$centers <- subs$mu
	if(is.null(subs$centers) & !is.null(subs$position)) subs$centers <- subs$position
	#jitter the duplicated score, since list name does not allow duplications
	if(any(duplicated(subs$score))) subs$score <- subs$score + rnorm(length(subs$score), 0, 0.0001)

	subs.chr <- split(subs, subs$chr)
	spaces <- names(subs.chr)
	nReg <- nrow(subs) #number of regions
	posP <- posM <- vector("list", nReg)
	names(posP) <- names(posM) <- subs$score

	for(sp in spaces)
	{
		tmpCts <- subs.chr[[sp]]
		tmpGD <- GD[[sp]]
		if(nrow(tmpCts)==0) next
		resP <- resM <- vector("list", nrow(tmpCts))
		for(i in 1:nrow(tmpCts))
		{
			idx <- as.character(tmpCts$score[i])
			if (tmpCts$strand[i] == "+") { tmpVecP <- tmpGD$"+"; tmpVecM <- tmpGD$"-"}
			if (tmpCts$strand[i] == "-") { tmpVecP <- tmpGD$"-"; tmpVecM <- tmpGD$"+"}
			tmpVecP <- tmpVecP[(tmpVecP + extend[2]) > (tmpCts$centers[i] - distance)]
			tmpVecP <- tmpVecP[(tmpVecP - extend[1]) < (tmpCts$centers[i] + distance)]
			direc <- (tmpCts$strand[i] == "+") - (tmpCts$strand[i] == "-") #direction of this TSS
			if (length(tmpVecP)>0)
			{
				if(shift) tmpVecP <- (tmpVecP - tmpCts$centers[i]) * direc + extra.shift
				posP[[idx]] <- IRanges(start = tmpVecP - extend[1], end = tmpVecP + extend[2])
			}else
			{
				posP[[idx]] <- IRanges()
			}
			tmpVecM <- tmpVecM[(tmpVecM + extend[1]) > (tmpCts$centers[i] - distance)]
			tmpVecM <- tmpVecM[(tmpVecM - extend[2]) < (tmpCts$centers[i] + distance)]
			if (length(tmpVecM)>0)
			{
				if(shift) tmpVecM <- (tmpVecM - tmpCts$centers[i]) * direc + extra.shift 
				posM[[idx]] <- IRanges(start = tmpVecM - extend[2], end = tmpVecM + extend[1])
			}else
			{
				posM[[idx]] <- IRanges()
			}
		}
	}
	return(list("+"=posP, "-"=posM, extra.shift=extra.shift))
}
 
 
extract.GD.reads <- function(GD, predictions, small=10^-6, calpha=1.5)
{
	browser()
	muF <- predictions$mu - predictions$delta/2
	muR <- predictions$mu + predictions$delta/2
	widthF <- sqrt(predictions$sigmaSqF) * calpha + small
	widthR <- sqrt(predictions$sigmaSqR) * calpha + small
	rd.F <- reduce(RangedData(IRanges(start=muF-widthF,end=muF+widthF), space=predictions$chr))
	rd.R <- reduce(RangedData(IRanges(start=muR-widthR,end=muR+widthR), space=predictions$chr))
	pred <- split(predictions, predictions$chr)
	for(sp in names(pred))
	{
		gd <- GD[[sp]]
		st <- list("+"=start(rd.F[sp]),"-"=start(rd.R[sp]))
		ed <- list("+"=end(rd.F[sp]), "-"=end(rd.R[sp]))
		nr.pred <- nrow(pred[[sp]])
		nReads <- sapply(gd, length)
		# combine reads and the boundaries
		tmp.F <- c(gd$"+", st$"+", ed$"+")
		tmp.R <- c(gd$"-", st$"-", ed$"-")
		# get ranks for the boundaries
		rk.F=rank(tmp.F, ties.method="random")[nReads["+"]+1:(nr.pred*2)]
		rk.R=rank(tmp.R, ties.method="random")[nReads["-"]+1:(nr.pred*2)]
		# take the difference of ranks
		
	
	as.vector(rbind(start()))

	}
}


# plot IRanges object
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),col = "black", sep = 0.5, ...)
{ 
 height <- 1 
 if (is(xlim, "Ranges")) xlim <- c(min(start(xlim)), max(end(xlim))) 
 bins <- disjointBins(IRanges(start(x), end(x) + 1)) 
 plot.new() 
 plot.window(xlim, c(0, max(bins) * (height + sep))) 
 ybottom <- bins * (sep + height) - height 
 rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, col = col, ...) 
 title(main)
 axis(1)
}



# ==========================================================
# = Find index and distance of matched pairs in two vector =
# = Input:
#  xx,yy: the two vectors to match
#  maxdist: distance less than this number is treated as matched
#  maxdup: max number of average match, used to allocate memory, say each elements in xx is expected to match up to 5 elements in yy
# = Output:
#  data frame with 3 columns, each line is a matched pair, first two cols are index in two vectors, and the third col is the distance between matched pair.
# ==========================================================
findMultimatch<- function(xx,yy,maxdist=50,maxdup=5)
{
	nx=length(xx)
	ny=length(yy)
	ox=order(xx)
	oy=order(yy)
	xx=xx[ox]
	yy=yy[oy]
	
	i=j=k=1
	nxy=maxdup*nx
	res=matrix(NA,nxy,3)
	while((i<=nx)&(j<=ny))
	{
		#print(paste("i=",i))
		#print(paste("j=",j))
		dist=xx[i]-yy[j]
		if(dist< -maxdist) 
		{
			i=i+1
		}else if(dist>maxdist)
		{
			j=j+1
		}else
		{
			jj=j; dist2=dist; ss=1
			while(ss==1)
			{
				res[k,]=c(ox[i],oy[jj],dist2)
				k=k+1
				if (jj==ny) 
				{
					ss=0
				}else
				{
					jj=jj+1
					dist2= xx[i]-yy[jj]
					if (abs(dist2)>=maxdist) ss=0
				}
			}
			i=i+1
		}
	}
	if(k>2) res=res[1:(k-1),]
	if(k==2) res=matrix(res[1,],1,3) 
	if(k<2) res=res[0,]
	return(res)
}

# ==========================================================
# = Find index and distance of matched pairs in two vector =
# = each element is only matched with nearest pair, if multiple match is available
# = Input:
#  xx,yy: the two vectors to match
#  maxdist: distance less than this number is treated as matched
#  maxdup: max number of average match, used to allocate memory, say each elements in xx is expected to match up to 5 elements in yy
# = Output:
#  data frame with 3 columns, each line is a matched pair, first two cols are index in two vectors, and the third col is the distance between matched pair.
# ==========================================================
findBestmatch<- function(xx,yy,maxdist=50,maxdup=5)
{
	nx=length(xx)
	ny=length(yy)
	ox=order(xx)
	oy=order(yy)
	xx=xx[ox]
	yy=yy[oy]
	
	i=j=k=1
	nxy=maxdup*nx
	res=matrix(NA,nxy,3)
	while((i<=nx)&(j<=ny))
	{
		#print(paste("i=",i))
		#print(paste("j=",j))
		dist=xx[i]-yy[j]
		if(dist< -maxdist) 
		{
			i=i+1
		}else if(dist>maxdist)
		{
			j=j+1
		}else
		{
			bestj=j; jj=j; bestd=dist; absbestd=abs(bestd); ss=1
			while(ss==1)
			{
				if ((jj==ny)|(bestd<0)) 
				{
					ss=0
				}else
				{
					jj=jj+1
					dist2= xx[i]-yy[jj]
					if (abs(dist2)>=maxdist) 
					{
						ss=0
					}else if(abs(dist2)<absbestd)
					{
						bestj=jj; bestd=dist2; absbestd=abs(bestd)
					}
				}
			}
			res[k,]=c(ox[i],oy[bestj],bestd)
			k=k+1
			i=i+1
		}
	}
	if(k>2) res=res[1:(k-1),]
	if(k==2) res=matrix(res[1,],1,3) 
	if(k<2) res=res[0,]
	return(res)
}


# ===================================================================
# = Usage: find matched/mismatched pairs of predictions, along with their distance
# = Input
# df1,df2: the data frames as results from two nuc positioning methods
# 		   these two dataframe must contain columns "chr"(chromosome name) and "center" (predicted nucleosome locations)
# topN: only topN nucs are used for matching
# maxdist: the max allowable distance between matched nucleosomes =
# maxdup: max number of average match, used to allocate memory, say each nuc in df1 can match up to 5 nucs in df2
# = Note
# the two df results must have same number of chrs and they are named using same rule (say both numeric or roman)
# ===================================================================
# matchResults<-function(df1=temp1$ping.df, df2=TF, topN=NULL, maxdist=50,maxdup=5,uniqueMatch=T)
matchResults<-function(df1, df2, topN=NULL, maxdist=50,maxdup=5,uniqueMatch=T)
{
	if(!is.null(topN))
	{
		if (nrow(df1)>topN) df1=df1[1:topN,]
		if (nrow(df2)>topN) df2=df2[1:topN,]
	}
	
	c1=split(df1,df1$chr)
	c2=split(df2,df2$chr)
	
	nn=length(c1)
	res=vector("list",nn); names(res)=names(c1)
	matched1=matched2=nomatched1=nomatched2=dists=NULL
	for(i in 1:nn)
	{
		if(uniqueMatch)  res[[i]]=findBestmatch(c1[[i]]$center,c2[[i]]$center,maxdist=maxdist,maxdup=maxdup)
		if(!uniqueMatch) res[[i]]=findMultimatch(c1[[i]]$center,c2[[i]]$center,maxdist=maxdist,maxdup=maxdup)
		matched1=rbind(matched1,c1[[i]][res[[i]][,1],])
		matched2=rbind(matched2,c2[[i]][res[[i]][,2],])
		nomatched1=rbind(nomatched1,c1[[i]][-res[[i]][,1],])
		nomatched2=rbind(nomatched2,c2[[i]][-res[[i]][,2],])
		dists=c(dists,res[[i]][,3])
	}

	return(list(matched1=matched1,matched2=matched2,nomatched1=nomatched1,nomatched2=nomatched2,
							distance=dists,n1=length(unique(matched1$rank)),n2=length(unique(matched2$rank))))
}	
#qq1=matchResults(df1=temp1$ping.df, df2=TF, topN=30000, maxdist=50,maxdup=5,uniqueMatch=T)	
#qq2=matchResults(df1=temp1$ping.df, df2=TF, topN=30000, maxdist=50,maxdup=5,uniqueMatch=F)


#area under curve of x vs y, truncate curve in area "x < cc"
trapezoid<- function(x,y,cc) 
{
	ordx=order(x)
	x=x[ordx]
	y=y[ordx]
	if (max(x)<cc)
	{
		#print("warning: cut point is outside of x range, all points are used")
		xcut=x
		ycut=y
	}else
	{
		xlcc=sum(x<cc)
		x1=x[xlcc]; x2=x[xlcc+1]
		y1=x[xlcc]; y2=y[xlcc+1]
		ycc=( y1*(cc-x1) + y2*(x2-cc) ) / (x2-x1)
		xcut=c(x[1:xlcc],cc)
		ycut=c(y[1:xlcc],ycc)
	}
	return(sum( diff(xcut) * ( ycut[-1] + head(ycut,-1) ) )/2)
}


############
# plot the scores of a list of PING or PINGN predictions, 
##########
## Input:
# res: a list of PING or PINGN predictions, each element of it is a dataframe with a "score" column 
# vpos, vcol, vlty: parameters for vertical "abline"
# logscore: indicator of using scores in log-scale
############
plotScores <- function(res, vpos=NULL, vcol=4, vlty=4, logscore=T, ...)
#plot.scores <- function(res, vpos=NULL, vcol=4, vlty=4, logscore=T, ...)
{
	l.res <- length(res)
	scores <- vector("list", l.res)
	names(scores) <- names(res)
	for(i in 1:l.res) 
	{
		scores[[i]] <- res[[i]]$score
		if(logscore) scores[[i]] <- log(scores[[i]])
	}
	ylim <- range(sapply(scores,range))
	xmax <- max(sapply(scores,length))
	if (logscore) 	ylab <- "log(score)" else ylab <- "Nucleosome score"
  plot(NA,xlab="ranks", ylab=ylab, ylim=ylim, xlim =c(0, xmax), ...)
	for(i in 1:l.res) lines(scores[[i]], lty=i, col=i, lwd=2)
	legend("topright", names(scores), lty=1:l.res,col=1:l.res, lwd=2, cex=3)
	abline(v=vpos,col=vcol,lty=vlty)
}


############
# check if parallel version of computing is ready to use
##########
## Output:
# result: a string indicate which version of parrall package to use
# nClust: number of CPUs to be used
############
check.parallel <- function()
{

	result <- "serial"
	nClust <- 1
  
  if(length(grep("parallel",loadedNamespaces()))>0)
  {
    cores<-getOption("cores")
    if(is.null(cores))
    {
      nClust<-parallel:::detectCores()
    }
    else
    {
      nClust<-cores
    }
    message("Using the parallel (parallel) version with ",nClust," cores")
		result <- "parallel"
  }
  
  if(length(grep("snowfall",loadedNamespaces()))==1)
  {
  	if(sfParallel())
  	{
			nClust<-sfCpus()
			message("Using the parallel (snowfall) version with ", nClust, " cpus or cores")
			result <- "snowfall"
	  }
	}

	if(result == "serial") message("Using the serial version")    
	return(list(result=result, nClust=nClust))
}