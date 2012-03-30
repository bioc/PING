########################################################
## Usage:
#   Post-process PING results
## Input:
# ping: PING result before post-process
# seg.boundary: a data.frame object gives the cutting points of segmentations (cut obtained in segmentation step, not the min/max of reads)
#								if seg.boundary is give, all predictions outside of boundary should be removed.
#								This data frame contains 3 columns "ID" for id of candidate regions and "seg.start"/"seg.end" for start/end points of segment
# min.dist=100, minimum distance of two adjacent nucs, smaller than that will be treated as duplicated prediction on boudary
#postPING <- function(ping, seg, rho=8, makePlot=F, datname="", seg.boundary=NULL, DupBound=NULL, IP=NULL, FragmentLenth=100,mart=NULL,sigmaB2=2500,rho1=0.8,alpha1=20,alpha2=100,beta2=100000,xi=160, PE=F, min.dist=100, lambda=-0.000064)
postPING <- function(ping, seg, rho=8, sigmaB2=2500,rho1=0.8,alpha1=20,alpha2=100,beta2=100000,xi=150,  min.dist=100, lambda=-0.000064)
{

	PE=F; makePlot=F; datname=""; seg.boundary=NULL; DupBound=NULL; IP=NULL; FragmentLenth=100; mart=NULL
	
# "Sonication":  {rho1=1.2; sigmaB2=6400;rho=15;alpha1=10; alpha2=98; beta2=200000}
# MNase:  {rho1=3; sigmaB2=4900; rho=8; alpha1=20; alpha2=100; beta2=100000}
	# produce the PING result dataframe and add rank info
	
	ps=as(ping,"data.frame")
	#ps=as.df(ping,seg)
	
	if(!is.null(seg.boundary)) # remove predictions outside of segmention boundary
	{
		temp <- merge(ps,seg.boundary, all.x=T, all.y=F, by="ID")
		idx <- ( (temp$mu >= temp$seg.start) & (temp$mu <= temp$seg.end) )
		ps <- ps[idx, ]
	}
	
	ps=ps[order(ps$score,decreasing=T),]
	ps$rank=1:nrow(ps)

	ps1=PostError(ps=ps, ping=ping, seg=seg, rho=rho, makePlot=makePlot, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=FragmentLenth,rho1=rho1,alpha1=alpha1,xi=xi, PE=PE, lambda=lambda)
	ps2=PostDelta(ps=ps1, ping=ping, seg=seg, rho=rho, makePlot=makePlot, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=FragmentLenth,sigmaB2=sigmaB2,rho1=rho1,alpha1=alpha1,xi=xi, PE=PE,lambda=lambda)
	ps3=PostSigma(ps=ps2, ping=ping, seg=seg, rho=rho, makePlot=makePlot, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=FragmentLenth,mart=mart,sigmaB2=sigmaB2,rho1=rho1,alpha1=alpha1,alpha2=alpha2,beta2=beta2,xi=xi, PE=PE,lambda=lambda)
	PS=PostDup(ps=ps3, ping=ping, seg=seg, rho=rho,rho1=rho1,alpha1=alpha1,xi=xi, PE=PE,min.dist=min.dist,lambda=lambda)

	return(PS)
}
#postPING(ping3, seg, rho=8, makePlot=F, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=100,mart=mart)

#PS$center=round(PS$mu)
#PS$start=PS$mu-PS$delta/2
#PS$end=PS$mu+PS$delta/2
#pdf(paste(datname,"_rho",rho,"_DupBound",DupBound,"_Estimates_Filters_AIC3_after_filterSigma.pdf",sep=""),width=8,height=8)	
#temp3=FilterPING(PS,detail=T,deltaB=c(80,250),sigmaB1=10000,sigmaB2=2500,seB=100)
#dev.off()

# #this function is same as postPING, but it return more details for inspection
# postPING2 <- function(ping, seg, rho=8, makePlot=F, datname="", DupBound=NULL, IP=NULL, FragmentLenth=100,mart,sigmaB2=2500,rho1=0.8,alpha1=20,alpha2=100,beta2=100000)
# {
# 	# produce the PING result dataframe and add rank info
# 	ps=as(ping,"data.frame")
# 	ps=ps[order(ps$score,decreasing=T),]
# 	ps$rank=1:nrow(ps)
# 
# 	ps1=PostError(ps=ps, ping=ping, seg=seg, rho=rho, makePlot=makePlot, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=FragmentLenth,rho1=rho1,alpha1=alpha1)
# 	ps2=PostDelta(ps=ps1, ping=ping, seg=seg, rho=rho, makePlot=makePlot, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=FragmentLenth,sigmaB2=sigmaB2,rho1=rho1,alpha1=alpha1)
# 	ps3=PostSigma(ps=ps2, ping=ping, seg=seg, rho=rho, makePlot=makePlot, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=FragmentLenth,mart=mart,sigmaB2=sigmaB2,rho1=rho1,alpha1=alpha1,alpha2=alpha2,beta2=beta2)
# 	PS=PostDup(ps=ps3, ping=ping, seg=seg, rho=rho,rho1=rho1,alpha1=alpha1)
# 
# 	return(list(before=ps,error=ps1,delta=ps2,sigma=ps3,after=PS))
# }
# #postPING(ping3, seg, rho=8, makePlot=F, datname=datname, DupBound=DupBound, IP=IP, FragmentLenth=100,mart=mart)



########################################################
## Usage:
#   Post-process PING results to solve the singularity problem of PING model fitting
## Input:
# ping: PING output before post-process
# seg:  PING segmentation results
# rho:  hyper-parameter in prior, which control strengh of prior on "delta", larger rho means more confidence on prior guess of "delta"
# ps:   dataframe converted from ping output, which might be already post-processed
# makePlot: indicator of whether or not make plots of pre- and post-processed PING results in pdf files.
## Following input parameters is only useful when we need to make plots in pdf files
# datname: name of datasets, which only used as plot file name
# DupBound: max # of allowed duplicated reads, which is only used in plot
# IP: the reads in IP data in "GenomeData" format
# FragmentLength: the length of XSET profile extention
## Output
# PS: the post-processed PING results used to replace input "ps"
########################################################
PostError <- function(ps, ping, seg, rho=8, makePlot=F, datname="", DupBound=NULL, IP=NULL, FragmentLenth=100,rho1,alpha1,xi, PE,lambda)
{
	idxE=which(code(ping)!="")
	if(length(idxE)==0)
	{
		print("No regions with pingerror")
		PS=ps
	}else
	{
		cat("\n The", length(idxE), "Regions with following IDs are reprocessed for singularity problem: \n")
		print(head(idxE))

		ssE=.summarySeg(seg)[idxE,]
		setParaPrior(xi=xi,rho=rho,alpha=alpha1,beta=20000,lambda=lambda,dMu=200)
		pingE=PING(seg[idxE])
		#change Prior back to default setting
		setParaPrior(xi=xi,rho=rho1,alpha=alpha1,beta=20000,lambda=lambda,dMu=200)
		
		PS1=as(pingE,"data.frame")
		#PS1=as.df(pingE,seg[idxE])
		
		PS1$ID=idxE[PS1$ID]
		ps$rank=NULL
		PS=rbind(ps,PS1)

		#modify rank info to the PING result dataframe
		PS=PS[order(PS$score,decreasing=T),]
		PS$rank=1:nrow(PS)

		
		if (makePlot) ### make plots
		{
			pdf(paste(datname,"_rho",rho,"_DupBound",DupBound,"_pingerror.pdf",sep=""),width=11,height=8.5)
			for(i in 1:nrow(ssE))
			{
				chr=ssE$chr[i]; minbase=ssE$min[i]; maxbase=ssE$max[i]
				Axis<-makeGenomeAxis(add53 = TRUE, add35 = TRUE, littleTicks = TRUE, dp = NULL)
				title = makeTitle(text =paste(datname, ", ", chr, ":", minbase,"-",maxbase, "(",round(maxbase-minbase),"bps), XSET extend ",FragmentLenth, " bps",sep=""), color = 'darkred')
				XSET.IP  = makeXSETtrack(IP, chr, m=minbase, M=maxbase, FragmentLenth=FragmentLenth)
				Reads.IP = new("RawRead",start=unlist(IP[[chr]]), end = unlist(IP[[chr]])+0.1, 
							   strand=rep(names(IP[[chr]]),lapply(IP[[chr]],length)),
							   dp = DisplayPars(size=4,lwd=3, color=c("red","blue", type="l")))
				gdPlot(list(title, XSET=XSET.IP,Axis,Reads=Reads.IP), minBase = minbase, maxBase =maxbase)
			}
			dev.off()
			
			pingE.df.f=FilterPING(PS1)$ping.df
			pdf(paste(datname,"_rho",rho,"_DupBound",DupBound,"_pingEreprocess.pdf",sep=""),width=11,height=8.5)
			for(i in 1:length(pingE))
			{
				plot(pingE[i],seg[idxE[i]])
				ff=pingE.df.f$mu[pingE.df.f$chr==seg[[idxE[i]]]@chr]
				abline(v=ff,col=2)
			}
			dev.off()
		}
	}
	
	#browser()

	
	return(PS)
}

########################################################
## Usage:
#   Post-process PING results to solve the problem of mismatched peaks (i.e. atypical delta)
## Input:
# ping: PING output before post-process
# seg:  PING segmentation results
# rho:  hyper-parameter in prior, which control strengh of prior on "delta", larger rho means more confidence on prior guess of "delta"
# ps:   dataframe converted from ping output, which might be already post-processed
# makePlot: indicator of whether or not make plots of pre- and post-processed PING results in pdf files.
## Following input parameters is only useful when we need to make plots in pdf files
# datname: name of datasets, which only used as plot file name
# DupBound: max # of allowed duplicated reads, which is only used in plot
# IP: the reads in IP data in "GenomeData" format
# FragmentLength: the length of XSET profile extention
## Output
# PS: the post-processed PING results used to replace input "ps"
########################################################

PostDelta <- function(ps, ping, seg, rho=8, makePlot=F, datname="", DupBound=NULL, IP=NULL, FragmentLenth=100, sigmaB2,rho1,alpha1,xi, PE,lambda)
{
	temp0=FilterPING(ps,detail=F,deltaB=c(80,250),sigmaB2=sigmaB2,sigmaB1=10000,seB=Inf,score=0)
	
	#find out who is filtered out by delta	
	idx=(ps$delta<temp0$myFilter$delta[1])|(ps$delta>temp0$myFilter$delta[2])

	if(sum(idx)==0)
	{
		print("No predictions with atypical delta")
		PS=ps
	}else
	{
		ff=ps[idx,]
		ff=ff[order(ff$rank),]
		idxFilt=c(unique(ff$ID))
		cat("\n The", length(idxFilt), "Regions with following IDs are reprocessed for atypical delta: \n")
		print(head(idxFilt))

		setParaPrior(xi=xi,rho=rho,alpha=alpha1,beta=20000,lambda=lambda,dMu=200)
		pingFilt=PING(seg[idxFilt])
		#change Prior back to default setting
		setParaPrior(xi=xi,rho=rho1,alpha=alpha1,beta=20000,lambda=lambda,dMu=200)

	
		tempPS1=as(pingFilt,"data.frame")
		#tempPS1=as.df(pingFilt,seg[idxFilt])
		
		tempPS1$ID=idxFilt[tempPS1$ID]
		tempPS2=subset(ps, !(ps$ID %in% idxFilt) )
		tempPS2$rank=NULL
		PS=rbind(tempPS2,tempPS1)
		PS=PS[order(PS$score,decreasing=T),]
		PS$rank=1:nrow(PS)
		
		if(makePlot)
		{
			pdf(paste(datname,"_rho",rho,"_DupBound",DupBound,"_pingFilt_delta.pdf",sep=""),width=11,height=8.5)
			for(i in 1:length(idxFilt))
			{
				plot(ping[idxFilt[i]],seg[idxFilt[i]])
				ff=temp0$ping.df$mu[temp0$ping.df$chr==seg[[idxFilt[i]]]@chr]
				abline(v=ff,col=2)
			}
			dev.off()
			
			ping.df.f=FilterPING(as(pingFilt,"data.frame"))$ping.df
			pdf(paste(datname,"_rho",rho,"_DupBound",DupBound,"_pingFilt_delta_Reprocess.pdf",sep=""),width=11,height=8.5)
			for(i in 1:length(idxFilt))
			{
					plot(pingFilt[i],seg[idxFilt[i]])
					ff=ping.df.f$mu[ping.df.f$chr==seg[[idxFilt[i]]]@chr]
					abline(v=ff,col=2)
			}
			dev.off()
		}
	}
		
	#browser()
	return(PS)
}

########################################################
## Usage:
#   Post-process PING results to solve the problem of wrongly merged peaks (i.e. too large sigma)
## Input:
# ping: PING output before post-process
# seg:  PING segmentation results
# ps:   dataframe converted from ping output, which might be already post-processed
# rho:  hyper-parameter in prior, which control strengh of prior on "delta", larger rho means more confidence on prior guess of "delta"
# makePlot: indicator of whether or not make plots of pre- and post-processed PING results in pdf files.
## Following input parameters is only useful when we need to make plots in pdf files
# datname: name of datasets, which only used as plot file name
# DupBound: max # of allowed duplicated reads, which is only used in plot
# IP: the reads in IP data in "GenomeData" format
# FragmentLength: the length of XSET profile extention
# mart: saved gene anotation info obtained from "ChIPpeakAnno"
## Output
# PS: the post-processed PING results used to replace input "ps"
########################################################
PostSigma <- function(ps, ping, seg, rho=8, makePlot=F, datname="", DupBound=NULL, IP=NULL, FragmentLenth=100,mart,sigmaB2,rho1,alpha1,alpha2,beta2,xi,PE,lambda)
{
	temp=FilterPING(ps,detail=F,deltaB=c(80,250),sigmaB2=sigmaB2,sigmaB1=10000,seB=Inf, score=0)
	
	#find out who is filtered out by SigmaSq2
	idx4=(ps$sigmaSqF>temp$myFilter$sigmaSq2[2])&(ps$sigmaSqR>temp$myFilter$sigmaSq2[2])

	if(sum(idx4)==0)
	{
		print("No predictions with atypical sigma")
		PS=ps
	}else
	{
			
		ff2=ps[idx4,]
		ff2=ff2[order(ff2$rank),]
		idxFiltSigma=ff2$rank
		#print("Peaks with following IDs are reprocessed for atypical sigma")
		cat("\n The", length(idxFiltSigma), "Peaks with following IDs are reprocessed for atypical sigma: \n")
		print(head(idxFiltSigma))

		newseg2=vector("list",length(idxFiltSigma))
		for(i in 1:length(idxFiltSigma))  {	newseg2[[i]]=processDup(paras=ps[idxFiltSigma[i],],seg=seg,nsigma=2,PE=PE) }
		segSigma=seg; segSigma@List=newseg2
		
		#change hyper-parameters for rho (to avoid atypical delta) and beta (to ask for smaller "sigma")
		setParaPrior(xi=xi,rho=rho,alpha=alpha2,beta=beta2,lambda=lambda,dMu=200)
		system.time(pingSigma<-PING(segSigma))
		#change Prior back to default setting
		setParaPrior(xi=xi,rho=rho1,alpha=alpha2,beta=200000,lambda=lambda,dMu=200)

		tempPS1=as(pingSigma,"data.frame")
		#tempPS1=as.df(pingSigma,segSigma)
		
		tempPS1$ID=ff2$ID[tempPS1$ID]+0.2 # I added 0.2 to the ID name to indicate these peaks are post-processed with large sigma
		tempPS2=ps[!idx4,]
		tempPS2$rank=NULL
		PS=rbind(tempPS2,tempPS1)
		PS=PS[order(PS$score,decreasing=T),]
		PS$rank=1:nrow(PS)
		
		if(makePlot)
		{
			tmp1=ps
			tmp1$center=round(tmp1$mu)
			tmp2=PS
			tmp2$center=round(tmp2$mu)
			temp$ping.df$center=round(temp$ping.df$mu)
		
			pdf(paste(datname,"_rho",rho,"_DupBound",DupBound,"_largeSigma.pdf",sep=""),width=8.5,height=11)	
			for(i in 1:nrow(ff2))
			{
				print(paste("plot figure",i,"of", nrow(ff2)))
				try(plotgene(chr=ff2$chr[i],predictions=list(PING_pre=tmp1,PING_post=tmp2), 
						 unfilter=list(fAIC3=temp$ping.df), unfilter.id=c(1), dif=NULL, 
						 datIP=IP, datCtl=NULL, datname=datname, FragmentLenth=FragmentLenth, genename=NULL, 
						 genetype="ensembl_gene_id", minbase=ff2$mu[i]-1000, maxbase=ff2$mu[i]+1000,mart=mart, seg=seg)	
						)
			}
			dev.off()
		}
	}

	#browser()
	return(PS)
}


########################################################
## Usage:
#   Post-process PING results to solve the problem of Duplicated predictions on segment boundaries
## Input:
# ping: PING result before post-process

PostDup <- function(ps, ping, seg, rho=8,rho1,alpha1,xi,PE,min.dist,lambda)
{
	ps=ps[order(ps$chr,ps$mu),]
	dups=which((diff(ps$ID)!=0)&(diff(ps$mu)<min.dist))
	dups=dups[ps$chr[dups]==ps$chr[dups+1]]
	ndup=length(dups)

	if(ndup==0)
	{
		print("No regions with Boudary problems")
		PS=ps
	}else
	{
		cat("\n The", length(dups), "regions with following IDs are reprocessed for Boundary problems: \n")
		print(head(dups))
	
		newseg=vector("list",ndup)
		for(i in 1:ndup) { newseg[[i]]=processDup(paras=ps[dups[i]+c(0,1),],seg=seg,PE=PE) }
		segDup=seg; segDup@List=newseg
		setParaEM(minK=1,maxK=2,tol=1e-4,B=100,mSelect="AIC3",mergePeaks=T,mapCorrect=T)
		setParaPrior(xi=xi,rho=rho,alpha=alpha1,beta=20000,lambda=lambda,dMu=200)
		system.time(pingDup<-PING(segDup))
		#change Prior back to default setting
		setParaPrior(xi=xi,rho=rho1,alpha=alpha1,beta=20000,lambda=lambda,dMu=200)
		setParaEM(minK=0,maxK=0,tol=1e-4,B=100,mSelect="AIC3",mergePeaks=T,mapCorrect=T)
		
		tempPS1=as(pingDup,"data.frame")
		#tempPS1=as.df(pingDup,segDup)
		
		tempPS1$ID=ps[dups[tempPS1$ID],"ID"]+0.5 # I add 0.5 to indicate these predictions are post-processed duplicated predictions on the boundary of two regions
		ps$rank=NULL
		PS=rbind(ps[-c(dups,dups+1),],tempPS1)
		PS=PS[order(PS$score,decreasing=T),]
		PS$rank=1:nrow(PS)

	}

	#browser()
	return(PS)
}



########################
##Purpose: 
# process duplicated nucs predicted in the overlaped segments
# for duplicated nucs, I extract the reads in their F/R peaks, refit the model
## INPUT:
# paras: From the data frame converted from PING result, we extract the rows of possible duplicated nucs predicted cross boundaries of segments.
# seg: PING segmentation results
## Output:
# a segReads object used to refit PING model
##################
processDup <- function(paras,seg,nsigma=1,PE)
{
	
	chr=as.character(paras$chr[1])
	muv=paras$mu
	deltav=paras$delta
	sigmaSqFv=paras$sigmaSqF
	sigmaSqRv=paras$sigmaSqR
	ID=paras$ID
	
	
	## find boundaries
	sigmaF=sqrt(sigmaSqFv)
	sigmaR=sqrt(sigmaSqRv)
	startF=min(muv-deltav/2-nsigma*sigmaF)
	endF  =max(muv-deltav/2+nsigma*sigmaF)
	startR=min(muv+deltav/2-nsigma*sigmaR)
	endR  =max(muv+deltav/2+nsigma*sigmaR)
	
	#process seg data
	IPF=IPR=CTLF=CTLR=numeric(0)
	map=matrix(as.integer(0),0,2)
	
	for(i in ID)
	{
		temp=seg[[i]]
		IPF=c(IPF,temp@yF)
		IPR=c(IPR,temp@yR)
		CTLF=c(CTLF,temp@cF)
		CTLR=c(CTLR,temp@cR)
		map=rbind(map,temp@map)	
	}
	
	IPF=subset(IPF,IPF>=startF)
	IPF=subset(IPF,IPF<=endF)
	IPR=subset(IPR,IPR>=startR)
	IPR=subset(IPR,IPR<=endR)
	CTLF=subset(CTLF,CTLF>=startF)
	CTLF=subset(CTLF,CTLF<=endF)
	CTLR=subset(CTLR,CTLR>=startR)
	CTLR=subset(CTLR,CTLR<=endR)
	
	if(PE)
	{
#		res=segReadsPE(as.numeric(IPF), as.numeric(IPR), numeric(0), numeric(0), 
#									 as.numeric(CTLF), as.numeric(CTLR), numeric(0), numeric(0), map, chr)		
	}else
	{
		res=segReads(as.numeric(IPF), as.numeric(IPR), as.numeric(CTLF), as.numeric(CTLR), map, chr)
	}
}

#processDup2 <- function(paras,IPdat,CTLdat=NULL, maps=NULL,jitter=F)
#{
#	
#	chr=as.character(paras$chr[1])
#	muv=paras$mu
#	deltav=paras$delta
#	sigmaSqFv=paras$sigmaSqF
#	sigmaSqRv=paras$sigmaSqR
#	
#	## find boundaries
#	sigmaF=sqrt(sigmaSqFv)
#	sigmaR=sqrt(sigmaSqRv)
#	startF=min(muv-deltav/2-sigmaF)
#	endF  =max(muv-deltav/2+sigmaF)
#	startR=min(muv+deltav/2-sigmaR)
#	endR  =max(muv+deltav/2+sigmaR)
#	
#	#process IP data
#	IP=IPdat[[chr]]
#	IPF=IP$"+"
#	IPF=subset(IPF,IPF>=startF)
#	IPF=subset(IPF,IPF<=endF)
#	IPR=IP$"-"
#	IPR=subset(IPR,IPR>=startR)
#	IPR=subset(IPR,IPR<=endR)
#	if (jitter) { IPF=IPF+rnorm(length(IPF),0,0.1);  IPR=IPR+rnorm(length(IPR),0,0.1)}
#	
#	#process control data	
#	if(!is.null(CTLdat))
#	{
#		CTL=CTLdat[[chr]]
#		CTLF=CTL$"+"
#		CTLF=subset(CTLF,CTLF>=startF)
#		CTLF=subset(CTLF,CTLF<=endF)
#		CTLR=CTL$"-"
#		CTLR=subset(CTLR,CTLR>=startR)
#		CTLR=subset(CTLR,CTLR<=endR)	
#		if (jitter) { CTLF=CTLF+rnorm(length(CTLF),0,0.1);  CTLR=CTLR+rnorm(length(CTLR),0,0.1)}
#	}else
#	{
#		CTLF=CTLR=numeric(0)
#	}
#	
#	#process mappability profile
#	if(!is.null(maps))
#	{
#		#map=maps[[chr]]
#		## TO BE ADDED
#	}else
#	{
#		map=matrix(0,0,2)
#	}
#	
#	res=segReads(as.numeric(IPF), as.numeric(IPR), as.numeric(CTLF), as.numeric(CTLR), map, chr)
#}