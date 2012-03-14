#make XSET track
makeXSETtrack<-function(dat, chr, m, M, FragmentLenth=200)
{
	extendReads <- function (reads, seqLen = 200, strand = c("+", "-")) 
	{
		if (is(reads, "GenomeDataList") || is(reads, "GenomeData")) 
			gdapply(reads, extendReads, seqLen = seqLen, strand = strand)
		else if (is(reads, "AlignedRead")) {
			rng <- IRanges(ifelse(strand(reads) == "+", position(reads), 
				position(reads) + width(reads) - seqLen), ifelse(strand(reads) == 
				"+", position(reads) + seqLen - 1L, position(reads) + 
				width(reads) - 1L))
			strandIdx <- strand(reads) %in% strand
			split(rng[strandIdx], chromosome(reads)[strandIdx])
		}
		else if (is.list(reads)) {
			seqLen <- as.integer(seqLen)
			if (all(c("+", "-") %in% names(reads))) {
				reads <- reads[strand]
				starts <- as.integer(c(reads[["+"]], reads[["-"]] - 
					seqLen + 1L))
				new("IRanges", start = starts, width = rep(seqLen, 
					length.out = length(starts)), NAMES = NULL)
			}
		}
		else stop("Invalid value for 'reads'")
	}

	chr=as.character(chr)

	EXT   = extendReads(dat[[chr]],seqLen=FragmentLenth,strand="+")
	EXT   = EXT[EXT@start>=m-FragmentLenth]
	EXT   = EXT[EXT@start<=M]
	XSET  = coverage(EXT)
	ends  = start(XSET)
	if (length(XSET)>0)
	{
		xsetP = makeSegmentation(value = head(runValue(XSET),-1), start = head(ends,-1), end = tail(ends,-1), dp = DisplayPars(color = 'blue', lwd = 2,lty = "solid"))
	}else
	{
		xsetP = makeSegmentation(value = 0, start = m, end = M, dp = DisplayPars(color = 'blue', lwd = 1))
	}
	
	EXT   = extendReads(dat[[chr]],seqLen=FragmentLenth,strand="-")
	EXT   = EXT[EXT@start>=m-FragmentLenth]
	EXT   = EXT[EXT@start<=M]
	XSET  = coverage(EXT)
	ends  = start(XSET)
	if (length(XSET)>0)
	{
		xsetM = makeSegmentation(value = head(runValue(XSET),-1), start = head(ends,-1), end = tail(ends,-1), dp = DisplayPars(color = 'red', lwd = 2,lty = "solid"))
	}else
	{
		xsetM = makeSegmentation(value = 0, start = m, end = M, dp = DisplayPars(color = 'red', lwd = 1))
	}

	if(FragmentLenth>147) #this will decide the ylim of plot, i.e. who is make segmentation and who is makebasetrack
	{
		EXT   = extendReads(dat[[chr]],seqLen=147)
		EXT   = EXT[EXT@start>=m-147]
		EXT   = EXT[EXT@start<=M]
		XSET  = coverage(EXT)
		ends  = start(XSET)
		if (length(XSET)>0)
		{
			xset147 = makeSegmentation(value = head(runValue(XSET),-1), start = head(ends,-1), end = tail(ends,-1), dp = DisplayPars(color = gray(0.6), lwd = 2,lty = "solid"))
		}else
		{
			xset147 = makeSegmentation(value = 0, start = m, end = M, dp = DisplayPars(color = gray(0.6), lwd = 2))
		}
	
	
		EXT   = extendReads(dat[[chr]],seqLen=FragmentLenth)
		EXT   = EXT[EXT@start>=m-FragmentLenth]
		EXT   = EXT[EXT@start<=M]
		XSET  = coverage(EXT)
		if (length(XSET)>0)
		{
			XSETs = makeBaseTrack(value=runValue(XSET), base = start(XSET), strand=".", dp = DisplayPars(size=8, lwd=2,color="black", type="l"),trackOverlay=list(xsetP,xsetM,xset147))
		}else
		{
			XSETs = makeBaseTrack(value=0, base = m, strand=".", dp = DisplayPars(size=8, lwd=2,color="black", type="l"), trackOverlay=list(xsetP,xsetM,xset147))
		}
	}else
	{
		EXT   = extendReads(dat[[chr]],seqLen=FragmentLenth)
		EXT   = EXT[EXT@start>=m-FragmentLenth]
		EXT   = EXT[EXT@start<=M]
		XSET  = coverage(EXT)
		ends  = start(XSET)
		if (length(XSET)>0)
		{
			xsetB = makeSegmentation(value = head(runValue(XSET),-1), start = head(ends,-1), end = tail(ends,-1), dp = DisplayPars(color = "black", lwd = 2,lty = "solid"))
		}else
		{
			xsetB = makeSegmentation(value = 0, start = m, end = M, dp = DisplayPars(color = "black", lwd = 2))
		}
	
	
		EXT   = extendReads(dat[[chr]],seqLen=147)
		EXT   = EXT[EXT@start>=m-147]
		EXT   = EXT[EXT@start<=M]
		XSET  = coverage(EXT)
		if (length(XSET)>0)
		{
			XSETs = makeBaseTrack(value=runValue(XSET), base = start(XSET), strand=".", dp = DisplayPars(size=8, lwd=2,color=gray(0.6), type="l"),trackOverlay=list(xsetP,xsetM,xsetB))
		}else
		{
			XSETs = makeBaseTrack(value=0, base = m, strand=".", dp = DisplayPars(size=8, lwd=2,color=gray(0.6), type="l"), trackOverlay=list(xsetP,xsetM,xsetB))
		}
		
	}
	
	return(XSETs)
}

#I overwrite "FragmentLenth" inside the function, the length will be calculated adaptively, as mode of density of predicted "delta" of nucleosomes in the region
# plotgene<-function(chr="chr12", seg, predictions, unfilter=NULL, unfilter.id=NULL,dif=res.WT$dif, datIP=WT, datCtl=KO, datname="WT vs KO", FragmentLenth=150, genename=NULL, genetype="ensembl_gene_id", minbase=NULL, maxbase=NULL, mart, nuc.rec.ID=1,colors=NULL,addGene=T, v.line=NULL)
plotgene<-function(chr="chr12", seg, predictions, unfilter=NULL, unfilter.id=NULL,dif, datIP, datCtl, datname="WT vs KO", FragmentLenth=150, genename=NULL, genetype="ensembl_gene_id", minbase=NULL, maxbase=NULL, mart, nuc.rec.ID=1,colors=NULL,addGene=T, v.line=NULL)
{
#chr="13"; ping=res.WT$ping; dif=res.WT$dif; datIP=WT; datCtl=KO; ipname="WT"; ctlname="KO"; FragmentLenth=200; genename="YML074C"; genetype="ensembl_gene_id"; minbase=NULL; maxbase=NULL

	if(is.null(colors)) colors=c(1,2:length(predictions)+1) #colors of nucleosomes from each predictions, if missing we use all colors except "red"
	chr=as.character(chr)
	ctl=!is.null(datCtl) #test if control data available
	nogenename=is.null(genename)


	Axis<-makeGenomeAxis(add53 = TRUE, add35 = TRUE, littleTicks = TRUE, dp = NULL)
	
	# Genes
	if (addGene)
	{
		if (!nogenename) 
		{
			myrange=range(dif[dif$feature ==genename,c("start","end","start_position","end_position")])
			if(is.null(minbase)) minbase<-round(myrange[1]-200) 
			if(is.null(maxbase)) maxbase<-round(myrange[2]+200)
			gene <- makeGene(id = genename, type = genetype, biomart = mart, dp = DisplayPars(plotId = TRUE, idRotation = 0,cex=0.8))
		}else
		{
			minbase<-round(minbase); maxbase<-round(maxbase)
			if(chr %in% c("M","X","Y")) 
			{
					chr.roman=chr
			}else if (chr =="MT") 
			{
					chr.roman="M"
			}else
			{
					chr.roman=as.roman(substring(chr,4))
			}
			genesplus <- makeGeneRegion(start = minbase, end = maxbase, strand = "+", chromosome =chr.roman, biomart = mart, dp = DisplayPars(plotId = TRUE, idRotation = 0, cex = .6)) 
			genesmin  <- makeGeneRegion(start = minbase, end = maxbase, strand = "-", chromosome =chr.roman, biomart = mart, dp = DisplayPars(plotId = TRUE, idRotation = 0, cex = .6)) 
		}	
	}
	
	nn=length(predictions)
	if (nn>0)
	{
		tt=vector("list",nn); names(tt)=names(predictions)
		for(i in 1:nn)
		{
			predictions[[i]]=predictions[[i]][predictions[[i]]$chr==chr,]
			predictions[[i]]=predictions[[i]][predictions[[i]]$center<=maxbase,]
			predictions[[i]]=predictions[[i]][predictions[[i]]$center>=minbase,]
			if (nrow(predictions[[i]])>0)
			{
				#tt[[i]] = makeOtherPrediction(start = predictions[[i]]$start, end = predictions[[i]]$end)
				if(is.null(predictions[[i]]$se)) predictions[[i]]$se=0
				predictions[[i]]$se=pmin(50,predictions[[i]]$se)
			}else
			{
				if(is.null(predictions[[i]]$se)) predictions[[i]]$se=double(0)
			}
			tt[[i]] = new("NucleosomeTrack",center=predictions[[i]]$center, se = predictions[[i]]$se, step=2,dp = DisplayPars(color=colors[i],size=2,lwd=2))
		}
	}
	
	
	if(nrow(predictions[[1]])==1) FragmentLenth=predictions[[1]]$delta
	if(nrow(predictions[[1]])==0) FragmentLenth=147
	if(nrow(predictions[[1]])>1)
	{
		aa=density(predictions[[1]]$delta)
		FragmentLenth=round(aa$x[which.max(aa$y)])
	}
	

	title = makeTitle(text =paste(datname, ", ", chr, ":", minbase,"-",maxbase, "(",round(maxbase-minbase),"bps), XSET extend ",FragmentLenth, " bps",sep=""), color = 'darkred')
	#transcript <- makeTranscript(id = genename, 	type = genetype, biomart = mart)
	#ideogram <- makeIdeogram(chromosome = chr)

	XSET.IP  = makeXSETtrack(datIP,  chr, m=minbase, M=maxbase, FragmentLenth=FragmentLenth)
	if (ctl) XSET.Ctl = makeXSETtrack(datCtl, chr, m=minbase, M=maxbase, FragmentLenth=FragmentLenth)
	
	#raw reads
	Reads.IP = new("RawRead",start=unlist(datIP[[chr]]), end = unlist(datIP[[chr]])+0.1, strand=rep(names(datIP[[chr]]),lapply(datIP[[chr]],length)),dp = DisplayPars(size=4,lwd=3, color=c("red","blue", type="l")))
	if (ctl) Reads.Ctl = new("RawRead",start=unlist(datCtl[[chr]]), end = unlist(datCtl[[chr]])+0.1, strand=rep(names(datCtl[[chr]]),lapply(datCtl[[chr]],length)),dp = DisplayPars(size=4, lwd=3, color=c("red","blue", type="l")))
	
	
	#MAPchr<-map[chr]
	#myMap<-makeRectangleOverlay(start = start(MAPchr), end = end(MAPchr), c(8,8), dp = DisplayPars(color = "grey", alpha = 0.1))
	#
	#gdPlot(list(Axis, Nucleosomes=Nucle,XSET=myXSET,"XSET+"=myXSETp,"XSET-"=myXSETm,AlignedReads=Reads), overlays=list(rectNucleosomes,myMap),minBase = minbase, maxBase =maxbase)

	#==
	#construct 'gdlist'
	#==
	nGeneTrack=(nogenename+1)*addGene
	gdlist=vector("list", 5+nn+ctl*2-1+nGeneTrack)
	kk=1
	gdlist[kk]=title; names(gdlist)[kk]=""; kk=kk+1

	if (addGene)
	{
		if(!nogenename) 
		{
			gdlist[kk+0:1]=list(gene, Axis)
			names(gdlist)[kk+0:1]=c("Sc1","")
			kk=kk+2
		} else  
		{
			gdlist[kk+0:2]=list(genesplus, Axis, genesmin)
			names(gdlist)[kk+0:2]=c("Sc1+","","Sc1-")
			kk=kk+3
		}
	}else
	{
			gdlist[kk]= Axis
			names(gdlist)[kk]=""
			kk=kk+1
	}
	
	if (ctl) 
	{
		gdlist[kk+0:3]=list( XSET.Ctl, Reads.Ctl, XSET.IP, Reads.IP)
		names(gdlist)[kk+0:3]=paste(rep(c("XSET.","Reads."),2),rep(c("Ctl","IP"),each=2),sep="")
		kk=kk+4
	} else
	{
		gdlist[kk+0:1]=list(XSET.IP, Reads.IP)
		names(gdlist)[kk+0:1]=c("XSET","Reads")
		kk=kk+2
	}
	
	if (nn>0) 
	{
		gdlist[kk-1+1:nn]=tt
		names(gdlist)[kk-1+1:nn]=names(tt)
		#kk=kk+nn
	}


	#==
	#construct 'rectangular overlays'
	#==	
	nfilter=length(unfilter)
	overlay.list=vector("list",nfilter+3)
	kk=1
	
	if (nrow(predictions[[nuc.rec.ID]])>0) #plot nucleosome overlay
	{
		overlay.list[[kk]]<-makeRectangleOverlay(start = predictions[[nuc.rec.ID]]$mu-73, end = predictions[[nuc.rec.ID]]$mu+73, 
				region =c(4-1+nGeneTrack,length(gdlist)), #region =c(4+nogenename,6+nogenename+ctl),  
				dp = DisplayPars(color = "black", alpha = 0.1))
		kk=kk+1
	}
	
	ss=.summarySeg(seg)
	ss=ss[ss$chr==chr,]
	idx=((ss$min>=minbase) & (ss$min<=maxbase)) | ((ss$max>=minbase) & (ss$max<=maxbase))
	ss=ss[idx,]
	if (nrow(ss)>0) #plot segmentation overlay
	{
		overlay.list[[kk]]<-makeRectangleOverlay(start = ss$min, end = ss$max, region = c(4-1+nGeneTrack,5-1+nGeneTrack+ctl),  dp = DisplayPars(color = "green", alpha = 0))
		kk=kk+1
	}

	if (nfilter>0) 
	{
		fOthers=vector("list",nfilter)
		for(i in 1:nfilter)
		{
			unfilter[[i]]=unfilter[[i]][unfilter[[i]]$chr==chr,]
			unfilter[[i]]=unfilter[[i]][unfilter[[i]]$center<=maxbase,]
			unfilter[[i]]=unfilter[[i]][unfilter[[i]]$center>=minbase,]

			fOthers[[i]]=predictions[[ unfilter.id[i] ]][!(predictions[[ unfilter.id[i] ]]$center %in% unfilter[[i]]$center),]
			if (nrow(fOthers[[i]])>0)
			{
				overlay.list[[kk]]=makeRectangleOverlay(start = fOthers[[i]]$center-73, end = fOthers[[i]]$center+73, region =c(unfilter.id[i]+5-1+nGeneTrack+ctl+c(-1,0)),  dp = DisplayPars(color = "red", alpha = 0))
				kk=kk+1
			}
		}
	}

	if(!is.null(v.line)) 
	{
		overlay.list[[kk]]=makeRectangleOverlay(start = v.line, end = v.line+10^(-5),  dp = DisplayPars(color = 6, alpha = 0))
		kk=kk+1	
	}
	overlay.list=overlay.list[1:(kk-1)]


	#browser()
	gdPlot(gdlist, overlays=overlay.list, minBase = minbase, maxBase =maxbase)
}
#plotgene(chr=ff2$chr[i],predictions=list(PING_N=temp2$ping.df,Post_PING=tmp2,TF=TF,NPS=NPS), 
#			 unfilter=list(fAIC3=temp2$ping.df), unfilter.id=c(1), dif=NULL, 
#			 datIP=IP, datCtl=NULL, datname=datname, FragmentLenth=100, genename=NULL, 
#			 genetype="ensembl_gene_id", minbase=ff2$mu[i]-1000, maxbase=ff2$mu[i]+1000,mart=mart, seg=seg)	
#			)


#  #Reads: a list of raw reads for the datasets in "GenomeData" fromat.
#  #seg: a list of PING segmentation results
#  #predictions: a list of PING results AFTER FILTER in data.frame format
#  #unfilter: a list of PING results BEFORE FILTER in data.frame format
#  #shortylab: indicator of whether or not use short version of label in y-axis (longer version might overlap with each other whan # of track is too many)
#  # plotcompare<-function(chr="chr12", seg, predictions, unfilter=NULL, dif=res.WT$dif, Reads=list(WT=WT,KO=KO), datname="WT vs KO", genename=NULL, genetype="ensembl_gene_id", minbase=NULL, maxbase=NULL, mart, colors=NULL,addGene=T, different=NULL,shortylab=F,dif.width=73)
#  plotcompare<-function(chr="chr12", seg, predictions, unfilter=NULL, dif, Reads=list(WT=WT,KO=KO), datname="WT vs KO", genename=NULL, genetype="ensembl_gene_id", minbase=NULL, maxbase=NULL, mart, colors=NULL,addGene=T, different=NULL,shortylab=F,dif.width=73)
#  {
#  #chr="13"; ping=res.WT$ping; dif=res.WT$dif; datIP=WT; datCtl=KO; ipname="WT"; ctlname="KO"; FragmentLenth=200; genename="YML074C"; genetype="ensembl_gene_id"; minbase=NULL; maxbase=NULL
#  
#  	if(is.null(colors)) colors=c(1,2:length(predictions)+1) #colors of nucleosomes from each predictions, if missing we use all colors except "red"
#  	chr=as.character(chr)
#  	nogenename=is.null(genename)
#  
#  	Axis<-makeGenomeAxis(add53 = TRUE, add35 = TRUE, littleTicks = TRUE, dp = NULL)
#  	
#  	# Genes
#  	if (addGene)
#  	{
#  		if (!nogenename) 
#  		{
#  			myrange=range(dif[dif$feature ==genename,c("start","end","start_position","end_position")])
#  			if(is.null(minbase)) minbase<-round(myrange[1]-200) 
#  			if(is.null(maxbase)) maxbase<-round(myrange[2]+200)
#  			gene <- makeGene(id = genename, type = genetype, biomart = mart, dp = DisplayPars(plotId = TRUE, idRotation = 0,cex=0.8))
#  		}else
#  		{
#  			minbase<-round(minbase); maxbase<-round(maxbase)
#  			if(chr %in% c("M","X","Y")) 
#  			{
#  					chr.roman=chr
#  			}else if (chr =="MT") 
#  			{
#  					chr.roman="M"
#  			}else
#  			{
#  					chr.roman=as.roman(substring(chr,4))
#  			}
#  			genesplus <- makeGeneRegion(start = minbase, end = maxbase, strand = "+", chromosome =chr.roman, biomart = mart, dp = DisplayPars(plotId = TRUE, idRotation = 0, cex = .6)) 
#  			genesmin  <- makeGeneRegion(start = minbase, end = maxbase, strand = "-", chromosome =chr.roman, biomart = mart, dp = DisplayPars(plotId = TRUE, idRotation = 0, cex = .6)) 
#  		}	
#  	}
#  	
#  	nn=length(predictions) # number of samples to be compared
#  	if (nn>0)
#  	{
#  		tt=vector("list",nn); names(tt)=names(predictions)
#  		for(i in 1:nn)
#  		{
#  			predictions[[i]]=predictions[[i]][predictions[[i]]$chr==chr,]
#  			predictions[[i]]=predictions[[i]][predictions[[i]]$center<=maxbase,]
#  			predictions[[i]]=predictions[[i]][predictions[[i]]$center>=minbase,]
#  			if (nrow(predictions[[i]])>0)
#  			{
#  				#tt[[i]] = makeOtherPrediction(start = predictions[[i]]$start, end = predictions[[i]]$end)
#  				if(is.null(predictions[[i]]$se)) predictions[[i]]$se=0
#  				predictions[[i]]$se=pmin(50,predictions[[i]]$se)
#  			}else
#  			{
#  				if(is.null(predictions[[i]]$se)) predictions[[i]]$se=double(0)
#  			}
#  			tt[[i]] = new("NucleosomeTrack",center=predictions[[i]]$center, se = predictions[[i]]$se, step=2,dp = DisplayPars(color=colors[i],size=2,lwd=2))
#  		}
#  	}
#  	
#  	title = makeTitle(text =paste(datname, ", ", chr, ":", minbase,"-",maxbase, "(",round(maxbase-minbase),"bps)",sep=""), color = 'darkred')
#  
#  	FragmentLenth=rep(147,nn)
#  	for(i in 1:nn)
#  	{
#  		if(nrow(predictions[[i]])==1) FragmentLenth[i]=predictions[[i]]$delta
#  		if(nrow(predictions[[i]])>1)
#  		{
#  			aa=density(predictions[[i]]$delta)
#  			FragmentLenth[i]=round(aa$x[which.max(aa$y)])
#  		}
#  	}
#  
#  	#transcript <- makeTranscript(id = genename, 	type = genetype, biomart = mart)
#  	#ideogram <- makeIdeogram(chromosome = chr)
#  
#  	XSET.IP=Reads.IP =vector("list",nn)
#  	for(i in 1:nn) 
#  	{
#  		#XSET
#  		XSET.IP[[i]]  = makeXSETtrack(Reads[[i]],  chr, m=minbase, M=maxbase, FragmentLenth=FragmentLenth[i])
#  		#raw reads
#  		Reads.IP[[i]] = new("RawRead",start=unlist(Reads[[i]][[chr]]), end = unlist(Reads[[i]][[chr]])+0.1, strand=rep(names(Reads[[i]][[chr]]),lapply(Reads[[i]][[chr]],length)),dp = DisplayPars(size=4,lwd=3, color=c("red","blue", type="l")))
#  	}
#  	
#  	#MAPchr<-map[chr]
#  	#myMap<-makeRectangleOverlay(start = start(MAPchr), end = end(MAPchr), c(8,8), dp = DisplayPars(color = "grey", alpha = 0.1))
#  	#
#  	#gdPlot(list(Axis, Nucleosomes=Nucle,XSET=myXSET,"XSET+"=myXSETp,"XSET-"=myXSETm,AlignedReads=Reads), overlays=list(rectNucleosomes,myMap),minBase = minbase, maxBase =maxbase)
#  
#  	#==
#  	#construct 'gdlist'
#  	#==
#  	nGeneTrack=(nogenename+1)*addGene
#  	gdlist=vector("list", 2+nn*3+nGeneTrack)
#  	
#  	
#  	# The title track
#  	kk=1; gdlist[kk]=title; names(gdlist)[kk]=""; kk=kk+1
#  
#  	# "nGeneTrack" tracks of gene name, and 1 track of axis
#  	if (addGene)  
#  	{
#  		if(!nogenename) 
#  		{
#  			gdlist[kk+0:1]=list(gene, Axis)
#  			names(gdlist)[kk+0:1]=c("Sc1","")
#  			kk=kk+2
#  		} else  
#  		{
#  			gdlist[kk+0:2]=list(genesplus, Axis, genesmin)
#  			if (shortylab) names(gdlist)[kk+0:2]=c("+","","-") else names(gdlist)[kk+0:2]=c("Sc1+","","Sc1-")
#  			kk=kk+3
#  		}
#  	}else
#  	{
#  			gdlist[kk]= Axis
#  			names(gdlist)[kk]=""
#  			kk=kk+1
#  	}
#  	
#  	# 3*nn tracks of XSET, Raw Reads, and predicted Nucs
#  	for(i in 1:nn) 
#  	{
#  		gdlist[kk+0:2]=list(XSET.IP[[i]], Reads.IP[[i]],tt[[i]])
#  		if (shortylab) names(gdlist)[kk+0:2]=c("",names(predictions)[i],"") else names(gdlist)[kk+0:2]=paste(names(predictions)[i],c("XSET","Reads","Nuc"),sep=".")
#  		kk=kk+3
#  	}
#  	#==
#  	# 'gdlist' constructed
#  	#==
#  	
#  
#  	#==
#  	#construct 'rectangular overlays'
#  	#==	
#  	overlay.list=vector("list",nn*2+1)
#  	kk=1
#  	
#  #	if (nrow(predictions[[nuc.rec.ID]])>0) #plot nucleosome overlay
#  #	{
#  #		overlay.list[[kk]]<-makeRectangleOverlay(start = predictions[[nuc.rec.ID]]$mu-73, end = predictions[[nuc.rec.ID]]$mu+73, 
#  #				region =c(4-1+nGeneTrack,length(gdlist)), #region =c(4+nogenename,6+nogenename+ctl),  
#  #				dp = DisplayPars(color = "black", alpha = 0.1))
#  #		kk=kk+1
#  #	}
#  	
#  	# nn tracks of overlay for segmentation results
#  	ss=vector("list",nn)
#  	for(i in 1:nn)
#  	{
#  		ss[[i]]=.summarySeg(seg[[i]])
#  		ss[[i]]=ss[[i]][ss[[i]]$chr==chr,]
#  		idx=((ss[[i]]$min>=minbase) & (ss[[i]]$min<=maxbase)) | ((ss[[i]]$max>=minbase) & (ss[[i]]$max<=maxbase))
#  		ss[[i]]=ss[[i]][idx,]
#  		if (nrow(ss[[i]])>0) #plot segmentation overlay
#  		{
#  			overlay.list[[kk]]<-makeRectangleOverlay(start = ss[[i]]$min, end = ss[[i]]$max, region = c(3+(i-1)*3+nGeneTrack+c(0,1)),  dp = DisplayPars(color = "green", alpha = 0))
#  			kk=kk+1
#  		}
#  	}
#  
#  	# nn tracks of overlay for filtered Nucs
#  	if(!is.null(unfilter))
#  	{
#  		fOthers=vector("list",nn)
#  		for(i in 1:nn)
#  		{
#  			unfilter[[i]]=unfilter[[i]][unfilter[[i]]$chr==chr,]
#  			unfilter[[i]]=unfilter[[i]][unfilter[[i]]$center<=maxbase,]
#  			unfilter[[i]]=unfilter[[i]][unfilter[[i]]$center>=minbase,]
#  	
#  			fOthers[[i]]=predictions[[i]][!(predictions[[i]]$center %in% unfilter[[i]]$center),]
#  			if (nrow(fOthers[[i]])>0)
#  			{
#  				overlay.list[[kk]]=makeRectangleOverlay(start = fOthers[[i]]$center-73, end = fOthers[[i]]$center+73, region =c(4+(i-1)*3+nGeneTrack+c(0,1)),  dp = DisplayPars(color = "red", alpha = 0))
#  				kk=kk+1
#  			}
#  		}
#  	}
#  
#      # mark differentially enriched nucleosomes
#      if (!is.null(different)) 
#      {
#      	overlay.list[[kk]]=makeRectangleOverlay(start = different-dif.width, end = different+dif.width, region =2+c(1,nn*3)+nGeneTrack,  dp = DisplayPars(color = "blue", alpha = 0))
#      	kk=kk+1
#      }
#  
#  	#remove empty overlays
#  	overlay.list=overlay.list[1:(kk-1)]
#  	#==
#  	# 'rectangular overlays' constructed
#  	#==	
#  
#  
#  	#browser()
#  	gdPlot(gdlist, overlays=overlay.list, minBase = minbase, maxBase =maxbase)
#  }
#plotgene(chr=ff2$chr[i],predictions=list(PING_N=temp2$ping.df,Post_PING=tmp2,TF=TF,NPS=NPS), 
#			 unfilter=list(fAIC3=temp2$ping.df), unfilter.id=c(1), dif=NULL, 
#			 datIP=IP, datCtl=NULL, datname=datname, FragmentLenth=100, genename=NULL, 
#			 genetype="ensembl_gene_id", minbase=ff2$mu[i]-1000, maxbase=ff2$mu[i]+1000,mart=mart, seg=seg)	
#			)