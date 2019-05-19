#Some functions (taken from Kevin's report on segmenting BAF data)
expit <- function(a) exp(a)/(1 + exp(a))
logit <- function(p) log(p/(1 - p))
ford <- function(y) log2(1/(4 * y * (1 - y)))
bacd <- function(w) (1 + sqrt(1 - 2^(-w)))/2


#General function to run segmentation algorithms:
runSegAlgs <- function(i, alg, data, respath,saveRes=TRUE,alpha=NULL,thresh=NULL){
  if(class(data)=='character'){
    simdat <- get(load(paste(data,'/simdat',i,'.rda',sep='')))
  }else{
    simdat <- data
  }
  lrr <- simdat$LRR
  baf <- simdat$BAF
  x <- logit(baf)
  #x[abs(x) > 3.5] <- NA
  y <- expit(x)
  yy <- ford(y)
  chrs <- simdat$chr
  chrom <- unique(chrs)
  df.lrr <- data.frame('LRR'=lrr)
  df.yy <- data.frame('YY'=yy)
  if(nrow(df.lrr) < ncol(df.lrr)){
    df.lrr <- as.data.frame(t(df.lrr))
    df.yy <- as.data.frame(t(df.yy))
  }
  colnames(df.lrr) <- colnames(df.yy) <- sapply(1:ncol(df.lrr),function(x){paste('V',x,sep='')})
  if(is.null(thresh)){
    merging <- "mergeLevels"
    mad.threshold <- 3
  }else{
    merging <- "MAD"
  }
  if(is.null(alpha)){
    alpha <- 0.01
  }
  t1 <- Sys.time()
  if(alg=='HMM'){
    res.lrr <- pSegmentHMM(df.lrr,chrs,mad.threshold = thresh,merging = merging)$outSmoothed
    res.baf <- pSegmentHMM(df.yy,chrs,mad.threshold = thresh,merging = merging)$outSmoothed
  }else if(alg=='GLAD'){
    res.lrr <- pSegmentGLAD(df.lrr,chrs)$outSmoothed
    res.baf <- pSegmentGLAD(df.yy,chrs)$outSmoothed
  }else if(alg=='Haar'){
    res.lrr <- pSegmentHaarSeg(df.lrr,chrs)$outSmoothed
    res.baf <- pSegmentHaarSeg(df.yy,chrs)$outSmoothed
  }else if(alg=='Wavelets'){
    res.lrr <- pSegmentWavelets(df.lrr,chrs)$outSmoothed
    res.baf <- pSegmentWavelets(df.yy,chrs)$outSmoothed
  }else if(alg=='DNAcopy'){
    res.lrr <- pSegmentDNAcopy(df.lrr,chrs,alpha=alpha)$outSmoothed
    res.baf <- pSegmentDNAcopy(df.yy,chrs,alpha=alpha)$outSmoothed
  }else if(alg=='CGH'){
    res.lrr <- pSegmentCGHseg(df.lrr,chrs)$outSmoothed
    res.baf <- pSegmentCGHseg(df.yy,chrs)$outSmoothed
  }
  t2 <- Sys.time()
  rt <- difftime(t2,t1,units = 'mins')
  combined <- data.frame('lrr'=res.lrr[,1],'baf'=res.baf[,1])
  diffdf <- data.frame('lrrDiff'=diff(combined$lrr),'bafDiff'=diff(combined$baf))
  ends <- c(which(diffdf$lrrDiff!=0 | diffdf$bafDiff!=0),nrow(combined))
  starts <- c(1,ends+1)[1:(length(ends))]
  nmarks <- ends - starts + 1
  lrrmed <- sapply(1:length(starts),function(j){median(lrr[starts[j]:ends[j]])})
  bafmed <- sapply(1:length(starts),function(j){median(baf[starts[j]:ends[j]])})
  segments <- data.frame('ID'=rep(simdat$id[1],length(starts)),'chrom'=rep(chrom,length(starts)),'loc.start'=sort(starts),
                         'loc.end'=sort(ends),'num.mark'=nmarks,'baf.median'=bafmed,'lrr.median'=lrrmed)
  segments <- segments[which(segments$num.mark>0),]
  obj <- list('segments'=segments,'runtime'=rt)
  if(saveRes==TRUE){
    dir <- paste(respath,'/',alg,sep='')
    if(!dir.exists(dir)){
      dir.create(dir)
    }
    filepath <- paste(dir,'/',alg,i,'.rda',sep='')
    save(obj,file=filepath)
    output <- NULL
  }else{
    output <- obj
  }
  output
}

#Run a set of many samples in parallel:
runSegAlgsOnSet <- function(algs, indices, cores=1, data, respath, alpha=NULL, thresh=NULL, save=FALSE){
  ncores <- min(detectCores(),cores)
  if(class(data)=='character'){
    indexLists <- lapply(1:length(algs),function(k){
      output <- indices
      respath <- paste(respath,'/',algs[k],sep='')
      files <- list.files(respath)
      numbers <- as.numeric(sapply(files,function(file){strsplit(strsplit(file,split=algs[k])[[1]][2],split='.rda')[[1]][1]}))
      alreadyDone <- which(indices %in% numbers)
      if(length(alreadyDone)>0){
        output <- indices[-alreadyDone]
      }
      output
    })
  }else{
    indexLists <- lapply(1:length(algs),function(k){1:length(data)})
  }
  algvec <- algs
  algPars <- data.frame('i'=unlist(indexLists),
      'alg'=unlist(lapply(1:length(indexLists),function(j){rep(algs[j],length(indexLists[[j]]))})))
  if(class(data)=='character'){
    d <- rep(data,nrow(algPars))
    rpath <- respath
  }else{
    d <- lapply(1:nrow(algPars),function(j){data[[as.numeric(algPars$i[j])]]})
    rpath <- NULL
  }
  if(cores>1){
    cl <- makeCluster(ncores)
    clusterExport(cl=cl,list('runSegAlgs','expit','logit','ford','bacd','dfdbd','algvec',
                             'algPars','d','rpath','thresh','alpha'),envir=environment())
    clusterEvalQ(cl, library("ADaCGH2","DNAcopy"))
    out <- parLapply(cl,1:nrow(algPars),function(j){
      runSegAlgs(i=algPars$i[j],alg=algPars$alg[j],data=d[[j]],respath=rpath,saveRes=TRUE)
    })
  }else{
    out <- lapply(1:nrow(algPars),function(j){
      runSegAlgs(i=algPars$i[j],alg=algPars$alg[j],data=d[[j]],respath=rpath,saveRes=FALSE)
    })
  }
  out
}

#Needed for DNAcopy:
dfdbd <- data(default.DNAcopy.bdry, package = "DNAcopy", 
              envir = environment())

#Matching inferred to true breaks:
computeR <- function(true,inferred){
  R <- rep(NA,length(true))
  if(length(inferred)>0){
    R <- sapply(1:length(true),function(i){min(abs(true[i]-inferred))})
  }
  R
}

#Assessment function:
assess <- function(alg, res, set, sims){
  if(class(res)=='character'){
    respath <- res
    objs <- lapply(set,function(i){get(load(paste(respath,'/',alg,'/',alg,i,'.rda',sep='')))})
  }else{
    objs <- res
  }
  if(class(sims)=='character'){
    chrSims <- get(load(sims))
  }else{
    chrSims <- sims
  }
  cnSegments <- Reduce(rbind,lapply(1:length(objs),function(i){
    if(is.null(objs[[i]]$segments)){
      out <- objs[[i]]
    }else{
      out <- objs[[i]]$segments
    }
    out
  }))
  runtimes <- lapply(1:length(objs),function(i){
    if(is.null(objs[[i]]$runtime)){
      out <- NA
    }else{
      out <- objs[[i]]$runtime
    }
    out
  })
  trueBreaks <- lapply(1:length(objs),function(i){
    tb <- chrSims[[i]]$clones[[1]]$cn$start[2:nrow(chrSims[[i]]$clones[[1]]$cn)]
    tb
  })
  inferredBreaks <- lapply(1:length(objs),function(i){
    breaks <- cnSegments[which(cnSegments$ID==paste('sample.',i,sep='')),]$loc.start[2:max(2,length(which(cnSegments$ID==paste('sample.',i,sep=''))))]
    breaks.loci <- loci[breaks]
    breaks.loci
  })
  falsePos <- sapply(1:length(objs),function(i){length(which(cnSegments$ID==paste('sample.',i,sep=''))) - (nrow(chrSims[[i]]$clones[[1]]$cn))})
  falseNeg <- -1*falsePos
  falseNeg[falseNeg<0] <- 0
  falsePos[falsePos<0] <- 0
  R <- lapply(1:length(inferredBreaks),function(i){computeR(inferred = inferredBreaks[[i]],true = trueBreaks[[i]])})
  list('R'=R,'falsePos'=falsePos, 'falseNeg'=falseNeg,'runtime'=runtimes,'alg'=alg)
}
