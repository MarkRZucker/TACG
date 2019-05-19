library(gtools)
#Chromosome lengths and SNP position from a SNP array are necessary for making realistic simulated data:
chlens <- get(load('chlens.rda'))
pos <- get(load('SNP_positions.rda'))

generateSimulationSet <- function(simPath, dataPath, nPerK, rounds=400, nu=0, pcnv=1, norm.contam=FALSE, dataPars=NULL){
  dataPars <- list('snps.seq'=1000000,'snps.cgh'=600000,'mu'=70,'sigma.reads'=25,'sigma0.lrr'=.15, 
                   'sigma0.baf'=.03,'density.sigma'=.1)
  threshold <- .04
  dirs <- list.dirs(simPath)
  dirs <- dirs[-1]
  if(length(dirs)>0){
    setIDs <- as.numeric(sapply(1:length(dirs),function(j){strsplit(dirs[j], split='sims-')[[1]][2]}))
    setID <- max(setIDs) + 1
  }else{
    setID <- 1
  }
  if(!dir.exists(simPath)){
    dir.create(simPath)
  }
  if(!dir.exists(dataPath)){
    dir.create(dataPath)
  }
  ks <- unlist(lapply(1:length(nPerK),function(x){rep(nPerK[x],x)}))
  psis <- get(load(file.path('psis.100.rda')))
  psis <- psis[which(sapply(1:nrow(psis), function(j) {
    length(which(psis[j, ]==0 | psis[j, ] > threshold))==ncol(psis)
  })), ]
  n <- nrow(psis)
  for(j in 1:nrow(psis)){
    go <- FALSE
    while(go==FALSE){
      tumor <- try(Tumor(psis[j,], rounds, nu, pcnv, norm.contam),silent=TRUE)
      if(class(tumor)!='try-error'){
        data <- generateTumorData(tumor, dataPars$snps.seq, dataPars$snps.cgh, dataPars$mu, dataPars$sigma.reads, 
                        dataPars$sigma0.lrr, dataPars$sigma0.baf, dataPars$density.sigma)
      }
      if(class(tumor)!='try-error' & class(data)!='try-error'){
        go <- TRUE
      }
    }
    filename.tumor <- paste(simPath, '/', 'sim-', j, '.rda', sep='')
    filename.data <- paste(dataPath, '/', 'dat-', j, '.rda', sep='')
    sim <- list('tumor'=tumor)
    save(sim, file=filename.tumor)
    dat <- list('dat'=data)
    save(dat, file=filename.data)
  }
  if(nu > 0 & pcnv > 0){
    alt.types <- 'both'
  }else if(nu > 0 & pcnv == 0){
    alt.types <- 'mut'
  }else if(nu == 0 & pcnv > 0){
    alt.types <- 'cnv'
  }else{
    alt.types <- 'none'
  }
  metadata <- list('tumor.params'=list(rounds, nu, pcnv, norm.contam), 'data.params'=dataPars, 'alt.types'=alt.types, 'date'=Sys.Date(),
                   'setID'=setID)
  save(metadata, file=paste(simPath, '/metadata', '.rda', sep=''))
}

#Generating one data set of each type: just CNVs; mutations and CNVs; and just mutaions:
#Just CNVs
set.seed(123)
generateSimulationSet(simPath='sims-cnv',dataPath='data-cnv', nPerK=rep(60,5), rounds=400, nu=0, pcnv=1, norm.contam=FALSE)

#Just mutations:
generateSimulationSet(simPath='sims-mut',dataPath='data-mut', nPerK=rep(60,5), rounds=400, nu=30, pcnv=0, norm.contam=FALSE)

#Both
generateSimulationSet(simPath='sims-both',dataPath='data-both', nPerK=rep(60,5), rounds=400, nu=30, pcnv=1, norm.contam=FALSE)



###Now code to generate a set of 300 artificial mixtures:
mix <- function(datList, psi, dataPath, mixPath, index, alter=FALSE){
  breaks <- lapply(1:22, function(j){
    brks <- sort(unique(unlist(sapply(1:length(datList), function(k){c(datList[[k]][datList[[k]]$chrom==j,]$loc.start, 
                                                                       datList[[k]][datList[[k]]$chrom==j,]$loc.end)}))))
    starts <- brks[1:(length(brks)-1)]
    ends <- c(starts[2:length(starts)]-1, brks[length(brks)])
    data.frame('start'=starts, 'end'=ends, 'chr'=rep(j, length(starts)))
  })
  coords <- Reduce(rbind, breaks)
  columns <- sample(1:2,length(psi),replace=TRUE)
  shifts <- sapply(1:length(datList),function(j){
    lrrs <- unlist(sapply(1:nrow(datList[[j]]),function(k){rep(datList[[j]]$seg.median[k],datList[[j]]$num.mark[k])}))
    dens <- density(na.omit(2*10^lrrs))
    dens$x[which.max(dens$y)] - 2
  })
  d2 <- lapply(1:length(datList),function(j){
    temp <- datList[[j]]
    res <- t(sapply(1:nrow(coords),function(k){
      chrdat <- temp[as.character(temp$chrom)==as.character(coords$chr[k]),]
      seg <- chrdat[chrdat$loc.end>=coords$start[k] & chrdat$loc.start<=coords$end[k],]
      baf <- seg$AvgBAF
      tcn <- 2*10^(seg$seg.median)
      tcn <- tcn - shifts[j]
      if(length(which(tcn<0))>0){
        tcn[which(tcn<0)] <- -tcn[which(tcn<0)]
      }
      X <- tcn*(1-baf)
      Y <- tcn*baf
      if(nrow(seg)>0){
        output <- c(X, Y)
      }else{
        output <- rep(NA, 2)
      }
      output
    }))
    colnames(res) <- c('X', 'Y')
    res
  })
  data.mixed <- as.data.frame(Reduce('+', lapply(1:length(d2),function(k){psi[k]*d2[[k]]})))
  notna <- which(sapply(1:nrow(data.mixed),function(x){length(which(is.na(data.mixed[x,])))==0}))
  for(j in 1:length(d2)){
    d2[[j]] <- d2[[j]][notna,]
  }
  data.mixed <- data.mixed[notna,]
  starts <- coords$start[notna]
  ends <- coords$end[notna]
  chrs <- coords$chr[notna]
  data.mixed$chr <- chrs
  data.mixed$seg <- 1:length(notna)
  data.mixed$LRR <- log10((data.mixed$X + data.mixed$Y)/2)
  data.mixed$BAF <- data.mixed$X/(data.mixed$X + data.mixed$Y)
  data.mixed$markers <- unlist(sapply(1:22,function(i){
    pos.chr <- pos[pos$Chr==i,]
    data.chr <- data.mixed[data.mixed$chr==i,]
    sapply(1:nrow(data.chr),function(j){
      length(which(pos.chr >= starts[which(data.mixed$chr==i)[j]] & pos.chr <= ends[which(data.mixed$chr==i)[j]]))
    })
  }))
  data.mixed$start <- starts
  data.mixed$end <- ends
  data.mixed <- data.mixed[,c(3,4,5,6,1,2,7,8,9)]
  clones <- lapply(1:length(d2),function(k){
    df <- data.frame('A'=rep(1,nrow(d2[[1]])),'B'=rep(1,nrow(d2[[1]])))
    colnames(df) <- c('A','B')
    df$chr <- data.mixed$chr
    df$start <- starts
    df$end <- ends
    df$seg <- 1:length(starts)
    df$parent.index <- rep(0, length(starts))
    df$markers <- data.mixed$markers
    df <- df[,c(3,4,5,1,2,6,7,8)]
    list('cn'=df, 'seq'=NULL)
  })
  for(k in 2:nrow(clones[[1]]$cn)){
    bin1 <- data.mixed$X[k]==data.mixed$X[k-1] & data.mixed$Y[k]==data.mixed$Y[k-1]
    bin2 <- clones[[1]]$cn$chr[k]==clones[[1]]$cn$chr[k-1]
    if(is.na(bin1)){
      bin1 <- FALSE
    }
    if(bin1==TRUE & bin2==TRUE){
      for(l in 1:length(clones)){
        clones[[l]]$cn$start[k] <- clones[[l]]$cn$start[k-1]
        clones[[l]]$cn$markers[k] <- clones[[l]]$cn$markers[k] + clones[[l]]$cn$markers[k-1]
        clones[[l]]$cn[k-1,] <- rep(NA, ncol(clones[[1]]$cn))
      }
      data.mixed$start[k] <- data.mixed$start[k-1]
      data.mixed$markers[k] <- data.mixed$markers[k-1]
      data.mixed[k-1,] <- rep(NA, ncol(data.mixed))
    }
  }
  notna <- unique(unlist(lapply(1:length(clones),function(j){which(!is.na(clones[[j]]$cn$seg))})))
  for(k in 1:length(clones)){
    clones[[k]]$cn <- clones[[k]]$cn[notna,]
    clones[[k]]$cn$seg <- rownames(clones[[k]]$cn) <- 1:nrow(clones[[k]]$cn)
    clones[[k]]$cn
  }
  data.mixed <- data.mixed[notna,]
  dat <- list('cn.data'=data.mixed, 'seq.data'=NULL)
  dat$cn.data <- na.omit(dat$cn.data)
  dat$cn.data$seg <- rownames(dat$cn.data) <- 1:nrow(dat$cn.data)
  if(alter){
    alt.pools <- lapply(1:length(psi[which(psi>0)]),function(j){
      x <- dat$cn.data$X
      y <- dat$cn.data$Y
      markers <- dat$cn.data$markers
      which(markers>1000 & x-shifts[j]/2>.95 & x-shifts[j]/2<1.05 & y-shifts[j]/2>.95 & y-shifts[j]/2<1.05)
    })
    altN <- sample(1:3,length(psi[which(psi>0)]),prob=c(.6,.3,.1),replace=TRUE)
    altered <- lapply(1:length(psi[which(psi>0)]),function(j){sample(alt.pools[[j]],altN[j])})
    change <- lapply(1:length(psi[which(psi>0)]),function(j){sample(c(-1,1),altN[j],replace=TRUE)})
    for(l in 1:length(psi[which(psi>0)])){
      if(length(altered[[l]])>0){
        for(m in 1:length(altered[[l]])){
          allele <- sample(c(1,2),1)
          clones[[l]]$cn[altered[[l]][[m]],allele+3] <- 
            clones[[l]]$cn[altered[[l]][[m]],allele+3] + change[[l]][[m]]
          dat$cn.data[altered[[l]][[m]],allele+4] <- dat$cn.data[altered[[l]][[m]],allele+4] + 
            change[[l]][[m]]*psi[l]
          if(dat$cn.data[altered[[l]][[m]],allele+4] < 0){
            dat$cn.data[altered[[l]][[m]],allele+4] <- -dat$cn.data[altered[[l]][[m]],allele+4]
          }
        }
      }
    }
  }
  tumor <- list('clones'=clones, 'psi'=psi, 'altered'=altered,'change'=change)
  save(dat, file=paste(dataPath, '/mixdat', '-', index, '.rda', sep=''))
  save(tumor, file=paste(mixPath, '/mixsim', '-', index, '.rda', sep=''))
}

generateMixtures <- function(dataPath, mixPath, nPerK, segmentedData, ID_pool){
  threshold <- .04
  psis <- get(load(file.path('psis.100.rda')))
  psis <- psis[which(sapply(1:nrow(psis), function(j) {
    length(which(psis[j, ]==0 | psis[j, ] > threshold))==ncol(psis)
  })), ]
  IDs <- lapply(1:nrow(psis), function(j){sample(ID_pool, length(which(psis[j,]>0)), replace=FALSE)})
  if(!dir.exists(mixPath)){
    dir.create(mixPath)
  }
  if(!dir.exists(dataPath)){
    dir.create(dataPath)
  }
  metadata <- list('IDs'=IDs, 'psis'=psis)
  save(metadata, file=paste(mixPath, '/metadata.rda', sep=''))
  sapply(1:nrow(psis),function(i){mix(lapply(1:length(IDs[[i]]),function(j){
    segmentedData[segmentedData$SamID==IDs[[i]][j],]}),psis[i,], dataPath, mixPath, i, alter=TRUE)})
}

segData <- get(load('rewind.hapmap.rda'))
IDset <- c('NA07019', 'NA12234', 'NA12249', 'NA12753', 'NA12761', 'NA18545', 'NA18975', 'NA18999', 'NA18517')

#Generating the data set:
generateMixtures(dataPath='mixdat', mixPath='mixsim', nPerK=rep(60,5), segmentedData=segData, ID_pool=IDset)
