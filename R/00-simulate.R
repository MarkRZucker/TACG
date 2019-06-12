### want to use Tumor objects as input
generateTumorData <- function(tumor, snps.seq, snps.ary, mu, sigma.reads, sigma0.lrr, sigma0.baf, density.sigma){
  if(!is.na(snps.ary)){
    cndat <- snpDataGen(tumor, snps.ary, sigma0.lrr, sigma0.baf, density.sigma)
  }else{
    cndat <- NA
  }
  if(!is.na(snps.seq)){
    seqdat <- seqDataGen(tumor, snps.seq, density.sigma, mu, sigma.reads)
  }else{
    seqdat <- NA
  }
  list('cn.data'=cndat, 'seq.data'=seqdat)
}

### TODO: Document the algorithm
snpDataGen <- function(tumor, snps.ary=600000, sigma0.lrr=.01, sigma0.baf=.01, segmented=TRUE, snp.rate=NULL, snp.loci=NULL, theta=1, density.sigma=.1, save=FALSE, fn, id){
  psi <- tumor@psi
  if(class(psi)=='WeightVector'){
    psi <- psi@psi
  }
  cn.clones <- lapply(1:length(tumor@clones), function(j){tumor@clones[[j]]$cn})
  chr <- cn.clones[[1]]$chr
  eta <- data.frame('A'=Reduce('+', lapply(1:length(which(psi>0)),function(j){psi[j]*cn.clones[[j]]$A})),
                    'B'=Reduce('+', lapply(1:length(which(psi>0)),function(j){psi[j]*cn.clones[[j]]$B})))
  eta.total <- eta$A + eta$B
  mu.baf <- eta$B/eta.total
  mu.lrr <- log10(eta.total/2)
  lens <- sapply(1:nrow(cn.clones[[1]]),function(j){cn.clones[[1]]$end[j] - cn.clones[[1]]$start[j] + 1})
  lens <- lens/sum(lens)
  markers <- round(lens*snps.ary)
  sigmas.lrr <- sigma0.lrr/(markers)^.5
  sigmas.baf <- sigma0.baf/(markers)^.5
  if(segmented){
    lrr <- rnorm(length(markers), mean=mu.lrr, sd=sigmas.lrr)
    baf <- rnorm(length(markers), mean=mu.baf, sd=sigmas.baf)
    baf[baf<0] <- -baf[baf<0]
    baf[baf>1] <- 1/baf[baf>1]
    total.backcomp <- 2*(10^(lrr))
    X <- total.backcomp*(baf)
    Y <- total.backcomp*(1 - baf)
    df <- data.frame('chr'=chr, seg=1:nrow(eta), 'LRR'=lrr, 'BAF'=baf,'X'=X, 'Y'=Y, 'markers'=markers)
  }else{
    if(is.null(snp.loci)){
      loci <- unlist(lapply(1:length(markers),function(j){sample(tumor@clones[[1]]$cn$start[j]:tumor@clones[[1]]$cn$end[j], markers[j], replace=FALSE)}))
    }else{
      loci <- snp.loci
      markers <- sapply(1:length(markers),function(j){length(which(loci >= tumor@clones[[1]]$cn$start[j] & loci <= tumor@clones[[1]]$cn$end[j]))})
    }
    loci <- sort(loci,decreasing=FALSE)
    lrr <- unlist(lapply(1:length(markers),function(j){rnorm(markers[j], mean=mu.lrr[j], sd=sigma0.lrr)}))
    baf <- unlist(lapply(1:length(markers),function(j){rnorm(markers[j], mean=mu.baf[j], sd=sigma0.baf)}))
    chrvec <- unlist(lapply(1:nrow(tumor@clones[[1]]$cn), function(j){rep(tumor@clones[[1]]$cn$chr[j],markers[j])}))
    hom <- sample(1:length(lrr),round((1-snp.rate)*length(lrr)),replace=FALSE)
    baf[hom] <- rnorm(length(hom),0,sigma0.baf/theta)
    baf[baf<0] <- -baf[baf<0]
    baf[baf>1] <- 1/baf[baf>1]
    invert <- sample(1:length(baf),round(.5*length(lrr)),replace=FALSE)
    baf[invert] <- 1 - baf[invert]
    df <- data.frame('id'=rep(paste('sample.',id,sep=''),length(loci)),'chr'=chrvec[1:length(loci)],
                     'LRR'=lrr[1:length(loci)], 'BAF'=baf[1:length(loci)], 'locus'=loci)
    df <- df[with(df,order(df$locus)),]
    rownames(df) <- 1:nrow(df)
  }
  if(save){
    save(df,file=fn)
    df <- NULL
  }
  df
}

#snp.info is an option for an already made snp df (like from snpDataGen) to generate
#read counts from those in addition to generating new ones.

#SeqDatGen requires auxiliary functions estBetaParams and rbeta2: rbeta but parametrized in terms of mean and sigma.
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  c('alpha' = abs(alpha), 'beta' = abs(beta)) # shouldn't negatives be an error ?! 
}

rbeta2 <- function(n, mu, sigma){
  params <- estBetaParams(mu, sigma^2)
  rbeta(n, params[1], params[2])
}

seqDataGen <- function(tumor, snps.seq=1000000, density.sigma, mu, sigma.reads){
  cn.clones <- lapply(1:length(tumor@clones), function(i){tumor@clones[[i]]$cn})
  chrs <- cn.clones[[1]]$chr
  psi <- as(tumor@psi, "numeric")
  eta <- data.frame('A'=Reduce('+', lapply(1:length(which(psi>0)), function(j){psi[j]*cn.clones[[j]]$A})), 
                    'B'=Reduce('+', lapply(1:length(which(psi>0)), function(j){psi[j]*cn.clones[[j]]$B})))
  eta.total <- eta$A + eta$B
  snp.dens <- rbeta2(1:nrow(eta), .5, density.sigma)
  markers <- as.vector(rmultinom(1, size=snps.seq, prob=snp.dens))
  starts <- cn.clones[[1]]$start
  ends <- cn.clones[[1]]$ends
  snpdf <- matrix(NA, nrow=sum(markers), ncol=7)
  colnames(snpdf) <- c('chr', 'seg', 'mut.id', 'refCounts', 'varCounts', 'VAF', 'totalCounts')
  for(j in 1:length(markers)){
    z <- sample(1:2, 2, replace=FALSE)
    pR <- eta[j, z[1]]/(eta[j, z[1]]+eta[j, z[2]])
    totalCounts <- round(rnorm(markers[j], 2*mu, sigma.reads))
    totalCounts[totalCounts<0] <- -totalCounts[totalCounts<0]
    refCounts <- rbinom(markers[j], totalCounts, pR)
    varCounts <- totalCounts - refCounts
    chr <- rep(chrs[j], markers[j])
    VAFs <- varCounts/(totalCounts)
    seg <- rep(j, markers[j])
    for(k in min(1,length(refCounts)):length(refCounts)){
      snpdf[sum(markers[1:j]) - markers[j] +k, ] <- c(chr[k], seg[k], NA, refCounts[k], varCounts[k], VAFs[k], totalCounts[k])
    }
  }
  snpdf <- as.data.frame(snpdf)
  snpdf$status <- rep('germline', nrow(snpdf))
  seq.clones <- lapply(1:length(tumor@clones), function(i){
    temp <- tumor@clones[[i]]$seq
    temp$uniqID <- sapply(1:nrow(temp), function(j){paste(temp$chr[j],'-',temp$start[j],sep='')})
    temp
  })
  cn.clones <- lapply(1:length(tumor@clones), function(i){tumor@clones[[i]]$cn})
  ## Here be monsters. Next part fails if norm.contam was set to TRUE earlier.
  joined <- unique(Reduce(rbind,lapply(1:length(which(psi > 0)),function(j){
    seq.clones[[j]]
  })))
  uniqID <- unique(sapply(1:nrow(joined),function(x){
    paste(joined$chr[x],'-',joined$start[x],sep='')
  }))
  genes <- joined$gene
  mutdf <- matrix(NA, nrow=length(uniqID), ncol=8)
  colnames(mutdf) <- c('chr', 'seg', 'start', 'refCounts', 'varCounts', 'VAF', 'totalCounts', 'genes')
  if(length(uniqID)>0){
    for(j in 1:length(uniqID)){
      indices <- which(sapply(1:length(which(psi>0)), function(k){uniqID[j] %in% seq.clones[[k]]$uniqID}))
      mutGene <- seq.clones[[indices[1]]]$gene[which(as.character(seq.clones[[indices[1]]]$uniqID)==as.character(uniqID[j]))]
      mutStart <- seq.clones[[indices[1]]]$start[which(as.character(seq.clones[[indices[1]]]$uniqID)==as.character(uniqID[j]))]
      indices.unmutated <- which(sapply(1:length(which(psi>0)), function(k){!uniqID[j] %in% seq.clones[[k]]$uniqID}))
      coefs <- psi[indices]
      mutated <- as.numeric(as.character(sapply(indices, function(k){seq.clones[[k]]$mutated.copies[which(as.character(seq.clones[[k]]$uniqID)==as.character(uniqID[j]))]})))
      normal <- as.numeric(as.character(sapply(indices, function(k){seq.clones[[k]]$normal.copies[which(as.character(seq.clones[[k]]$uniqID)==as.character(uniqID[j]))]})))
      seg <- seq.clones[[indices[1]]]$seg[which(as.character(seq.clones[[indices[1]]]$uniqID)==as.character(uniqID[j]))]
      if(length(indices.unmutated)>0){
        coefs.unmutated <- psi[indices.unmutated]
        eta.unmutated <- sum(psi[indices.unmutated]*sapply(indices.unmutated, function(i){
          cn.clones[[i]]$A[seg]+cn.clones[[i]]$B[seg]}))
      }else{
        eta.unmutated <- 0
      }
      eta <- c(sum(coefs*mutated), sum(coefs*normal) + eta.unmutated)
      pR <- eta[2]/(eta[2]+eta[1])
      totalCount <- round(rnorm(1, 2*mu, sigma.reads))
      totalCount[totalCount<0] <- -totalCount[totalCount<0]
      refCount <- rbinom(1, totalCount, pR)
      varCount <- totalCount - refCount
      VAF <- varCount/(varCount+refCount)
      VAF <- round(VAF,digits = 4)
      if(is.nan(VAF)){
        VAF <- NA
      }
      chr <- seq.clones[[indices[1]]]$chr[which(as.character(seq.clones[[indices[1]]]$uniqID)==as.character(uniqID[j]))]
      mutdf[j, ] <- c(as.character(chr), as.character(seg), as.character(mutStart), refCount, varCount, 
                      VAF, totalCount, as.character(mutGene))
    }
  }
  mutdf <- as.data.frame(mutdf)
  mutdf$status <- rep('somatic', nrow(mutdf))
  vardf <- rbind(snpdf, mutdf)
  vardf <- vardf[with(vardf, order(seg)), ]
  rownames(vardf) <- NULL
  vardf
}


#A plot function to visualize data and verify that the simulation is working.
plotTumorData <- function(tumor, data){
  snpdata <- data$cn.data
  seqdata <- data$seq.data
  cn.clones <- lapply(1:length(tumor@clones), function(i){tumor@clones[[i]]$cn})
  psi <- as(tumor@psi, "numeric")
  markers <- snpdata$markers
  starts <- sapply(1:length(markers), function(j){sum(markers[1:j])-markers[j]})
  ends <- sapply(1:length(markers), function(j){sum(markers[1:j])})
  eta <- data.frame('A'=Reduce('+', lapply(1:length(psi), function(j){psi[j]*cn.clones[[j]]$A})), 
                    'B'=Reduce('+', lapply(1:length(psi), function(j){psi[j]*cn.clones[[j]]$B})))
  lesser <- sapply(1:nrow(eta), function(j){min(eta$A[j], eta$B[j])})
  lesser.data <- sapply(1:nrow(eta), function(j){min(snpdata$X[j], snpdata$Y[j])})
  greater <- sapply(1:nrow(eta), function(j){max(eta$A[j], eta$B[j])})
  greater.data <- sapply(1:nrow(eta), function(j){max(snpdata$X[j], snpdata$Y[j])})
  errors <- c(lesser.data - lesser, greater.data - greater)
  opar <- par(mfrow=c(3, 1))
  on.exit(par(opar))
  plot(0, type='n', xlim=c(0, sum(markers)), ylim=c(0, 5), main='Number of Lesser Allele Copies', 
       xlab='Marker Index', ylab='Expected Copy Number')
  sapply(1:nrow(eta), function(j){
    segments(x0=starts[j], x1=ends[j], y0=lesser[j], y1=lesser[j], col='black')
    segments(x0=starts[j], x1=ends[j], y0=lesser.data[j], y1=lesser.data[j], col='red')
  })
  plot(0, type='n', xlim=c(0, sum(markers)), ylim=c(0, 5), main='Number of Greater Allele Copies', 
       xlab='Marker Index', ylab='Expected Copy Number')
  sapply(1:nrow(eta), function(j){
    segments(x0=starts[j], x1=ends[j], y0=greater[j], y1=greater[j], col='black')
    segments(x0=starts[j], x1=ends[j], y0=greater.data[j], y1=greater.data[j], col='blue')
  })
  plot(density(errors), xlab='Error Index', main='Density of Segment Errors')
  invisible(data)
}
