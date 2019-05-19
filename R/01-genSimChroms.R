#Generates simulations and data:
genSimChroms <- function(N, chr=17, loci, minLen=200000, maxLen=10000000, datapath=NULL, chlens, save=FALSE){
  kset <- sample(1:2,N,replace=TRUE)
  nestedness <- sapply(kset,function(k){sample(0:(k-1),1)})
  psiset <- lenset <- matrix(NA,nrow=N,ncol=2)
  devmat <- startmat <- matrix(NA,nrow=N,ncol=2)
  for(i in 1:nrow(psiset)){
    lenset[i,1] <- round(runif(1,minLen,maxLen))
    if(kset[i]>1){
      psiset[i,1:kset[i]] <- as.vector(rdirichlet(1,alpha=rep(1,kset[i]+1)))[1:kset[i]]
      if(nestedness[i]==0){
        lenset[i,2] <- round(runif(1,minLen,maxLen))
      }else{
        lenset[i,2] <- round(runif(1,.5*lenset[i,1],.8*lenset[i,1]))
      }
    }else{
      psiset[i,1:kset[i]] <- runif(1,0,1)
    }
    if(nestedness[i]==1){
      mid <- round(runif(1,2,chlens[chr]-lenset[i,1]+1))
      start1 <- sample(c(1,mid,chlens[chr]-lenset[i,1]+1),1)
      start2 <- round(runif(1,start1+1,start1+lenset[i,1]-lenset[i,2]))
      startmat[i,] <- c(start1,start2)
      devmat[i,] <- rep(sample(c(1,-1),1),2)
    }else{
      mid <- round(runif(1,max(na.omit(lenset[i,]))+1,chlens[chr]-sum(na.omit(lenset[i,]))-1))
      posits <- sample(c('left','mid','right'),kset[i],replace=FALSE)
      startmat[i,which(posits=='left')] <- 1
      startmat[i,which(posits=='mid')] <- mid
      startmat[i,which(posits=='right')] <- chlens[chr] - lenset[i,which(posits=='right')] + 1
      devmat[i,1:kset[i]] <- sample(c(1,-1),kset[i],replace=TRUE)
    }
  }
  parset <- as.data.frame(cbind(psiset,startmat,lenset,devmat,nestedness))
  colnames(parset) <- c('psi1','psi2','start1','start2','len1','len2','dev1','dev2','nested')
  ###Slight correction: for nested cases, we want psi for the the smaller, interior CNV
  #to be the larger psi, i.e., psi1; I'm just gonna fix it instead of redoing it:
  for(i in which(parset$nested==1)){
    currPsi1 <- parset$psi1[i]
    parset$psi1[i] <- parset$psi2[i]
    parset$psi2[i] <- currPsi1
  }
  sigmaSet <- c(.03,.05)
  kset <- sapply(1:nrow(parset),function(i){length(na.omit(c(parset$psi1[i],parset$psi2[i])))})
  nestedness <- parset$nested
  temp <- lapply(1:nrow(parset),function(q){
    if(kset[q]==1){
      lens <- c(parset$len1[q])
      starts <- c(parset$start1[q])
      A <- matrix(1,ncol=1,nrow=1)
      B <- matrix(1+parset$dev1[q],ncol=1,nrow=1)
      desc <- 0
    }else{
      if(nestedness[q]==0){
        firstStart <- min(c(parset$start1[q],parset$start2[q]))
        secStart <- max(c(parset$start1[q],parset$start2[q]))
        firstLen <- c(parset$len1[q],parset$len2[q])[which.min(c(parset$start1[q],parset$start2[q]))]
        secLen <- c(parset$len1[q],parset$len2[q])[which.max(c(parset$start1[q],parset$start2[q]))]
        lens <- c(parset$len1[q],secStart-(firstStart+firstLen),parset$len2[q])
        starts <- c(parset$start1[q],firstStart+firstLen,parset$start2[q])
        A <- matrix(rep(1,6),nrow=2,ncol=3)
        if(parset$start1[q]<parset$start2[q]){
          B <- matrix(rbind(c(1+parset$dev1[q],1,1),c(1,1,1+parset$dev2[q])),nrow=2,ncol=3)
        }else{
          B <- matrix(rbind(c(1,1,1+parset$dev1[q]),c(1+parset$dev2[q],1,1)),nrow=2,ncol=3)
        }
        desc <- rep(0,0)
      }else{
        starts <- c(parset$start1[q],parset$start2[q],parset$start2[q]+parset$len2[q])
        lens <- c(starts[2]-starts[1],parset$len2[q],starts[1]+parset$len1[q]-starts[3])
        A <- matrix(rep(1,6),nrow=2,ncol=3)
        B <- rbind(c(1,1,1)+parset$dev1[q],c(1,1+parset$dev2[q],1))
        desc <- c(0,1)
      }
    }
    order <- sort(starts,index.return=TRUE)$ix
    starts <- starts[order]
    lens <- lens[order]
    cols <- ncol(A)
    rows <- nrow(A)
    A <- A[,order]
    B <- B[,order]
    if(cols==1 | rows==1){
      A <- matrix(A,nrow=rows,ncol=cols)
      B <- matrix(B,nrow=rows,ncol=cols)
    }
    psi <- na.omit(c(parset$psi1[q],parset$psi2[q]))
    if(starts[1]>1){
      starts <- c(1,starts)
      lens <- c(starts[2]-1,lens)
      A <- cbind(rep(1,nrow(A)),A)
      B <- cbind(rep(1,nrow(B)),B)
    }
    if(starts[length(starts)]+lens[length(starts)]-1<chlens[chr]){
      starts <- c(starts,starts[length(starts)]+lens[length(lens)])
      lens <- c(lens,chlens[chr]-(starts[length(starts)])+1)
      A <- cbind(A,rep(1,nrow(A)))
      B <- cbind(B,rep(1,nrow(B)))
    }
    A <- rbind(A,rep(1,ncol(A)))
    B <- rbind(B,rep(1,ncol(B)))
    desc <- c(desc,0)
    clones <- lapply(1:nrow(A),function(j){
      list('cn'=data.frame('chr'=rep(chr,ncol(A)),'start'=starts,'end'=starts+lens-1,'A'=A[j,],'B'=B[j,],
                           'parent.index'=desc[j],'seg'=1:ncol(A)),'seq'=NULL)
    })
    sim <- new("Tumor",psi=new("WeightVector",psi=psi),clones=clones)
    if(save){
      snpDataGen(sim,snp.loci=loci,snps.cgh=length(loci),sigma0.baf=sigmaSet[1],sigma0.lrr=sigmaSet[1]*5,segmented=FALSE,snp.rate=.33,save=TRUE,
                 fn=paste(datpath,'/simdat',q,'.rda',sep=''),id=q)
      out <- sim
    }else{
      dat <- snpDataGen(sim,snp.loci=loci,snps.cgh=length(loci),sigma0.baf=sigmaSet[1],sigma0.lrr=sigmaSet[1]*5,segmented=FALSE,snp.rate=.33,save=FALSE,
                 fn=NULL,id=q)
      out <- list('dat'=dat,'sim'=sim)
    }
    out
  })
  if(save){
    chrSims <- lapply(1:length(temp),function(q){temp[[q]]$sim})
    save(chrSims,file=paste(datapath,'chrSims.rda',sep=''))
  }else{
    out <- temp
  }
  out
}
