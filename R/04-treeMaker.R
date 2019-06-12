#Mutual exclusivities: B and C:
#transitionMat[3,3] <- transitionMat[4,2] <- 0
#for(i in 1:nrow(transitionMat)){
#  transitionMat[i,] <- transitionMat[i,]/sum(transitionMat[i,])
#}

#Other possible inputs: gene names?
#genes <- c('A','B','C','D','E')
setClass("mutationProbs", slots=c(mutation.probs = "matrix", genenames = "character"))

#Let's make a function:
generateMutationProbs <- function(loci, chrs, genenames, mutually.exclusive, probs=NULL){
  rownameset <- sapply(1:length(loci),function(i){paste(chrs[i],'-',loci[i],sep='')})
  transitionMat <- matrix(0,nrow=length(loci)+1,ncol=length(loci))
  namevec <- c(0, rownameset[1:length(loci)])
  rownames(transitionMat) <- namevec
  colnames(transitionMat) <- rownameset[1:length(loci)]
  if(is.null(probs)){
    probset <- sort(runif(length(loci),0,1),decreasing=TRUE)
    probset <- probset/sum(probset)
  }
  transitionMat[1,] <- probset
  for(i in 2:nrow(transitionMat)){
    rowProbs <- probset
    rowProbs[i-1] <- 0
    transitionMat[i,] <- rowProbs
  }
  
  #Mutual exclusivities:
  for(x in 1:length(mutually.exclusive)){
    a <- mutually.exclusive[[x]][1]
    b <- mutually.exclusive[[x]][2]
    transitionMat[3,3] <- transitionMat[4,2] <- 0
    for(i in 1:nrow(transitionMat)){
      transitionMat[i,] <- transitionMat[i,]/sum(transitionMat[i,])
    } 
  }
  list('mutation.probs' = transitionMat, 'genenames' = genenames)
  new("mutationProbs", mutation.probs = transitionMat, genenames = genenames)
}

###To do:
#Add minimum transition prob for (more or less) mutually exclusive mutations.



###Second order probabilities:
#ngenotypes <- 1 + nmuts + factorial(nmuts)/(factorial(2)*factorial(nmuts-2))
##Add poisson parameter for irrelevant mutation acquisition for each genotype?
#transitionMat <- matrix(0,nrow=ngenotypes,ncol=nmuts)
#combos <- combn(loci[1:nmuts],2)
#namevec <- c('Normal', loci[1:nmuts],sapply(1:ncol(combos),function(i){paste(combos[1,i],combos[2,i],sep='')}))
#rownames(transitionMat) <- namevec
#colnames(transitionMat) <- loci[1:nmuts]

#For first and second order probabilities are essentially the same, just removing whichever mutations
#are already present, since you can't get the same mutations twice:
#for(i in 1:nrow(transitionMat)){
#  present <- strsplit(rownames(transitionMat)[i],split='')[[1]]
#  absent  <- which(!colnames(transitionMat) %in% present)
#  transitionMat[i,absent] <- transitionMat[1,absent]
#}

#Except: let's posit two mutual exclusivities: C and D rarely or never co-occur, and BC rarely or never co-occurs
#with E.
#transitionMat[c(3,8,12,16:18),4] <- 0
#transitionMat[c(3,8,12,16:18),4] <- 0

