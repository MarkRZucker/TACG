chdata <- read.table("chlens.txt", sep="\t") # chromnumber, length
chnames <- paste("chr", c(1:22, "X", "Y"), sep="")
chdata$Chrom <- factor(chnames[chdata$V1], levels = chnames)
chlens <- as.vector(tapply(chdata$V2, list(chdata$Chrom), max))
names(chlens) <- chnames

save(chlens, file = "../../Data/sysdata.rda")
