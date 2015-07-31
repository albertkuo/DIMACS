library(ggplot2)
library(reshape2)
## === Exploratory Analysis with eRNA ===

## === enhancer prediction with eRNA
d = read.delim("~/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt")
selected <- names(d) %in% c("Id","CNhs12331", "CNhs12334", "CNhs12325", "CNhs12328") 
subd <- d[selected]

# Plot histograms
t <- melt(subd)
# TPM > 0
t0 = subset(t, value>=0)
ggplot(t0,aes(x = value)) + 
  facet_wrap(~variable,scales = "free") + 
  geom_histogram()

# TPM >1
t1 = subset(t, value>1)
ggplot(t1,aes(x = value)) + 
  facet_wrap(~variable,scales = "free") + 
  geom_histogram() + 
  ggtitle("TPM > 1")

# TPM >5
t2 = subset(t, value>5)
ggplot(t2,aes(x = value)) + 
  facet_wrap(~variable,scales = "free") + 
  geom_histogram() + 
  ggtitle("TPM > 5")

# change "ID" to interval range
splitId <- function(str){
  ls = strsplit(str,":")[[1]]
  chr = ls[1]
  ls = strsplit(ls[2],"-")[[1]]
  lower = strtoi(ls[1])
  upper = strtoi(ls[2])
  return(list(chr,lower,upper))
}

subd$Idlist <- lapply(subd$Id, function(x) splitId(as.character(x)))
subd$chr <- lapply(subd$Idlist, function(x) unlist(x[1]))
subd$lower <- lapply(subd$Idlist, function(x) unlist(x[2]))
subd$upper <- lapply(subd$Idlist, function(x) unlist(x[3]))
subd$Idlist <- NULL
View(subd)

## === enhancer prediction with chromatin marks
# File Names
s_filenames <- list("GM12878_Spectacle_strong_enhancers.bed.gz", 
                    "HepG2_Spectacle_strong_enhancers.bed.gz", 
                    "K562_Spectacle_strong_enhancers.bed.gz",
                    "HeLa-S3_Spectacle_strong_enhancers.bed.gz")

w_filenames <- list("GM12878_Spectacle_weak_enhancers.bed.gz", 
                    "HepG2_Spectacle_weak_enhancers.bed.gz", 
                    "K562_Spectacle_weak_enhancers.bed.gz",
                    "HeLa-S3_Spectacle_weak_enhancers.bed.gz")

# Get data given filenames
read_enhancers <- function(filenames){
  enhancers <- list()
  for (i in 1:length(filenames)){
    enhancers[[i]] <- read.table(gzfile(filenames[[i]]), header=T, stringsAsFactors=FALSE)
    enhancers[[i]] <- enhancers[[i]][c(1,2,3)]
    colnames(enhancers[[i]]) = c("chr","lower","upper")
  }
  return(enhancers)
}

merge.all <- function(x, y) {
  merge(x, y, all=T, by="Id") 
}

strong <- read_enhancers(s_filenames)
weak <- read_enhancers(w_filenames)

## === Check range overlap (match strong/weak enhancers eRNA data)
library(IRanges)
rangesA <- split(IRanges(rangesA$start, rangesA$stop), rangesA$chrom)
rangesB <- split(IRanges(rangesB$start, rangesB$stop), rangesB$chrom)
overlap <- function(chr_subd, lower_subd, upper_subd, chr_e, lower_e, upper_e){
    query <- IRanges(lower_subd, upper_subd)
    qpartition <- factor(chr_subd)
    qlist <- split(query,qpartition)
    
    subject <- IRanges(lower_e, upper_e)
    spartition <- factor(chr_e)
    slist <- split(subject, spartition)
    return(unlist(overlapsAny(qlist,slist)))
  }

merged <- subd
for (i in 1:4){
  s = strong[[i]]
  w = weak[[i]]
  merged[,length(merged)+1] <- overlap(unlist(subd$chr),unlist(subd$lower),unlist(subd$upper),unlist(s$chr),unlist(s$lower),unlist(s$upper))
  merged[,length(merged)+1] <- overlap(unlist(subd$chr),unlist(subd$lower),unlist(subd$upper),unlist(w$chr),unlist(w$lower),unlist(w$upper))
}
# Rename columns
colnames(merged)[c(9,11,13,15)] = c("hs12331_s","hs12334_s","hs12325_s","hs12328_s")   
colnames(merged)[c(10,12,14,16)] = c("hs12331_w","hs12334_w","hs12325_w","hs12328_w") 

# Label with Types: 0 is none, 1 is Weak, 2 is Strong, 3 is Strong and Weak (both types overlap with interval)
for(i in 17:20){
  merged[,i] <- 2*merged[,2*i-25]+merged[,2*i-24]
}
names(merged)[17:20] = c("hs12331_type","hs12334_type","hs12325_type","hs12328_type") 
condensed <- merged[,c(1,2,17,3,18,4,19,5,20)] 

# Plot data
# 4 plots, one for each type
ggplot(condensed,aes(x = CNhs12331)) + 
  facet_wrap(~hs12331_type,scales = "free") + 
  geom_histogram() +
  labs(title = "Histogram for Enhancer Types")

# Layered histogram
library(scales)
c = subset(condensed, CNhs12331>5)
attach(c)
ggplot(c, aes(x = CNhs12331)) + 
  geom_bar(data=subset(c, hs12331_type == '0'),aes(y = (..count..)/sum(..count..)),fill = "yellow", alpha = 0.2) + 
  geom_bar(data=subset(c, hs12331_type == '1'),aes(y = (..count..)/sum(..count..)),fill = "red", alpha = 0.2) + 
  geom_bar(data=subset(c, hs12331_type == '2'),aes(y = (..count..)/sum(..count..)),fill = "purple", alpha = 0.2) + 
  geom_bar(data=subset(c, hs12331_type == '3'),aes(y = (..count..)/sum(..count..)),fill = "black", alpha = 0.2) + 
  scale_y_continuous(labels = percent) +
  labs(ylab = "Percent") + 
  labs(xlab = "TPM")
detach(c)

# Scatter plot
attach(condensed)
plot(CNhs12331[hs12331_type==2])
plot(CNhs12331,col=c("yellow","red","purple","black")[hs12331_type+1],ylab="CNhs12331 TPM")
legend("topleft", pch = c(1,1), col=c("yellow","red","purple","black"), c("None","Weak","Strong","Weak+Strong"), bty="o",  box.col="darkgreen", cex=.8)
aggregate(CNhs12331, list(hs12331_type), mean)
aggregate(CNhs12334, list(hs12334_type), summary)
aggregate(CNhs12325, list(hs12325_type), summary)
aggregate(CNhs12328, list(hs12325_type), summary)
detach(condensed)
