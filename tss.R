# Plots distances to TSS
library(ggplot2)
tss <- read.delim("~/DIMACS/GM12878_tss.txt")
tss_sample <- subset(tss, tss$state %in% c('E1','E2','E3','E4','E5','E12','E13','E14','E6','E20'))
tss_sample$state <- factor(tss_sample$state, levels = c('E2','E1','E20','E3','E4','E5','E12','E6','E13','E14'))
attach(tss_sample)

# Grouping
get_group <- function(state){
  if(state=='E2'){
    return(1)
  } else if(state=='E1' || state=='E20'){
    return(2)
  } else if(state=='E3' || state=='E4' || state=='E5'){
    return(3)
  } else{
    return(4)
  }
}

tss_sample$group <- sapply(tss_sample$state, get_group)

#Wilcoxon Test
pairwise.wilcox.test(distance, state, paired=FALSE, p.adjust="holm")

#Plotting 
ggplot(data = tss_sample, aes(x = distance, color = factor(state))) + geom_histogram()
#Box plots
ggplot(data = tss_sample, aes(x = state, y=distance, fill=factor(group))) + geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=5, size=3)


ord <- c(20,1,2,3,4,5,12,13,14,6)
# Means
means <- as.data.frame(aggregate(distance, list(state), mean))
means <- means[c(6,1,5,7,8,9,2,3,4,10),]
plot(y=means$x, x=1:nrow(means), xaxt='n', xlab='state', ylab='Mean Distance to TSS', col = 'red')
axis(1, at=1:nrow(means), labels=as.character(ord))
# Medians
medians <- as.data.frame(aggregate(distance, list(state), median))
medians <- means[c(6,1,5,7,8,9,2,3,4,10),]
points(y=medians$x, x=1:nrow(medians), xaxt='n', xlab='state', col = 'blue') #legend
