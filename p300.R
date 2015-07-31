# Plots p300 overlap 

library(ggplot2)
library(data.table)
library(plyr)

fileList <- list.files(path="~/DIMACS/data/p300_data", pattern="*.txt")
data = lapply(paste("~/DIMACS/data/p300_data/",fileList, sep=""), read.delim)
data_c = lapply(data, function(df) cbind(state = rep(df$state[1],2),as.data.frame(table(df$p300))))
d <- rbindlist(data_c)
setnames(d,1:3, c('state','p300','freq'))
d$state <- factor(d$state, levels = c('E2','E1','E20','E3','E4','E5','E12','E6','E13','E14'))

# Plot
ggplot(d, aes(factor(state))) + geom_bar()
qplot(state, data=d, geom="bar", fill=factor(p300))
# same as qplot
ggplot(d, aes(state, fill=factor(p300))) + geom_bar() 
  + scale_fill_discrete(name="P300 Overlap") # rename legend

# Plot stacked bar chart with labels 
ggplot(d, aes(x = state, fill = factor(p300))) +
  geom_bar(stat = "identity", aes(y = freq), position = "dodge") + 
  geom_text(aes(x = state, y = freq, ymax = freq, label = freq, hjust = ifelse(as.numeric(p300)==1,1,0)),
            , position = "dodge", size = 3.5) + 
  scale_fill_discrete(name="P300 Overlap")

# Plot percentage overlap
po <- data.frame(state=character(),percentage=integer())
state <- d$state[c(F,T)]
percentage <- d$freq[c(F,T)]/(d$freq[c(T,F)]+d$freq[c(F,T)])
po <- data.frame(state=state,percentage=percentage)
po$group <- c(2, 3, 4, 4, 1, 2, 3, 3, 3, 4)

ggplot(po, aes(x = state)) +
  geom_bar(stat = "identity", aes(y = percentage, fill=factor(group))) +
  scale_fill_discrete(name="Enhancer\nGroup", 
                      breaks=c("1","2","3","4"),
                      labels=c("1 High-High", "2 Low-High", "3 High-Low", "4 Moderate-Low"))
