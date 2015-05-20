################
#PRE-PROCESSING OF FIA DATA FOR PREDICTION IN STAN
################

#First, we will subset the TREE table such that:
#1. Only the most recent complete sampling cycle (2009-2013) is included
#2. Only live trees (STATUSCD==1) are included.
#3. Only trees >5 cm DBH are included.
tree_tab <- read.csv("MN_TREE.csv")
recent <- subset(tree_tab, INVYR >= 2009 & INVYR <= 2013) #Extract most recent FIA cycle
live.trees <- subset(recent, STATUSCD==1) #Extract only live trees

live.trees$DIA <- live.trees$DIA *2.54  #Convert diameter measurements to cm
live.trees$TPA_UNADJ <- live.trees$TPA_UNADJ * 2.47 #Convert the adjustment factors to trees/hectare

#trees <- subset(live.trees, DIA >= 5) #Remove all trees < 5 cm DBH from table
#write.csv(trees, "MN_2009_2013_5cm_dbh.csv") #Save this file
#data <- read.csv("MN_2009_2013_5cm_dbh.csv")

#We'll need the Jenkins species groups for this analysis, so we'll merge the REF_SPECIES table
#on to our subset of the TREES table.
#data <- read.csv("MN_2009_2013_5cm_dbh.csv")
ref.data <- read.csv("REF_SPECIES.csv") #Reference table with Jenkins et al. parameters
merge.dat <- merge(live.trees, ref.data, by="SPCD")
write.csv(merge.dat, "MN_2009_2013_jenkins.csv")  #This file will be passed to the Stan model for prediction in FIA

