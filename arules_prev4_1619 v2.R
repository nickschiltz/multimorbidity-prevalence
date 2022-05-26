#==============================================================
# Title: Prevalence of Multimorbidity, 2016-2019
# Funder: KL2 and PhRMA
# Date: October 5, 2018
# Last modified: January 11, 2022
# Author: Nicholas K. Schiltz, PhD
# Description: This analysis identifies the most commonly
#   occuring combinations of chronic disease among US adults
#   in the MEPS database, 2016-19
#=============================================================

#load in packages
library("survey")
library("arules")
library("stringr")
library("tidyr")
library("svMisc")
library("data.table")
library("arulesViz")

#set working directory
setwd("C:/Users/nks8/Box/General_Reference/M/MEPS KL2 project/assoc_rules")
setwd("C:/Users/Nicholas/Box/General_Reference/M/MEPS KL2 project/assoc_rules")

#read in data from csv
datain <- read.csv("tx_chrcond_2016_2019_v1.csv")

#make copy - faster if need to start over
mydata <- datain 
#mydata <- subset(mydata, AGE < 45)
#mydata <- subset(mydata, AGE > 44 & AGE < 65)
#mydata <- subset(mydata, AGE > 64)

#summary of data
summary(mydata)
names(mydata)

#drop the n=1 cases with missing self-rated health
mydata <- mydata[complete.cases(mydata), ]

#convert categorical to factor
mydata[,c(5,9,15)] <- lapply(mydata[,c(5,9,15)],as.factor)
str(mydata)

#add a spline at 56 and older and 55 and under (based on MARS model)
mydata$AGE1859 <- ifelse(mydata$AGE<=59,59-mydata$AGE,0) 
mydata$AGE5985 <- ifelse(mydata$AGE>=59,mydata$AGE-59,0) 

## Remove "ccb_" from each variable's name
names(mydata) <- gsub("^.*?_","",names(mydata))
names(mydata)

#subset only the chronic conditions for arules
ccdata <- mydata[,c(17:36)]
ccdata <- mydata[,c(18:34,36)]

#sort them by order - will make items appear in consistent meaningful order for arules()
ccdata <- ccdata[,c(order(-colSums(ccdata)))]

#convert to transaction format
cctrans <- as(as.matrix(ccdata),"transactions")
inspect(head(cctrans, n=20))

#-------- Frequent Pattern Mining --------

#set parameters (make less stringent - will prune after applying survey weights)
params <- new("APparameter", support = 0.0025, ext=T, minlen=1, maxlen=10,target="frequent itemsets")
conts <- new("APcontrol", sort=-1, verbose = FALSE)

#run the frequent pattern mining
fpm1 <- apriori(cctrans, parameter=params, control=conts)

#add some measures and clean output
quality(fpm1) <- cbind(quality(fpm1),lift = interestMeasure(fpm1, measure="lift", transactions = cctrans))
quality(fpm1) <- round(quality(fpm1), digits=4)
summary(fpm1) #longest itemset is 6 conditions
inspect(head(fpm1,n=60))

#extract output to data frame
sto1 <- as(fpm1, "data.frame");
#clean item lists
sto1$items <- as.character(sto1$items)
sto1$items2 <- gsub('^.|.$','', sto1$items)
#make separate columns for each condition
sto2 <- separate(sto1, col=items2, into = c("CC1","CC2","CC3","CC4","CC5","CC6"), 
                  extra = "drop", fill = "right", sep = ",", remove=F)
#create variable counting number of conditions in each itemset
sto2$NumItems <- str_count(sto2$items2, ",") + 1
#create variable with combinations displayed between ampersands (nicer for publication)
sto2$Combination <- gsub(',',' & ', sto2$items2)

#sort by count
sto3 <- sto2[order(-sto2$count),]

#Begin process of flagging combinations in original data (mydata)

#make a dummy variable in mydata that is all 1 (needed for svytotal)
mydata$DUMMY <- 1
#make variables for sum of number of conditions and 0 condition flag (needed later)
mydata$numbercc <- rowSums(ccdata)
mydata$nocc <- ifelse(mydata$numbercc == 0,1,0)

#specify survey design for weighting
mepsdsgn = svydesign(id = ~VARPSU, strata = ~VARSTR, weights = ~AVGPERSWT,
                     data = mydata, nest = TRUE)

options(scipen=12)
a <- svyby(~DUMMY, ~numbercc, mepsdsgn, svytotal)
summary(a)

b <- svyby(~TOTEXP2019, ~HYPERTEN, mepsdsgn, quantile=c(0.25,0.5,0.75), svyquantile)
svyquantile(~TOTEXP2019, mepsdsgn, quantile=c(0.25,0.5,0.75), ci=TRUE)


#initialize values
Flag <- matrix(data = NA,nrow = nrow(mydata), ncol = nrow(sto2))
Label <- character(ncol(Flag))
WeightN <- numeric(ncol(Flag))
Wt_LCI_N <- numeric(ncol(Flag))
Wt_UCI_N <- numeric(ncol(Flag))
WeightN_LCI <- numeric(ncol(Flag))
WeightN_UCI <- numeric(ncol(Flag))
WeightPct <- numeric(ncol(Flag))
WeightPct_LCI <- numeric(ncol(Flag))
WeightPct_UCI <- numeric(ncol(Flag))
TotExpWt <- numeric(ncol(Flag))
TotExpWt_LCI <- numeric(ncol(Flag))
TotExpWt_UCI <- numeric(ncol(Flag))
AvgExpWt <- numeric(ncol(Flag))
AvgExpWt_LCI <- numeric(ncol(Flag))
AvgExpWt_UCI <- numeric(ncol(Flag))
MedExpWt <- numeric(ncol(Flag))
Q1ExpWt <- numeric(ncol(Flag))
Q3ExpWt <- numeric(ncol(Flag))
PoorHlthPct <- numeric(ncol(Flag))
PoorHlthPct_LCI <- numeric(ncol(Flag))
PoorHlthPct_UCI <- numeric(ncol(Flag))
# MedExpWt_LCI <- numeric(ncol(Flag))
# MedExpWt_UCI <- numeric(ncol(Flag))
# Q1ExpWt_LCI <- numeric(ncol(Flag))
# Q1ExpWt_UCI <- numeric(ncol(Flag))
# Q3ExpWt_LCI <- numeric(ncol(Flag))
# Q3ExpWt_UCI <- numeric(ncol(Flag)


#Loop for combinations  (this takes 12+ minutes)
start_time <- Sys.time()
for (i in 1:nrow(sto2))
{
  library(survey)
  if (sto2[i,"NumItems"] == 1){
     Flag[,i] <- ifelse(mydata[,sto2[i,"CC1"]]==1,1,0)
     Label[i] <- paste0(sto2[i,"CC1"]) 
  } else if (sto2[i,"NumItems"] == 2){
     Flag[,i] <- ifelse(mydata[,sto2[i,"CC1"]]==1 & mydata[,sto2[i,"CC2"]]==1,1,0)
     Label[i] <- paste0(sto2[i,"CC1"]," & ",sto2[i,"CC2"])       
  } else if (sto2[i,"NumItems"] == 3){
    Flag[,i] <- ifelse(mydata[,sto2[i,"CC1"]]==1 & mydata[,sto2[i,"CC2"]]==1 & mydata[,sto2[i,"CC3"]]==1,1,0)
    Label[i] <- paste0(sto2[i,"CC1"]," & ",sto2[i,"CC2"]," & ",sto2[i,"CC3"])     
  } else if (sto2[i,"NumItems"] == 4){
    Flag[,i] <- ifelse(mydata[,sto2[i,"CC1"]]==1 & mydata[,sto2[i,"CC2"]]==1 & mydata[,sto2[i,"CC3"]]==1 
                       & mydata[,sto2[i,"CC4"]]==1,1,0)
    Label[i] <- paste0(sto2[i,"CC1"]," & ",sto2[i,"CC2"]," & ",sto2[i,"CC3"]," & ",sto2[i,"CC4"])        
  } else if (sto2[i,"NumItems"] == 5){
    Flag[,i] <- ifelse(mydata[,sto2[i,"CC1"]]==1 & mydata[,sto2[i,"CC2"]]==1 & mydata[,sto2[i,"CC3"]]==1 
                       & mydata[,sto2[i,"CC4"]]==1 & mydata[,sto2[i,"CC5"]]==1,1,0)
    Label[i] <- paste0(sto2[i,"CC1"]," & ",sto2[i,"CC2"]," & ",sto2[i,"CC3"]," & ",sto2[i,"CC4"],
                       " & ",sto2[i,"CC5"])       
  } else if (sto2[i,"NumItems"] == 6){
    Flag[,i] <- ifelse(mydata[,sto2[i,"CC1"]]==1 & mydata[,sto2[i,"CC2"]]==1 & mydata[,sto2[i,"CC3"]]==1 
                       & mydata[,sto2[i,"CC4"]]==1 & mydata[,sto2[i,"CC5"]]==1 & mydata[,sto2[i,"CC6"]]==1,1,0)
    Label[i] <- paste0(sto2[i,"CC1"]," & ",sto2[i,"CC2"]," & ",sto2[i,"CC3"]," & ",sto2[i,"CC4"],
                       " & ",sto2[i,"CC5"]," & ",sto2[i,"CC6"])       
  }
  #mydata$Flag <- Indicator[,i]
  WeightN[i] <- coef(svytotal(~Flag[,i], mepsdsgn))
  WeightN_LCI[i] <- confint(svytotal(~Flag[,i], mepsdsgn, level=0.95))[,1]
  WeightN_UCI[i] <- confint(svytotal(~Flag[,i], mepsdsgn, level=0.95))[,2]
  weightpct <- svyciprop(~I(Flag[,i]==1), mepsdsgn, method="be")
  WeightPct[i] <- weightpct[1]
  WeightPct_LCI[i] <- confint(weightpct)[1]
  WeightPct_UCI[i] <- confint(weightpct)[2]
  totalexp <- svyby(~TOTEXP2019, ~Flag[,i], mepsdsgn, svytotal)
  TotExpWt[i] <- coef(totalexp[2,])
  TotExpWt_LCI[i] <- confint(totalexp[2,])[1]
  TotExpWt_UCI[i] <- confint(totalexp[2,])[2]
  means <- svyby(~TOTEXP2019, ~Flag[,i], mepsdsgn, svymean) 
  AvgExpWt[i] <- coef(means[2,])
  AvgExpWt_LCI[i] <- confint(means[2,])[1]
  AvgExpWt_UCI[i] <- confint(means[2,])[2]
  medians <- svyby(~TOTEXP2019, ~Flag[,i], mepsdsgn, quantile=c(0.25,0.5,0.75), svyquantile)
  MedExpWt[i] <- medians[2,3]
  #MedExpWt_LCI[i] <- medians[2,6]
  #MedExpWt_UCI[i] <- medians[2,9]
  Q1ExpWt[i] <- medians[2,2]
  #Q1ExpWt_LCI[i] <- medians[2,5]
  #Q1ExpWt_UCI[i] <- medians[2,8]
  Q3ExpWt[i] <- medians[2,4]
  #Q3ExpWt_LCI[i] <- medians[2,7]
  #Q3ExpWt_UCI[i] <- medians[2,10]
  poorhlth <- svyby(~POORHLTH, ~Flag[,i], mepsdsgn, svymean)
  PoorHlthPct[i] <- coef(poorhlth[2,])
  PoorHlthPct_LCI[i] <- confint(poorhlth[2,])[1]
  PoorHlthPct_UCI[i] <- confint(poorhlth[2,])[2]
  Sys.sleep(0.1)
  print(i)
}
print(Sys.time()-start_time)
rm(i)


#create a data frame
out1 <- data.frame(cbind(Label,WeightN,WeightN_LCI, WeightN_UCI, WeightPct, WeightPct_LCI, WeightPct_UCI,
                         TotExpWt, TotExpWt_LCI, TotExpWt_UCI, AvgExpWt, AvgExpWt_LCI, AvgExpWt_UCI,
                         MedExpWt, Q1ExpWt, Q3ExpWt, PoorHlthPct, PoorHlthPct_LCI, PoorHlthPct_UCI))
#rename variable
names(out1)[1] <- "Combination"

#make a copy and clean variables
out2 <- out1
out2$WeightN <- round(as.numeric(as.character(out2$WeightN)),0)
out2$WeightN_LCI <- round(as.numeric(as.character(out2$WeightN_LCI)),0)
out2$WeightN_UCI <- round(as.numeric(as.character(out2$WeightN_UCI)),0)
out2$WeightPct <- round(as.numeric(as.character(out2$WeightPct)),5)*100
out2$WeightPct_LCI <- round(as.numeric(as.character(out2$WeightPct_LCI)),5)*100
out2$WeightPct_UCI <- round(as.numeric(as.character(out2$WeightPct_UCI)),5)*100
out2$TotExpWt <- round(as.numeric(as.character(out2$TotExpWt)),0)
out2$TotExpWt_LCI <- round(as.numeric(as.character(out2$TotExpWt_LCI)),0)
out2$TotExpWt_UCI <- round(as.numeric(as.character(out2$TotExpWt_UCI)),0)
out2$AvgExpWt <- round(as.numeric(as.character(out2$AvgExpWt)),0)
out2$AvgExpWt_LCI <- round(as.numeric(as.character(out2$AvgExpWt_LCI)),0)
out2$AvgExpWt_UCI <- round(as.numeric(as.character(out2$AvgExpWt_UCI)),0)
out2$MedExpWt <- round(as.numeric(as.character(out2$MedExpWt)),0)
out2$Q1ExpWt <- round(as.numeric(as.character(out2$Q1ExpWt)),0)
out2$Q3ExpWt <- round(as.numeric(as.character(out2$Q3ExpWt)),0)
out2$WeightNmil <- round(out2$WeightN/1000000, 2)
out2$WeightNmil_LCI <- round(out2$WeightN_LCI/1000000, 2)
out2$WeightNmil_UCI <- round(out2$WeightN_UCI/1000000, 2)
out2$TotExpWtbil <- round(out2$TotExpWt/1000000000, 2)
out2$TotExpWtbil_LCI <- round(out2$TotExpWt_LCI/1000000000, 2)
out2$TotExpWtbil_UCI <- round(out2$TotExpWt_UCI/1000000000, 2)
out2$PoorHlthPct <- round(as.numeric(as.character(out2$PoorHlthPct)),5)*100
out2$PoorHlthPct_LCI <- round(as.numeric(as.character(out2$PoorHlthPct_LCI)),5)*100
out2$PoorHlthPct_UCI <- round(as.numeric(as.character(out2$PoorHlthPct_UCI)),5)*100

#merge weighted output with original arules output
out3 <- merge(sto2,out2, by.x = "Combination", sort=F)

#calculate 95% CI for lift
out3$lift_SE <- sqrt((1/out3$count) + (1/((out3$support*out3$lift)*nrow(mydata))) - (2/nrow(mydata)))
out3$Lift_LCI <- exp(log(out3$lift) - 1.96*out3$lift_SE)
out3$Lift_UCI <- exp(log(out3$lift) + 1.96*out3$lift_SE)

#calculate relative standard error - RSE over 0.30 might be unstable, RSE > 0.50 should not be reported 
out3$WeightNSE <- (out3$WeightN - out3$WeightN_LCI)/1.96
out3$WeightNRSE <- out3$WeightNSE/out3$WeightN
out3$WeightPSE <- (out3$WeightPct - out3$WeightPct_LCI)/1.96
out3$WeightPRSE <- out3$WeightPSE/out3$WeightPct
out3$TotExpSE <- (out3$TotExpWt - out3$TotExpWt_LCI)/1.96
out3$TotExpRSE <- out3$TotExpSE/out3$TotExpWt
out3$AvgExpSE <- (out3$AvgExpWt - out3$AvgExpWt_LCI)/1.96
out3$AvgExpRSE <- out3$AvgExpSE/out3$AvgExpWt

#no RSE were close to over 0.30 
# printing the highest value for each RSE variable to demonstrate
max(out3$WeightPRSE)
max(out3$TotExpRSE)
max(out3$AvgExpRSE)

#clean up names
out3$CleanName <- out3$Combination
out3$CleanName <- gsub("HYPERTEN","Hypertension",out3$CleanName)
out3$CleanName <- gsub("HYPERLIP","Hyperlipidemia",out3$CleanName)
out3$CleanName <- gsub("MUSCSKELETAL","Musculoskeletal",out3$CleanName)
out3$CleanName <- gsub("DIABETES","Diabetes",out3$CleanName)
out3$CleanName <- gsub("DEPR_ANXIETY","Depression/Anxiety",out3$CleanName)
out3$CleanName <- gsub("ARTHRITIS","Arthritis",out3$CleanName)
out3$CleanName <- gsub("ASTHCOPD","Asthma/COPD",out3$CleanName)
out3$CleanName <- gsub("URINARY","Urinary problems",out3$CleanName)
out3$CleanName <- gsub("STOMACH","Stomach problems",out3$CleanName)
out3$CleanName <- gsub("THYROID","Thyroid disorder",out3$CleanName)
out3$CleanName <- gsub("CARDIO","Cardiovascular disease",out3$CleanName)
out3$CleanName <- gsub("CANCER","Cancer",out3$CleanName)
out3$CleanName <- gsub("STROKE","Stroke",out3$CleanName)
out3$CleanName <- gsub("HEARTFAIL","Heart Failure",out3$CleanName)
out3$CleanName <- gsub("OSTEOPOROSIS","Osteoporosis",out3$CleanName)
out3$CleanName <- gsub("HEPATITIS","Hepatitis",out3$CleanName)
out3$CleanName <- gsub("DEMENTIA","Dementia",out3$CleanName)
out3$CleanName <- gsub("OBESITY","Obesity",out3$CleanName)
out3$CleanName <- gsub("COLON","Colon problems",out3$CleanName)
out3$CleanName <- gsub("KIDNEY","Kidney disease",out3$CleanName)

#calculate "improvement" for average expenses
nocc <- coef(svymean(~TOTEXP2019, mepsdsgn))
MaxExpSubset <- numeric(nrow(out3))
ImpExp_Pct <- numeric(nrow(out3))
for (i in 1:nrow(out3)){
  if (out3[i,"NumItems"] == 1){
    MaxExpSubset[i] = nocc
    ImpExp_Pct[i] <- (round((out3[i,"AvgExpWt"] - MaxExpSubset[i])/(MaxExpSubset[i]),4))*100
  }
  else if (out3[i,"NumItems"] > 1){
    a <- as(fpm1[is.superset(fpm1[i],fpm1, sparse=F)], "data.frame")
    a <- as.character(a$items)[1:(nrow(a)-1)]
    b <- subset(out3, items %in% a)
    MaxExpSubset[i] <- max(b$AvgExpWt)
    ImpExp_Pct[i] <- (round((out3[i,"AvgExpWt"] - MaxExpSubset[i])/(MaxExpSubset[i]),6))*100
  }
}
out4 <- cbind(out3, MaxExpSubset, ImpExp_Pct)
summary(out4$ImpExp_Pct)

#output all - by combo length
out5 <- out4[order(out4$NumItems, -out4$WeightN),]
write.table(out5, file="./treatedprev1619/arules_prev_byNumItem_v2.csv", row.names=FALSE, sep=",")

#output all - by total weighted N
out5 <- out4[order(-out4$WeightN),]
write.table(out5, file="./treatedprev1619/arules_prev_byWtN_v2.csv", row.names=FALSE, sep=",")

#total number of combos by num items with over 1 million
out5 <- subset(out5, WeightN >= 1000000 & NumItems > 0)
out5 <- out5[order(out5$NumItems, -out5$WeightN),]
table(out5$NumItems)
write.table(out5, file="./treatedprev1619/arules_prev_byWtN_1mil.csv", row.names=FALSE, sep=",")

#output only 1 million plus and only numitmes >= 2
out5 <- subset(out5, WeightN >= 1000000 & NumItems > 1)
out5 <- out5[order(-out5$WeightN),]
write.table(out5, file="./treatedprev1619/arules_prev_byWtN_1mil_2L.csv", row.names=FALSE, sep=",")

#output table 2 (cleaner and only NumItems >1)
out5t <- out5[,c("CleanName","WeightPct","WeightPct_LCI","WeightPct_UCI","WeightN",
                       "WeightN_LCI","WeightN_UCI","count","lift","Lift_LCI","Lift_UCI",
                       "NumItems","TotExpWtbil","TotExpWtbil_LCI","TotExpWtbil_UCI")]
write.table(out5t, file="./treatedprev1619/table1.csv", row.names=FALSE, sep=",")

#total number of combos by num items with over 1 million
table(out5t$NumItems)

###-not used
#output table 3 - total expenditures
out6 <- subset(out4, WeightN >= 1000000)
out6 <- out6[order(-out6$TotExpWt),]
out6 <- out6[,c("CleanName" , "NumItems" , "WeightN" , "WeightPct", "count" , "TotExpWt" , 
                "TotExpWt_LCI" , "TotExpWt_UCI", "TotExpWtbil" , "TotExpWtbil_LCI" , 
                "TotExpWtbil_UCI" , "AvgExpWt" ,  "MedExpWt" )]
write.table(out6, file="./treatedprev1619/table3a.csv", row.names=FALSE, sep=",")
out6.2 <- head(subset(out6, out5$NumItems==2), n=10)
out6.3 <- head(subset(out6, out5$NumItems==3), n=10)
out6.4 <- head(subset(out6, out5$NumItems==4), n=10)
out6.t <- rbind(out6.2, out6.3, out6.4)
write.table(out6.t, file="./treatedprev1619/table3a_byItems.csv", row.names=FALSE, sep=",")

#output table 3b - avg expenditures
out7 <- subset(out4, WeightN >= 1000000 & NumItems > 1)
out7 <- out7[order(-out7$AvgExpWt),]
out7 <- out7[,c("CleanName" , "NumItems" , "WeightN" , "WeightPct", "count" , "AvgExpWt" , 
                "AvgExpWt_LCI" , "AvgExpWt_UCI", "TotExpWtbil" , "TotExpWt" , 
                "MedExpWt" , "Q1ExpWt" ,  "Q3ExpWt", "MaxExpSubset", "ImpExp_Pct" )]
write.table(out7, file="./treatedprev1619/table3b.csv", row.names=FALSE, sep=",")
#### end not used

#output only those that add to mean by 10% or more over simpler itemset
#use this one for final paper
out8 <- subset(out7, ImpExp_Pct >= 10 & NumItems > 1)
write.table(out8, file="./treatedprev1619/table3.csv", row.names=FALSE, sep=",")

#output all - by combo length and Weight N
out9 <- subset(out4, out4$WeightN >= 1000000)
out9 <- out9[order(out9$NumItems, -out9$WeightN),]
write.table(out9, file="./treatedprev1619/master_alldata.csv", row.names=FALSE, sep=",")



##### - Association Rule Mining #########



#make a copy of the data
mydata2 <- mydata

#subset only the chronic conditions and poor health for arules
ccdata2 <- mydata2[,c(16:36)]
head(ccdata2)

#convert to transaction format
cctrans2 <- as(as.matrix(ccdata2),"transactions")

#sort columns by item frequency - when arules are made, highest prev items will be listed first
ccdata2 <- ccdata2[,c(order(-itemFrequency(cctrans2)))]

#re-convert to transaction format
cctrans2 <- as(as.matrix(ccdata2),"transactions")
inspect(head(cctrans2, n=20))

#set parameters (make less stringent - will prune after applying survey weights)
params <- new("APparameter", support = 0.0025, confidence = .10, ext=T, minlen=1, maxlen=10,target="rules", 
              ext=T, originalSupport=F)
appear <- list(rhs = c("POORHLTH"))
conts <- new("APcontrol", sort=-1, verbose = FALSE)

#run the associaiton rule mining
arm1 <- apriori(cctrans2, parameter=params, appearance = appear, control=conts)
#arm1 <- apriori(cctrans2, parameter=params, control=conts)

#add some measures and clean output
Measures <- interestMeasure(arm1, c("improvement","boost","oddsRatio","phi","hyperLift","hyperConfidence","coverage",
                                     "conviction","cosine","doc","gini","leverage","RLD","addedValue",
                                    "chiSquared","certainty","collectiveStrength","confirmedConfidence",
                                    "casualConfidence","casualSupport","counterexample",
                                    "imbalance","implicationIndex","importance","jaccard", "jMeasure",
                                    "kappa","lambda","laplace","leastContradiction",
                                    "lerman","maxConfidence","mutualInformation","rulePowerFactor",
                                    "sebag","varyingLiaison","yuleQ","yuleY"), transactions=cctrans2)

quality(arm1) <- cbind(quality(arm1), Measures)
quality(arm1) <- round(quality(arm1), digits=4)
summary(arm1) #longest itemset is 6 conditions
inspect(head(arm1,n=60))
#extract output to data frame
sto3 <- DATAFRAME(arm1, separate=T, setStart= '', itemSep = ' & ', setEnd = '')
#sto3 <- as(arm1, "data.frame");

#clean item lists
sto3$LHSchar <- as.character(sto3$LHS)
sto4 <- separate(sto3, col=LHSchar, into = c("CC1","CC2","CC3","CC4","CC5","CC6"), 
                 extra = "drop", fill = "right", sep = " & ", remove=F)
#create variable counting number of conditions in each itemset
sto4$NumItems <- str_count(sto4$LHSchar, "&") + 1
#remove the NULL set (row 1)
sto4 <- sto4[c(-1),]

#explore association rules
arm2 <- apriori(cctrans2, parameter=params, appearance = appear, control=conts)
#ruleExplorer(arm2)

b <- rules2groupedMatrix(arm2, k=5)
b$m
table(b$clustering_rules)
inspect(arm2[b$clustering_rules==1])
plot(arm2, method = "grouped matrix", k = 5)

g <- associations2igraph(arm2)
g

plot(g)

#explore frequent patterns
#fpm2 <- apriori(cctrans, parameter=params, control=conts)
#ruleExplorer(cctrans)
#inspectDT(fpm1)

#plot(arm1[2:545], measure="support", shading = "oddsRatio", interactive = T)

#Begin process of flagging combinations in original data (mydata)

#make a dummy variable in mydata that is all 1 (needed for svytotal)
mydata2$DUMMY <- 1

#specify survey design for weighting
mepsdsgn2 = svydesign(id = ~VARPSU, strata = ~VARSTR, weights = ~AVGPERSWT,
                     data = mydata2, nest = TRUE)

svytotal(~DUMMY, mepsdsgn)

#initialize values
Flag <- matrix(data = NA,nrow = nrow(mydata2), ncol = nrow(sto4))
Label <- character(ncol(Flag))
WeightN <- numeric(ncol(Flag))
WeightN_LCI <- numeric(ncol(Flag))
WeightN_UCI <- numeric(ncol(Flag))
WeightPct <- numeric(ncol(Flag))
WeightPct_LCI <- numeric(ncol(Flag))
WeightPct_UCI <- numeric(ncol(Flag))
aOR <-  numeric(ncol(Flag))
aOR_LCI <- numeric(ncol(Flag))
aOR_UCI <- numeric(ncol(Flag))
aOR_pval <- numeric(ncol(Flag))
OR <-  numeric(ncol(Flag))
OR_LCI <- numeric(ncol(Flag))
OR_UCI <- numeric(ncol(Flag))
OR_pval <- numeric(ncol(Flag))


#Loop for combinations  (this takes 30+ minutes)
start_time <- Sys.time()
#for (i in 1:5)
for (i in 1:nrow(sto4))
{
  if (sto4[i,"NumItems"] == 1){
    Flag[,i] <- ifelse(mydata2[,sto4[i,"CC1"]]==1,1,0)
    Label[i] <- paste0(sto4[i,"CC1"]) 
  } else if (sto4[i,"NumItems"] == 2){
    Flag[,i] <- ifelse(mydata2[,sto4[i,"CC1"]]==1 & mydata2[,sto4[i,"CC2"]]==1,1,0)
    Label[i] <- paste0(sto4[i,"CC1"]," & ",sto4[i,"CC2"])       
  } else if (sto4[i,"NumItems"] == 3){
    Flag[,i] <- ifelse(mydata2[,sto4[i,"CC1"]]==1 & mydata2[,sto4[i,"CC2"]]==1 & mydata2[,sto4[i,"CC3"]]==1,1,0)
    Label[i] <- paste0(sto4[i,"CC1"]," & ",sto4[i,"CC2"]," & ",sto4[i,"CC3"])     
  } else if (sto4[i,"NumItems"] == 4){
    Flag[,i] <- ifelse(mydata2[,sto4[i,"CC1"]]==1 & mydata2[,sto4[i,"CC2"]]==1 & mydata2[,sto4[i,"CC3"]]==1 
                       & mydata2[,sto4[i,"CC4"]]==1,1,0)
    Label[i] <- paste0(sto4[i,"CC1"]," & ",sto4[i,"CC2"]," & ",sto4[i,"CC3"]," & ",sto4[i,"CC4"])        
  } else if (sto4[i,"NumItems"] == 5){
    Flag[,i] <- ifelse(mydata2[,sto4[i,"CC1"]]==1 & mydata2[,sto4[i,"CC2"]]==1 & mydata2[,sto4[i,"CC3"]]==1 
                       & mydata2[,sto4[i,"CC4"]]==1 & mydata2[,sto4[i,"CC5"]]==1,1,0)
    Label[i] <- paste0(sto4[i,"CC1"]," & ",sto4[i,"CC2"]," & ",sto4[i,"CC3"]," & ",sto4[i,"CC4"],
                       " & ",sto4[i,"CC5"])       
  } else if (sto4[i,"NumItems"] == 6){
    Flag[,i] <- ifelse(mydata2[,sto4[i,"CC1"]]==1 & mydata2[,sto4[i,"CC2"]]==1 & mydata2[,sto4[i,"CC3"]]==1 
                       & mydata2[,sto4[i,"CC4"]]==1 & mydata2[,sto4[i,"CC5"]]==1 & mydata2[,sto4[i,"CC6"]]==1,1,0)
    Label[i] <- paste0(sto4[i,"CC1"]," & ",sto4[i,"CC2"]," & ",sto4[i,"CC3"]," & ",sto4[i,"CC4"],
                       " & ",sto4[i,"CC5"]," & ",sto4[i,"CC6"])       
  }
  mydata2$indicator <- Flag[,i] 
  mepsdsgn2 = svydesign(id = ~VARPSU, strata = ~VARSTR, weights = ~AVGPERSWT,
                       data = mydata2, nest = TRUE)
  Model1 <- svyglm(POORHLTH ~ indicator + AGE1859 + AGE5985 + FEMALE + RACETHX, design=mepsdsgn2,
                   family=quasibinomial())
  Summ1 <- summary(Model1)
  aOR[i] <-  exp(Summ1$coefficients[2])
  aOR_LCI[i] <- exp(Summ1$coefficients[2] - 1.96*Summ1$coefficients[2,2])
  aOR_UCI[i] <- exp(Summ1$coefficients[2] + 1.96*Summ1$coefficients[2,2])
  aOR_pval[i] <- Summ1$coefficients[2,4]
  Model2 <- svyglm(POORHLTH ~ indicator, design=mepsdsgn2, family=quasibinomial())
  Summ2 <- summary(Model2)
  OR[i] <-  exp(Summ2$coefficients[2])
  OR_LCI[i] <- exp(Summ2$coefficients[2] - 1.96*Summ2$coefficients[2,2])
  OR_UCI[i] <- exp(Summ2$coefficients[2] + 1.96*Summ2$coefficients[2,2])
  OR_pval[i] <- Summ2$coefficients[2,4]
  WeightN[i] <- coef(svytotal(~Flag[,i], mepsdsgn2))
  WeightN_LCI[i] <- confint(svytotal(~Flag[,i], mepsdsgn2, level=0.95))[,1]
  WeightN_UCI[i] <- confint(svytotal(~Flag[,i], mepsdsgn2, level=0.95))[,2]
  weightpct <- svyciprop(~I(Flag[,i]==1), mepsdsgn2, method="be")
  WeightPct[i] <- weightpct[1]
  WeightPct_LCI[i] <- confint(weightpct)[1]
  WeightPct_UCI[i] <- confint(weightpct)[2]
  Sys.sleep(0.1)
  print(i)
}

end_time <- Sys.time()

end_time - start_time

#create a data frame
out11 <- data.frame(cbind(Label,WeightN,WeightN_LCI, WeightN_UCI, WeightPct, WeightPct_LCI, WeightPct_UCI,
                      OR, OR_LCI, OR_UCI, OR_pval, aOR, aOR_LCI, aOR_UCI, aOR_pval))
write.table(out11, file="./treatedprev1619/out11.csv", row.names=FALSE, sep=",")

#rename variable
names(out11)[1] <- "LHSchar"

#make a copy and clean variables
out12 <- out11
out12$WeightN <- round(as.numeric(as.character(out12$WeightN)),0)
out12$WeightN_LCI <- round(as.numeric(as.character(out12$WeightN_LCI)),0)
out12$WeightN_UCI <- round(as.numeric(as.character(out12$WeightN_UCI)),0)
out12$WeightPct <- round(as.numeric(as.character(out12$WeightPct)),5)*100
out12$WeightPct_LCI <- round(as.numeric(as.character(out12$WeightPct_LCI)),5)*100
out12$WeightPct_UCI <- round(as.numeric(as.character(out12$WeightPct_UCI)),5)*100

out12$aOR <- round(as.numeric(as.character(out12$aOR)),4)
out12$aOR_LCI <- round(as.numeric(as.character(out12$aOR_LCI)),4)
out12$aOR_UCI <- round(as.numeric(as.character(out12$aOR_UCI)),4)
out12$aOR_pval <- round(as.numeric(as.character(out12$aOR_pval)),6)
out12$OR <- round(as.numeric(as.character(out12$OR)),4)
out12$OR_LCI <- round(as.numeric(as.character(out12$OR_LCI)),4)
out12$OR_UCI <- round(as.numeric(as.character(out12$OR_UCI)),4)
out12$OR_pval <- round(as.numeric(as.character(out12$OR_pval)),6)

#merge weighted output with original
out13 <- merge(sto4,out12, by.x = "LHSchar")

#calculate improvement percentage
out13$improve_pct <- out13$improvement / out13$confidence

#calculate 95% CI for lift
out13$lift_SE <- sqrt((1/out13$count) + (1/((out13$support*out13$lift)*nrow(mydata2))) - (2/nrow(mydata2)))
out13$Lift_LCI <- exp(log(out13$lift) - 1.96*out13$lift_SE)
out13$Lift_UCI <- exp(log(out13$lift) + 1.96*out13$lift_SE)

#calculate relative standard error - RSE over 0.30 might be unstable, RSE > 0.50 should not be reported 
out13$aORSE <- (out13$aOR - out13$aOR_LCI)/1.96
out13$aORRSE <- out13$aORSE/out13$aOR
out13$ORSE <- (out13$OR - out13$OR_LCI)/1.96
out13$ORRSE <- out13$ORSE/out13$OR
#no RSE were close to over 0.30 
max(out13$aORRSE)
max(out13$ORRSE)

#clean up names
out13$CleanName <- out13$LHSchar
out13$CleanName <- gsub("HYPERTEN","Hypertension",out13$CleanName)
out13$CleanName <- gsub("HYPERLIP","Hyperlipidemia",out13$CleanName)
out13$CleanName <- gsub("MUSCSKELETAL","Musculoskeletal",out13$CleanName)
out13$CleanName <- gsub("DIABETES","Diabetes",out13$CleanName)
out13$CleanName <- gsub("DEPR_ANXIETY","Depression/Anxiety",out13$CleanName)
out13$CleanName <- gsub("ARTHRITIS","Arthritis",out13$CleanName)
out13$CleanName <- gsub("ASTHCOPD","Asthma/COPD",out13$CleanName)
out13$CleanName <- gsub("URINARY","Urinary problems",out13$CleanName)
out13$CleanName <- gsub("STOMACH","Stomach problems",out13$CleanName)
out13$CleanName <- gsub("THYROID","Thyroid disorder",out13$CleanName)
out13$CleanName <- gsub("CARDIO","Cardiovascular disease",out13$CleanName)
out13$CleanName <- gsub("CANCER","Cancer",out13$CleanName)
out13$CleanName <- gsub("STROKE","Stroke",out13$CleanName)
out13$CleanName <- gsub("HEARTFAIL","Heart Failure",out13$CleanName)
out13$CleanName <- gsub("OSTEOPOROSIS","Osteoporosis",out13$CleanName)
out13$CleanName <- gsub("HEPATITIS","Hepatitis",out13$CleanName)
out13$CleanName <- gsub("DEMENTIA","Dementia",out13$CleanName)
out13$CleanName <- gsub("OBESITY","Obesity",out13$CleanName)
out13$CleanName <- gsub("COLON","Colon problems",out13$CleanName)
out13$CleanName <- gsub("KIDNEY","Kidney disease",out13$CleanName)

#output all - by combo length
out14 <- out13[order(out13$NumItems, -out13$aOR),]
write.table(out14, file="./treatedprev1619/table4_alldata.csv", row.names=FALSE, sep=",")

#output only over 1 million, only with improve_pct >= 0.10, and sorted by aOR
out15 <- subset(out13, out13$WeightN >= 1000000 & NumItems > 1) #drops 58 rules
out15 <- subset(out15, out15$improve_pct >= 0.10) #drops 197 rules
out15 <- out15[order(-out15$aOR),]
out15 <- out15[,c(76, 54:56, 8, 5, 68, 9, 10, 60:67, 11, 53)]
write.table(out15, file="./treatedprev1619/table4.csv", row.names=FALSE, sep=",")

#output all - by combo length and Weight N
out16 <- subset(out13, out13$WeightN >= 1000000)
out17 <- out16[order(out16$NumItems, -out16$WeightN),]
write.table(out17, file="./treatedprev1619/master_aordata.csv", row.names=FALSE, sep=",")

#save image
save.image()

#### - Population Characteristics - Weighted #########

varlist <- c("Agecat1","FEMALE","RACETHX", "SRHEALTH", "numcccat1")

for (i in 1:dim(varlist)){
  
  
}
table(mydata$DUMMY)
coef(svytotal(~DUMMY, mepsdsgn))

table(mydata$Agecat1)
a <- coef(svytotal(~Agecat1, mepsdsgn))[1]
svyciprop(~I(Flag[,i]==1), mepsdsgn2, method="be")
svyciprop(~I(FEMALE==0), mepsdsgn, method="be")



coef(svytotal(~Agecat1, mepsdsgn))[1]
svyciprop(~I(Agecat1=="[18,40)"), mepsdsgn, method="be", level=0.95)
confint(svyciprop(~I(Agecat1=="[18,40)"), mepsdsgn, method="be", level=0.95))[,1]
confint(svyciprop(~I(Agecat1=="[18,40)"), mepsdsgn, method="be", level=0.95))[,2]


weightpct <- svyciprop(~I(Flag[,i]==1), mepsdsgn, method="be")
WeightPct[i] <- weightpct[1]
WeightPct_LCI[i] <- confint(weightpct)[1]
WeightPct_UCI[i] <- confint(weightpct)[2]


library(tableone)
#specify survey design for weighting
mydata3 <- mydata
mydata3$Agecat1<-cut(mydata$AGE, c(18,40,65,86), right = F)
mydata3$numcccat1<-cut(mydata$numbercc, c(0,1,2,3,4,5,6,30), right = F)
mydata3$mmcat1 <- cut(mydata$numbercc, c(0,1,2,30), right=F)
mydata3$MM0_1 <- cut(mydata$numbercc, c(0,1,30), right=F)
mydata3$MM01_2 <- cut(mydata$numbercc, c(0,2,30), right=F)
mydata3$MM02_3 <- cut(mydata$numbercc, c(0,3,30), right=F)
mydata3$MM03_4 <- cut(mydata$numbercc, c(0,4,30), right=F)
mydata3$MM04_5 <- cut(mydata$numbercc, c(0,5,30), right=F)
mydata3$TOPQUART <- ifelse(mydata3$TOTEXP2019 >= quantile(mydata3$TOTEXP2019, 0.75, na.rm=T), 1, 0)
#quantile(mydata$TOTEXP2015, prob=(c(0,.75,.90,.95,.99,.999,.9999,.99999,1)))
#mean(mydata$TOTEXP2015)

library(gtsummary)

mepsdsgn = svydesign(id = ~VARPSU, strata = ~VARSTR, weights = ~AVGPERSWT,
                     data = mydata3, nest = TRUE)
tab1 <- svyCreateTableOne(vars = c("Agecat1","FEMALE","RACETHX", "POORHLTH", "TOPQUART",
                                   "mmcat1", "MM0_1","MM01_2","MM02_3","MM03_4","MM04_5"),
                  strata = "MM0_1", data = mepsdsgn, 
                  factorVars = c("FEMALE","POORHLTH","TOPQUART","RACETHX"))

tab1d <- print(tab1, showAllLevels=T, quote=T, format="pf")
write.table(tab1d, file="./treatedprev1619/eTab1.csv", row.names=FALSE)

tbl_svysummary(data = mepsdsgn, by = MM0_1, percent = "row", digits= list(c("Agecat1","FEMALE","RACETHX", "POORHLTH") ~ c(0,1)), 
               include = c("Agecat1","FEMALE","RACETHX", "POORHLTH", "TOTEXP2019"), type=list(FEMALE ~ "categorical"))

tbl_svysummary(data = mepsdsgn, by = MM01_2, percent = "row", digits= list(c("Agecat1","FEMALE","RACETHX", "POORHLTH") ~ c(0,1)), 
               include = c("Agecat1","FEMALE","RACETHX", "POORHLTH", "TOTEXP2019"), type=list(FEMALE ~ "categorical"))

tbl_svysummary(data = mepsdsgn, by = MM02_3, percent = "row", digits= list(c("Agecat1","FEMALE","RACETHX", "POORHLTH") ~ c(0,1)), 
               include = c("Agecat1","FEMALE","RACETHX", "POORHLTH", "TOTEXP2019"), type=list(FEMALE ~ "categorical"))

tbl_svysummary(data = mepsdsgn, by = MM03_4, percent = "row", digits= list(c("Agecat1","FEMALE","RACETHX", "POORHLTH") ~ c(0,1)), 
               include = c("Agecat1","FEMALE","RACETHX", "POORHLTH", "TOTEXP2019"), type=list(FEMALE ~ "categorical"))

table(mydata3$DUMMY)
table(mydata3$Agecat1)
table(mydata3$FEMALE)
table(mydata3$RACETHX)
table(mydata3$SRHEALTH)
table(mydata3$numcccat1)
table(mydata3$MM0_1)
table(mydata3$MM01_2)
table(mydata3$MM02_3)
table(mydata3$MM03_4)
table(mydata3$MM04_5)


svyby(~TOTEXP2019, ~numcccat1, mepsdsgn, svymean) 
svyby(~mmcat1, ~Agecat1, mepsdsgn, svytotal)

svyby(~TOTEXP2019, ~MM02_3, mepsdsgn, svymean) 
svyby(~POORHLTH, ~MM02_3, mepsdsgn, svymean) 




########### END PROGRAM ##########

#--------------Upset Plot---------------------------

library(scales)
library(gridExtra)
library(grid)

combosdata <- read.csv(file="./treatedprev1619/arules_prev_byWtN_1mil.csv")
ordata <- read.csv( file="./treatedprev1619/master_aordata.csv")

combosdata2 <- combosdata %>%
          select(CleanName, )

combosdata2 <- left_join(combosdata, ordata, by="CleanName")

ordata2 <- read.csv( file="./treatedprev1619/table4.csv")


#version 2
p1 <- combosdata %>%
  filter(NumItems > 1 & WeightN > 5000000) %>%
  #slice_max(WeightN, n=30) %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN, .fun='median' )) %>%
  ggplot(aes(x = WeightNmil, y = CleanName)) +
  geom_point(size = 3.5, color = "#FF1F5B") +
  geom_errorbarh(aes(xmax = WeightNmil_UCI, xmin = WeightNmil_LCI), size = .5, height = 
                   .2, color = "gray50") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  ylab("") +
  xlab("Prevalence (Millions)")

p2 <- combosdata %>%
  filter(NumItems > 1 & WeightN > 5000000) %>%
  #slice_max(WeightN, n=30) %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN, .fun='median' )) %>%
  ggplot(aes(x = TotExpWtbil, y = CleanName)) +
  geom_point(size = 3.5, color = "#009ADE") +
  geom_errorbarh(aes(xmax = TotExpWtbil_UCI, xmin = TotExpWtbil_LCI), size = .5, height = 
                   .2, color = "gray50") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
  ylab("") +
  xlab("Cost (Billions USD)")

p3 <- combosdata %>%
    filter(NumItems > 1 & WeightN > 5000000) %>%
    #slice_max(WeightN, n=30) %>%
    mutate(CleanName = fct_reorder(CleanName, WeightN, .fun='median' )) %>%
    ggplot(aes(x = PoorHlthPct, y = CleanName)) +
    geom_point(size = 3.5, color = "#FFC61E") +
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = PoorHlthPct_UCI, xmin = PoorHlthPct_LCI), size = .5, height = 
                     .2, color = "gray50") +
    #scale_x_log10(breaks = pretty_breaks(n = 5)) +
    theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
    ylab("") +
    xlab("Percent with poor self-rated health")


grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))

pdf("figures/Figure1_topPrev.pdf", width=11,height=6)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))
dev.off()


# Those in top 20 of poor health or average expense
topdat1 <- combosdata2 %>%
  filter(NumItems.x > 1 & improve_pct > 0.10) %>%
  slice_max(PoorHlthPct, n=20)

topdat2 <- combosdata2 %>%
  filter(NumItems.x > 1 & improve_pct > 0.10) %>%
  slice_max(AvgExpWt, n=20) %>%
  select(CleanName)

topdat3 <- full_join(topdat1, topdat2, by="CleanName") %>%
            select(CleanName)

topdat4 <- left_join(topdat3, combosdata2, by="CleanName")

#not used - just to see how many are in both criteria (12)
#topdat5 <- inner_join(topdat1, topdat2, by="CleanName")

p4 <- topdat4 %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN.x, .fun='median' )) %>%
  ggplot(aes(x = WeightNmil, y = CleanName)) +
  geom_point(size = 3.5, color = "#FF1F5B") +
  geom_errorbarh(aes(xmax = WeightNmil_UCI, xmin = WeightNmil_LCI), size = .5, height = 
                   .2, color = "gray25") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  ylab("") +
  xlab("Prevalence (Millions)")

p5 <- topdat4 %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN.x, .fun='median' )) %>%
  ggplot(aes(x = AvgExpWt, y = CleanName)) +
  geom_point(size = 3.5, color = "#009ADE") +
  geom_errorbarh(aes(xmax = AvgExpWt_UCI, xmin = AvgExpWt_LCI), size = .5, height = 
                   .2, color = "gray25") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
  ylab("") +
  xlab("Per Capita Cost (USD)")

p6 <- topdat4 %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN.x, .fun='median' )) %>%
  ggplot(aes(x = aOR, y = CleanName)) +
  geom_point(size = 3.5, color = "#FFC61E") +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = aOR_UCI, xmin = aOR_LCI), size = .5, height = 
                   .2, color = "gray25") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
  ylab("") +
  xlab("Adjusted odds of poor health")

grid.newpage()
grid.draw(cbind(ggplotGrob(p4), ggplotGrob(p5), ggplotGrob(p6), size = "last"))

pdf("figures/Figure2_topImpact.pdf", width=11,height=6)
grid.newpage()
grid.draw(cbind(ggplotGrob(p4), ggplotGrob(p5), ggplotGrob(p6), size = "last"))
dev.off()


library(stringr)
nohyper <- combosdata2 %>%
            filter(!str_detect(CleanName, 'Hypertension')) %>%
            filter(!str_detect(CleanName, 'Hyperlipidemia'))

#version 2
p1 <- nohyper %>%
  filter(NumItems.x > 1) %>%
  slice_max(WeightN.x, n=25) %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN.x, .fun='median' )) %>%
  ggplot(aes(x = WeightNmil, y = CleanName)) +
  geom_point(size = 3.5, color = "#FF1F5B") +
  geom_errorbarh(aes(xmax = WeightNmil_UCI, xmin = WeightNmil_LCI), size = .5, height = 
                   .2, color = "gray50") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  ylab("") +
  xlab("Prevalence (Millions)")

p2 <- nohyper %>%
  filter(NumItems.x > 1) %>%
  slice_max(WeightN.x, n=25) %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN.x, .fun='median' )) %>%
  ggplot(aes(x = TotExpWtbil, y = CleanName)) +
  geom_point(size = 3.5, color = "#009ADE") +
  geom_errorbarh(aes(xmax = TotExpWtbil_UCI, xmin = TotExpWtbil_LCI), size = .5, height = 
                   .2, color = "gray50") +
  scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
  ylab("") +
  xlab("Cost (Billions USD)")

p3 <- nohyper %>%
  filter(NumItems.x > 1) %>%
  slice_max(WeightN.x, n=25) %>%
  mutate(CleanName = fct_reorder(CleanName, WeightN.x, .fun='median' )) %>%
  ggplot(aes(x = PoorHlthPct, y = CleanName)) +
  geom_point(size = 3.5, color = "#FFC61E") +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = PoorHlthPct_UCI, xmin = PoorHlthPct_LCI), size = .5, height = 
                   .2, color = "gray50") +
  #scale_x_log10(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
  ylab("") +
  xlab("Percent with poor self-rated health")

grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))

pdf("figures/Figure_noHyper_topPrev.pdf", width=11,height=6)
grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))
dev.off()




#-----Not Used----
# Those in top 25 of high improvement over single 

impdat1 <- combosdata2 %>%
  filter(NumItems.x == 2 & improve_pct > 0.10) %>%
  slice_max(improve_pct, n=25)

impdat2 <- combosdata2 %>%
  mutate(ImpExp_Tot = AvgExpWt - MaxExpSubset) %>%
  filter(NumItems.x ==2) %>%
  slice_max(ImpExp_Tot, n=25) %>%
  select(CleanName, ImpExp_Tot, ImpExp_Pct, AvgExpWt, MaxExpSubset)

impdat3 <- inner_join(impdat1, impdat2, by="CleanName")

p7 <- impdat3 %>%
  mutate(CleanName = fct_reorder(CleanName, ImpExp_Tot, .fun='median' )) %>%
  ggplot(aes(x = ImpExp_Tot, y = CleanName)) +
  geom_bar(stat="identity", fill = "#009ADE") +
  #scale_x(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Additional Cost (USD)")

p8 <- impdat3 %>%
  mutate(CleanName = fct_reorder(CleanName, ImpExp_Tot, .fun='median' )) %>%
  ggplot(aes(x = improve_pct, y = CleanName)) +
  geom_bar(stat="identity", fill = "#FFC61E") +
  #scale_x(breaks = pretty_breaks(n = 5)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text.y=element_blank()) +
  ylab("") +
  xlab("Percent poor health")

grid.newpage()
grid.draw(cbind(ggplotGrob(p7), ggplotGrob(p8), size = "last"))





library(gghighlight)





grid.arrange(p,p2, ncol=2)

plot(p2)


plot(p1)



#version 1
ordata2 %>%
  slice_max(aOR, n=30) %>%
  mutate(CleanName = fct_reorder(CleanName, aOR, .fun='median' )) %>%
  ggplot(aes(x = aOR, y = CleanName)) +
  geom_point(size = 3.5, color = "orange") +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = aOR_UCI, xmin = aOR_LCI), size = .5, height = 
                   .2, color = "gray50") +
  scale_x_log10(breaks = pretty_breaks(n = 10)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Odds ratio (log scale)")

  
    axis_combmatrix(sep = "&") +
  coord_flip()


ggplot(mtcars, aes(x=combined)) +
  geom_bar() +
  axis_combmatrix(sep = "_")



library(ComplexUpset)


#specify list of SDOH indicators
cc_vars = colnames(mydata)[17:36]


upset(
  mydata,
  cc_vars,
  mode='inclusive_intersection',
  min_size=250,
  n_intersections=30,
  width_ratio=0.2
)

library(ggplot2)
mtcars$combined <- paste0("Cyl: ", mtcars$cyl, "_Gears: ", mtcars$gear)
head(mtcars)
ggplot(mtcars, aes(x=combined)) +
  geom_bar() +
  axis_combmatrix(sep = "_")



ccdata %>%
  distinct(title, year, length, .keep_all=TRUE) %>%
  ggplot(aes(x=Genres)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20)

tidy_movies











#test if age should be discretized
library(earth)
agesplit <- earth(POORHLTH ~ AGE + RACETHX + FEMALE, data=mydata2, thresh=.001)
summary(agesplit, style="max")
#show cut points
cuts <- agesplit$cuts[agesplit$selected.terms, ]
print(cuts)
evimp(agesplit)
plot(agesplit)
boxplot(AGE~POORHLTH, data=mydata2)




