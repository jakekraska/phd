# Data Analysis for the Item Calibration Study (ICS) CHC Test
# Original code developed by Mr Jake Kraksa (Monash University)
# Built on R 3.6.3 with R Studio 1.2.5033
# set the working directory to source file directory manually
# restart R session prior to beginning

#---                                                                      ---#
################################## Libraries ##################################
#---                                                                      ---#

# Missing Data

library(psych) # version 1.9.12
library(mice) # version 3.8.0
library(tibble) # version 2.1.3
library(reshape2) # version 1.4.3
library(VIM) # version 5.1.1
library(BaylorEdPsych) # version 0.5

# CFA

library(lavaan) # version 0.6-5
library(semTools) # version 0.5-2

# Mokken

library(mokken) # version 2.8.11

# IRT

library(mirt) # version 1.3.1
library(WrightMap) # version 1.2.2

# DIF

library(difR) # version 5.0

# Visualise Results

library(ggplot2) # version 3.3.0
library(knitr) # version 1.28
library(latticeExtra)

# Data cleaning

library(dplyr) # version 0.8.5
library(stringr) # version 1.4.0
library(tidyr) # version 1.0.2

#---                                                                      ---#
################################## Set Options ##################################
#---                                                                      ---#

set.seed(145)
options(max.print=10000)

#---                                                                              ---#
################################## Create functions ##################################
#---                                                                              ---#

percentage.missing <- function(x) {
  sum(is.na(x))/length(x)*100
} # evaluates a list and determines the percentage of NA

ggplot_addline <- function(x,...) {
  gsub('\\s','\n',x)
}

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

cfa.parameters <- function(chcability, cfaobject, cfanumber, imp, print = TRUE) {
  cfa.par <- parameterEstimates(cfaobject, standardized = TRUE) %>% 
    dplyr::filter(op == "=~") %>%
    dplyr::select(Item = rhs, B = est, SE = se, Z = z, "p-value" = pvalue, "Std. Beta" = std.all) 
  
  write.table(cfa.par, 
              file = paste0("ics_tables/",chcability,"_cfa_",cfanumber,"_imp_",imp,".csv"), 
              sep = ",", quote = FALSE, 
              row.names = F)
  print("--------------")
  print(cfa.par %>%  
          kable(digits = 3, format="pandoc", caption=paste("Factor Loadings for ", deparse(substitute(object)), sep = "")))
  print("--------------")
  return(cfa.par)
}

mokken.outcomes <- function(chcability, mokkenobject, mokkennumber, imp, print = TRUE) {
  x <- coefH(mokkenobject, nice.output = FALSE)
  write.table(x$Hi, file = paste0("ics_tables/",chcability,"_mokken_",mokkennumber,"_imp_",imp,".csv"), sep = ",", quote = FALSE, row.names = T)
  print("--------------")
  print("Item Scalability Loevinger's H")
  print(x$Hi)
  print("Item Scalability SE")
  print(x$se.Hi)
  print("--------------")
  print("Scale Scalability Loevinger's H")
  print(x$H)
  print("Scale Scalability SE")
  print(x$se.H)
  print("--------------")
  return(x)
}

aisp.outcomes <- function(chcability, aispobject, aispnumber,imp, print = TRUE) {
  aisp <- aisp(aispobject)
  write.table(aisp, file = paste0("ics_tables/",chcability,"_aisp_",aispnumber,"_imp_",imp,".csv"), sep = ",", quote = FALSE, row.names = T)
  print("--------------")
  print(aisp)
  print("--------------")
  return(aisp)
}

rasch.outcomes <- function(chcability, raschobject, raschnumber, imp, print = TRUE) {
  rasch <- mirt(raschobject, model = 1, itemtype = "Rasch", SE = TRUE)
  rasch.item.fit <- itemfit(rasch)
  rasch.scale.fit <- M2(rasch)
  print("--------------")
  print(paste0("Rasch ", raschnumber, ": Imputation ", imp))
  print("--------------")
  print("Rasch Item Fit")
  print(rasch.item.fit)
  write.table(rasch.item.fit, file = paste0("ics_tables/",chcability,"_rasch_item_fit_",raschnumber,"_imp_",imp,".csv"), sep = ",", quote = FALSE, row.names = T)
  print("--------------")
  print("M2 stats")
  print(rasch.scale.fit)
  write.table(rasch.scale.fit, file = paste0("ics_tables/",chcability,"_rasch_scale_fit_",raschnumber,"_imp_",imp,".csv"), sep = ",", quote = FALSE, row.names = T)
  print("--------------")
  print("Marginal Rxx")
  print(marginal_rxx(rasch))
  print("--------------")
  print(plot(rasch, type= "trace", facet_items = FALSE, main =""))
  print(plot(rasch, type = "infoSE", xlim=c(-6,6), main = "") + layer(panel.abline(h=.5, col='red')))
  print(plot(rasch, type = "rxx", xlim=c(-6,6), main = ""))
  print("--------------")
  return(list(rasch,rasch.item.fit))
}

rasch.poor.fit <- function(dataset, raschobject, imp, print = TRUE) {
  nan <- colnames(select(dataset, as.vector(filter(raschobject, p.S_X2 == "NaN")$item)))
  na <- colnames(select(dataset, as.vector(filter(raschobject, is.na(p.S_X2))$item)))
  below.01 <- colnames(select(dataset, as.vector(filter(raschobject, p.S_X2 < 0.01)$item)))
  print("----------")
  print(paste0("Items with NaN for Imputation ",imp))
  print(nan)
  print("----------")
  print(paste0("Items with NA for Imputation ",imp))
  print(na)
  print("----------")
  print(paste0("Items Below .01 for Imputation ",imp))
  print(below.01)
  return(list("nan" = nan, "na" = na, "below.01" = below.01))
}

remove.items.from.sets <- function(chcability,items) {
  # Get the item set objects for items
  p2a.items <- get(paste0(chcability,".p2a.items"))
  p2u.items <- get(paste0(chcability,".p2u.items"))
  seta.items <- get(paste0(chcability,".seta.items"))
  setb.items <- get(paste0(chcability,".setb.items"))
  anchor.items <- get(paste0(chcability,".anchor.items"))
  
  # Get the item set objects for times
  p2a.times <- get(paste0(chcability,".p2a.times"))
  p2u.times <- get(paste0(chcability,".p2u.times"))
  seta.times <- get(paste0(chcability,".seta.times"))
  setb.times <- get(paste0(chcability,".setb.times"))
  anchor.times <- get(paste0(chcability,".anchor.times"))
  
  # remove the items from the item sets
  p2a.items <- p2a.items[!p2a.items %in% items]
  p2u.items <- p2u.items[!p2u.items %in% items]
  seta.items <- seta.items[!seta.items %in% items]
  setb.items <- setb.items[!setb.items %in% items]
  anchor.items <- anchor.items[!anchor.items %in% items]
  
  # remove the times from the item sets
  p2a.times <- p2a.times[!p2a.times %in% items]
  p2u.times <- p2u.times[!p2u.times %in% items]
  seta.times <- seta.times[!seta.times %in% items]
  setb.times <- setb.times[!setb.times %in% items]
  anchor.times <- anchor.times[!anchor.times %in% items]
  
  # assign the item sets back to the objects
  assign(paste0(chcability,".p2a.items"),p2a.items)
  assign(paste0(chcability,".p2u.items"),p2u.items)
  assign(paste0(chcability,".seta.items"),seta.items)
  assign(paste0(chcability,".setb.items"),setb.items)
  assign(paste0(chcability,".anchor.items"),anchor.items)
  
  # assign the item sets back to the objects
  assign(paste0(chcability,".p2a.times"),p2a.times)
  assign(paste0(chcability,".p2u.times"),p2u.times)
  assign(paste0(chcability,".seta.times"),seta.times)
  assign(paste0(chcability,".setb.times"),setb.times)
  assign(paste0(chcability,".anchor.times"),anchor.times)
}

#---                                                                      ---#
################################## Load Data ##################################
#---                                                                      ---#

data <- read.csv("data.csv", stringsAsFactors = FALSE)
data$phase <- recode(data$phase, itos = "ITOS", p2a = "ICS-A", p2u = "ICS-U")
data$phase <- as.factor(data$phase)
data$gender <- as.factor(data$gender)
data$nationality <- as.factor(data$nationality)
data$device <- as.factor(data$device)


#---                                                                                            ---#
################################## Demographics and Data Cleaning ##################################
#---                                                                                            ---#

### Gender Demographics Prior to Data Cleaning

participants.by.gender <- table(data$phase, data$gender)
participants.by.gender <- rbind(participants.by.gender,colSums(participants.by.gender))
participants.by.gender

### Age Demographics Prior to Data Cleaning

data$age.group <- factor(cut(data$age, breaks = c(0, 6, 18, 30, 40, 50, 60, 70, 80, 91, Inf), 
                             labels = c("0-5","6-17","18-29","30-39","40-49","50-59","60-69","70-79","80-90","91+"), right = FALSE),
                         exclude = NULL)

participants.by.age <- table(data$age.group, data$phase)
participants.by.age <- cbind(participants.by.age,rowSums(participants.by.age))
participants.by.age

ggplot(data = data, aes(x = age.group, fill = phase)) +
  geom_bar(position = position_dodge(width = .7)) +
  geom_text(stat = "count", aes(label = ..count..), position = position_dodge(width = .7)) +
  xlab("Age Group") +
  scale_y_continuous(name = "Number of Participants", breaks = seq.int(0,1010,125), limits = c(0,1010)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(fill = guide_legend(title = "Phase"))

### Nationality Demographics Prior to Data Cleaning

participants.by.nationality <- table(data$nationality, data$phase)
participants.by.nationality <- cbind(participants.by.nationality, rowSums(participants.by.nationality))
participants.by.nationality

### Remove participants that did not complete any of the tests

table(data$phase) # number in each phase before removal

non.completers <- select(filter(data, gcscore == 0, gfscore == 0, gvscore == 0, gwmscore == 0), 
                         id, age, gender, phase, gcscore, gfscore, gvscore, gwmscore) # df of those that did not complete

table(non.completers$phase) # number in each phase to be removed

data <- anti_join(x = data, y = non.completers, by = "id") # remove the non completers from the data set
rm(non.completers) # remove non completers df from environment

table(data$phase) # number in each phase after removal

### Participants by CHC Factor

outcomes <- c("gcscore","gfscore","gvscore","gwmscore")
phases <- c("ITOS","ICS-A","ICS-U")
participants.by.chc <- 
  sapply(outcomes, function(x) { 
    sapply(phases, function(y) {
      length(which(data$phase == y & data[x] != 0))
    })
  })
participants.by.chc <- rbind(participants.by.chc,colSums(participants.by.chc))
participants.by.chc

### Gender Demographics After Non-Responder Cleaning

participants.by.gender <- table(data$phase, data$gender)
participants.by.gender <- rbind(participants.by.gender,colSums(participants.by.gender))
participants.by.gender

### Age Demographics After Non-Responder Cleaning

itos <- filter(data, phase == "ITOS") # separate out ITOS data
icsa <- filter(data, phase == "ICS-A") # separate out ICS-A data
icsu <- filter(data, phase == "ICS-U") # separate out ICS-U data

itos$age <- replace(itos$age, itos$age < 18 | itos$age > 90, NA) # change invalid data to NA
icsa$age <- replace(icsa$age, icsa$age < 18 | icsa$age > 90, NA) # change invalid data to NA
icsu$age <- replace(icsu$age, icsu$age < 6 | icsu$age > 18, NA) # change invalid data to NA

data <- rbind(itos,icsa,icsu) # combine the data back together
rm(itos,icsa,icsu) # remove individual dfs

data$age.group <- factor(cut(data$age, breaks = c(0, 6, 18, 30, 40, 50, 60, 70, 80, 91, Inf), 
                             labels = c("0-5","6-17","18-29","30-39","40-49","50-59","60-69","70-79","80-90","91+"), right = FALSE),
                         exclude = NULL)

participants.by.age <- table(data$age.group, data$phase)
participants.by.age <- cbind(participants.by.age,rowSums(participants.by.age))
participants.by.age

# Age distribution
ggplot(data = data, aes(x = age.group, fill = phase)) +
  geom_bar(position = position_dodge(width = .7)) +
  geom_text(stat = "count", aes(label = ..count..), position = position_dodge(width = .7)) +
  xlab("Age Group") +
  scale_y_continuous(name = "Number of Participants", breaks = seq.int(0,1010,125), limits = c(0,1010)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(fill = guide_legend(title = "Phase"))

# Boxplot based on phase
ggplot(data) +
  geom_boxplot() +
  aes(x = phase, y = age, fill = phase) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous("Age (years old)", breaks = c(6,18,30,40,50,60,70,80,90), limits = c(6,90)) +
  theme(legend.position = "none")

# Boxplot for entire sample
ggplot(data) +
  geom_boxplot() +
  aes(y = age) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous("Age (years old)", breaks = c(6,18,30,40,50,60,70,80,90), limits = c(6,90))

### Nationality demographics after non-responder cleaning

participants.by.nationality <- table(data$nationality, data$phase)
participants.by.nationality <- cbind(participants.by.nationality, rowSums(participants.by.nationality))
participants.by.nationality

### Demographics for final data set

data[,c("gcscore","gfscore","gvscore","gwmscore")][data[, c("gcscore","gfscore","gvscore","gwmscore")] == 0] <- NA

psych::describe(data[,c("age","gender","nationality","device","gcscore","gfscore","gvscore","gwmscore",
                        "gctime","gftime","gvtime","gwmtime")])

#---                                                                                      ---#
################################## Subset Data for Analysis ##################################
#---                                                                                      ---#

### Names of Items by Set

gc.seta.items <- sapply(c(4:19,21:34), function(x) {paste0("score_gc",x)}) # set the names of the items for set a
gc.seta.times <- sapply(c(4:19,21:34), function(x) {paste0("timeTaken_gc",x)}) # set the names of the items for set a
gc.setb.items <- sapply(c(56:60,62:67,69:107), function(x) {paste0("score_gc",x)}) # set the names for the items for set b
gc.setb.times <- sapply(c(56:60,62:67,69:107), function(x) {paste0("timeTaken_gc",x)}) # set the names for the items for set b
gc.anchor.items <- sapply(c(36:46,48:51,53:55), function(x) {paste0("score_gc",x)}) # set the names for the items for anchor items
gc.anchor.times <- sapply(c(36:46,48:51,53:55), function(x) {paste0("timeTaken_gc",x)}) # set the names for the items for anchor items

gf.30.seta.items <- sapply(c(1:10), function(x) {paste0("score_gf",x)}) # set the names of the items for set a
gf.30.seta.times <- sapply(c(1:10), function(x) {paste0("timeTaken_gf",x)}) # set the names of the items for set a
gf.30.setb.items <- sapply(c(37,40:42,45:49,52:54,57:59,62:64,67:69,72:78), function(x) {paste0("score_gf",x)}) # set the names for the items for set b
gf.30.setb.times <- sapply(c(37,40:42,45:49,52:54,57:59,62:64,67:69,72:78), function(x) {paste0("timeTaken_gf",x)}) # set the names for the items for set b
gf.30.anchor.items <- sapply(c(11:15), function(x) {paste0("score_gf",x)}) # set the names for the items for anchor items
gf.30.anchor.times <- sapply(c(11:15), function(x) {paste0("timeTaken_gf",x)}) # set the names for the items for anchor items

gf.60.seta.items <- sapply(c(16:28), function(x) {paste0("score_gf",x)}) # set the names of the items for set a
gf.60.seta.times <- sapply(c(16:28), function(x) {paste0("timeTaken_gf",x)}) # set the names of the items for set a
gf.60.setb.items <- sapply(c(34,38:39,43:44,50:51,55:56,60:61,65:66,70:71,79:80), function(x) {paste0("score_gf",x)}) # set the names for the items for set b
gf.60.setb.times <- sapply(c(34,38:39,43:44,50:51,55:56,60:61,65:66,70:71,79:80), function(x) {paste0("timeTaken_gf",x)}) # set the names for the items for set b
gf.60.anchor.items <- sapply(c(29:33), function(x) {paste0("score_gf",x)}) # set the names for the items for anchor items
gf.60.anchor.times <- sapply(c(29:33), function(x) {paste0("timeTaken_gf",x)}) # set the names for the items for anchor items

gf.seta.items <- c(gf.30.seta.items, gf.60.seta.items)
gf.seta.times <- c(gf.30.seta.times, gf.60.seta.times)
gf.anchor.items <- c(gf.30.anchor.items, gf.60.anchor.items)
gf.anchor.times <- c(gf.30.anchor.times, gf.60.anchor.times)
gf.setb.items <- c(gf.30.setb.items, gf.60.setb.items)
gf.setb.times <- c(gf.30.setb.times, gf.60.setb.times)

gv.seta.items <- sapply(c(1:21), function(x) {paste0("score_gv",x)}) # set the names of the items for set a
gv.seta.times <- sapply(c(1:21), function(x) {paste0("timeTaken_gv",x)}) # set the names of the items for set a
gv.setb.items <- sapply(c(53:72), function(x) {paste0("score_gv",x)}) # set the names for the items for set b
gv.setb.times <- sapply(c(53:72), function(x) {paste0("timeTaken_gv",x)}) # set the names for the items for set b
gv.anchor.items <- sapply(c(22:52), function(x) {paste0("score_gv",x)}) # set the names for the items for anchor items
gv.anchor.times <- sapply(c(22:52), function(x) {paste0("timeTaken_gv",x)}) # set the names for the items for anchor items

gwm.seta.items <- sapply(c(11:21), function(x) {paste0("score_gwm",x)}) # set the names of the items for set a
gwm.seta.times <- sapply(c(11:21), function(x) {paste0("timeTaken_gwm",x)}) # set the names of the items for set a
gwm.setb.items <- sapply(c(33:44), function(x) {paste0("score_gwm",x)}) # set the names for the items for set b
gwm.setb.times <- sapply(c(33:44), function(x) {paste0("timeTaken_gwm",x)}) # set the names for the items for set b
gwm.anchor.items <- sapply(c(22:32), function(x) {paste0("score_gwm",x)}) # set the names for the items for anchor items
gwm.anchor.times <- sapply(c(22:32), function(x) {paste0("timeTaken_gwm",x)}) # set the names for the items for anchor items

### Names of Items by Phase

gc.itos.items <- sapply(1:55, function(x) {paste0("score_gc",x)})
gc.itos.times <- sapply(1:55, function(x) {paste0("timeTaken_gc",x)})
gc.icsa.items <- sapply(c(36:46,48:51,53:55,56:60,62:67,69:107,4:19,21:34), function(x) {paste0("score_gc",x)})
gc.icsa.times <- sapply(c(36:46,48:51,53:55,56:60,62:67,69:107,4:19,21:34), function(x) {paste0("timeTaken_gc",x)})
gc.icsu.items <- sapply(c(36:46,48:51,53:55,4:19,21:34,56:60,62:67,69:107), function(x) {paste0("score_gc",x)})
gc.icsu.times <- sapply(c(36:46,48:51,53:55,4:19,21:34,56:60,62:67,69:107), function(x) {paste0("timeTaken_gc",x)})
gc.items.by.phase <- list(itos = gc.itos.items, icsa = gc.icsa.items, icsu = gc.icsu.items)

gf.itos.items <- sapply(1:33, function(x) {paste0("score_gf",x)})
gf.itos.times <- sapply(1:33, function(x) {paste0("timeTaken_gf",x)})
gf.icsa.items <- sapply(c(11:15, # anchor 30 seconds
                         37,40:42,45:49,52:54,57:59,62:64,67:69,72:78, # set b 30 seconds
                         1:10, # set a 30 seconds
                         29:33, # anchor 60 seconds
                         34,38:39,43:44,50:51,55:56,60:61,65:66,70:71,79:80, # set b 60 seconds
                         16:28),  # set a 60 seconds
                       function(x) {paste0("score_gf",x)})
gf.icsa.times <- sapply(c(11:15, # anchor 30 seconds
                         37,40:42,45:49,52:54,57:59,62:64,67:69,72:78, # set b 30 seconds
                         1:10, # set a 30 seconds
                         29:33, # anchor 60 seconds
                         34,38:39,43:44,50:51,55:56,60:61,65:66,70:71,79:80, # set b 60 seconds
                         16:28),  # set a 60 seconds
                       function(x) {paste0("timeTaken_gf",x)})
gf.icsu.items <- sapply(c(11:15, # anchor 30 seconds
                         1:10, # set a 30 seconds
                         37,40:42,45:49,52:54,57:59,62:64,67:69,72:78, # set b 30 seconds
                         29:33, # anchor 60 seconds
                         16:28,  # set a 60 seconds
                         34,38:39,43:44,50:51,55:56,60:61,65:66,70:71,79:80), # set b 60 seconds
                       function(x) {paste0("score_gf",x)})
gf.icsu.times <- sapply(c(11:15, # anchor 30 seconds
                         1:10, # set a 30 seconds
                         37,40:42,45:49,52:54,57:59,62:64,67:69,72:78, # set b 30 seconds
                         29:33, # anchor 60 seconds
                         16:28,  # set a 60 seconds
                         34,38:39,43:44,50:51,55:56,60:61,65:66,70:71,79:80), # set b 60 seconds
                       function(x) {paste0("timeTaken_gf",x)})
gf.items.by.phase <- list(itos = gf.itos.items, icsa = gf.icsa.items, icsu = gf.icsu.items)

gv.itos.items <- sapply(1:52, function(x) {paste0("score_gv",x)})
gv.itos.times <- sapply(1:52, function(x) {paste0("timeTaken_gv",x)})
gv.icsa.items <- sapply(c(22:52,53:72,1:21), function(x) {paste0("score_gv",x)})
gv.icsa.times <- sapply(c(22:52,53:72,1:21), function(x) {paste0("timeTaken_gv",x)})
gv.icsu.items <- sapply(c(22:52,1:21,53:72), function(x) {paste0("score_gv",x)})
gv.icsu.times <- sapply(c(22:52,1:21,53:72), function(x) {paste0("timeTaken_gv",x)})
gv.items.by.phase <- list(itos = gv.itos.items, icsa = gv.icsa.items, icsu = gv.icsu.items)

gwm.itos.items <- sapply(1:38, function(x) {paste0("score_gwm",x)})
gwm.itos.times <- sapply(1:38, function(x) {paste0("timeTaken_gwm",x)})
gwm.icsa.items <- sapply(11:44, function(x) {paste0("score_gwm",x)})
gwm.icsa.times <- sapply(11:44, function(x) {paste0("timeTaken_gwm",x)})
gwm.icsu.items <- sapply(11:44, function(x) {paste0("score_gwm",x)})
gwm.icsu.times <- sapply(11:44, function(x) {paste0("timeTaken_gwm",x)})
gwm.items.by.phase <- list(itos = gwm.itos.items, icsa = gwm.icsa.items, icsu = gwm.icsu.items)

### Create columns with indicator of missingness, 1 = some data missing, 0 = full data
data <- mutate(data, 
               gc.set.a = ceiling(rowSums(is.na(data[gc.seta.items]))/length(gc.seta.items)),
               gc.set.b = ceiling(rowSums(is.na(data[gc.setb.items]))/length(gc.setb.items)),
               gc.anchor = ceiling(rowSums(is.na(data[gc.anchor.items]))/length(gc.anchor.items)),
               gf.30.set.a = ceiling(rowSums(is.na(data[gf.30.seta.items]))/length(gf.30.seta.items)),
               gf.30.set.b = ceiling(rowSums(is.na(data[gf.30.setb.items]))/length(gf.30.setb.items)),
               gf.30.anchor = ceiling(rowSums(is.na(data[gf.30.anchor.items]))/length(gf.30.anchor.items)),
               gf.60.set.a = ceiling(rowSums(is.na(data[gf.60.seta.items]))/length(gf.60.seta.items)),
               gf.60.set.b = ceiling(rowSums(is.na(data[gf.60.setb.items]))/length(gf.60.setb.items)),
               gf.60.anchor = ceiling(rowSums(is.na(data[gf.60.anchor.items]))/length(gf.60.anchor.items)),
               gv.set.a = ceiling(rowSums(is.na(data[gv.seta.items]))/length(gv.seta.items)),
               gv.set.b = ceiling(rowSums(is.na(data[gv.setb.items]))/length(gv.setb.items)),
               gv.anchor = ceiling(rowSums(is.na(data[gv.anchor.items]))/length(gv.anchor.items)),
               gwm.set.a = ceiling(rowSums(is.na(data[gwm.seta.items]))/length(gwm.seta.items)),
               gwm.set.b = ceiling(rowSums(is.na(data[gwm.setb.items]))/length(gwm.setb.items)),
               gwm.anchor = ceiling(rowSums(is.na(data[gwm.anchor.items]))/length(gwm.anchor.items))) 

### Recode missingness # NA = missing data in set, 1 = full data
data$gc.set.a <- data$gc.set.a %>% na_if(1) %>% recode(`0` = 1)
data$gc.set.b <- data$gc.set.b %>% na_if(1) %>% recode(`0` = 1)
data$gc.anchor <- data$gc.anchor %>% na_if(1) %>% recode(`0` = 1)
data$gf.30.set.a <- data$gf.30.set.a %>% na_if(1) %>% recode(`0` = 1)
data$gf.30.set.b <- data$gf.30.set.b %>% na_if(1) %>% recode(`0` = 1)
data$gf.30.anchor <- data$gf.30.anchor %>% na_if(1) %>% recode(`0` = 1)
data$gf.60.set.a <- data$gf.60.set.a %>% na_if(1) %>% recode(`0` = 1)
data$gf.60.set.b <- data$gf.60.set.b %>% na_if(1) %>% recode(`0` = 1)
data$gf.60.anchor <- data$gf.60.anchor %>% na_if(1) %>% recode(`0` = 1)
data$gv.set.a <- data$gv.set.a %>% na_if(1) %>% recode(`0` = 1)
data$gv.set.b <- data$gv.set.b %>% na_if(1) %>% recode(`0` = 1)
data$gv.anchor <- data$gv.anchor %>% na_if(1) %>% recode(`0` = 1)
data$gwm.set.a <- data$gwm.set.a %>% na_if(1) %>% recode(`0` = 1)
data$gwm.set.b <- data$gwm.set.b %>% na_if(1) %>% recode(`0` = 1)
data$gwm.anchor <- data$gwm.anchor %>% na_if(1) %>% recode(`0` = 1)

### Data by phase

itos <- filter(data, phase == "ITOS") # separate out ITOS data
icsa <- filter(data, phase == "ICS-A") # separate out ICS-A data
icsu <- filter(data, phase == "ICS-U") # separate out ICS-U data

### CHC factors with all variables

gc <- select(filter(data, gcscore > 0),id,phase,age,age.group,gender,nationality,device,
             starts_with("gc"),starts_with("score_gc"),starts_with("timeTaken_gc"),starts_with("wisc."))
gf <- select(filter(data, gfscore > 0),id,phase,age,age.group,gender,nationality,device,
             starts_with("gf"),starts_with("score_gf"),starts_with("timeTaken_gf"),starts_with("wisc."))
gv <- select(filter(data, gvscore > 0),id,phase,age,age.group,gender,nationality,device,
             starts_with("gv"),starts_with("score_gv"),starts_with("timeTaken_gv"),starts_with("wisc."))
gwm <- select(filter(data, gwmscore > 0),id,phase,age,age.group,gender,nationality,device,
              starts_with("gwm"),starts_with("score_gwm"),starts_with("timeTaken_gwm"),starts_with("wisc."))

#---                                                                              ---#
################################## Data Exploration ##################################
#---                                                                              ---#

#### Gc Raw Data ####

# Descriptive Statistics

psych::describe(gc$gcscore)
psych::describeBy(gc$gcscore, group = gc$phase, digits = 2, mat = TRUE)

psych::describe(gc$gctime)
psych::describeBy(gc$gctime, group = gc$phase, digits = 2, mat = TRUE)

# Histogram of scores
ggplot(data, aes(x=gcscore)) + 
  geom_histogram(binwidth = 5, aes(fill = phase)) +
  facet_grid(phase ~ .) +
  scale_x_continuous(breaks = c(seq(0,100,5)), name = "Raw Score") +
  scale_y_continuous(breaks = c(seq(0,700,100)), name = "Frequency") +
  theme(legend.position = "none")

# ANOVA of differences in score by phase
summary(aov(formula = gcscore ~ phase, data = gc))
confint(aov(formula = gcscore ~ phase, data = gc))
TukeyHSD(aov(formula = gcscore ~ phase, data = gc))

# Boxplot of total score
ggplot(gc) + 
  geom_boxplot(mapping = aes(x = phase, y = gcscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,100,5)), limits = c(0,100), name = "Total Score") +
  theme(legend.position = "none")

# ANOVA of differences in time by phase
summary(aov(formula = gctime ~ phase, data = gc))
confint(aov(formula = gctime ~ phase, data = gc))
TukeyHSD(aov(formula = gctime ~ phase, data = gc))

# Boxplot of time
ggplot(gc) + 
  geom_boxplot(mapping = aes(x = phase, y = gctime/60, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  theme(legend.position = "none")

# get the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
gc.correct <- as.data.frame(
  sapply(gc[unique(c(gc.itos.items,gc.icsa.items,gc.icsu.items))], function(x) {
    (1 - (table(x)[[1]]/length(na.omit(x)))) * 100
  })
) 

# Turn df into usable data for ggplot
row.names(gc.correct) <- str_remove(row.names(gc.correct), "score_gc")
names(gc.correct) <- "Correct"
gc.correct <- rownames_to_column(gc.correct, var = "Item")

# plot the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
ggplot(data = gc.correct, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gc.correct$Item), by = 2)) + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

# Time by Score
ggplot(gc, aes(x = (gctime/60), y = gcscore)) +
  geom_point() +
  facet_grid(phase ~ .) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  scale_y_continuous(breaks = c(seq(0,100,10)), limits = c(0,100), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# Gc Gender Differences
ggplot(gc) +
  geom_boxplot(mapping = aes(x = gender, y = gcscore, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,100,5)), limits = c(0,100), name = "Total Score") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

summary(aov(formula = gcscore ~ gender, data = gc))
confint(aov(formula = gcscore ~ gender, data = gc))
TukeyHSD(aov(formula = gcscore ~ gender, data = gc))

ggplot(gc) +
  geom_boxplot(mapping = aes(x = gender, y = gctime/60, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

summary(aov(formula = gctime ~ gender, data = gc))
confint(aov(formula = gctime ~ gender, data = gc))
TukeyHSD(aov(formula = gctime ~ gender, data = gc))

#### Gf Raw Data ####

# Descriptive Statistics

psych::describe(gf$gfscore)
psych::describeBy(gf$gfscore, group = gf$phase, digits = 2, mat = TRUE)

psych::describe(gf$gftime)
psych::describeBy(gf$gftime, group = gf$phase, digits = 2, mat = TRUE)

# Histogram of scores
ggplot(data, aes(x=gfscore)) + 
  geom_histogram(binwidth = 5, aes(fill = phase)) +
  facet_grid(phase ~ .) +
  scale_x_continuous(breaks = c(seq(0,100,5)), name = "Raw Score") +
  scale_y_continuous(breaks = c(seq(0,700,100)), name = "Frequency") +
  theme(legend.position = "none")

# ANOVA of differences in score by phase
summary(aov(formula = gfscore ~ phase, data = gf))
confint(aov(formula = gfscore ~ phase, data = gf))
TukeyHSD(aov(formula = gfscore ~ phase, data = gf))

# Boxplot of total score
ggplot(gf) + 
  geom_boxplot(mapping = aes(x = phase, y = gfscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,50,5)), limits = c(0,50), name = "Total Score") +
  theme(legend.position = "none")

# ANOVA of differences in time by phase
summary(aov(formula = gftime ~ phase, data = gf))
confint(aov(formula = gftime ~ phase, data = gf))
TukeyHSD(aov(formula = gftime ~ phase, data = gf))

# Boxplot of time
ggplot(gf) + 
  geom_boxplot(mapping = aes(x = phase, y = gftime/60, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  theme(legend.position = "none")

# get the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
gf.correct <- as.data.frame(
  sapply(gf[unique(c(gf.itos.items,gf.icsa.items,gf.icsu.items))], function(x) {
    (1 - (table(x)[[1]]/length(na.omit(x)))) * 100
  })
) 

# Turn df into usable data for ggplot
row.names(gf.correct) <- str_remove(row.names(gf.correct), "score_gf")
names(gf.correct) <- "Correct"
gf.correct <- rownames_to_column(gf.correct, var = "Item")

# plot the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
ggplot(data = gf.correct, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gf.correct$Item), by = 2)) + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

# Time by Score
ggplot(gf, aes(x = (gftime/60), y = gfscore)) +
  geom_point() +
  facet_grid(phase ~ .) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  scale_y_continuous(breaks = c(seq(0,50,5)), limits = c(0,50), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# Gf Gender Differences by score
ggplot(gf) +
  geom_boxplot(mapping = aes(x = gender, y = gfscore, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,50,5)), limits = c(0,50), name = "Total Score") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

# ANOVA of score between gender
summary(aov(formula = gfscore ~ gender, data = gf))
confint(aov(formula = gfscore ~ gender, data = gf))
TukeyHSD(aov(formula = gfscore ~ gender, data = gf))

# Gf Gender differences by time
ggplot(gf) +
  geom_boxplot(mapping = aes(x = gender, y = gftime/60, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

# ANOVA of time between gender
summary(aov(formula = gftime ~ gender, data = gf))
confint(aov(formula = gftime ~ gender, data = gf))
TukeyHSD(aov(formula = gftime ~ gender, data = gf))

#### Gv Raw Data ####

# Descriptive Statistics

psych::describe(gv$gvscore)
psych::describeBy(gv$gvscore, group = gv$phase, digits = 2, mat = TRUE)

psych::describe(gv$gvtime)
psych::describeBy(gv$gvtime, group = gv$phase, digits = 2, mat = TRUE)

# Histogram of scores
ggplot(data, aes(x=gvscore)) + 
  geom_histogram(binwidth = 5, aes(fill = phase)) +
  facet_grid(phase ~ .) +
  scale_x_continuous(breaks = c(seq(0,100,5)), name = "Raw Score") +
  scale_y_continuous(breaks = c(seq(0,700,100)), name = "Frequency") +
  theme(legend.position = "none")

# ANOVA of differences in score by phase
summary(aov(formula = gvscore ~ phase, data = gv))
confint(aov(formula = gvscore ~ phase, data = gv))
TukeyHSD(aov(formula = gvscore ~ phase, data = gv))

# Boxplot of total score
ggplot(gv) + 
  geom_boxplot(mapping = aes(x = phase, y = gvscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,70,5)), limits = c(0,70), name = "Total Score") +
  theme(legend.position = "none")

# ANOVA of differences in time by phase
summary(aov(formula = gvtime ~ phase, data = gv))
confint(aov(formula = gvtime ~ phase, data = gv))
TukeyHSD(aov(formula = gvtime ~ phase, data = gv))

# Boxplot of time
ggplot(gv) + 
  geom_boxplot(mapping = aes(x = phase, y = gvtime/60, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  theme(legend.position = "none")

# get the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
gv.correct <- as.data.frame(
  sapply(gv[unique(c(gv.itos.items,gv.icsa.items,gv.icsu.items))], function(x) {
    (1 - (table(x)[[1]]/length(na.omit(x)))) * 100
  })
) 

# Turn df into usable data for ggplot
row.names(gv.correct) <- str_remove(row.names(gv.correct), "score_gv")
names(gv.correct) <- "Correct"
gv.correct <- rownames_to_column(gv.correct, var = "Item")

# plot the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
ggplot(data = gv.correct, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gv.correct$Item), by = 2)) + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

# Gv Time by Score
ggplot(gv, aes(x = (gvtime/60), y = gvscore)) +
  geom_point() +
  facet_grid(phase ~ .) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  scale_y_continuous(breaks = c(seq(0,100,10)), limits = c(0,100), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# Gv Gender Differences by Score
ggplot(gv) +
  geom_boxplot(mapping = aes(x = gender, y = gvscore, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,70,5)), limits = c(0,70), name = "Total Score") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

# ANOVA of Gv score by gender
summary(aov(formula = gvscore ~ gender, data = gv))
confint(aov(formula = gvscore ~ gender, data = gv))
TukeyHSD(aov(formula = gvscore ~ gender, data = gv))

# Gv Gender Differences by Time
ggplot(gv) +
  geom_boxplot(mapping = aes(x = gender, y = gvtime/60, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

# ANOVA of Gv time by gender
summary(aov(formula = gvtime ~ gender, data = gv))
confint(aov(formula = gvtime ~ gender, data = gv))
TukeyHSD(aov(formula = gvtime ~ gender, data = gv))

#### Gwm Raw Data ####

# Descriptive Statistics

psych::describe(gwm$gwmscore)
psych::describeBy(gwm$gwmscore, group = gwm$phase, digits = 2, mat = TRUE)

psych::describe(gwm$gwmtime)
psych::describeBy(gwm$gwmtime, group = gwm$phase, digits = 2, mat = TRUE)

# Histogram of scores
ggplot(data, aes(x=gwmscore)) + 
  geom_histogram(binwidth = 5, aes(fill = phase)) +
  facet_grid(phase ~ .) +
  scale_x_continuous(breaks = c(seq(0,40,5)), name = "Raw Score") +
  scale_y_continuous(breaks = c(seq(0,200,25)), name = "Frequency") +
  theme(legend.position = "none")

# ANOVA of differences in score by phase
summary(aov(formula = gwmscore ~ phase, data = gwm))
confint(aov(formula = gwmscore ~ phase, data = gwm))
TukeyHSD(aov(formula = gwmscore ~ phase, data = gwm))

# Boxplot of total score
ggplot(gwm) + 
  geom_boxplot(mapping = aes(x = phase, y = gwmscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,40,5)), limits = c(0,40), name = "Total Score") +
  theme(legend.position = "none")

# ANOVA of differences in time by phase
summary(aov(formula = gwmtime ~ phase, data = gwm))
confint(aov(formula = gwmtime ~ phase, data = gwm))
TukeyHSD(aov(formula = gwmtime ~ phase, data = gwm))

# Boxplot of time
ggplot(gwm) + 
  geom_boxplot(mapping = aes(x = phase, y = gwmtime/60, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  theme(legend.position = "none")

# get the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
gwm.correct <- as.data.frame(
  sapply(gwm[unique(c(gwm.itos.items,gwm.icsa.items,gwm.icsu.items))], function(x) {
    (1 - (table(x)[[1]]/length(na.omit(x)))) * 100
  })
) 

# Turn df into usable data for ggplot
row.names(gwm.correct) <- str_remove(row.names(gwm.correct), "score_gwm")
names(gwm.correct) <- "Correct"
gwm.correct <- rownames_to_column(gwm.correct, var = "Item")

# plot the proportion of participants that obtained a correct answer 
# based on the number of participants that attempted the question
ggplot(data = gwm.correct, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gwm.correct$Item), by = 1)) + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

# Time by Score
ggplot(gwm, aes(x = (gwmtime/60), y = gwmscore)) +
  geom_point() +
  facet_grid(phase ~ .) +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken (mins)") +
  scale_y_continuous(breaks = c(seq(0,40,5)), limits = c(0,40), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# Gwm Gender Differences by score
ggplot(gwm) +
  geom_boxplot(mapping = aes(x = gender, y = gwmscore, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,40,5)), limits = c(0,40), name = "Total Score") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

# ANOVA of score on gender
summary(aov(formula = gwmscore ~ gender, data = gwm))
confint(aov(formula = gwmscore ~ gender, data = gwm))
TukeyHSD(aov(formula = gwmscore ~ gender, data = gwm))

# Gwm gender differences by time
ggplot(gwm) +
  geom_boxplot(mapping = aes(x = gender, y = gwmtime/60, fill = gender), width = .25) +
  scale_y_continuous(breaks = c(seq(0,30,1)), limits = c(0,30), name = "Time Taken") +
  scale_x_discrete(labels = c("Female", "Male", "Other", "PNTS"), name = "Gender") +
  theme(legend.position = "none")

# ANOVA of time on gender
summary(aov(formula = gwmtime ~ gender, data = gwm))
confint(aov(formula = gwmtime ~ gender, data = gwm))
TukeyHSD(aov(formula = gwmtime ~ gender, data = gwm))

#---                                                                                      ---#
################################## Missing Data Exploration ##################################
#---                                                                                      ---#

#### Gc Missing Data ####
#### Use the individual CHC factor data.frame because this excludes
#### total scores of 0

psych::describe(apply(select(gc, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item

gc.missing <- sapply(1:3, function(x) {
  apply(filter(select(gc, gc.items.by.phase[[x]]), gc$phase == phases[[x]]), 2, percentage.missing)
}) # create list of each phase item data missing

sapply(1:3, function(x) {
  round(psych::describe(gc.missing[[x]]),2)
}) # summary of gc missing data by phase

gc.missing <- bind_rows(lapply(gc.missing, as.data.frame.list)) # combine list into dataframe
row.names(gc.missing) <- phases # name rows
gc.missing <- rownames_to_column(gc.missing, var = "phase")

ggplot(reshape2::melt(gc.missing, "phase"), aes(x=as.numeric(variable), y = value, colour = phase)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_vline(xintercept = 34) +
  geom_vline(xintercept = 56) +
  annotate("text", x = 16, y = 40, label = "Set A Items") +
  annotate("text", x = 46, y = 40, label = "Anchor Items") +
  annotate("text", x = 85, y = 40, label = "Set B Items")

gc.itos <- filter(gc, phase == "ITOS")
gc.icsa <- filter(gc, phase == "ICS-A")
gc.icsu <- filter(gc, phase == "ICS-U")

gc.itos$gcmissing <- apply(select(gc.itos, gc.items.by.phase[[1]]), 1, percentage.missing) # missing data by items administered
gc.icsa$gcmissing <- apply(select(gc.icsa, gc.items.by.phase[[2]]), 1, percentage.missing) # missing data by items administered
gc.icsu$gcmissing <- apply(select(gc.icsu, gc.items.by.phase[[3]]), 1, percentage.missing)
gc <- rbind(gc.itos,gc.icsa,gc.icsu) # combine the data back together

psych::describeBy(gc$gcmissing, group = gc$phase)

aggr(gc[c("gc.set.a","gc.set.b","gc.anchor")], combined = TRUE,
     labels = c("Set A", "Set B", "Anchor")) # Visualise Missing Data Patterns

#### Gf Missing Data ####
#### Use the individual CHC factor data.frame because this excludes
#### scores of 0

psych::describe(apply(select(gf, starts_with("score_gf")), 2, percentage.missing))

gf.missing <- sapply(1:3, function(x) {
  apply(filter(select(gf, gf.items.by.phase[[x]]), gf$phase == phases[[x]]), 2, percentage.missing)
}) # create list of each phase item data missing

sapply(1:3, function(x) {
  psych::describe(gf.missing[[x]])
}) # summary of gf missing data by phase♪

gf.missing <- bind_rows(lapply(gf.missing, as.data.frame.list)) # combine list into dataframe
row.names(gf.missing) <- phases # name rows
gf.missing <- rownames_to_column(gf.missing, var = "phase")

ggplot(reshape2::melt(gf.missing, "phase"), aes(x=as.numeric(variable), y = value, colour = phase)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

gf.itos <- filter(gf, phase == "ITOS")
gf.icsa <- filter(gf, phase == "ICS-A")
gf.icsu <- filter(gf, phase == "ICS-U")

gf.itos$gfmissing <- apply(select(gf.itos, gf.items.by.phase[[1]]), 1, percentage.missing)
gf.icsa$gfmissing <- apply(select(gf.icsa, gf.items.by.phase[[2]]), 1, percentage.missing)
gf.icsu$gfmissing <- apply(select(gf.icsu, gf.items.by.phase[[3]]), 1, percentage.missing)
gf <- rbind(gf.itos,gf.icsa,gf.icsu) # combine the data back together

psych::describeBy(gf$gfmissing, group = gf$phase)

aggr(gf[c("gf.30.set.a","gf.30.set.b","gf.30.anchor","gf.60.set.a","gf.60.set.b","gf.60.anchor")], combined = TRUE,
     labels = c("30s Set A", "30s Set B", "30s Anc", "60s Set A", "60s Set B", "60s Anc")) # Visualise Missing Data Patterns

#### Gv Missing Data ####
#### Use the individual CHC factor data.frame because this excludes
#### scores of 0

psych::describe(apply(select(gv, starts_with("score_gv")), 2, percentage.missing))

gv.missing <- sapply(1:3, function(x) {
  apply(filter(select(gv, gv.items.by.phase[[x]]), gv$phase == phases[[x]]), 2, percentage.missing)
}) # create list of each phase item data missing

sapply(1:3, function(x) {
  round(psych::describe(gv.missing[[x]]),2)
}) # summary of gv missing data by phase

gv.missing <- bind_rows(lapply(gv.missing, as.data.frame.list)) # combine list into dataframe
row.names(gv.missing) <- phases # name rows
gv.missing <- rownames_to_column(gv.missing, var = "phase")

ggplot(reshape2::melt(gv.missing, "phase"), aes(x=as.numeric(variable), y = value, colour = phase)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_vline(xintercept = 21) +
  geom_vline(xintercept = 52) +
  annotate("text", x = 10, y = 40, label = "Set A Items") +
  annotate("text", x = 37, y = 40, label = "Anchor Items") +
  annotate("text", x = 64, y = 40, label = "Set B Items")

gv.itos <- filter(gv, phase == "ITOS")
gv.icsa <- filter(gv, phase == "ICS-A")
gv.icsu <- filter(gv, phase == "ICS-U")

gv.itos$gvmissing <- apply(select(gv.itos, gv.items.by.phase[[1]]), 1, percentage.missing)
gv.icsa$gvmissing <- apply(select(gv.icsa, gv.items.by.phase[[2]]), 1, percentage.missing)
gv.icsu$gvmissing <- apply(select(gv.icsu, gv.items.by.phase[[3]]), 1, percentage.missing)
gv <- rbind(gv.itos,gv.icsa,gv.icsu) # combine the data back together

psych::describeBy(gv$gvmissing, group = gv$phase)

aggr(gv[c("gv.set.a","gv.set.b","gv.anchor")], combined = TRUE,
     labels = c("Set A", "Set B", "Anchor")) # Visualise Missing Data Patterns

#### Gwm Missing Data ####
#### Use the individual CHC factor data.frame because this excludes
#### scores of 0

psych::describe(apply(select(gwm, starts_with("score_gwm")), 2, percentage.missing))

gwm.missing <- sapply(1:3, function(x) {
  apply(filter(select(gwm, gwm.items.by.phase[[x]]), gwm$phase == phases[[x]]), 2, percentage.missing)
}) # create list of each phase item data missing

sapply(1:3, function(x) {
  round(psych::describe(gwm.missing[[x]]),2)
}) # summary of gwm missing data by phase

gwm.missing <- bind_rows(lapply(gwm.missing, as.data.frame.list)) # combine list into dataframe
row.names(gwm.missing) <- phases # name rows
gwm.missing <- rownames_to_column(gwm.missing, var = "phase")

ggplot(reshape2::melt(gwm.missing, "phase"), aes(x=as.numeric(variable), y = value, colour = phase)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_vline(xintercept = 21) +
  geom_vline(xintercept = 32) +
  annotate("text", x = 10, y = 40, label = "Set A Items") +
  annotate("text", x = 26, y = 40, label = "Anchor Items") +
  annotate("text", x = 38, y = 40, label = "Set B Items")

gwm.itos <- filter(gwm, phase == "ITOS")
gwm.icsa <- filter(gwm, phase == "ICS-A")
gwm.icsu <- filter(gwm, phase == "ICS-U")

gwm.itos$gwmmissing <- apply(select(gwm.itos, gwm.items.by.phase[[1]]), 1, percentage.missing)
gwm.icsa$gwmmissing <- apply(select(gwm.icsa, gwm.items.by.phase[[2]]), 1, percentage.missing)
gwm.icsu$gwmmissing <- apply(select(gwm.icsu, gwm.items.by.phase[[3]]), 1, percentage.missing)
gwm <- rbind(gwm.itos,gwm.icsa,gwm.icsu) # combine the data back together

psych::describeBy(gwm$gwmmissing, group = gwm$phase)

aggr(gwm[c("gwm.set.a","gwm.set.b","gwm.anchor")], combined = TRUE,
     labels = c("Set A", "Set B", "Anchor")) # Visualise Missing Data Patterns

#---                                                                                           ---#
################################## Missing Data Plots and Summary ##################################
#---                                                                                           ---#

data <- merge(data, gc[,c("id","gcmissing")], by = "id", all.x = TRUE)
data <- merge(data, gf[,c("id","gfmissing")], by = "id", all.x = TRUE)
data <- merge(data, gv[,c("id","gvmissing")], by = "id", all.x = TRUE)
data <- merge(data, gwm[,c("id","gwmmissing")], by = "id", all.x = TRUE)

ggplot(reshape2::melt(data = data[,c("phase","gcmissing","gfmissing","gvmissing","gwmmissing")], id.vars = "phase", 
                      measure.vars = c("gcmissing","gfmissing","gvmissing","gwmmissing")), 
       aes(x = value, fill = variable)) +
  geom_histogram(binwidth = 25) +
  facet_grid(phase~variable) + 
  theme(axis.line=element_line()) +
  scale_y_continuous(name = "Count", breaks = seq(0,1750,250)) +
  scale_x_continuous(name = "Percentage of Items Missing", breaks = seq(0,100,25)) +
  scale_fill_discrete(name = "CHC Factor", labels = c("Gc:VL", "Gf:I", "Gv:Vz", "Gwm:Wv"))

ggplot(reshape2::melt(data = data[,c("phase","gcmissing","gfmissing","gvmissing","gwmmissing")], id.vars = "phase", 
                      measure.vars = c("gcmissing","gfmissing","gvmissing","gwmmissing")), 
       aes(x = value, fill = variable)) +
  geom_histogram(binwidth = 10) +
  facet_grid(variable~.) + 
  theme(axis.line=element_line()) +
  scale_y_continuous(name = "Count", breaks = seq(0,1750,250)) +
  scale_x_continuous(name = "Percentage of Items Missing", breaks = seq(0,100,25)) +
  scale_fill_discrete(name = "CHC Factor", labels = c("Gc:VL", "Gf:I", "Gv:Vz", "Gwm:Wv"))

ggplot(reshape2::melt(data[,c("gcmissing","gfmissing","gvmissing","gwmmissing")]), 
       aes(x = variable, y = value, colour = variable)) +
  geom_boxplot() +
  ylab("Percentage of Items Missing") +
  scale_x_discrete(name = "CHC Factor", labels = c("Gc:VL","Gf:I","Gv:Vz","Gwm:Wv")) +
  theme(legend.position = "none")

item.sets <- c("gc.set.a","gc.set.b","gc.anchor", "gf.30.set.a","gf.30.set.b","gf.30.anchor",
               "gf.60.set.a","gf.60.set.b","gf.60.anchor","gv.set.a","gv.set.b","gv.anchor",
               "gwm.set.a","gwm.set.b","gwm.anchor")

sapply(item.sets, function(x) {
  round(table(data[x])[[1]]/nrow(data[x])*100,2)
})

### Save data for simulation

write.csv(data, "simulation_data.csv", row.names = FALSE)

#---                                                                               ---#
################################## Gc Data Imputation ##################################
#---                                                                               ---#

# initial
gc.cleaned <- gc
gc.cleaned$gcmissing.initial <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.initial) # obtain missing data by participant
gc.missing1 <- select(gc.cleaned, id, gcmissing.initial)

# remove people who dont have a full set of items for anchor items
gc.cleaned <- filter(gc.cleaned, !is.na(gc.anchor)) # set new data sets
gc.cleaned$gcmissing.fullanchor <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.fullanchor)
gc.missing2 <- select(gc.cleaned, id, gcmissing.fullanchor)

# remove items that were removed as part of itos
gc.itos.removed <- c("score_gc1","score_gc2","score_gc3","score_gc20","score_gc35","score_gc47","score_gc52",
                     "timeTaken_gc1","timeTaken_gc2","timeTaken_gc3","timeTaken_gc20","timeTaken_gc35","timeTaken_gc47","timeTaken_gc52")
gc.cleaned <- select(gc.cleaned, -all_of(gc.itos.removed)) # remove unused items
gc.cleaned$gcmissing.removeditems <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.removeditems)
gc.missing3 <- select(gc.cleaned, id, gcmissing.removeditems)

# remove items that have low rates of response
gc.low.response <- apply(select(gc.cleaned, starts_with("score_gc"), starts_with("timeTaken")), 2, percentage.missing) > 75
gc.low.response <- names(gc.low.response[gc.low.response == TRUE])
gc.low.response

gc.icsa.items <- gc.icsa.items[!gc.icsa.items %in% gc.low.response]
gc.icsa.times <- gc.icsa.times[!gc.icsa.times %in% gc.low.response]
gc.icsu.items <- gc.icsu.items[!gc.icsu.items %in% gc.low.response]
gc.icsu.times <- gc.icsu.times[!gc.icsu.times %in% gc.low.response]
gc.seta.items <- gc.seta.items[!gc.seta.items %in% gc.low.response]
gc.seta.times <- gc.seta.times[!gc.seta.times %in% gc.low.response]
gc.setb.items <- gc.setb.items[!gc.setb.items %in% gc.low.response]
gc.setb.times <- gc.setb.times[!gc.setb.times %in% gc.low.response]
gc.anchor.items <- gc.anchor.items[!gc.anchor.items %in% gc.low.response]
gc.anchor.times <- gc.anchor.times[!gc.anchor.times %in% gc.low.response]

gc.cleaned <- select(gc.cleaned, -all_of(gc.low.response))
gc.cleaned$gcmissing.lowresponse <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.lowresponse)
gc.missing4 <- select(gc.cleaned, id, gcmissing.lowresponse)

# remove participants with high rates of missingness
gc.cleaned <- filter(gc.cleaned, gcmissing.lowresponse < 33)
gc.cleaned$gcmissing.removedhighmissingness <- apply(select(gc.cleaned, c(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.removedhighmissingness)
gc.missing5 <- select(gc.cleaned, id, gcmissing.removedhighmissingness)

# remove participants with missing age
gc.cleaned <- filter(gc.cleaned, !is.na(age))
gc.cleaned$gcmissing.missingage <- apply(select(gc.cleaned, c(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.missingage)
gc.missing6 <- select(gc.cleaned, id, gcmissing.missingage)

# recalculate times and scores for cleaned data
gc.cleaned <- mutate(gc.cleaned, gcscore = rowSums(select(gc.cleaned, starts_with("score_gc")), na.rm = TRUE), 
                     gctime = rowSums(select(gc.cleaned, starts_with("timeTaken_gc")), na.rm = TRUE))

# remove univariate time outliers by 1st and 99th percentile
gc.cleaned <- filter(gc.cleaned, gctime > quantile(gc.cleaned$gctime, 0.05) & gctime < quantile(gc.cleaned$gctime, 0.95))
gc.cleaned$gcmissing.univariatetime <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.univariatetime) # obtain missing data by participant
gc.missing7 <- select(gc.cleaned, id, gcmissing.univariatetime)

# remove univariate score outliers by 1st and 99th percentile
gc.cleaned <- filter(gc.cleaned, gcscore > quantile(gc.cleaned$gcscore, 0.01) & gcscore < quantile(gc.cleaned$gcscore, 0.99))
gc.cleaned$gcmissing.univariatescore <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.univariatescore) # obtain missing data by participant
gc.missing8 <- select(gc.cleaned, id, gcmissing.univariatescore)

# remove multivariate outliers
gc.lm.model <- paste("gcscore ~ age + gctime + gender + phase")
gc.lm <- lm(gc.lm.model, data = gc.cleaned)
gc.cd <- cooks.distance(gc.lm)
plot(gc.cd, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean(gc.cd, na.rm = T), col = "red")
text(x=1:length(gc.cd)+1, y=gc.cd, labels=ifelse(gc.cd>4*mean(gc.cd, na.rm = TRUE ), names(gc.cd),""), col = "red")
gc.influential <- as.numeric(names(gc.cd)[(gc.cd > 4*mean(gc.cd, na.rm = TRUE))])
gc.influential <- gc.influential[!is.na(gc.influential)]
gc.cleaned <- gc.cleaned[-gc.influential,]
gc.cleaned$gcmissing.multivariate <- apply(select(gc.cleaned, all_of(gc.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gc.cleaned$gcmissing.multivariate) # obtain missing data by participant
gc.missing9 <- select(gc.cleaned, id, gcmissing.multivariate)

# create set column
gc.cleaned$set <- "cleaned"

# create df of item data missing
gc.missing.prior <- round(as.data.frame(apply(select(gc, starts_with("score_gc")), 2, percentage.missing))) # select original Gc data
names(gc.missing.prior)[1] <- "missing"
gc.missing.prior$stage <- "Original"
gc.missing.prior <- rownames_to_column(gc.missing.prior, "item")
gc.missing.after <- round(as.data.frame(apply(select(gc.cleaned, starts_with("score_gc")), 2, percentage.missing))) # select cleaned data
names(gc.missing.after)[1] <- "missing"
gc.missing.after$stage <- "After Cleaning"
gc.missing.after <- rownames_to_column(gc.missing.after, "item")
gc.missing.items <- rbind(gc.missing.prior,gc.missing.after)
gc.missing.items$item <- as.numeric(str_extract(gc.missing.items$item, "[0-9]+"))
rm(gc.missing.after,gc.missing.prior)

# Visualise missing data by item
ggplot(gc.missing.items, aes(x = item, y = missing, colour = stage)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom") +
  labs(colour = "Stage") +
  geom_vline(xintercept = 34) +
  geom_vline(xintercept = 56) +
  annotate("text", x = 16, y = 50, label = "Set A Items") +
  annotate("text", x = 46, y = 50, label = "Anchor Items") +
  annotate("text", x = 85, y = 50, label = "Set B Items")

# combine missing data columns
gc.missing.participants <- full_join(gc.missing1, gc.missing2)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing3)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing4)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing5)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing6)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing7)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing8)
gc.missing.participants <- full_join(gc.missing.participants, gc.missing9)
rm(gc.missing1,gc.missing2,gc.missing3,gc.missing4,gc.missing5,gc.missing6,gc.missing7,gc.missing8,gc.missing9)

# visualise the changes in missing data across steps
ggplot(reshape2::melt(gc.missing.participants[,c("gcmissing.initial","gcmissing.fullanchor","gcmissing.removeditems","gcmissing.lowresponse",
                                                 "gcmissing.removedhighmissingness", "gcmissing.missingage","gcmissing.univariatetime",
                                                 "gcmissing.univariatescore","gcmissing.multivariate")]), 
       aes(x = variable, y = value, colour = variable)) +
  geom_boxplot(width = .25) +
  ylab("Percentage of Items Missing") +
  scale_x_discrete(name = "", labels = ggplot_addline(c("Full Data","Full Anchor","Items Removed","Low Item Response Rate",
                                                        "High Missingness Removed","Missing Age Removed",
                                                        "Outliers by Time","Outliers by Score",
                                                        "Multivariate Outliers"))) +
  theme(legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text")

# Recalculate the set missingness
gc.cleaned <- mutate(gc.cleaned, 
                     gc.set.a = ceiling(rowSums(is.na(gc.cleaned[gc.seta.items]))/length(gc.seta.items)),
                     gc.set.b = ceiling(rowSums(is.na(gc.cleaned[gc.setb.items]))/length(gc.setb.items)),
                     gc.anchor = ceiling(rowSums(is.na(gc.cleaned[gc.anchor.items]))/length(gc.anchor.items)))

# Recode missingness # NA = missing data in set, 1 = full data
gc.cleaned$gc.set.a <- gc.cleaned$gc.set.a %>% na_if(1) %>% recode(`0` = 1)
gc.cleaned$gc.set.b <- gc.cleaned$gc.set.b %>% na_if(1) %>% recode(`0` = 1)
gc.cleaned$gc.anchor <- gc.cleaned$gc.anchor %>% na_if(1) %>% recode(`0` = 1)

# Visualise Missing Data Patterns
aggr(gc.cleaned[c("gc.set.a","gc.set.b","gc.anchor")], combined = TRUE,
     labels = c("Set A", "Set B", "Anchor")) 

# Set Up Data Set
gc.imputed <- select(gc.cleaned, id, phase, age, age.group, gender, nationality, device, all_of(gc.icsu.items), all_of(gc.icsu.times)) # create data set
gc.imputed[gc.icsu.items] <- lapply(gc.imputed[gc.icsu.items], factor)

# prepare prediction model
gc.qp <- quickpred(gc.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# flag logged events items
gc.ini <- mice(gc.imputed, maxit = 1, pred = gc.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)
gc.outlist.logged <- as.character(gc.ini$loggedEvents[, "out"])
gc.outlist.logged

# flag outflux problem items
gc.flux <- flux(gc.imputed)
gc.fluxplot <- fluxplot(gc.imputed, main = NULL)
gc.outlist.outflux <- row.names(gc.flux)[gc.flux$outflux < 0.3]
gc.outlist.outflux

# remove problem items for mi
if (exists("gc.outlist.logged") & exists("gc.outlist.outflux")) {
  gc.outlist <- unique(c(gc.outlist.logged, gc.outlist.outflux))
} else if (exists("gc.outlist.logged")) {
  gc.outlist <- gc.outlist.logged
} else if (exists("gc.outlist.outflux")) {
  gc.outlist <- gc.outlist.outflux
} else {
  print("gc outlist not created")
}

if (exists("gc.outlist")) {
  gc.outlist <- as.numeric(str_extract(gc.outlist, "[0-9]+"))
  gc.outlist <- c(paste0("score_gc",gc.outlist),paste0("timeTaken_gc",gc.outlist))
  gc.outlist <- unique(gc.outlist)
  gc.imputed <- select(gc.imputed, -all_of(gc.outlist))
  gc.icsa.items <- gc.icsa.items[!gc.icsa.items %in% gc.outlist]
  gc.icsa.times <- gc.icsa.times[!gc.icsa.times %in% gc.outlist]
  gc.icsu.items <- gc.icsu.items[!gc.icsu.items %in% gc.outlist]
  gc.icsu.times <- gc.icsu.times[!gc.icsu.times %in% gc.outlist]
  gc.seta.items <- gc.seta.items[!gc.seta.items %in% gc.outlist]
  gc.seta.times <- gc.seta.times[!gc.seta.times %in% gc.outlist]
  gc.setb.items <- gc.setb.items[!gc.setb.items %in% gc.outlist]
  gc.setb.times <- gc.setb.times[!gc.setb.times %in% gc.outlist]
  gc.anchor.items <- gc.anchor.items[!gc.anchor.items %in% gc.outlist]
  gc.anchor.times <- gc.anchor.times[!gc.anchor.times %in% gc.outlist]
}

# outflux after removal
if (exists("gc.outlist.outflux")) {
  gc.fluxplot.2 <- fluxplot(gc.imputed, main = NULL)
}

# prepare prediction model
gc.qp <- quickpred(gc.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# conduct mi - impute data using random forests
gc.mice <- mice(gc.imputed, pred = gc.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)

# create data sets for use in analyses
gc.imputed <- NULL # create empty variable
for (i in 1:5) gc.imputed[[i]] <- mice::complete(gc.mice, action = i, inc = FALSE) # put imputed data into list
for (i in 1:5) gc.imputed[[i]][gc.icsu.items] <- lapply(gc.imputed[[i]][gc.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric
for (i in 1:5) gc.imputed[[i]]$set <-  paste0("imp",i)
gc.imputed <- lapply(gc.imputed, function(x) {
  x <- mutate(x, gcscore = rowSums(select(x, starts_with("score_gc")), na.rm = TRUE))
  x <- mutate(x, gctime = rowSums(select(x, starts_with("timeTaken_gc")), na.rm = TRUE))
})
gc.backup <- gc.imputed

# recalculate times and scores for cleaned data
gc.cleaned <- mutate(gc.cleaned, gcscore = rowSums(select(gc.cleaned, starts_with("score_gc")), na.rm = TRUE), 
                     gctime = rowSums(select(gc.cleaned, starts_with("timeTaken_gc")), na.rm = TRUE))

# Combine the results together
gc.comparison <- rbind(gc.imputed[[1]][c("set","gcscore","gctime")],gc.cleaned[c("set","gcscore","gctime")])
gc.comparison <- rbind(gc.imputed[[2]][c("set","gcscore","gctime")],gc.comparison)
gc.comparison <- rbind(gc.imputed[[3]][c("set","gcscore","gctime")],gc.comparison)
gc.comparison <- rbind(gc.imputed[[4]][c("set","gcscore","gctime")],gc.comparison)
gc.comparison <- rbind(gc.imputed[[5]][c("set","gcscore","gctime")],gc.comparison)

# Boxplot of total score
ggplot(gc.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gcscore, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(30,length(gf.icsu.items),10)), limits = c(30,length(gf.icsu.items)), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# Boxplot of time
ggplot(gc.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gctime/60, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(5,max(gc.cleaned$gctime)/60,5)), limits = c(5,max(gc.cleaned$gctime)/60), name = "Total Time") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# visualise proportion correct
gc.cleaned.correct <- sapply(gc.cleaned[gc.icsu.items], function(x) {(1 - (table(x)[[1]]/length(na.omit(x)))) * 100 })
gc.imputed.correct <- lapply(gc.imputed, function(x) { 
  sapply(x[gc.icsu.items], function(y) {
    (1 - (table(y)[[1]]/length(na.omit(y)))) * 100 
  })
})
gc.combined.correct <- data.frame(cleaned = gc.cleaned.correct, imp1 = gc.imputed.correct[[1]], imp2 = gc.imputed.correct[[2]],
                                  imp3 = gc.imputed.correct[[3]], imp4 = gc.imputed.correct[[4]], imp5 = gc.imputed.correct[[5]])
row.names(gc.combined.correct) <- str_remove(row.names(gc.combined.correct), "score_gc")
gc.combined.correct <- rownames_to_column(gc.combined.correct, var = "item")
gc.combined.correct <- pivot_longer(data = gc.combined.correct, cols = -item, names_to = "set", values_to = "correct")
gc.combined.correct$item <- factor(gc.combined.correct$item, levels = c(sort(unique(as.numeric(gc.combined.correct$item)))))
gc.combined.correct$set <- factor(gc.combined.correct$set)
ggplot(gc.combined.correct, aes(x = item, y = correct, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gc.combined.correct$item), by = 2)) +
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Percentage Correct") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# visualise performance by age
gc.cleaned.age <- select(gc.cleaned, age, gcscore)
gc.imputed.age <- lapply(gc.imputed, function(x) { select(x, gcscore)})
gc.combined.age <- data.frame(age = gc.cleaned.age$age, cleaned = gc.cleaned.age$gcscore, imp1 = gc.imputed.age[[1]]$gcscore, 
                              imp2 = gc.imputed.age[[2]]$gcscore, imp3 = gc.imputed.age[[3]]$gcscore, imp4 = gc.imputed.age[[4]]$gcscore, 
                              imp5 = gc.imputed.age[[5]]$gcscore)
gc.combined.age <- pivot_longer(data = gc.combined.age, cols = -age, names_to = "set", values_to = "score")
gc.combined.age$set <- factor(gc.combined.age$set)
ggplot(gc.combined.age, aes(x = age, y = score, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(6,90,2), limits = c(6,90), name = "Age") +
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Total Score") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  geom_smooth(method = "lm")

#---                                                                           ---#
################################## Gc Reliability ##################################
#---                                                                          ---#

gc.reliability <- NULL
for (i in 1:5) gc.reliability[[i]] <- psych::alpha(gc.imputed[[i]][gc.icsu.items], warnings = FALSE)$total
gc.reliability <- sapply(gc.reliability, function(x) { x[,"raw_alpha"] })
mean(gc.reliability)

#---                                                                              ---#
################################## Gc Rasch Analysis ##################################
#---                                                                              ---#

# run rasch
gc.rasch <- list()
for (i in 1:5) gc.rasch[[i]] <- rasch.outcomes(chcability = "gc", raschobject = gc.imputed[[i]][gc.icsu.items], raschnumber = 1, imp = i)

# person fit
gc.rasch.person <- which(personfit(gc.rasch[[1]][[1]])$Zh < -2)
gc.rasch.person <- append(gc.rasch.person, which(personfit(gc.rasch[[1]][[1]])$Zh > 2))
gc.rasch.person
gc.imputed[[1]][gc.rasch.person,]$age

# wright map
WrightMap::wrightMap(fscores(gc.rasch[[1]][[1]]), main.title = NA)

# determine poor item fits
gc.poor.fit <- NULL
for (i in 1:5) gc.poor.fit[[i]] <- rasch.poor.fit(dataset = gc.imputed[[i]], raschobject = gc.rasch[[i]][[2]], imp = i)

# analyse the NA items
na.performance <- lapply(gc.poor.fit[[1]]$na,function(x) {table(gc.imputed[[1]][[x]], gc.imputed[[1]]$age.group)})
names(na.performance) <- gc.poor.fit[[1]]$na
na.performance
gc.rasch.removal <- c(7, 9, 11, 4, 5, 8) # this is a qualitative decision
gc.rasch.removal <- c(paste0("score_gc",gc.rasch.removal),paste0("timeTaken_gc",gc.rasch.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gc.imputed[[i]] <- select(gc.imputed[[i]], -all_of(gc.rasch.removal))
# remove.items.from.sets(chcability = "gc", items = gc.rasch.removal)
gc.icsa.items <- gc.icsa.items[!gc.icsa.items %in% gc.rasch.removal]
gc.icsa.times <- gc.icsa.times[!gc.icsa.times %in% gc.rasch.removal]
gc.icsu.items <- gc.icsu.items[!gc.icsu.items %in% gc.rasch.removal]
gc.icsu.times <- gc.icsu.times[!gc.icsu.times %in% gc.rasch.removal]
gc.seta.items <- gc.seta.items[!gc.seta.items %in% gc.rasch.removal]
gc.seta.times <- gc.seta.times[!gc.seta.times %in% gc.rasch.removal]
gc.setb.items <- gc.setb.items[!gc.setb.items %in% gc.rasch.removal]
gc.setb.times <- gc.setb.times[!gc.setb.times %in% gc.rasch.removal]
gc.anchor.items <- gc.anchor.items[!gc.anchor.items %in% gc.rasch.removal]
gc.anchor.times <- gc.anchor.times[!gc.anchor.times %in% gc.rasch.removal]

# run rasch 2
gc.rasch.2 <- list()
for (i in 1:5) gc.rasch.2[[i]] <- rasch.outcomes(chcability = "gc", raschobject = gc.imputed[[i]][gc.icsu.items], raschnumber = 2, imp = i)

# determine poor item fits 2
gc.poor.fit.2 <- NULL
for (i in 1:5) gc.poor.fit.2[[i]] <- rasch.poor.fit(dataset = gc.imputed[[i]], raschobject = gc.rasch.2[[i]][[2]], imp = i)

# analyse the poor fitting items outside of NA
below.01.performance <- lapply(gc.poor.fit[[1]]$below.01,function(x) {table(gc.imputed[[1]][[x]], gc.imputed[[1]]$age.group)})
names(below.01.performance) <- gc.poor.fit[[1]]$below.01
below.01.performance
gc.rasch.removal.2 <- c(54, 55, 42, 75, 79, 82, 83, 84, 88, 90, 98, 65) # this is a qualitative decision
gc.rasch.removal.2 <- c(paste0("score_gc",gc.rasch.removal.2),paste0("timeTaken_gc",gc.rasch.removal.2))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gc.imputed[[i]] <- select(gc.imputed[[i]], -all_of(gc.rasch.removal.2))
# remove.items.from.sets(chcability = "gc", items = gc.rasch.removal)
gc.icsa.items <- gc.icsa.items[!gc.icsa.items %in% gc.rasch.removal.2]
gc.icsa.times <- gc.icsa.times[!gc.icsa.times %in% gc.rasch.removal.2]
gc.icsu.items <- gc.icsu.items[!gc.icsu.items %in% gc.rasch.removal.2]
gc.icsu.times <- gc.icsu.times[!gc.icsu.times %in% gc.rasch.removal.2]
gc.seta.items <- gc.seta.items[!gc.seta.items %in% gc.rasch.removal.2]
gc.seta.times <- gc.seta.times[!gc.seta.times %in% gc.rasch.removal.2]
gc.setb.items <- gc.setb.items[!gc.setb.items %in% gc.rasch.removal.2]
gc.setb.times <- gc.setb.times[!gc.setb.times %in% gc.rasch.removal.2]
gc.anchor.items <- gc.anchor.items[!gc.anchor.items %in% gc.rasch.removal.2]
gc.anchor.times <- gc.anchor.times[!gc.anchor.times %in% gc.rasch.removal.2]

# run rasch 3
gc.rasch.3 <- list()
for (i in 1:5) gc.rasch.3[[i]] <- rasch.outcomes(chcability = "gc", raschobject = gc.imputed[[i]][gc.icsu.items], raschnumber = 3, imp = i)

# determine poor item fits 2
gc.poor.fit.3 <- NULL
for (i in 1:5) gc.poor.fit.3[[i]] <- rasch.poor.fit(dataset = gc.imputed[[i]], raschobject = gc.rasch.3[[i]][[2]], imp = i)

#---                                                                                ---#
################################## Gc Mokken Analysis ##################################
#---                                                                               ---#

# conduct first mokken
gc.mokken <- list()
for (i in 1:5) gc.mokken[[i]] <- mokken.outcomes(chcability = "gc", mokkenobject = gc.imputed[[i]][gc.icsu.items], mokkennumber = 1, imp = i)

# determine poorly ordered items and put them in a vector for removal
gc.mokken.prob <- NULL
for (i in 1:5) gc.mokken.prob[[i]] <- rownames(gc.mokken[[i]]$Hi)[which(gc.mokken[[i]]$Hi < .3)]
gc.mokken.prob
gc.mokken.prob <- sort(unique(unlist(gc.mokken.prob)))
gc.mokken.prob

gc.mokken.removal <- c(14, 48, 63, 71, 76, 80) # this is a qualitative decision
gc.mokken.removal <- c(paste0("score_gc",gc.mokken.removal),paste0("timeTaken_gc",gc.mokken.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gc.imputed[[i]] <- select(gc.imputed[[i]], -all_of(gc.mokken.removal))
# remove.items.from.sets(chcability = "gc", items = gc.mokken.removal)
gc.icsa.items <- gc.icsa.items[!gc.icsa.items %in% gc.mokken.removal]
gc.icsa.times <- gc.icsa.times[!gc.icsa.times %in% gc.mokken.removal]
gc.icsu.items <- gc.icsu.items[!gc.icsu.items %in% gc.mokken.removal]
gc.icsu.times <- gc.icsu.times[!gc.icsu.times %in% gc.mokken.removal]
gc.seta.items <- gc.seta.items[!gc.seta.items %in% gc.mokken.removal]
gc.seta.times <- gc.seta.times[!gc.seta.times %in% gc.mokken.removal]
gc.setb.items <- gc.setb.items[!gc.setb.items %in% gc.mokken.removal]
gc.setb.times <- gc.setb.times[!gc.setb.times %in% gc.mokken.removal]
gc.anchor.items <- gc.anchor.items[!gc.anchor.items %in% gc.mokken.removal]
gc.anchor.times <- gc.anchor.times[!gc.anchor.times %in% gc.mokken.removal]

# conduct second mokken
gc.mokken.2 <- list()
for (i in 1:5) gc.mokken.2[[i]] <- mokken.outcomes(chcability = "gc", mokkenobject = gc.imputed[[i]][gc.icsu.items], mokkennumber = 2, imp = i)

#---                                                                                 ---#
################################## Gc Local Dependency ##################################
#---                                                                                 ---#

gc.local.dependency <- list()
for (i in 1:5) gc.local.dependency[[i]] <- residuals(gc.rasch[[i]][[1]], type = "Q3", digits = 2, suppress = +.2)

# determine items with high local dependency
gc.local.dependency.removal <- "score_gc37"

# remove the items that do not fit well to the rasch model
for (i in 1:5) gc.imputed[[i]] <- select(gc.imputed[[i]], -all_of(gc.local.dependency.removal))
# remove.items.from.sets(chcability = "gc", items = gc.local.dependency.removal)
gc.icsa.items <- gc.icsa.items[!gc.icsa.items %in% gc.local.dependency.removal]
gc.icsa.times <- gc.icsa.times[!gc.icsa.times %in% gc.local.dependency.removal]
gc.icsu.items <- gc.icsu.items[!gc.icsu.items %in% gc.local.dependency.removal]
gc.icsu.times <- gc.icsu.times[!gc.icsu.times %in% gc.local.dependency.removal]
gc.seta.items <- gc.seta.items[!gc.seta.items %in% gc.local.dependency.removal]
gc.seta.times <- gc.seta.times[!gc.seta.times %in% gc.local.dependency.removal]
gc.setb.items <- gc.setb.items[!gc.setb.items %in% gc.local.dependency.removal]
gc.setb.times <- gc.setb.times[!gc.setb.times %in% gc.local.dependency.removal]
gc.anchor.items <- gc.anchor.items[!gc.anchor.items %in% gc.local.dependency.removal]
gc.anchor.times <- gc.anchor.times[!gc.anchor.times %in% gc.local.dependency.removal]

#---                                                                                              ---#
################################## Gc Differential Item Functioning ##################################
#---                                                                                              ---#

gc.dif.gender <- list()
for (i in 1:5) gc.dif.gender[[i]] <- difLogistic(Data = gc.imputed[[i]][gc.icsu.items], group = gc.imputed[[i]]$gender, focal.name = "m", alpha = .01)
for (i in 1:5) plot(gc.dif.gender[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gc.dif.nationality <- list()
for (i in 1:5) gc.dif.nationality[[i]] <- difLogistic(Data = gc.imputed[[i]][gc.icsu.items], group = gc.imputed[[i]]$nationality, focal.name = "au", alpha = .01)
for (i in 1:5) plot(gc.dif.nationality[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gc.dif.device <- list()
for (i in 1:5) gc.dif.device[[i]] <- difLogistic(Data = gc.imputed[[i]][gc.icsu.items], group = gc.imputed[[i]]$nationality, focal.name = "au", alpha = .01)
for (i in 1:5) plot(gc.dif.device[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

# determine dif items
gc.dif.items <- NULL
for (i in 1:5) gc.dif.items[[i]] <- { 
  gender <- gc.dif.gender[[i]]$DIFitems
  nationality <- gc.dif.nationality[[i]]$DIFitems
  device <- gc.dif.device[[i]]$DIFitems
  list(gender,nationality,device)
}

gc.dif.items <- sort(as.numeric(unique(unlist(gc.dif.items))), na.last = NA)
gc.dif.items
gc.dif.items <- gc.icsu.items[gc.dif.items]
gc.dif.items

for (i in gc.dif.items) { plot(gc.dif.gender[[1]], plot = "itemCurve", item = i) }

#---                                                                                             ---#
################################## Gc Confirmatory Factor Analysis ##################################
#---                                                                                            ---#

# run the CFA on the multivariate imputation data sets
for (i in 1:5) gc.imputed[[i]][gc.icsu.items] <- lapply(gc.imputed[[i]][gc.icsu.items], function(x) { x <- ordered(x, levels = c(0,1)) }) # Change all items to ordered factor
gc.model <- paste("gc =~ ", paste(gc.icsu.items, collapse = "+"), sep = "") # establish model
gc.cfa.mi <- runMI(gc.model, data = gc.imputed, fun = "cfa", meanstructure = TRUE) # run cfa --> takes a very long time
gc.cfa.mi.summary <- summary(gc.cfa.mi)
gc.cfa.mi.fit <- fitMeasures(gc.cfa.mi, fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))

# determine what cfas were successful
gc.cfa.success <- NULL
for (i in 1:5) gc.cfa.success[i] <- !is.null(gc.cfa.mi@miList[[i]])

# based on the cfas that were successful, run them 1 by 1 for diagnostic and reporting purposes, fit the cfa
gc.cfa <- list()
for (i in 1:5) gc.cfa[[i]] <- if(gc.cfa.success[i] == TRUE) {cfa(gc.model, data = gc.imputed[[i]])}

# store the cfa summary for item selection later
gc.cfa.summary <- NULL
for (i in 1:5) gc.cfa.summary[[i]] <- if(gc.cfa.success[i] == TRUE) { summary(object = gc.cfa[[i]], standardized = TRUE) }

# print and save the parameters using predefined function
gc.cfa.parameters <- NULL
for (i in 1:5) gc.cfa.parameters[[i]] <- if(gc.cfa.success[i] == TRUE) {
  cfa.parameters(chcability = "gc", cfaobject = gc.cfa[[i]], cfanumber = 1, imp = i)}

# print and store the overall fit
gc.cfa.fit <- NULL
for (i in 1:5) gc.cfa.fit[[i]] <- if(gc.cfa.success[i] == TRUE) {
  fitMeasures(gc.cfa[[i]], fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
}
gc.cfa.fit

# determine low loading items and put them in a vector for removal
gc.cfa.poor.fit <- list()
for (i in 1:5) gc.cfa.poor.fit[[i]] <- if(gc.cfa.success[i] == TRUE) {
  select(filter(gc.cfa.summary[[i]]$PE, lhs == "gc", std.lv < .3, rhs != ""), rhs)$rhs } else {NA}
gc.cfa.poor.fit
gc.cfa.poor.fit <- sort(unique(unlist(gc.cfa.poor.fit)))
gc.cfa.poor.fit

#---                                                                                     ---#
################################## Gc Final Rasch Analysis ##################################
#---                                                                                    ---#

for (i in 1:5) gc.imputed[[i]][gc.icsu.items] <- lapply(gc.imputed[[i]][gc.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric

# run rasch 4
gc.rasch.4 <- list()
for (i in 1:5) gc.rasch.4[[i]] <- rasch.outcomes(chcability = "gc", raschobject = gc.imputed[[i]][gc.icsu.items], raschnumber = 4, imp = i)

# parameters
gc.classic.parameters <- list()
for (i in 1:5) gc.classic.parameters[[i]] <- coef(gc.rasch.4[[i]][[1]], IRTpars = TRUE, simplify = TRUE)$items
gc.classic.parameters

gc.parameters <- list()
for (i in 1:5) gc.parameters[[i]] <- coef(gc.rasch.4[[i]][[1]], simplify = TRUE)$items
gc.parameters

# pool the parameters
gc.rasch.average <- list()
for (i in 1:5) gc.rasch.average[[i]] <- as.data.frame(coef(gc.rasch.4[[i]][[1]], printSE = TRUE, as.data.frame = TRUE))
for (i in 1:5) gc.rasch.average[[i]] <- rownames_to_column(gc.rasch.average[[i]], var = "item")
for (i in 1:5) gc.rasch.average[[i]] <- filter(gc.rasch.average[[i]], par > -10, par < 10, par != 1)
gc.item.list <- lapply(gc.rasch.average, function(x) { x$item})
gc.par.list <- lapply(gc.rasch.average, function(x) { x$par })
gc.separ.list <- lapply(gc.rasch.average, function(x) { x$SE })
gc.rasch.average <- averageMI(gc.par.list,gc.separ.list,as.data.frame = TRUE)
gc.rasch.average$item <- str_remove(gc.item.list[[1]], "score_")
gc.rasch.average$item <- str_remove(gc.rasch.average$item, ".d")
gc.rasch.average$a1 <- 1
gc.rasch.average$g <- 0
gc.rasch.average$u <- 1
colnames(gc.rasch.average)[1] <- "d"
gc.rasch.average <- select(gc.rasch.average, item, a1, d, g, u)
gc.rasch.average <- filter(gc.rasch.average, grepl(pattern = "gc", item))

# get f scores
gc.fscores <- gc.imputed[[1]]
gc.fscores <- mutate(gc.fscores, gcscore = rowSums(select(gc.fscores, starts_with("score_gc"))))
gc.fscores$fscores <- fscores(object = gc.rasch.4[[1]][[1]])
gc.fscores <- select(gc.fscores, -starts_with("score_gc"), -starts_with("timeTaken"), -set)
gc.wisc <- select(filter(data, gcscore > 0, wisc.fsiq.comp > 0), id, starts_with("wisc"))
gc.fscores <- left_join(gc.wisc,gc.fscores)

# save parameters, remaining items and f scores
write.table(gc.rasch.average, file = "gc_rasch_parameters.csv", sep = ",", quote = FALSE, row.names = F)
write(gc.icsu.items, file = "gcitems.txt")
write.table(gc.fscores, file = "gc_fscores.csv", sep = ",", quote = FALSE, row.names = F)

#---                                                                               ---#
################################## Gf Data Imputation ##################################
#---                                                                               ---#

# initial
gf.cleaned <- gf
gf.cleaned$gfmissing.initial <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.initial) # obtain missing data by participant
gf.missing1 <- select(gf.cleaned, id, gfmissing.initial)

# remove people who dont have a full set of items for anchor items
gf.cleaned <- filter(gf.cleaned, !is.na(gf.30.anchor) & !is.na(gf.60.anchor)) # set new data sets
gf.cleaned$gfmissing.fullanchor <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.fullanchor)
gf.missing2 <- select(gf.cleaned, id, gfmissing.fullanchor)

# remove items that were removed as part of itos
gf.itos.removed <- c("score_gf22", "timeTaken_gf22", "score_gf24", "timeTaken_gf24", "score_gf25", "timeTaken_gf25",
                     "score_gf30", "timeTaken_gf30", "score_gf31", "timeTaken_gf31", "score_gf33", "timeTaken_gf33")

gf.icsa.items <- gf.icsa.items[!gf.icsa.items %in% gf.itos.removed]
gf.icsa.times <- gf.icsa.times[!gf.icsa.times %in% gf.itos.removed]
gf.icsu.items <- gf.icsu.items[!gf.icsu.items %in% gf.itos.removed]
gf.icsu.times <- gf.icsu.times[!gf.icsu.times %in% gf.itos.removed]
gf.seta.items <- gf.seta.items[!gf.seta.items %in% gf.itos.removed]
gf.seta.times <- gf.seta.times[!gf.seta.times %in% gf.itos.removed]
gf.setb.items <- gf.setb.items[!gf.setb.items %in% gf.itos.removed]
gf.setb.times <- gf.setb.times[!gf.setb.times %in% gf.itos.removed]
gf.anchor.items <- gf.anchor.items[!gf.anchor.items %in% gf.itos.removed]
gf.anchor.times <- gf.anchor.times[!gf.anchor.times %in% gf.itos.removed]

gf.cleaned <- select(gf.cleaned, -all_of(gf.itos.removed)) # remove unused items
gf.cleaned$gfmissing.removeditems <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.removeditems)
gf.missing3 <- select(gf.cleaned, id, gfmissing.removeditems)

# remove items that have low rates of response
gf.low.response <- apply(select(gf.cleaned, starts_with("score_gf"), starts_with("timeTaken")), 2, percentage.missing) > 75
gf.low.response <- names(gf.low.response[gf.low.response == TRUE])
gf.low.response

gf.icsa.items <- gf.icsa.items[!gf.icsa.items %in% gf.low.response]
gf.icsa.times <- gf.icsa.times[!gf.icsa.times %in% gf.low.response]
gf.icsu.items <- gf.icsu.items[!gf.icsu.items %in% gf.low.response]
gf.icsu.times <- gf.icsu.times[!gf.icsu.times %in% gf.low.response]
gf.seta.items <- gf.seta.items[!gf.seta.items %in% gf.low.response]
gf.seta.times <- gf.seta.times[!gf.seta.times %in% gf.low.response]
gf.setb.items <- gf.setb.items[!gf.setb.items %in% gf.low.response]
gf.setb.times <- gf.setb.times[!gf.setb.times %in% gf.low.response]
gf.anchor.items <- gf.anchor.items[!gf.anchor.items %in% gf.low.response]
gf.anchor.times <- gf.anchor.times[!gf.anchor.times %in% gf.low.response]

gf.cleaned <- select(gf.cleaned, -all_of(gf.low.response))
gf.cleaned$gfmissing.lowresponse <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.lowresponse)
gf.missing4 <- select(gf.cleaned, id, gfmissing.lowresponse)

# remove participants with high rates of missingness
gf.cleaned <- filter(gf.cleaned, gfmissing.lowresponse < 33)
gf.cleaned$gfmissing.removedhighmissingness <- apply(select(gf.cleaned, c(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.removedhighmissingness)
gf.missing5 <- select(gf.cleaned, id, gfmissing.removedhighmissingness)

# remove participants with missing age
gf.cleaned <- filter(gf.cleaned, !is.na(age))
gf.cleaned$gfmissing.missingage <- apply(select(gf.cleaned, c(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.missingage)
gf.missing6 <- select(gf.cleaned, id, gfmissing.missingage)

# recalculate times and scores for cleaned data
gf.cleaned <- mutate(gf.cleaned, gfscore = rowSums(select(gf.cleaned, starts_with("score_gf")), na.rm = TRUE), 
                     gftime = rowSums(select(gf.cleaned, starts_with("timeTaken_gf")), na.rm = TRUE))

# remove univariate time outliers by 1st and 99th percentile
gf.cleaned <- filter(gf.cleaned, gftime > quantile(gf.cleaned$gftime, 0.05) & gftime < quantile(gf.cleaned$gftime, 0.95))
gf.cleaned$gfmissing.univariatetime <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.univariatetime) # obtain missing data by participant
gf.missing7 <- select(gf.cleaned, id, gfmissing.univariatetime)

# remove univariate score outliers by 1st and 99th percentile
gf.cleaned <- filter(gf.cleaned, gfscore > quantile(gf.cleaned$gfscore, 0.01) & gfscore < quantile(gf.cleaned$gfscore, 0.99))
gf.cleaned$gfmissing.univariatescore <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.univariatescore) # obtain missing data by participant
gf.missing8 <- select(gf.cleaned, id, gfmissing.univariatescore)

# remove multivariate outliers
gf.lm.model <- paste("gfscore ~ age + gftime + gender + phase")
gf.lm <- lm(gf.lm.model, data = gf.cleaned)
gf.cd <- cooks.distance(gf.lm)
plot(gf.cd, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean(gf.cd, na.rm = T), col = "red")
text(x=1:length(gf.cd)+1, y=gf.cd, labels=ifelse(gf.cd>4*mean(gf.cd, na.rm = TRUE ), names(gf.cd),""), col = "red")
gf.influential <- as.numeric(names(gf.cd)[(gf.cd > 4*mean(gf.cd, na.rm = TRUE))])
gf.influential <- gf.influential[!is.na(gf.influential)]
gf.cleaned <- gf.cleaned[-gf.influential,]
gf.cleaned$gfmissing.multivariate <- apply(select(gf.cleaned, all_of(gf.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gf.cleaned$gfmissing.multivariate) # obtain missing data by participant
gf.missing9 <- select(gf.cleaned, id, gfmissing.multivariate)

# create set column
gf.cleaned$set <- "cleaned"

# create df of item data missing
gf.missing.prior <- round(as.data.frame(apply(select(gf, starts_with("score_gf")), 2, percentage.missing))) # select original gf data
names(gf.missing.prior)[1] <- "missing"
gf.missing.prior$stage <- "Original"
gf.missing.prior <- rownames_to_column(gf.missing.prior, "item")
gf.missing.after <- round(as.data.frame(apply(select(gf.cleaned, starts_with("score_gf")), 2, percentage.missing))) # select cleaned data
names(gf.missing.after)[1] <- "missing"
gf.missing.after$stage <- "After Cleaning"
gf.missing.after <- rownames_to_column(gf.missing.after, "item")
gf.missing.items <- rbind(gf.missing.prior,gf.missing.after)
gf.missing.items$item <- as.numeric(str_extract(gf.missing.items$item, "[0-9]+"))
rm(gf.missing.after,gf.missing.prior)

# Visualise missing data by item
ggplot(gf.missing.items, aes(x = item, y = missing, colour = stage)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,88,2))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom") +
  labs(colour = "Stage")

# combine missing data columns
gf.missing.participants <- full_join(gf.missing1, gf.missing2)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing3)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing4)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing5)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing6)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing7)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing8)
gf.missing.participants <- full_join(gf.missing.participants, gf.missing9)
rm(gf.missing1,gf.missing2,gf.missing3,gf.missing4,gf.missing5,gf.missing6,gf.missing7,gf.missing8,gf.missing9)

# visualise the changes in missing data across steps
ggplot(reshape2::melt(gf.missing.participants[,c("gfmissing.initial","gfmissing.fullanchor","gfmissing.removeditems","gfmissing.lowresponse",
                                                 "gfmissing.removedhighmissingness", "gfmissing.missingage","gfmissing.univariatetime",
                                                 "gfmissing.univariatescore","gfmissing.multivariate")]), 
       aes(x = variable, y = value, colour = variable)) +
  geom_boxplot(width = .25) +
  ylab("Percentage of Items Missing") +
  scale_x_discrete(name = "Stage of Data Cleaning", labels = ggplot_addline(c("Full Data","Full Anchor","Items Removed","Low Response Items",
                                                                              "High Missingness Removed","Missing Age Removed",
                                                                              "Outliers by Time","Outliers by Score",
                                                                              "Multivariate Outliers"))) +
  theme(legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text")

# Recalculate the set missingness
gf.cleaned <- mutate(gf.cleaned, 
                     gf.set.a = ceiling(rowSums(is.na(gf.cleaned[gf.seta.items]))/length(gf.seta.items)),
                     gf.set.b = ceiling(rowSums(is.na(gf.cleaned[gf.setb.items]))/length(gf.setb.items)),
                     gf.anchor = ceiling(rowSums(is.na(gf.cleaned[gf.anchor.items]))/length(gf.anchor.items)))

# Recode missingness # NA = missing data in set, 1 = full data
gf.cleaned$gf.set.a <- gf.cleaned$gf.set.a %>% na_if(1) %>% recode(`0` = 1)
gf.cleaned$gf.set.b <- gf.cleaned$gf.set.b %>% na_if(1) %>% recode(`0` = 1)
gf.cleaned$gf.anchor <- gf.cleaned$gf.anchor %>% na_if(1) %>% recode(`0` = 1)

# Visualise Missing Data Patterns
aggr(gf.cleaned[c("gf.30.set.a","gf.30.set.b","gf.30.anchor","gf.60.set.a","gf.60.set.b","gf.60.anchor")], combined = TRUE,
     labels = c("30s SetA", "30s SetB", "30s Anc", "60s SetA", "60s SetB", "60s Anc")) 

# Set Up Data Set
gf.imputed <- select(gf.cleaned, id, phase, age, age.group, gender, nationality, device, all_of(gf.icsu.items), all_of(gf.icsu.times)) # create data set
gf.imputed[gf.icsu.items] <- lapply(gf.imputed[gf.icsu.items], factor)

# prepare prediction model
gf.qp <- quickpred(gf.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# flag logged events items
gf.ini <- mice(gf.imputed, maxit = 1, pred = gf.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)
gf.outlist.logged <- as.character(gf.ini$loggedEvents[, "out"])

# flag outflux problem items
gf.flux <- flux(gf.imputed)
gf.fluxplot <- fluxplot(gf.imputed, main = NULL)
gf.outlist.outflux <- row.names(gf.flux)[gf.flux$outflux < 0.3]

# remove problem items for mi
if (exists("gf.outlist.logged") & exists("gf.outlist.outflux")) {
  gf.outlist <- unique(c(gf.outlist.logged, gf.outlist.outflux))
} else if (exists("gf.outlist.logged")) {
  gf.outlist <- gf.outlist.logged
} else if (exists("gf.outlist.outflux")) {
  gf.outlist <- gf.outlist.outflux
} else {
  print("gf outlist not created")
}

if (exists("gf.outlist")) {
  gf.outlist <- as.numeric(str_extract(gf.outlist, "[0-9]+"))
  gf.outlist <- c(paste0("score_gf",gf.outlist),paste0("timeTaken_gf",gf.outlist))
  gf.outlist <- unique(gf.outlist)
  gf.imputed <- select(gf.imputed, -all_of(gf.outlist))
  gf.cleaned <- select(gf.cleaned, -all_of(gf.outlist))
  gf.icsa.items <- gf.icsa.items[!gf.icsa.items %in% gf.outlist]
  gf.icsa.times <- gf.icsa.times[!gf.icsa.times %in% gf.outlist]
  gf.icsu.items <- gf.icsu.items[!gf.icsu.items %in% gf.outlist]
  gf.icsu.times <- gf.icsu.times[!gf.icsu.times %in% gf.outlist]
  gf.seta.items <- gf.seta.items[!gf.seta.items %in% gf.outlist]
  gf.seta.times <- gf.seta.times[!gf.seta.times %in% gf.outlist]
  gf.setb.items <- gf.setb.items[!gf.setb.items %in% gf.outlist]
  gf.setb.times <- gf.setb.times[!gf.setb.times %in% gf.outlist]
  gf.anchor.items <- gf.anchor.items[!gf.anchor.items %in% gf.outlist]
  gf.anchor.times <- gf.anchor.times[!gf.anchor.times %in% gf.outlist]
}

# outflux after removal
if (exists("gf.outlist.outflux")) {
  gf.fluxplot.2 <- fluxplot(gf.imputed, main = NULL)
}

# prepare prediction model
gf.qp <- quickpred(gf.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# conduct mi - impute data using random forests
gf.mice <- mice(gf.imputed, pred = gf.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)

# create data sets for use in analyses
gf.imputed <- NULL # create empty variable
for (i in 1:5) gf.imputed[[i]] <- mice::complete(gf.mice, action = i, inc = FALSE) # put imputed data into list
for (i in 1:5) gf.imputed[[i]][gf.icsu.items] <- lapply(gf.imputed[[i]][gf.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric
for (i in 1:5) gf.imputed[[i]]$set <-  paste0("imp",i)
gf.imputed <- lapply(gf.imputed, function(x) {
  x <- mutate(x, gfscore = rowSums(select(x, starts_with("score_gf")), na.rm = TRUE))
  x <- mutate(x, gftime = rowSums(select(x, starts_with("timeTaken_gf")), na.rm = TRUE))
})
gf.backup <- gf.imputed

# recalculate times and scores for cleaned data
gf.cleaned <- mutate(gf.cleaned, gfscore = rowSums(select(gf.cleaned, starts_with("score_gf")), na.rm = TRUE), 
                     gftime = rowSums(select(gf.cleaned, starts_with("timeTaken_gf")), na.rm = TRUE))

# Combine the results together
gf.comparison <- rbind(gf.imputed[[1]][c("set","gfscore","gftime")],gf.cleaned[c("set","gfscore","gftime")])
gf.comparison <- rbind(gf.imputed[[2]][c("set","gfscore","gftime")],gf.comparison)
gf.comparison <- rbind(gf.imputed[[3]][c("set","gfscore","gftime")],gf.comparison)
gf.comparison <- rbind(gf.imputed[[4]][c("set","gfscore","gftime")],gf.comparison)
gf.comparison <- rbind(gf.imputed[[5]][c("set","gfscore","gftime")],gf.comparison)

# Boxplot of total score
ggplot(gf.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gfscore, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(0,length(gf.icsu.items),10)), limits = c(0,length(gf.icsu.items)), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# Boxplot of time
ggplot(gf.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gftime/60, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(0,max(gf.cleaned$gftime)/60,5)), limits = c(0,max(gf.cleaned$gftime)/60), name = "Total Time") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# visualise proportion correct
gf.cleaned.correct <- sapply(gf.cleaned[gf.icsu.items], function(x) {(1 - (table(x)[[1]]/length(na.omit(x)))) * 100 })
gf.imputed.correct <- lapply(gf.imputed, function(x) { sapply(x[gf.icsu.items], function(y) { (1 - (table(y)[[1]]/length(na.omit(y)))) * 100}) })
gf.combined.correct <- data.frame(cleaned = gf.cleaned.correct, imp1 = gf.imputed.correct[[1]], imp2 = gf.imputed.correct[[2]],
                                  imp3 = gf.imputed.correct[[3]], imp4 = gf.imputed.correct[[4]], imp5 = gf.imputed.correct[[5]])
row.names(gf.combined.correct) <- str_remove(row.names(gf.combined.correct), "score_gf")
gf.combined.correct <- rownames_to_column(gf.combined.correct, var = "item")
gf.combined.correct <- pivot_longer(data = gf.combined.correct, cols = -item, names_to = "set", values_to = "correct")
gf.combined.correct$item <- factor(gf.combined.correct$item, levels = c(sort(unique(as.numeric(gf.combined.correct$item)))))
gf.combined.correct$set <- factor(gf.combined.correct$set)
ggplot(gf.combined.correct, aes(x = item, y = correct, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gf.combined.correct$item), by = 2)) +
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Percentage Correct") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# visualise performance by age
gf.cleaned.age <- select(gf.cleaned, age, gfscore)
gf.imputed.age <- lapply(gf.imputed, function(x) { select(x, gfscore)})
gf.combined.age <- data.frame(age = gf.cleaned.age$age, cleaned = gf.cleaned.age$gfscore, imp1 = gf.imputed.age[[1]]$gfscore, 
                              imp2 = gf.imputed.age[[2]]$gfscore, imp3 = gf.imputed.age[[3]]$gfscore, imp4 = gf.imputed.age[[4]]$gfscore, 
                              imp5 = gf.imputed.age[[5]]$gfscore)
gf.combined.age <- pivot_longer(data = gf.combined.age, cols = -age, names_to = "set", values_to = "score")
gf.combined.age$set <- factor(gf.combined.age$set)
ggplot(gf.combined.age, aes(x = age, y = score, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(6,90,2), limits = c(6,90), name = "Age") +
  scale_y_continuous(breaks = seq(0,30,10), limits = c(0,30), name = "Total Score") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  geom_smooth(method = "lm")

#---                                                                           ---#
################################## Gf Reliability ##################################
#---                                                                          ---#

gf.reliability <- NULL
for (i in 1:5) gf.reliability[[i]] <- psych::alpha(gf.imputed[[i]][gf.icsu.items], warnings = FALSE)$total
gf.reliability <- sapply(gf.reliability, function(x) { x[,"raw_alpha"] })
mean(gf.reliability)

#---                                                                              ---#
################################## Gf Rasch Analysis ##################################
#---                                                                              ---#

# run rasch
gf.rasch <- list()
for (i in 1:5) gf.rasch[[i]] <- rasch.outcomes(chcability = "gf", raschobject = gf.imputed[[i]][gf.icsu.items], raschnumber = 1, imp = i)

# person fit
gf.rasch.person <- which(personfit(gf.rasch[[1]][[1]])$Zh < -2)
gf.rasch.person <- append(gf.rasch.person, which(personfit(gf.rasch[[1]][[1]])$Zh > 2))
gf.rasch.person
gf.imputed[[1]][gf.rasch.person,]$age

# wright map
WrightMap::wrightMap(fscores(gf.rasch[[1]][[1]]), main.title = NA)

# determine poor item fits
gf.poor.fit <- NULL
for (i in 1:5) gf.poor.fit[[i]] <- rasch.poor.fit(dataset = gf.imputed[[i]], raschobject = gf.rasch[[i]][[2]], imp = i)

# analyse the poor fit
poor.fit <- list()
for (i in 1:5) poor.fit[[i]] <- lapply(gf.poor.fit[[i]][[3]],function(x) {table(gf.imputed[[1]][[x]], gf.imputed[[1]]$age.group)})
for (i in 1:5) names(poor.fit[[i]]) <- gf.poor.fit[[i]][[3]]
poor.fit
gf.rasch.removal <- c(17, 52, 49) # this is a qualitative decision
gf.rasch.removal <- c(paste0("score_gf",gf.rasch.removal),paste0("timeTaken_gf",gf.rasch.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gf.imputed[[i]] <- select(gf.imputed[[i]], -all_of(gf.rasch.removal))
# remove.items.from.sets(chcability = "gf", items = gf.rasch.removal)
gf.icsa.items <- gf.icsa.items[!gf.icsa.items %in% gf.rasch.removal]
gf.icsa.times <- gf.icsa.times[!gf.icsa.times %in% gf.rasch.removal]
gf.icsu.items <- gf.icsu.items[!gf.icsu.items %in% gf.rasch.removal]
gf.icsu.times <- gf.icsu.times[!gf.icsu.times %in% gf.rasch.removal]
gf.seta.items <- gf.seta.items[!gf.seta.items %in% gf.rasch.removal]
gf.seta.times <- gf.seta.times[!gf.seta.times %in% gf.rasch.removal]
gf.setb.items <- gf.setb.items[!gf.setb.items %in% gf.rasch.removal]
gf.setb.times <- gf.setb.times[!gf.setb.times %in% gf.rasch.removal]
gf.anchor.items <- gf.anchor.items[!gf.anchor.items %in% gf.rasch.removal]
gf.anchor.times <- gf.anchor.times[!gf.anchor.times %in% gf.rasch.removal]

# run rasch 2
gf.rasch.2 <- list()
for (i in 1:5) gf.rasch.2[[i]] <- rasch.outcomes(chcability = "gf", raschobject = gf.imputed[[i]][gf.icsu.items], raschnumber = 2, imp = i)

# determine poor item fits 2
gf.poor.fit.2 <- NULL
for (i in 1:5) gf.poor.fit.2[[i]] <- rasch.poor.fit(dataset = gf.imputed[[i]], raschobject = gf.rasch.2[[i]][[2]], imp = i)

#---                                                                                ---#
################################## Gf Mokken Analysis ##################################
#---                                                                               ---#

# conduct first mokken
gf.mokken <- list()
for (i in 1:5) gf.mokken[[i]] <- mokken.outcomes(chcability = "gf", mokkenobject = gf.imputed[[i]][gf.icsu.items], mokkennumber = 1, imp = i)

# determine poorly ordered items and put them in a vector for removal
gf.mokken.prob <- NULL
for (i in 1:5) gf.mokken.prob[[i]] <- rownames(gf.mokken[[i]]$Hi)[which(gf.mokken[[i]]$Hi < .3)]
gf.mokken.prob
gf.mokken.prob <- sort(unique(unlist(gf.mokken.prob)))
gf.mokken.prob

gf.mokken.removal <- c(13, 15, 18, 29, 40, 41, 42, 48, 54, 58, 63) # this is a qualitative decision
gf.mokken.removal <- c(paste0("score_gf",gf.mokken.removal),paste0("timeTaken_gf",gf.mokken.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gf.imputed[[i]] <- select(gf.imputed[[i]], -all_of(gf.mokken.removal))
# remove.items.from.sets(chcability = "gf", items = gf.mokken.removal)
gf.icsa.items <- gf.icsa.items[!gf.icsa.items %in% gf.mokken.removal]
gf.icsa.times <- gf.icsa.times[!gf.icsa.times %in% gf.mokken.removal]
gf.icsu.items <- gf.icsu.items[!gf.icsu.items %in% gf.mokken.removal]
gf.icsu.times <- gf.icsu.times[!gf.icsu.times %in% gf.mokken.removal]
gf.seta.items <- gf.seta.items[!gf.seta.items %in% gf.mokken.removal]
gf.seta.times <- gf.seta.times[!gf.seta.times %in% gf.mokken.removal]
gf.setb.items <- gf.setb.items[!gf.setb.items %in% gf.mokken.removal]
gf.setb.times <- gf.setb.times[!gf.setb.times %in% gf.mokken.removal]
gf.anchor.items <- gf.anchor.items[!gf.anchor.items %in% gf.mokken.removal]
gf.anchor.times <- gf.anchor.times[!gf.anchor.times %in% gf.mokken.removal]

# conduct second mokken
gf.mokken.2 <- list()
for (i in 1:5) gf.mokken.2[[i]] <- mokken.outcomes(chcability = "gf", mokkenobject = gf.imputed[[i]][gf.icsu.items], mokkennumber = 2, imp = i)

#---                                                                                 ---#
################################## Gf Local Dependency ##################################
#---                                                                                 ---#

gf.local.dependency <- list()
for (i in 1:5) gf.local.dependency[[i]] <- residuals(gf.rasch[[i]][[1]], type = "Q3", digits = 2, suppress = +.2)

# determine items with high local dependency
gf.local.dependency.removal <- ""

# remove the items that do not fit well to the rasch model
for (i in 1:5) gf.imputed[[i]] <- select(gf.imputed[[i]], -all_of(gf.local.dependency.removal))
# remove.items.from.sets(chcability = "gf", items = gf.local.dependency.removal)
gf.icsa.items <- gf.icsa.items[!gf.icsa.items %in% gf.local.dependency.removal]
gf.icsa.times <- gf.icsa.times[!gf.icsa.times %in% gf.local.dependency.removal]
gf.icsu.items <- gf.icsu.items[!gf.icsu.items %in% gf.local.dependency.removal]
gf.icsu.times <- gf.icsu.times[!gf.icsu.times %in% gf.local.dependency.removal]
gf.seta.items <- gf.seta.items[!gf.seta.items %in% gf.local.dependency.removal]
gf.seta.times <- gf.seta.times[!gf.seta.times %in% gf.local.dependency.removal]
gf.setb.items <- gf.setb.items[!gf.setb.items %in% gf.local.dependency.removal]
gf.setb.times <- gf.setb.times[!gf.setb.times %in% gf.local.dependency.removal]
gf.anchor.items <- gf.anchor.items[!gf.anchor.items %in% gf.local.dependency.removal]
gf.anchor.times <- gf.anchor.times[!gf.anchor.times %in% gf.local.dependency.removal]

#---                                                                                              ---#
################################## Gf Differential Item Functioning ##################################
#---                                                                                              ---#

gf.dif.gender <- list()
for (i in 1:5) gf.dif.gender[[i]] <- difLogistic(Data = gf.imputed[[i]][gf.icsu.items], group = gf.imputed[[i]]$gender, focal.name = "m", alpha = .01)
for (i in 1:5) plot(gf.dif.gender[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gf.dif.nationality <- list()
for (i in 1:5) gf.dif.nationality[[i]] <- difLogistic(Data = gf.imputed[[i]][gf.icsu.items], group = gf.imputed[[i]]$nationality, focal.name = "au", alpha = .01)
for (i in 1:5) plot(gf.dif.nationality[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gf.dif.device <- list()
for (i in 1:5) gf.dif.device[[i]] <- difLogistic(Data = gf.imputed[[i]][gf.icsu.items], group = gf.imputed[[i]]$nationality, focal.name = "au", alpha = .01)
for (i in 1:5) plot(gf.dif.device[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

# determine dif items
gf.dif.items <- NULL
for (i in 1:5) gf.dif.items[[i]] <- { 
  gender <- gf.dif.gender[[i]]$DIFitems
  nationality <- gf.dif.nationality[[i]]$DIFitems
  device <- gf.dif.device[[i]]$DIFitems
  list(gender,nationality,device)
}

gf.dif.items <- sort(as.numeric(unique(unlist(gf.dif.items))), na.last = NA)
gf.dif.items
gf.dif.items <- gf.icsu.items[gf.dif.items]
gf.dif.items

for (x in 1:5) for (i in gf.dif.items) { plot(gf.dif.gender[[x]], plot = "itemCurve", item = i) }
for (x in 1:5) for (i in gf.dif.items) { plot(gf.dif.nationality[[x]], plot = "itemCurve", item = i) }
for (x in 1:5) for (i in gf.dif.items) { plot(gf.dif.device[[x]], plot = "itemCurve", item = i) }

#---                                                                                             ---#
################################## Gf Confirmatory Factor Analysis ##################################
#---                                                                                            ---#

# run the CFA on the multivariate imputation data sets
for (i in 1:5) gf.imputed[[i]][gf.icsu.items] <- lapply(gf.imputed[[i]][gf.icsu.items], function(x) { x <- ordered(x, levels = c(0,1)) }) # Change all items to ordered factor
gf.model <- paste("gf =~ ", paste(gf.icsu.items, collapse = "+"), sep = "") # establish model
gf.cfa.mi <- runMI(gf.model, data = gf.imputed, fun = "cfa", meanstructure = TRUE) # run cfa --> takes a very long time
gf.cfa.mi.summary <- summary(gf.cfa.mi)
gf.cfa.mi.fit <- fitMeasures(gf.cfa.mi, test = "D2", pool.robust = TRUE)

# determine what cfas were successful
gf.cfa.success <- NULL
for (i in 1:5) gf.cfa.success[i] <- !is.null(gf.cfa.mi@miList[[i]])

# based on the cfas that were successful, run them 1 by 1 for diagnostic and reporting purposes, fit the cfa
gf.cfa <- list()
for (i in 1:5) gf.cfa[[i]] <- if(gf.cfa.success[i] == TRUE) {cfa(gf.model, data = gf.imputed[[i]])}

# store the cfa summary for item selection later
gf.cfa.summary <- NULL
for (i in 1:5) gf.cfa.summary[[i]] <- if(gf.cfa.success[i] == TRUE) { summary(object = gf.cfa[[i]], standardized = TRUE) }

# print and save the parameters using predefined function
gf.cfa.parameters <- NULL
for (i in 1:5) gf.cfa.parameters[[i]] <- if(gf.cfa.success[i] == TRUE) {
  cfa.parameters(chcability = "gf", cfaobject = gf.cfa[[i]], cfanumber = 1, imp = i)}

# print and store the overall fit
gf.cfa.fit <- NULL
for (i in 1:5) gf.cfa.fit[[i]] <- if(gf.cfa.success[i] == TRUE) {
  fitMeasures(gf.cfa[[i]], fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
}
gf.cfa.fit

# determine low loading items and put them in a vector for removal
gf.cfa.poor.fit <- list()
for (i in 1:5) gf.cfa.poor.fit[[i]] <- if(gf.cfa.success[i] == TRUE) {
  select(filter(gf.cfa.summary[[i]]$PE, lhs == "gf", std.lv < .3, rhs != ""), rhs)$rhs } else {NA}
gf.cfa.poor.fit
gf.cfa.poor.fit <- sort(unique(unlist(gf.cfa.poor.fit)))
gf.cfa.poor.fit

#---                                                                                     ---#
################################## Gf Final Rasch Analysis ##################################
#---                                                                                    ---#

for (i in 1:5) gf.imputed[[i]][gf.icsu.items] <- lapply(gf.imputed[[i]][gf.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric

# run rasch 3
gf.rasch.3 <- list()
for (i in 1:5) gf.rasch.3[[i]] <- rasch.outcomes(chcability = "gf", raschobject = gf.imputed[[i]][gf.icsu.items], raschnumber = 3, imp = i)

# parameters
gf.classic.parameters <- list()
for (i in 1:5) gf.classic.parameters[[i]] <- coef(gf.rasch.3[[i]][[1]], IRTpars = TRUE, simplify = TRUE)$items
gf.classic.parameters

gf.parameters <- list()
for (i in 1:5) gf.parameters[[i]] <- coef(gf.rasch.3[[i]][[1]], simplify = TRUE)$items
gf.parameters

# pool the parameters
gf.rasch.average <- list()
for (i in 1:5) gf.rasch.average[[i]] <- as.data.frame(coef(gf.rasch.3[[i]][[1]], printSE = TRUE, as.data.frame = TRUE))
for (i in 1:5) gf.rasch.average[[i]] <- rownames_to_column(gf.rasch.average[[i]], var = "item")
for (i in 1:5) gf.rasch.average[[i]] <- filter(gf.rasch.average[[i]], par > -10, par < 10, par != 1)
gf.item.list <- lapply(gf.rasch.average, function(x) { x$item})
gf.par.list <- lapply(gf.rasch.average, function(x) { x$par })
gf.separ.list <- lapply(gf.rasch.average, function(x) { x$SE })
gf.rasch.average <- averageMI(gf.par.list,gf.separ.list,as.data.frame = TRUE)
gf.rasch.average$item <- str_remove(gf.item.list[[1]], "score_")
gf.rasch.average$item <- str_remove(gf.rasch.average$item, ".d")
gf.rasch.average$a1 <- 1
gf.rasch.average$g <- 0
gf.rasch.average$u <- 1
colnames(gf.rasch.average)[1] <- "d"
gf.rasch.average <- select(gf.rasch.average, item, a1, d, g, u)
gf.rasch.average <- filter(gf.rasch.average, grepl(pattern = "gf", item))

# get f scores
gf.fscores <- gf.imputed[[1]]
gf.fscores <- mutate(gf.fscores, gfscore = rowSums(select(gf.fscores, starts_with("score_gf"))))
gf.fscores$fscores <- fscores(object = gf.rasch.3[[1]][[1]])
gf.fscores <- select(gf.fscores, -starts_with("score_gf"), -starts_with("timeTaken"), -set)
gf.wisc <- select(filter(data, gfscore > 0, wisc.fsiq.comp > 0), id, starts_with("wisc"))
gf.fscores <- left_join(gf.wisc,gf.fscores)

# save parameters, remaining items and f scores
write.table(gf.rasch.average, file = "gf_rasch_parameters.csv", sep = ",", quote = FALSE, row.names = F)
write(gf.icsu.items, file = "gfitems.txt")
write.table(gf.fscores, file = "gf_fscores.csv", sep = ",", quote = FALSE, row.names = F)

#---                                                                               ---#
################################## Gv Data Imputation ##################################
#---                                                                               ---#

gv.itos.removed <- ""

# initial
gv.cleaned <- gv
gv.cleaned$gvmissing.initial <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.initial) # obtain missing data by participant
gv.missing1 <- select(gv.cleaned, id, gvmissing.initial)

# remove people who dont have a full set of items for anchor items
gv.cleaned <- filter(gv.cleaned, !is.na(gv.anchor)) # set new data sets
gv.cleaned$gvmissing.fullanchor <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.fullanchor)
gv.missing2 <- select(gv.cleaned, id, gvmissing.fullanchor)

# remove items that were removed as part of itos
gv.cleaned <- select(gv.cleaned, -all_of(gv.itos.removed)) # remove unused items
gv.cleaned$gvmissing.removeditems <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.removeditems)
gv.missing3 <- select(gv.cleaned, id, gvmissing.removeditems)

# remove items that have low rates of response
gv.low.response <- apply(select(gv.cleaned, starts_with("score_gv"), starts_with("timeTaken")), 2, percentage.missing) > 75
gv.low.response <- names(gv.low.response[gv.low.response == TRUE])
gv.low.response

gv.icsa.items <- gv.icsa.items[!gv.icsa.items %in% gv.low.response]
gv.icsa.times <- gv.icsa.times[!gv.icsa.times %in% gv.low.response]
gv.icsu.items <- gv.icsu.items[!gv.icsu.items %in% gv.low.response]
gv.icsu.times <- gv.icsu.times[!gv.icsu.times %in% gv.low.response]
gv.seta.items <- gv.seta.items[!gv.seta.items %in% gv.low.response]
gv.seta.times <- gv.seta.times[!gv.seta.times %in% gv.low.response]
gv.setb.items <- gv.setb.items[!gv.setb.items %in% gv.low.response]
gv.setb.times <- gv.setb.times[!gv.setb.times %in% gv.low.response]
gv.anchor.items <- gv.anchor.items[!gv.anchor.items %in% gv.low.response]
gv.anchor.times <- gv.anchor.times[!gv.anchor.times %in% gv.low.response]

gv.cleaned <- select(gv.cleaned, -all_of(gv.low.response))
gv.cleaned$gvmissing.lowresponse <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.lowresponse)
gv.missing4 <- select(gv.cleaned, id, gvmissing.lowresponse)

# remove participants with high rates of missingness
gv.cleaned <- filter(gv.cleaned, gvmissing.lowresponse < 33)
gv.cleaned$gvmissing.removedhighmissingness <- apply(select(gv.cleaned, c(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.removedhighmissingness)
gv.missing5 <- select(gv.cleaned, id, gvmissing.removedhighmissingness)

# remove participants with missing age
gv.cleaned <- filter(gv.cleaned, !is.na(age))
gv.cleaned$gvmissing.missingage <- apply(select(gv.cleaned, c(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.missingage)
gv.missing6 <- select(gv.cleaned, id, gvmissing.missingage)

# recalculate times and scores for cleaned data
gv.cleaned <- mutate(gv.cleaned, gvscore = rowSums(select(gv.cleaned, starts_with("score_gv")), na.rm = TRUE), 
                     gvtime = rowSums(select(gv.cleaned, starts_with("timeTaken_gv")), na.rm = TRUE))

# remove univariate time outliers by 1st and 99th percentile
gv.cleaned <- filter(gv.cleaned, gvtime > quantile(gv.cleaned$gvtime, 0.05) & gvtime < quantile(gv.cleaned$gvtime, 0.95))
gv.cleaned$gvmissing.univariatetime <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.univariatetime) # obtain missing data by participant
gv.missing7 <- select(gv.cleaned, id, gvmissing.univariatetime)

# remove univariate score outliers by 1st and 99th percentile
gv.cleaned <- filter(gv.cleaned, gvscore > quantile(gv.cleaned$gvscore, 0.01) & gvscore < quantile(gv.cleaned$gvscore, 0.99))
gv.cleaned$gvmissing.univariatescore <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.univariatescore) # obtain missing data by participant
gv.missing8 <- select(gv.cleaned, id, gvmissing.univariatescore)

# remove multivariate outliers
gv.lm.model <- paste("gvscore ~ age + gvtime + gender + phase")
gv.lm <- lm(gv.lm.model, data = gv.cleaned)
gv.cd <- cooks.distance(gv.lm)
plot(gv.cd, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean(gv.cd, na.rm = T), col = "red")
text(x=1:length(gv.cd)+1, y=gv.cd, labels=ifelse(gv.cd>4*mean(gv.cd, na.rm = TRUE ), names(gv.cd),""), col = "red")
gv.influential <- as.numeric(names(gv.cd)[(gv.cd > 4*mean(gv.cd, na.rm = TRUE))])
gv.influential <- gv.influential[!is.na(gv.influential)]
gv.cleaned <- gv.cleaned[-gv.influential,]
gv.cleaned$gvmissing.multivariate <- apply(select(gv.cleaned, all_of(gv.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gv.cleaned$gvmissing.multivariate) # obtain missing data by participant
gv.missing9 <- select(gv.cleaned, id, gvmissing.multivariate)

# create set column
gv.cleaned$set <- "cleaned"

# create df of item data missing
gv.missing.prior <- round(as.data.frame(apply(select(gv, starts_with("score_gv")), 2, percentage.missing))) # select original gv data
names(gv.missing.prior)[1] <- "missing"
gv.missing.prior$stage <- "Original"
gv.missing.prior <- rownames_to_column(gv.missing.prior, "item")
gv.missing.after <- round(as.data.frame(apply(select(gv.cleaned, starts_with("score_gv")), 2, percentage.missing))) # select cleaned data
names(gv.missing.after)[1] <- "missing"
gv.missing.after$stage <- "After Cleaning"
gv.missing.after <- rownames_to_column(gv.missing.after, "item")
gv.missing.items <- rbind(gv.missing.prior,gv.missing.after)
gv.missing.items$item <- as.numeric(str_extract(gv.missing.items$item, "[0-9]+"))
rm(gv.missing.after,gv.missing.prior)

# Visualise missing data by item
ggplot(gv.missing.items, aes(x = item, y = missing, colour = stage)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom") +
  labs(colour = "Stage") +
  geom_vline(xintercept = 34) +
  geom_vline(xintercept = 56) +
  annotate("text", x = 16, y = 50, label = "Set A Items") +
  annotate("text", x = 46, y = 50, label = "Anchor Items") +
  annotate("text", x = 85, y = 50, label = "Set B Items")

# combine missing data columns
gv.missing.participants <- full_join(gv.missing1, gv.missing2)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing3)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing4)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing5)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing6)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing7)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing8)
gv.missing.participants <- full_join(gv.missing.participants, gv.missing9)
rm(gv.missing1,gv.missing2,gv.missing3,gv.missing4,gv.missing5,gv.missing6,gv.missing7,gv.missing8,gv.missing9)

# visualise the changes in missing data across steps
ggplot(reshape2::melt(gv.missing.participants[,c("gvmissing.initial","gvmissing.fullanchor","gvmissing.removeditems","gvmissing.lowresponse",
                                                 "gvmissing.removedhighmissingness", "gvmissing.missingage","gvmissing.univariatetime",
                                                 "gvmissing.univariatescore","gvmissing.multivariate")]), 
       aes(x = variable, y = value, colour = variable)) +
  geom_boxplot(width = .25) +
  ylab("Percentage of Items Missing") +
  scale_x_discrete(name = "Stage of Data Cleaning", labels = ggplot_addline(c("Full Data","Full Anchor","Items Removed","Low Response Items",
                                                                              "High Missingness Removed","Missing Age Removed",
                                                                              "Outliers by Time","Outliers by Score",
                                                                              "Multivariate Outliers"))) +
  theme(legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text")

# Recalculate the set missingness
gv.cleaned <- mutate(gv.cleaned, 
                     gv.set.a = ceiling(rowSums(is.na(gv.cleaned[gv.seta.items]))/length(gv.seta.items)),
                     gv.set.b = ceiling(rowSums(is.na(gv.cleaned[gv.setb.items]))/length(gv.setb.items)),
                     gv.anchor = ceiling(rowSums(is.na(gv.cleaned[gv.anchor.items]))/length(gv.anchor.items)))

# Recode missingness # NA = missing data in set, 1 = full data
gv.cleaned$gv.set.a <- gv.cleaned$gv.set.a %>% na_if(1) %>% recode(`0` = 1)
gv.cleaned$gv.set.b <- gv.cleaned$gv.set.b %>% na_if(1) %>% recode(`0` = 1)
gv.cleaned$gv.anchor <- gv.cleaned$gv.anchor %>% na_if(1) %>% recode(`0` = 1)

# Visualise Missing Data Patterns
aggr(gv.cleaned[c("gv.set.a","gv.set.b","gv.anchor")], combined = TRUE, labels = c("Set A", "Set B", "Anchor")) 

# Set Up Data Set
gv.imputed <- select(gv.cleaned, id, phase, age, age.group, gender, nationality, device, all_of(gv.icsu.items), all_of(gv.icsu.times)) # create data set
gv.imputed[gv.icsu.items] <- lapply(gv.imputed[gv.icsu.items], factor)

# prepare prediction model
gv.qp <- quickpred(gv.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# flag outflux problem items
gv.flux <- flux(gv.imputed)
gv.fluxplot <- fluxplot(gv.imputed, main = NULL)
gv.outlist.outflux <- row.names(gv.flux)[gv.flux$outflux < 0.3]

# remove problem items for mi
if (exists("gv.outlist.logged") & exists("gv.outlist.outflux")) {
  gv.outlist <- unique(c(gv.outlist.logged, gv.outlist.outflux))
} else if (exists("gv.outlist.logged")) {
  gv.outlist <- gv.outlist.logged
} else if (exists("gv.outlist.outflux")) {
  gv.outlist <- gv.outlist.outflux
} else {
  print("gv outlist not created")
}

if (exists("gv.outlist")) {
  gv.outlist <- as.numeric(str_extract(gv.outlist, "[0-9]+"))
  gv.outlist <- c(paste0("score_gv",gv.outlist),paste0("timeTaken_gv",gv.outlist))
  gv.outlist <- unique(gv.outlist)
  gv.imputed <- select(gv.imputed, -all_of(gv.outlist))
  gv.cleaned <- select(gv.cleaned, -all_of(gv.outlist))
  gv.icsa.items <- gv.icsa.items[!gv.icsa.items %in% gv.outlist]
  gv.icsa.times <- gv.icsa.times[!gv.icsa.times %in% gv.outlist]
  gv.icsu.items <- gv.icsu.items[!gv.icsu.items %in% gv.outlist]
  gv.icsu.times <- gv.icsu.times[!gv.icsu.times %in% gv.outlist]
  gv.seta.items <- gv.seta.items[!gv.seta.items %in% gv.outlist]
  gv.seta.times <- gv.seta.times[!gv.seta.times %in% gv.outlist]
  gv.setb.items <- gv.setb.items[!gv.setb.items %in% gv.outlist]
  gv.setb.times <- gv.setb.times[!gv.setb.times %in% gv.outlist]
  gv.anchor.items <- gv.anchor.items[!gv.anchor.items %in% gv.outlist]
  gv.anchor.times <- gv.anchor.times[!gv.anchor.times %in% gv.outlist]
}

# outflux after removal
if (exists("gv.outlist.outflux")) { gv.fluxplot.2 <- fluxplot(gv.imputed, main = NULL) }

# prepare prediction model
gv.qp <- quickpred(gv.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# conduct mi - impute data using random forests
gv.mice <- mice(gv.imputed, pred = gv.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)

# create data sets for use in analyses
gv.imputed <- NULL # create empty variable
for (i in 1:5) gv.imputed[[i]] <- mice::complete(gv.mice, action = i, inc = FALSE) # put imputed data into list
for (i in 1:5) gv.imputed[[i]][gv.icsu.items] <- lapply(gv.imputed[[i]][gv.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric
for (i in 1:5) gv.imputed[[i]]$set <-  paste0("imp",i)
gv.imputed <- lapply(gv.imputed, function(x) {
  x <- mutate(x, gvscore = rowSums(select(x, starts_with("score_gv")), na.rm = TRUE))
  x <- mutate(x, gvtime = rowSums(select(x, starts_with("timeTaken_gv")), na.rm = TRUE))
})
gv.backup <- gv.imputed

# recalculate times and scores for cleaned data
gv.cleaned <- mutate(gv.cleaned, gvscore = rowSums(select(gv.cleaned, starts_with("score_gv")), na.rm = TRUE), 
                     gvtime = rowSums(select(gv.cleaned, starts_with("timeTaken_gv")), na.rm = TRUE))

# Combine the results together
gv.comparison <- rbind(gv.imputed[[1]][c("set","gvscore","gvtime")],gv.cleaned[c("set","gvscore","gvtime")])
gv.comparison <- rbind(gv.imputed[[2]][c("set","gvscore","gvtime")],gv.comparison)
gv.comparison <- rbind(gv.imputed[[3]][c("set","gvscore","gvtime")],gv.comparison)
gv.comparison <- rbind(gv.imputed[[4]][c("set","gvscore","gvtime")],gv.comparison)
gv.comparison <- rbind(gv.imputed[[5]][c("set","gvscore","gvtime")],gv.comparison)

# Boxplot of total score
ggplot(gv.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gvscore, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(20,60,10)), limits = c(20,60), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# Boxplot of time
ggplot(gv.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gvtime/60, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(0,20,5)), limits = c(0,20), name = "Total Time") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# visualise proportion correct
gv.cleaned.correct <- sapply(gv.cleaned[gv.icsu.items], function(x) {(1 - (table(x)[[1]]/length(na.omit(x)))) * 100 })
gv.imputed.correct <- lapply(gv.imputed, function(x) {sapply(x[gv.icsu.items], function(y) {(1 - (table(y)[[1]]/length(na.omit(y)))) * 100 }) })
gv.combined.correct <- data.frame(cleaned = gv.cleaned.correct, imp1 = gv.imputed.correct[[1]], imp2 = gv.imputed.correct[[2]],
                                  imp3 = gv.imputed.correct[[3]], imp4 = gv.imputed.correct[[4]], imp5 = gv.imputed.correct[[5]])
row.names(gv.combined.correct) <- str_remove(row.names(gv.combined.correct), "score_gv")
gv.combined.correct <- rownames_to_column(gv.combined.correct, var = "item")
gv.combined.correct <- pivot_longer(data = gv.combined.correct, cols = -item, names_to = "set", values_to = "correct")
gv.combined.correct$item <- factor(gv.combined.correct$item, levels = c(sort(unique(as.numeric(gv.combined.correct$item)))))
gv.combined.correct$set <- factor(gv.combined.correct$set)
ggplot(gv.combined.correct, aes(x = item, y = correct, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gv.combined.correct$item), by = 2)) +
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Percentage Correct") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# visualise performance by age
gv.cleaned.age <- select(gv.cleaned, age, gvscore)
gv.imputed.age <- lapply(gv.imputed, function(x) { select(x, gvscore)})
gv.combined.age <- data.frame(age = gv.cleaned.age$age, cleaned = gv.cleaned.age$gvscore, imp1 = gv.imputed.age[[1]]$gvscore, 
                              imp2 = gv.imputed.age[[2]]$gvscore, imp3 = gv.imputed.age[[3]]$gvscore, imp4 = gv.imputed.age[[4]]$gvscore, 
                              imp5 = gv.imputed.age[[5]]$gvscore)
gv.combined.age <- pivot_longer(data = gv.combined.age, cols = -age, names_to = "set", values_to = "score")
gv.combined.age$set <- factor(gv.combined.age$set)
ggplot(gv.combined.age, aes(x = age, y = score, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(6,90,2), limits = c(6,90), name = "Age") +
  scale_y_continuous(breaks = seq(20,50,10), limits = c(20,50), name = "Total Score") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  geom_smooth(method = "lm")

#---                                                                           ---#
################################## Gv Reliability ##################################
#---                                                                          ---#

gv.reliability <- NULL
for (i in 1:5) gv.reliability[[i]] <- psych::alpha(gv.imputed[[i]][gv.icsu.items], warnings = FALSE)$total
gv.reliability <- sapply(gv.reliability, function(x) { x[,"raw_alpha"] })
mean(gv.reliability)

#---                                                                              ---#
################################## Gv Rasch Analysis ##################################
#---                                                                              ---#

# run rasch
gv.rasch <- list()
for (i in 1:5) gv.rasch[[i]] <- rasch.outcomes(chcability = "gv", raschobject = gv.imputed[[i]][gv.icsu.items], raschnumber = 1, imp = i)

# person fit
gv.rasch.person <- which(personfit(gv.rasch[[1]][[1]])$Zh < -2)
gv.rasch.person <- append(gv.rasch.person, which(personfit(gv.rasch[[1]][[1]])$Zh > 2))
gv.rasch.person
gv.imputed[[1]][gv.rasch.person,]$age

# wright map
WrightMap::wrightMap(fscores(gv.rasch[[1]][[1]]), main.title = NA)

# determine poor item fits
gv.poor.fit <- NULL
for (i in 1:5) gv.poor.fit[[i]] <- rasch.poor.fit(dataset = gv.imputed[[i]], raschobject = gv.rasch[[i]][[2]], imp = i)
gv.rasch.removal <- c(23, 44, 52, 10, 13, 17, 35, 43, 48, 50) # this is a qualitative decision
gv.rasch.removal <- c(paste0("score_gv",gv.rasch.removal),paste0("timeTaken_gv",gv.rasch.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gv.imputed[[i]] <- select(gv.imputed[[i]], -all_of(gv.rasch.removal))
# remove.items.from.sets(chcability = "gv", items = gv.rasch.removal)
gv.icsa.items <- gv.icsa.items[!gv.icsa.items %in% gv.rasch.removal]
gv.icsa.times <- gv.icsa.times[!gv.icsa.times %in% gv.rasch.removal]
gv.icsu.items <- gv.icsu.items[!gv.icsu.items %in% gv.rasch.removal]
gv.icsu.times <- gv.icsu.times[!gv.icsu.times %in% gv.rasch.removal]
gv.seta.items <- gv.seta.items[!gv.seta.items %in% gv.rasch.removal]
gv.seta.times <- gv.seta.times[!gv.seta.times %in% gv.rasch.removal]
gv.setb.items <- gv.setb.items[!gv.setb.items %in% gv.rasch.removal]
gv.setb.times <- gv.setb.times[!gv.setb.times %in% gv.rasch.removal]
gv.anchor.items <- gv.anchor.items[!gv.anchor.items %in% gv.rasch.removal]
gv.anchor.times <- gv.anchor.times[!gv.anchor.times %in% gv.rasch.removal]

# run rasch 2
gv.rasch.2 <- list()
for (i in 1:5) gv.rasch.2[[i]] <- rasch.outcomes(chcability = "gv", raschobject = gv.imputed[[i]][gv.icsu.items], raschnumber = 2, imp = i)

#---                                                                                ---#
################################## Gv Mokken Analysis ##################################
#---                                                                               ---#

# conduct first mokken
gv.mokken <- list()
for (i in 1:5) gv.mokken[[i]] <- mokken.outcomes(chcability = "gv", mokkenobject = gv.imputed[[i]][gv.icsu.items], mokkennumber = 1, imp = i)

# determine poorly ordered items and put them in a vector for removal
gv.mokken.prob <- NULL
for (i in 1:5) gv.mokken.prob[[i]] <- rownames(gv.mokken[[i]]$Hi)[which(gv.mokken[[i]]$Hi < .3)]
gv.mokken.prob
gv.mokken.prob <- sort(unique(unlist(gv.mokken.prob)))
gv.mokken.prob

gv.mokken.removal <- c(1, 11, 16, 30, 32, 37, 38, 4, 42, 7, 8) # this is a qualitative decision
gv.mokken.removal <- c(paste0("score_gv",gv.mokken.removal),paste0("timeTaken_gv",gv.mokken.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gv.imputed[[i]] <- select(gv.imputed[[i]], -all_of(gv.mokken.removal))
# remove.items.from.sets(chcability = "gv", items = gv.mokken.removal)
gv.icsa.items <- gv.icsa.items[!gv.icsa.items %in% gv.mokken.removal]
gv.icsa.times <- gv.icsa.times[!gv.icsa.times %in% gv.mokken.removal]
gv.icsu.items <- gv.icsu.items[!gv.icsu.items %in% gv.mokken.removal]
gv.icsu.times <- gv.icsu.times[!gv.icsu.times %in% gv.mokken.removal]
gv.seta.items <- gv.seta.items[!gv.seta.items %in% gv.mokken.removal]
gv.seta.times <- gv.seta.times[!gv.seta.times %in% gv.mokken.removal]
gv.setb.items <- gv.setb.items[!gv.setb.items %in% gv.mokken.removal]
gv.setb.times <- gv.setb.times[!gv.setb.times %in% gv.mokken.removal]
gv.anchor.items <- gv.anchor.items[!gv.anchor.items %in% gv.mokken.removal]
gv.anchor.times <- gv.anchor.times[!gv.anchor.times %in% gv.mokken.removal]

# conduct second mokken
gv.mokken.2 <- list()
for (i in 1:5) gv.mokken.2[[i]] <- mokken.outcomes(chcability = "gv", mokkenobject = gv.imputed[[i]][gv.icsu.items], mokkennumber = 2, imp = i)

#---                                                                                 ---#
################################## Gv Local Dependency ##################################
#---                                                                                 ---#

gv.local.dependency <- list()
for (i in 1:5) gv.local.dependency[[i]] <- residuals(gv.rasch[[i]][[1]], type = "Q3", digits = 2, suppress = +.2)

# determine items with high local dependency
gv.local.dependency.removal <- ""

# remove the items that do not fit well to the rasch model
for (i in 1:5) gv.imputed[[i]] <- select(gv.imputed[[i]], -all_of(gv.local.dependency.removal))
# remove.items.from.sets(chcability = "gv", items = gv.local.dependency.removal)
gv.icsa.items <- gv.icsa.items[!gv.icsa.items %in% gv.local.dependency.removal]
gv.icsa.times <- gv.icsa.times[!gv.icsa.times %in% gv.local.dependency.removal]
gv.icsu.items <- gv.icsu.items[!gv.icsu.items %in% gv.local.dependency.removal]
gv.icsu.times <- gv.icsu.times[!gv.icsu.times %in% gv.local.dependency.removal]
gv.seta.items <- gv.seta.items[!gv.seta.items %in% gv.local.dependency.removal]
gv.seta.times <- gv.seta.times[!gv.seta.times %in% gv.local.dependency.removal]
gv.setb.items <- gv.setb.items[!gv.setb.items %in% gv.local.dependency.removal]
gv.setb.times <- gv.setb.times[!gv.setb.times %in% gv.local.dependency.removal]
gv.anchor.items <- gv.anchor.items[!gv.anchor.items %in% gv.local.dependency.removal]
gv.anchor.times <- gv.anchor.times[!gv.anchor.times %in% gv.local.dependency.removal]

#---                                                                                              ---#
################################## Gv Differential Item Functioning ##################################
#---                                                                                              ---#

gv.dif.gender <- list()
for (i in 1:5) gv.dif.gender[[i]] <- difLogistic(Data = gv.imputed[[i]][gv.icsu.items], group = gv.imputed[[i]]$gender, focal.name = "m", alpha = .01)
for (i in 1:5) plot(gv.dif.gender[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gv.dif.nationality <- list()
for (i in 1:5) gv.dif.nationality[[i]] <- difLogistic(Data = gv.imputed[[i]][gv.icsu.items], group = gv.imputed[[i]]$nationality, focal.name = "au", alpha = .01)
for (i in 1:5) plot(gv.dif.nationality[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gv.dif.device <- list()
for (i in 1:5) gv.dif.device[[i]] <- difLogistic(Data = gv.imputed[[i]][gv.icsu.items], group = gv.imputed[[i]]$nationality, focal.name = "ipad", alpha = .01)
for (i in 1:5) plot(gv.dif.device[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

# determine dif items
gv.dif.items <- NULL
for (i in 1:5) gv.dif.items[[i]] <- { 
  gender <- gv.dif.gender[[i]]$DIFitems
  nationality <- gv.dif.nationality[[i]]$DIFitems
  device <- gv.dif.device[[i]]$DIFitems
  list(gender,nationality,device)
}

gv.dif.items <- sort(as.numeric(unique(unlist(gv.dif.items))), na.last = NA)
gv.dif.items
gv.dif.items <- gv.icsu.items[gv.dif.items]
gv.dif.items

for (y in 1:5) for (i in gv.dif.items) { plot(gv.dif.gender[[y]], plot = "itemCurve", item = i) }
for (y in 1:5) for (i in gv.dif.items) { plot(gv.dif.nationality[[y]], plot = "itemCurve", item = i) }
for (y in 1:5) for (i in gv.dif.items) { plot(gv.dif.device[[y]], plot = "itemCurve", item = i) }

#---                                                                                             ---#
################################## Gv Confirmatory Factor Analysis ##################################
#---                                                                                            ---#

# run the CFA on the multivariate imputation data sets
for (i in 1:5) gv.imputed[[i]][gv.icsu.items] <- lapply(gv.imputed[[i]][gv.icsu.items], function(x) { x <- ordered(x, levels = c(0,1)) }) # Change all items to ordered factor
gv.model <- paste("gv =~ ", paste(gv.icsu.items, collapse = "+"), sep = "") # establish model
gv.cfa.mi <- runMI(gv.model, data = gv.imputed, fun = "cfa", meanstructure = TRUE) # run cfa --> takes a very long time
gv.cfa.mi.summary <- summary(gv.cfa.mi)
gv.cfa.mi.fit <- fitMeasures(gv.cfa.mi, fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))

# determine what cfas were successful
gv.cfa.success <- NULL
for (i in 1:5) gv.cfa.success[i] <- !is.null(gv.cfa.mi@miList[[i]])

# based on the cfas that were successful, run them 1 by 1 for diagnostic and reporting purposes, fit the cfa
gv.cfa <- list()
for (i in 1:5) gv.cfa[[i]] <- if(gv.cfa.success[i] == TRUE) {cfa(gv.model, data = gv.imputed[[i]])}

# store the cfa summary for item selection later
gv.cfa.summary <- NULL
for (i in 1:5) gv.cfa.summary[[i]] <- if(gv.cfa.success[i] == TRUE) { summary(object = gv.cfa[[i]], standardized = TRUE) }

# print and save the parameters using predefined function
gv.cfa.parameters <- NULL
for (i in 1:5) gv.cfa.parameters[[i]] <- if(gv.cfa.success[i] == TRUE) {
  cfa.parameters(chcability = "gv", cfaobject = gv.cfa[[i]], cfanumber = 1, imp = i)}

# print and store the overall fit
gv.cfa.fit <- NULL
for (i in 1:5) gv.cfa.fit[[i]] <- if(gv.cfa.success[i] == TRUE) {
  fitMeasures(gv.cfa[[i]], fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
}
gv.cfa.fit

# determine low loading items and put them in a vector for removal
gv.cfa.poor.fit <- list()
for (i in 1:5) gv.cfa.poor.fit[[i]] <- if(gv.cfa.success[i] == TRUE) {
  select(filter(gv.cfa.summary[[i]]$PE, lhs == "gv", std.lv < .3, rhs != ""), rhs)$rhs } else {NA}
gv.cfa.poor.fit
gv.cfa.poor.fit <- sort(unique(unlist(gv.cfa.poor.fit)))
gv.cfa.poor.fit

# determine items with high local dependency
gv.cfa.removal <- "score_gv25"

# remove the items that do not fit well to the rasch model
for (i in 1:5) gv.imputed[[i]] <- select(gv.imputed[[i]], -all_of(gv.cfa.removal))
# remove.items.from.sets(chcability = "gv", items = gv.cfa.removal)
gv.icsa.items <- gv.icsa.items[!gv.icsa.items %in% gv.cfa.removal]
gv.icsa.times <- gv.icsa.times[!gv.icsa.times %in% gv.cfa.removal]
gv.icsu.items <- gv.icsu.items[!gv.icsu.items %in% gv.cfa.removal]
gv.icsu.times <- gv.icsu.times[!gv.icsu.times %in% gv.cfa.removal]
gv.seta.items <- gv.seta.items[!gv.seta.items %in% gv.cfa.removal]
gv.seta.times <- gv.seta.times[!gv.seta.times %in% gv.cfa.removal]
gv.setb.items <- gv.setb.items[!gv.setb.items %in% gv.cfa.removal]
gv.setb.times <- gv.setb.times[!gv.setb.times %in% gv.cfa.removal]
gv.anchor.items <- gv.anchor.items[!gv.anchor.items %in% gv.cfa.removal]
gv.anchor.times <- gv.anchor.times[!gv.anchor.times %in% gv.cfa.removal]

#---                                                                                     ---#
################################## Gv Final Rasch Analysis ##################################
#---                                                                                    ---#

for (i in 1:5) gv.imputed[[i]][gv.icsu.items] <- lapply(gv.imputed[[i]][gv.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric

# run rasch 3
gv.rasch.3 <- list()
for (i in 1:5) gv.rasch.3[[i]] <- rasch.outcomes(chcability = "gv", raschobject = gv.imputed[[i]][gv.icsu.items], raschnumber = 3, imp = i)

# parameters
gv.classic.parameters <- list()
for (i in 1:5) gv.classic.parameters[[i]] <- coef(gv.rasch.3[[i]][[1]], IRTpars = TRUE, simplify = TRUE)$items
gv.classic.parameters

gv.parameters <- list()
for (i in 1:5) gv.parameters[[i]] <- coef(gv.rasch.3[[i]][[1]], simplify = TRUE)$items
gv.parameters

# pool the parameters
gv.rasch.average <- list()
for (i in 1:5) gv.rasch.average[[i]] <- as.data.frame(coef(gv.rasch.3[[i]][[1]], printSE = TRUE, as.data.frame = TRUE))
for (i in 1:5) gv.rasch.average[[i]] <- rownames_to_column(gv.rasch.average[[i]], var = "item")
for (i in 1:5) gv.rasch.average[[i]] <- filter(gv.rasch.average[[i]], par > -10, par < 10, par != 1)
gv.item.list <- lapply(gv.rasch.average, function(x) { x$item})
gv.par.list <- lapply(gv.rasch.average, function(x) { x$par })
gv.separ.list <- lapply(gv.rasch.average, function(x) { x$SE })
gv.rasch.average <- averageMI(gv.par.list,gv.separ.list,as.data.frame = TRUE)
gv.rasch.average$item <- str_remove(gv.item.list[[1]], "score_")
gv.rasch.average$item <- str_remove(gv.rasch.average$item, ".d")
gv.rasch.average$a1 <- 1
gv.rasch.average$g <- 0
gv.rasch.average$u <- 1
colnames(gv.rasch.average)[1] <- "d"
gv.rasch.average <- select(gv.rasch.average, item, a1, d, g, u)
gv.rasch.average <- filter(gv.rasch.average, grepl(pattern = "gv", item))

# get f scores
gv.fscores <- gv.imputed[[1]]
gv.fscores <- mutate(gv.fscores, gvscore = rowSums(select(gv.fscores, starts_with("score_gv"))))
gv.fscores$fscores <- fscores(object = gv.rasch.3[[1]][[1]])
gv.fscores <- select(gv.fscores, -starts_with("score_gv"), -starts_with("timeTaken"), -set)
gv.wisc <- select(filter(data, gvscore > 0, wisc.fsiq.comp > 0), id, starts_with("wisc"))
gv.fscores <- left_join(gv.wisc,gv.fscores)

# save parameters, remaining items and f scores
write.table(gv.rasch.average, file = "gv_rasch_parameters.csv", sep = ",", quote = FALSE, row.names = F)
write(gv.icsu.items, file = "gvitems.txt")
write.table(gv.fscores, file = "gv_fscores.csv", sep = ",", quote = FALSE, row.names = F)

#---                                                                               ---#
################################## Gwm Data Imputation ##################################
#---                                                                               ---#

# initial
gwm.cleaned <- filter(gwm, phase != "ITOS")
gwm.cleaned$gwmmissing.initial <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.initial) # obtain missing data by participant
gwm.missing1 <- select(gwm.cleaned, id, gwmmissing.initial)

# remove people who dont have a full set of items for anchor items
gwm.cleaned <- filter(gwm.cleaned, !is.na(gwm.anchor)) # set new data sets
gwm.cleaned$gwmmissing.fullanchor <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.fullanchor)
gwm.missing2 <- select(gwm.cleaned, id, gwmmissing.fullanchor)

# remove items that were removed as part of itos
gwm.itos.removed <- c(1,2,3,4,5,6,7,8,9,10)
gwm.itos.removed <- c(paste0("score_gwm",gwm.itos.removed),paste0("timeTaken_gwm",gwm.itos.removed))
gwm.cleaned <- select(gwm.cleaned, -all_of(gwm.itos.removed)) # remove unused items
gwm.cleaned$gwmmissing.removeditems <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.removeditems)
gwm.missing3 <- select(gwm.cleaned, id, gwmmissing.removeditems)

# remove items that have low rates of response
gwm.low.response <- apply(select(gwm.cleaned, starts_with("score_gwm"), starts_with("timeTaken")), 2, percentage.missing) > 75
gwm.low.response <- names(gwm.low.response[gwm.low.response == TRUE])
gwm.low.response

gwm.icsa.items <- gwm.icsa.items[!gwm.icsa.items %in% gwm.low.response]
gwm.icsa.times <- gwm.icsa.times[!gwm.icsa.times %in% gwm.low.response]
gwm.icsu.items <- gwm.icsu.items[!gwm.icsu.items %in% gwm.low.response]
gwm.icsu.times <- gwm.icsu.times[!gwm.icsu.times %in% gwm.low.response]
gwm.seta.items <- gwm.seta.items[!gwm.seta.items %in% gwm.low.response]
gwm.seta.times <- gwm.seta.times[!gwm.seta.times %in% gwm.low.response]
gwm.setb.items <- gwm.setb.items[!gwm.setb.items %in% gwm.low.response]
gwm.setb.times <- gwm.setb.times[!gwm.setb.times %in% gwm.low.response]
gwm.anchor.items <- gwm.anchor.items[!gwm.anchor.items %in% gwm.low.response]
gwm.anchor.times <- gwm.anchor.times[!gwm.anchor.times %in% gwm.low.response]

gwm.cleaned <- select(gwm.cleaned, -all_of(gwm.low.response))
gwm.cleaned$gwmmissing.lowresponse <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.lowresponse)
gwm.missing4 <- select(gwm.cleaned, id, gwmmissing.lowresponse)

# remove participants with high rates of missingness
gwm.cleaned <- filter(gwm.cleaned, gwmmissing.lowresponse < 33)
gwm.cleaned$gwmmissing.removedhighmissingness <- apply(select(gwm.cleaned, c(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.removedhighmissingness)
gwm.missing5 <- select(gwm.cleaned, id, gwmmissing.removedhighmissingness)

# remove participants with missing age
gwm.cleaned <- filter(gwm.cleaned, !is.na(age))
gwm.cleaned$gwmmissing.missingage <- apply(select(gwm.cleaned, c(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.missingage)
gwm.missing6 <- select(gwm.cleaned, id, gwmmissing.missingage)

# recalculate times and scores for cleaned data
gwm.cleaned <- mutate(gwm.cleaned, gwmscore = rowSums(select(gwm.cleaned, starts_with("score_gwm")), na.rm = TRUE), 
                     gwmtime = rowSums(select(gwm.cleaned, starts_with("timeTaken_gwm")), na.rm = TRUE))

# remove univariate time outliers by 1st and 99th percentile
gwm.cleaned <- filter(gwm.cleaned, gwmtime > quantile(gwm.cleaned$gwmtime, 0.05) & gwmtime < quantile(gwm.cleaned$gwmtime, 0.95))
gwm.cleaned$gwmmissing.univariatetime <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.univariatetime) # obtain missing data by participant
gwm.missing7 <- select(gwm.cleaned, id, gwmmissing.univariatetime)

# remove univariate score outliers by 1st and 99th percentile
gwm.cleaned <- filter(gwm.cleaned, gwmscore > quantile(gwm.cleaned$gwmscore, 0.01) & gwmscore < quantile(gwm.cleaned$gwmscore, 0.99))
gwm.cleaned$gwmmissing.univariatescore <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.univariatescore) # obtain missing data by participant
gwm.missing8 <- select(gwm.cleaned, id, gwmmissing.univariatescore)

# remove multivariate outliers
gwm.lm.model <- paste("gwmscore ~ age + gwmtime + gender + phase")
gwm.lm <- lm(gwm.lm.model, data = gwm.cleaned)
gwm.cd <- cooks.distance(gwm.lm)
plot(gwm.cd, pch = "*", cex = 2, main = "Influential Obs by Cooks Distance")
abline(h = 4*mean(gwm.cd, na.rm = T), col = "red")
text(x=1:length(gwm.cd)+1, y=gwm.cd, labels=ifelse(gwm.cd>4*mean(gwm.cd, na.rm = TRUE ), names(gwm.cd),""), col = "red")
gwm.influential <- as.numeric(names(gwm.cd)[(gwm.cd > 4*mean(gwm.cd, na.rm = TRUE))])
gwm.influential <- gwm.influential[!is.na(gwm.influential)]
gwm.cleaned <- gwm.cleaned[-gwm.influential,]
gwm.cleaned$gwmmissing.multivariate <- apply(select(gwm.cleaned, all_of(gwm.icsu.items)), 1, percentage.missing)
psych::describe(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing)) # obtain missing data by column/item
psych::describe(gwm.cleaned$gwmmissing.multivariate) # obtain missing data by participant
gwm.missing9 <- select(gwm.cleaned, id, gwmmissing.multivariate)

# create set column
gwm.cleaned$set <- "cleaned"

# create df of item data missing
gwm.missing.prior <- round(as.data.frame(apply(select(gwm, starts_with("score_gwm")), 2, percentage.missing))) # select original Gwm data
names(gwm.missing.prior)[1] <- "missing"
gwm.missing.prior$stage <- "Original"
gwm.missing.prior <- rownames_to_column(gwm.missing.prior, "item")
gwm.missing.after <- round(as.data.frame(apply(select(gwm.cleaned, starts_with("score_gwm")), 2, percentage.missing))) # select cleaned data
names(gwm.missing.after)[1] <- "missing"
gwm.missing.after$stage <- "After Cleaning"
gwm.missing.after <- rownames_to_column(gwm.missing.after, "item")
gwm.missing.items <- rbind(gwm.missing.prior,gwm.missing.after)
gwm.missing.items$item <- as.numeric(str_extract(gwm.missing.items$item, "[0-9]+"))
rm(gwm.missing.after,gwm.missing.prior)

# Visualise missing data by item
ggplot(gwm.missing.items, aes(x = item, y = missing, colour = stage)) +
  geom_point() +
  ylab("Percentage Missing") +
  xlab("Item Number") +
  scale_x_continuous(breaks = c(seq.int(1,107,3))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom") +
  labs(colour = "Stage") +
  geom_vline(xintercept = 22) +
  geom_vline(xintercept = 33) +
  annotate("text", x = 11, y = 50, label = "Set A Items") +
  annotate("text", x = 27, y = 50, label = "Anchor Items") +
  annotate("text", x = 39, y = 50, label = "Set B Items")

# combine missing data columns
gwm.missing.participants <- full_join(gwm.missing1, gwm.missing2)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing3)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing4)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing5)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing6)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing7)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing8)
gwm.missing.participants <- full_join(gwm.missing.participants, gwm.missing9)
rm(gwm.missing1,gwm.missing2,gwm.missing3,gwm.missing4,gwm.missing5,gwm.missing6,gwm.missing7,gwm.missing8,gwm.missing9)

# visualise the changes in missing data across steps
ggplot(reshape2::melt(gwm.missing.participants[,c("gwmmissing.initial","gwmmissing.fullanchor","gwmmissing.removeditems","gwmmissing.lowresponse",
                                                 "gwmmissing.removedhighmissingness", "gwmmissing.missingage","gwmmissing.univariatetime",
                                                 "gwmmissing.univariatescore","gwmmissing.multivariate")]), 
       aes(x = variable, y = value, colour = variable)) +
  geom_boxplot(width = .25) +
  ylab("Percentage of Items Missing") +
  scale_x_discrete(name = "", labels = ggplot_addline(c("Full Data","Full Anchor","Items Removed","Low Item Response Rate",
                                                        "High Missingness Removed","Missing Age Removed",
                                                        "Outliers by Time","Outliers by Score",
                                                        "Multivariate Outliers"))) +
  theme(legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text")

# Recalculate the set missingness
gwm.cleaned <- mutate(gwm.cleaned, 
                     gwm.set.a = ceiling(rowSums(is.na(gwm.cleaned[gwm.seta.items]))/length(gwm.seta.items)),
                     gwm.set.b = ceiling(rowSums(is.na(gwm.cleaned[gwm.setb.items]))/length(gwm.setb.items)),
                     gwm.anchor = ceiling(rowSums(is.na(gwm.cleaned[gwm.anchor.items]))/length(gwm.anchor.items)))

# Recode missingness # NA = missing data in set, 1 = full data
gwm.cleaned$gwm.set.a <- gwm.cleaned$gwm.set.a %>% na_if(1) %>% recode(`0` = 1)
gwm.cleaned$gwm.set.b <- gwm.cleaned$gwm.set.b %>% na_if(1) %>% recode(`0` = 1)
gwm.cleaned$gwm.anchor <- gwm.cleaned$gwm.anchor %>% na_if(1) %>% recode(`0` = 1)

# Visualise Missing Data Patterns
aggr(gwm.cleaned[c("gwm.set.a","gwm.set.b","gwm.anchor")], combined = TRUE,
     labels = c("Set A", "Set B", "Anchor")) 

# Set Up Data Set
gwm.imputed <- select(gwm.cleaned, id, phase, age, age.group, gender, nationality, device, all_of(gwm.icsu.items), all_of(gwm.icsu.times)) # create data set
gwm.imputed[gwm.icsu.items] <- lapply(gwm.imputed[gwm.icsu.items], factor)

# prepare prediction model
gwm.qp <- quickpred(gwm.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# flag logged events items
gwm.ini <- mice(gwm.imputed, maxit = 1, pred = gwm.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)
gwm.outlist.logged <- as.character(gwm.ini$loggedEvents[, "out"])
gwm.outlist.logged

# flag outflux problem items
gwm.flux <- flux(gwm.imputed)
gwm.fluxplot <- fluxplot(gwm.imputed, main = NULL)
gwm.outlist.outflux <- row.names(gwm.flux)[gwm.flux$outflux < 0.3]
gwm.outlist.outflux

# prepare prediction model
gwm.qp <- quickpred(gwm.imputed, mincor = 0.2, include = c("age"), exclude = c("gender","nationality","phase","id","age.group","device"))

# conduct mi - impute data using random forests
gwm.mice <- mice(gwm.imputed, pred = gwm.qp, m = 5, method = 'rf', visitSequence = "monotone", print = FALSE)

# create data sets for use in analyses
gwm.imputed <- NULL # create empty variable
for (i in 1:5) gwm.imputed[[i]] <- mice::complete(gwm.mice, action = i, inc = FALSE) # put imputed data into list
for (i in 1:5) gwm.imputed[[i]][gwm.icsu.items] <- lapply(gwm.imputed[[i]][gwm.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric
for (i in 1:5) gwm.imputed[[i]]$set <-  paste0("imp",i)
gwm.imputed <- lapply(gwm.imputed, function(x) {
  x <- mutate(x, gwmscore = rowSums(select(x, starts_with("score_gwm")), na.rm = TRUE))
  x <- mutate(x, gwmtime = rowSums(select(x, starts_with("timeTaken_gwm")), na.rm = TRUE))
})
gwm.backup <- gwm.imputed

# recalculate times and scores for cleaned data
gwm.cleaned <- mutate(gwm.cleaned, gwmscore = rowSums(select(gwm.cleaned, starts_with("score_gwm")), na.rm = TRUE), 
                     gwmtime = rowSums(select(gwm.cleaned, starts_with("timeTaken_gwm")), na.rm = TRUE))

# Combine the results together
gwm.comparison <- rbind(gwm.imputed[[1]][c("set","gwmscore","gwmtime")],gwm.cleaned[c("set","gwmscore","gwmtime")])
gwm.comparison <- rbind(gwm.imputed[[2]][c("set","gwmscore","gwmtime")],gwm.comparison)
gwm.comparison <- rbind(gwm.imputed[[3]][c("set","gwmscore","gwmtime")],gwm.comparison)
gwm.comparison <- rbind(gwm.imputed[[4]][c("set","gwmscore","gwmtime")],gwm.comparison)
gwm.comparison <- rbind(gwm.imputed[[5]][c("set","gwmscore","gwmtime")],gwm.comparison)

# Boxplot of total score
ggplot(gwm.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gwmscore, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(0,30,5)), limits = c(0,30), name = "Total Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# Boxplot of time
ggplot(gwm.comparison) + 
  geom_boxplot(mapping = aes(x = set, y = gwmtime/60, fill = set), width = 0.25) +
  scale_x_discrete(name = "Data Set", labels = c("Cleaned", 
                                                 "Imputation 1", 
                                                 "Imputation 2", 
                                                 "Imputation 3", 
                                                 "Imputation 4", 
                                                 "Imputation 5")) +
  scale_y_continuous(breaks = c(seq(5,max(gwm.cleaned$gwmtime)/60,5)), limits = c(5,max(gwm.cleaned$gwmtime)/60), name = "Total Time") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5), legend.position = "none")

# visualise proportion correct
gwm.cleaned.correct <- sapply(gwm.cleaned[gwm.icsu.items], function(x) {(1 - (table(x)[[1]]/length(na.omit(x)))) * 100 })
gwm.imputed.correct <- lapply(gwm.imputed, function(x) { 
  sapply(x[gwm.icsu.items], function(y) {
    (1 - (table(y)[[1]]/length(na.omit(y)))) * 100 
  })
})
gwm.combined.correct <- data.frame(cleaned = gwm.cleaned.correct, imp1 = gwm.imputed.correct[[1]], imp2 = gwm.imputed.correct[[2]],
                                  imp3 = gwm.imputed.correct[[3]], imp4 = gwm.imputed.correct[[4]], imp5 = gwm.imputed.correct[[5]])
row.names(gwm.combined.correct) <- str_remove(row.names(gwm.combined.correct), "score_gwm")
gwm.combined.correct <- rownames_to_column(gwm.combined.correct, var = "item")
gwm.combined.correct <- pivot_longer(data = gwm.combined.correct, cols = -item, names_to = "set", values_to = "correct")
gwm.combined.correct$item <- factor(gwm.combined.correct$item, levels = c(sort(unique(as.numeric(gwm.combined.correct$item)))))
gwm.combined.correct$set <- factor(gwm.combined.correct$set)
ggplot(gwm.combined.correct, aes(x = item, y = correct, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(name = "Item Number", breaks = seq(1, length(gwm.combined.correct$item), by = 2)) +
  scale_y_continuous(breaks = seq(0,100,10), limits = c(0,100), name = "Percentage Correct") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))

# visualise performance by age
gwm.cleaned.age <- select(gwm.cleaned, age, gwmscore)
gwm.imputed.age <- lapply(gwm.imputed, function(x) { select(x, gwmscore)})
gwm.combined.age <- data.frame(age = gwm.cleaned.age$age, cleaned = gwm.cleaned.age$gwmscore, imp1 = gwm.imputed.age[[1]]$gwmscore, 
                              imp2 = gwm.imputed.age[[2]]$gwmscore, imp3 = gwm.imputed.age[[3]]$gwmscore, imp4 = gwm.imputed.age[[4]]$gwmscore, 
                              imp5 = gwm.imputed.age[[5]]$gwmscore)
gwm.combined.age <- pivot_longer(data = gwm.combined.age, cols = -age, names_to = "set", values_to = "score")
gwm.combined.age$set <- factor(gwm.combined.age$set)
ggplot(gwm.combined.age, aes(x = age, y = score, group = set, colour = set, linetype = set)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(6,90,2), limits = c(6,90), name = "Age") +
  scale_y_continuous(breaks = seq(0,30,5), limits = c(0,30), name = "Total Score") +
  scale_colour_discrete(name = "Data Set", labels = c("Cleaned", 
                                                      "Imputation 1", 
                                                      "Imputation 2", 
                                                      "Imputation 3", 
                                                      "Imputation 4", 
                                                      "Imputation 5")) +
  scale_linetype_discrete(name = "Data Set", labels = c("Cleaned", 
                                                        "Imputation 1", 
                                                        "Imputation 2", 
                                                        "Imputation 3", 
                                                        "Imputation 4", 
                                                        "Imputation 5")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
  geom_smooth(method = "lm")

#---                                                                           ---#
################################## Gwm Reliability ##################################
#---                                                                          ---#

gwm.reliability <- NULL
for (i in 1:5) gwm.reliability[[i]] <- psych::alpha(gwm.imputed[[i]][gwm.icsu.items], warnings = FALSE)$total
gwm.reliability <- sapply(gwm.reliability, function(x) { x[,"raw_alpha"] })
mean(gwm.reliability)

#---                                                                              ---#
################################## Gwm Rasch Analysis ##################################
#---                                                                              ---#

# run rasch
gwm.rasch <- list()
for (i in 1:5) gwm.rasch[[i]] <- rasch.outcomes(chcability = "gwm", raschobject = gwm.imputed[[i]][gwm.icsu.items], raschnumber = 1, imp = i)

# person fit
gwm.rasch.person <- which(personfit(gwm.rasch[[1]][[1]])$Zh < -2)
gwm.rasch.person <- append(gwm.rasch.person, which(personfit(gwm.rasch[[1]][[1]])$Zh > 2))
gwm.rasch.person
gwm.imputed[[1]][gwm.rasch.person,]$age

# wright map
WrightMap::wrightMap(fscores(gwm.rasch[[1]][[1]]), main.title = NA)

# determine poor item fits
gwm.poor.fit <- NULL
for (i in 1:5) gwm.poor.fit[[i]] <- rasch.poor.fit(dataset = gwm.imputed[[i]], raschobject = gwm.rasch[[i]][[2]], imp = i)

#---                                                                                ---#
################################## Gwm Mokken Analysis ##################################
#---                                                                               ---#

# conduct first mokken
gwm.mokken <- list()
for (i in 1:5) gwm.mokken[[i]] <- mokken.outcomes(chcability = "gwm", mokkenobject = gwm.imputed[[i]][gwm.icsu.items], mokkennumber = 1, imp = i)

# determine poorly ordered items and put them in a vector for removal
gwm.mokken.prob <- NULL
for (i in 1:5) gwm.mokken.prob[[i]] <- rownames(gwm.mokken[[i]]$Hi)[which(gwm.mokken[[i]]$Hi < .3)]
gwm.mokken.prob
gwm.mokken.prob <- sort(unique(unlist(gwm.mokken.prob)))
gwm.mokken.prob

gwm.mokken.removal <- c(37:44) # this is a qualitative decision
gwm.mokken.removal <- c(paste0("score_gwm",gwm.mokken.removal),paste0("timeTaken_gwm",gwm.mokken.removal))

# remove the items that do not fit well to the rasch model
for (i in 1:5) gwm.imputed[[i]] <- select(gwm.imputed[[i]], -all_of(gwm.mokken.removal))
# remove.items.from.sets(chcability = "gwm", items = gwm.mokken.removal)
gwm.icsa.items <- gwm.icsa.items[!gwm.icsa.items %in% gwm.mokken.removal]
gwm.icsa.times <- gwm.icsa.times[!gwm.icsa.times %in% gwm.mokken.removal]
gwm.icsu.items <- gwm.icsu.items[!gwm.icsu.items %in% gwm.mokken.removal]
gwm.icsu.times <- gwm.icsu.times[!gwm.icsu.times %in% gwm.mokken.removal]
gwm.seta.items <- gwm.seta.items[!gwm.seta.items %in% gwm.mokken.removal]
gwm.seta.times <- gwm.seta.times[!gwm.seta.times %in% gwm.mokken.removal]
gwm.setb.items <- gwm.setb.items[!gwm.setb.items %in% gwm.mokken.removal]
gwm.setb.times <- gwm.setb.times[!gwm.setb.times %in% gwm.mokken.removal]
gwm.anchor.items <- gwm.anchor.items[!gwm.anchor.items %in% gwm.mokken.removal]
gwm.anchor.times <- gwm.anchor.times[!gwm.anchor.times %in% gwm.mokken.removal]

# conduct second mokken
gwm.mokken.2 <- list()
for (i in 1:5) gwm.mokken.2[[i]] <- mokken.outcomes(chcability = "gwm", mokkenobject = gwm.imputed[[i]][gwm.icsu.items], mokkennumber = 2, imp = i)

#---                                                                                 ---#
################################## Gwm Local Dependency ##################################
#---                                                                                 ---#

gwm.local.dependency <- list()
for (i in 1:5) gwm.local.dependency[[i]] <- residuals(gwm.rasch[[i]][[1]], type = "Q3", digits = 2, suppress = +.2)

# determine items with high local dependency
gwm.local.dependency.removal <- ""

# remove the items that do not fit well to the rasch model
for (i in 1:5) gwm.imputed[[i]] <- select(gwm.imputed[[i]], -all_of(gwm.local.dependency.removal))
# remove.items.from.sets(chcability = "gwm", items = gwm.local.dependency.removal)
gwm.icsa.items <- gwm.icsa.items[!gwm.icsa.items %in% gwm.local.dependency.removal]
gwm.icsa.times <- gwm.icsa.times[!gwm.icsa.times %in% gwm.local.dependency.removal]
gwm.icsu.items <- gwm.icsu.items[!gwm.icsu.items %in% gwm.local.dependency.removal]
gwm.icsu.times <- gwm.icsu.times[!gwm.icsu.times %in% gwm.local.dependency.removal]
gwm.seta.items <- gwm.seta.items[!gwm.seta.items %in% gwm.local.dependency.removal]
gwm.seta.times <- gwm.seta.times[!gwm.seta.times %in% gwm.local.dependency.removal]
gwm.setb.items <- gwm.setb.items[!gwm.setb.items %in% gwm.local.dependency.removal]
gwm.setb.times <- gwm.setb.times[!gwm.setb.times %in% gwm.local.dependency.removal]
gwm.anchor.items <- gwm.anchor.items[!gwm.anchor.items %in% gwm.local.dependency.removal]
gwm.anchor.times <- gwm.anchor.times[!gwm.anchor.times %in% gwm.local.dependency.removal]

#---                                                                                              ---#
################################## Gwm Differential Item Functioning ##################################
#---                                                                                              ---#

gwm.dif.gender <- list()
for (i in 1:5) gwm.dif.gender[[i]] <- difLogistic(Data = gwm.imputed[[i]][gwm.icsu.items], group = gwm.imputed[[i]]$gender, focal.name = "m", alpha = .01)
for (i in 1:5) plot(gwm.dif.gender[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gwm.dif.nationality <- list()
for (i in 1:5) gwm.dif.nationality[[i]] <- difLogistic(Data = gwm.imputed[[i]][gwm.icsu.items], group = gwm.imputed[[i]]$nationality, focal.name = "au", alpha = .01)
for (i in 1:5) plot(gwm.dif.nationality[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

gwm.dif.device <- list()
for (i in 1:5) gwm.dif.device[[i]] <- difLogistic(Data = gwm.imputed[[i]][gwm.icsu.items], group = gwm.imputed[[i]]$device, focal.name = "iphone", alpha = .01)
for (i in 1:5) plot(gwm.dif.device[[i]], plot = "lrStat") # igore the "the plot was not captured" warning, this is for saving

# determine dif items
gwm.dif.items <- NULL
for (i in 1:5) gwm.dif.items[[i]] <- { 
  gender <- gwm.dif.gender[[i]]$DIFitems
  nationality <- gwm.dif.nationality[[i]]$DIFitems
  device <- gwm.dif.device[[i]]$DIFitems
  list(gender,nationality,device)
}

gwm.dif.items <- sort(as.numeric(unique(unlist(gwm.dif.items))), na.last = NA)
gwm.dif.items
gwm.dif.items <- gwm.icsu.items[gwm.dif.items]
gwm.dif.items

for (y in 1:5) for (i in gwm.dif.items) { plot(gwm.dif.gender[[y]], plot = "itemCurve", item = i) }
for (y in 1:5) for (i in gwm.dif.items) { plot(gwm.dif.nationality[[y]], plot = "itemCurve", item = i) }
for (y in 1:5) for (i in gwm.dif.items) { plot(gwm.dif.device[[y]], plot = "itemCurve", item = i) }

#---                                                                                             ---#
################################## Gwm Confirmatory Factor Analysis ##################################
#---                                                                                            ---#

# run the CFA on the multivariate imputation data sets
for (i in 1:5) gwm.imputed[[i]][gwm.icsu.items] <- lapply(gwm.imputed[[i]][gwm.icsu.items], function(x) { x <- ordered(x, levels = c(0,1)) }) # Change all items to ordered factor
gwm.model <- paste("gwm =~ ", paste(gwm.icsu.items, collapse = "+"), sep = "") # establish model
gwm.cfa.mi <- runMI(gwm.model, data = gwm.imputed, fun = "cfa", meanstructure = TRUE) # run cfa --> takes a very long time
gwm.cfa.mi.summary <- summary(gwm.cfa.mi)
gwm.cfa.mi.fit <- fitMeasures(gwm.cfa.mi, fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))

# determine what cfas were successful
gwm.cfa.success <- NULL
for (i in 1:5) gwm.cfa.success[i] <- !is.null(gwm.cfa.mi@miList[[i]])

# based on the cfas that were successful, run them 1 by 1 for diagnostic and reporting purposes, fit the cfa
gwm.cfa <- list()
for (i in 1:5) gwm.cfa[[i]] <- if(gwm.cfa.success[i] == TRUE) {cfa(gwm.model, data = gwm.imputed[[i]])}

# store the cfa summary for item selection later
gwm.cfa.summary <- NULL
for (i in 1:5) gwm.cfa.summary[[i]] <- if(gwm.cfa.success[i] == TRUE) { summary(object = gwm.cfa[[i]], standardized = TRUE) }

# print and save the parameters using predefined function
gwm.cfa.parameters <- NULL
for (i in 1:5) gwm.cfa.parameters[[i]] <- if(gwm.cfa.success[i] == TRUE) {
  cfa.parameters(chcability = "gwm", cfaobject = gwm.cfa[[i]], cfanumber = 1, imp = i)}

# print and store the overall fit
gwm.cfa.fit <- NULL
for (i in 1:5) gwm.cfa.fit[[i]] <- if(gwm.cfa.success[i] == TRUE) {
  fitMeasures(gwm.cfa[[i]], fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
}
gwm.cfa.fit

# determine low loading items and put them in a vector for removal
gwm.cfa.poor.fit <- list()
for (i in 1:5) gwm.cfa.poor.fit[[i]] <- if(gwm.cfa.success[i] == TRUE) {
  select(filter(gwm.cfa.summary[[i]]$PE, lhs == "gwm", std.lv < .3, rhs != ""), rhs)$rhs } else {NA}
gwm.cfa.poor.fit
gwm.cfa.poor.fit <- sort(unique(unlist(gwm.cfa.poor.fit)))
gwm.cfa.poor.fit

# determine items to remove
gwm.cfa.removal <- "score_gwm32"

# remove the items that do not fit well to the rasch model
for (i in 1:5) gwm.imputed[[i]] <- select(gwm.imputed[[i]], -all_of(gwm.cfa.removal))
# remove.items.from.sets(chcability = "gwm", items = gwm.cfa.removal)
gwm.icsa.items <- gwm.icsa.items[!gwm.icsa.items %in% gwm.cfa.removal]
gwm.icsa.times <- gwm.icsa.times[!gwm.icsa.times %in% gwm.cfa.removal]
gwm.icsu.items <- gwm.icsu.items[!gwm.icsu.items %in% gwm.cfa.removal]
gwm.icsu.times <- gwm.icsu.times[!gwm.icsu.times %in% gwm.cfa.removal]
gwm.seta.items <- gwm.seta.items[!gwm.seta.items %in% gwm.cfa.removal]
gwm.seta.times <- gwm.seta.times[!gwm.seta.times %in% gwm.cfa.removal]
gwm.setb.items <- gwm.setb.items[!gwm.setb.items %in% gwm.cfa.removal]
gwm.setb.times <- gwm.setb.times[!gwm.setb.times %in% gwm.cfa.removal]
gwm.anchor.items <- gwm.anchor.items[!gwm.anchor.items %in% gwm.cfa.removal]
gwm.anchor.times <- gwm.anchor.times[!gwm.anchor.times %in% gwm.cfa.removal]

#---                                                                                     ---#
################################## Gwm Final Rasch Analysis ##################################
#---                                                                                    ---#

for (i in 1:5) gwm.imputed[[i]][gwm.icsu.items] <- lapply(gwm.imputed[[i]][gwm.icsu.items], function(x) { x <- as.numeric(levels(x)[x]) }) # Change all items to numeric

# run rasch 2
gwm.rasch.2 <- list()
for (i in 1:5) gwm.rasch.2[[i]] <- rasch.outcomes(chcability = "gwm", raschobject = gwm.imputed[[i]][gwm.icsu.items], raschnumber = 2, imp = i)

# parameters
gwm.classic.parameters <- list()
for (i in 1:5) gwm.classic.parameters[[i]] <- coef(gwm.rasch.2[[i]][[1]], IRTpars = TRUE, simplify = TRUE)$items
gwm.classic.parameters

gwm.parameters <- list()
for (i in 1:5) gwm.parameters[[i]] <- coef(gwm.rasch.2[[i]][[1]], simplify = TRUE)$items
gwm.parameters

# pool the parameters
gwm.rasch.average <- list()
for (i in 1:5) gwm.rasch.average[[i]] <- as.data.frame(coef(gwm.rasch.2[[i]][[1]], printSE = TRUE, as.data.frame = TRUE))
for (i in 1:5) gwm.rasch.average[[i]] <- rownames_to_column(gwm.rasch.average[[i]], var = "item")
for (i in 1:5) gwm.rasch.average[[i]] <- filter(gwm.rasch.average[[i]], par > -10, par < 10, par != 1)
gwm.item.list <- lapply(gwm.rasch.average, function(x) { x$item})
gwm.par.list <- lapply(gwm.rasch.average, function(x) { x$par })
gwm.separ.list <- lapply(gwm.rasch.average, function(x) { x$SE })
gwm.rasch.average <- averageMI(gwm.par.list,gwm.separ.list,as.data.frame = TRUE)
gwm.rasch.average$item <- str_remove(gwm.item.list[[1]], "score_")
gwm.rasch.average$item <- str_remove(gwm.rasch.average$item, ".d")
gwm.rasch.average$a1 <- 1
gwm.rasch.average$g <- 0
gwm.rasch.average$u <- 1
colnames(gwm.rasch.average)[1] <- "d"
gwm.rasch.average <- select(gwm.rasch.average, item, a1, d, g, u)
gwm.rasch.average <- filter(gwm.rasch.average, grepl(pattern = "gwm", item))

# get f scores
gwm.fscores <- gwm.imputed[[1]]
gwm.fscores <- mutate(gwm.fscores, gwmscore = rowSums(select(gwm.fscores, starts_with("score_gwm"))))
gwm.fscores$fscores <- fscores(object = gwm.rasch.2[[1]][[1]])
gwm.fscores <- select(gwm.fscores, -starts_with("score_gwm"), -starts_with("timeTaken"), -set)
gwm.wisc <- select(filter(data, gwmscore > 0, wisc.fsiq.comp > 0), id, starts_with("wisc"))
gwm.fscores <- left_join(gwm.wisc,gwm.fscores)

# save parameters, remaining items and f scores
write.table(gwm.rasch.average, file = "gwm_rasch_parameters.csv", sep = ",", quote = FALSE, row.names = F)
write(gwm.icsu.items, file = "gwmitems.txt")
write.table(gwm.fscores, file = "gwm_fscores.csv", sep = ",", quote = FALSE, row.names = F)
