# Data Analysis for the Simulation of the CHC-CAT Test
# Original code developed by Mr Jake Kraksa (Monash University)
# Built on R 3.6.3 with R Studio 1.2.5033
# set the working directory to source file directory manually
# restart R session prior to beginning

#---                                                                      ---#
################################## Libraries ##################################
#---                                                                      ---#

# Data exploration

library(psych) # version 1.9.12

# CAT

library(mirtCAT) # version 1.9.3

# Visualise Results

library(ggplot2) # version 3.3.0
library(knitr) # version 1.28
library(latticeExtra)

# Data cleaning

library(plyr) # version 1.8.6
library(dplyr) # version 0.8.5
library(stringr) # version 1.4.0
library(tidyr) # version 1.0.2
library(tibble)

# Parellel computation

library(parallel) # version 3.6.3

#---                                                                      ---#
################################## Set Options ##################################
#---                                                                      ---#

### Set Seed for replication
set.seed(145)

### Set print options
options(max.print = 10000)

#---                                                                             ---#
################################## Load Functions ##################################
#---                                                                             ---#

percentage_missing <- function(x) { sum(is.na(x))/length(x)*100 } # evaluates a list and determines the percentage of NA

calculate_sem <- function(rxx,sd) {sd * sqrt(1-rxx)} # calculates sem

cat_simulation <- function(chc_mo, sem) {
  Theta <- matrix(rnorm(n = 5000, mean = 0, sd = 1)) # Simulate theta for 5000 people
  pattern <- generate_pattern(mo = chc_mo, Theta = Theta) ### Generate response pattern
  cat_design <- list(min_SEM = sem) ### the design factors of the simulation
  cl <- makeCluster(detectCores()) ### for parallel processing
  cat <- mirtCAT(mo = chc_mo, local_pattern = pattern, start_item = "MI", method = "EAP", criteria = "MI", design = cat_design, cl = cl) # run
  stopCluster(cl) # stop the cluster
  est_Theta <- laply(cat, function(y) y$thetas) # list of theta
  average_items <- mean(laply(cat, function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(Theta[,1], est_Theta) ### correlation between true theta and estimated theta
  bias <- mean(Theta[,1] - est_Theta) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((Theta[,1] - est_Theta)^2))
  results <- data.frame("sem" = sem, "averageItems" = average_items, "correlation" = correlation,
                        "bias" = bias, "rmsd" = rmsd) # put the results into a dataframe
  return(results)
}

#---                                                                      ---#
################################## Load Data ##################################
#---                                                                      ---#

### Load main data
data <- read.csv(file = "simulation_data.csv", stringsAsFactors = FALSE)
data$phase <- recode(data$phase, itos = "ITOS", p2a = "ICS-A", p2u = "ICS-U")
data$phase <- as.factor(data$phase)
data$gender <- as.factor(data$gender)
data$nationality <- as.factor(data$nationality)
data$device <- as.factor(data$device)
colnames(data) <- str_remove(colnames(data), "score_")

#---                                                                            ---#
################################## Simulation Data ##################################
#---                                                                            ---#

reliability <- c(seq(from = 1.0, to = .50, by = -.1))
sem <- calculate_sem(reliability,1)
Theta_groups <- seq(-3,3,.6)

#---                                                                               ---#
################################## Lexical Knowledge  ##################################
#---                                                                               ---#

### --- Get Item List
gc_items <- scan(file = "gcitems.txt", character(), quote = "") # retrieve retained items from ICS
gc_items <- str_remove(gc_items, "score_") # standardise item names consistent with item parameters file

### --- Get ICS Parameters
gc_pars <- read.csv(file = "gc_rasch_parameters.csv") ### Parameters from ICS
rownames(gc_pars) <- gc_pars[, 1] ### put items as row names
gc_pars <- gc_pars[, -1] ### remove item names from df

### --- Clean Data
gc_data_all <- select(filter(data, gcscore > 0), phase, age, gender, nationality, device, all_of(gc_items)) ### raw data for all participants 
gc_data_all <- gc_data_all[c(which(apply(gc_data_all, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gc_data_all <- mutate(gc_data_all, gcscore = rowSums(select(gc_data_all, starts_with("gc")))) ### calculate total score

gc_data_underage <- select(filter(data, gcscore > 0, age < 19), phase, age, gender, nationality, device, all_of(gc_items)) ### raw data for underage participants
gc_data_underage <- gc_data_underage[c(which(apply(gc_data_underage, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gc_data_underage <- mutate(gc_data_underage, gcscore = rowSums(select(gc_data_underage, starts_with("gc")))) ### calculate total score

### --- Explore Data for all data
psych::describe(select(gc_data_all, age, gcscore))

ggplot(data = gc_data_all, mapping = aes(x = age, y = gcscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gc_data_all) + 
  geom_boxplot(mapping = aes(x = phase, y = gcscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gc_data_all$gcscore),5)), limits = c(0,max(gc_data_all$gcscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gc_data_all) + 
  geom_boxplot(mapping = aes(x = nationality, y = gcscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gc_data_all$gcscore),5)), limits = c(0,max(gc_data_all$gcscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gc_data_all) + 
  geom_boxplot(mapping = aes(x = device, y = gcscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gc_data_all$gcscore),5)), limits = c(0,max(gc_data_all$gcscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gc_correct_all <- as.data.frame(sapply(gc_data_all[gc_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gc_correct_all) <- str_remove(row.names(gc_correct_all), "gc")
names(gc_correct_all) <- "Correct"
gc_correct_all <- rownames_to_column(gc_correct_all, var = "Item")
gc_correct_all <- gc_correct_all[order(-gc_correct_all$Correct),]
ggplot(data = gc_correct_all, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Explore data for underage
psych::describe(select(gc_data_underage, age, gcscore))

ggplot(data = gc_data_underage, mapping = aes(x = age, y = gcscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gc_data_underage) + 
  geom_boxplot(mapping = aes(x = phase, y = gcscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gc_data_underage$gcscore),5)), limits = c(0,max(gc_data_underage$gcscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gc_data_underage) + 
  geom_boxplot(mapping = aes(x = nationality, y = gcscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gc_data_underage$gcscore),5)), limits = c(0,max(gc_data_underage$gcscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gc_data_underage) + 
  geom_boxplot(mapping = aes(x = device, y = gcscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gc_data_underage$gcscore),5)), limits = c(0,max(gc_data_underage$gcscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gc_correct_underage <- as.data.frame(sapply(gc_data_underage[gc_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gc_correct_underage) <- str_remove(row.names(gc_correct_underage), "gc")
names(gc_correct_underage) <- "Correct"
gc_correct_underage <- rownames_to_column(gc_correct_underage, var = "Item")
gc_correct_underage <- gc_correct_underage[order(-gc_correct_underage$Correct),]
ggplot(data = gc_correct_underage, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Generate MIRT object using ICS Parameters
gc_mirt_ics <- generate.mirt_object(parameters = gc_pars, itemtype = "Rasch")
plot(gc_mirt_ics, type = "infoSE")
coef(gc_mirt_ics, simplify = TRUE)$item

### -- Simulation for ICS parameters
gc_sim_ics <- list()
for (i in 1:length(sem)) { gc_sim_ics[[i]] <- cat_simulation(chc_mo = gc_mirt_ics, sem = sem[i]) }
gc_sim_ics <- bind_rows(gc_sim_ics)
gc_sim_ics$reliability <- reliability
gc_sim_ics
write.table(x = gc_sim_ics, file = paste0("sim_tables/","gc_sim_ics.csv"), sep = ",", row.names = F)

### --- Generate MIRT object using all age data
gc_mirt_all <- mirt(data = gc_data_all[gc_items], model = 1, itemtype = "Rasch")
plot(gc_mirt_all, type = "infoSE")
coef(gc_mirt_all, simplify = TRUE)$item

### -- Simulation using all age data
gc_sim_all <- list()
for (i in 1:length(sem)) { gc_sim_all[[i]] <- cat_simulation(chc_mo = gc_mirt_all, sem = sem[i]) }
gc_sim_all <- bind_rows(gc_sim_all)
gc_sim_all$reliability <- reliability
gc_sim_all
write.table(x = gc_sim_all, file = paste0("sim_tables/","gc_sim_all.csv"), sep = ",", row.names = F)

### --- Generate MIRT object for underage participants, no multiple imputation required
gc_mirt_underage <- mirt(data = gc_data_underage[gc_items], model = 1, itemtype = "Rasch")
plot(gc_mirt_underage, type = "infoSE")
coef(gc_mirt_underage, simplify = TRUE)$item

### -- Simulation for underage data
gc_sim_underage <- list()
for (i in 1:length(sem)) { gc_sim_underage[[i]] <- cat_simulation(chc_mo = gc_mirt_underage, sem = sem[i]) }
gc_sim_underage <- bind_rows(gc_sim_underage)
gc_sim_underage$reliability <- reliability
gc_sim_underage
write.table(x = gc_sim_underage, file = paste0("sim_tables/","gc_sim_underage.csv"), sep = ",", row.names = F)

### --- Performance of participants by group
gc_Theta <- matrix(rnorm(n = 5000, mean = 0, sd = 1)) # Simulate theta for 5000 people
gc_pattern <- generate_pattern(mo = gc_mirt_underage, Theta = gc_Theta) ### Generate response pattern
gc_cat_design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
gc_cat <- mirtCAT(mo = gc_mirt_underage, local_pattern = gc_pattern, start_item = "MI", method = "EAP", criteria = "MI", design = gc_cat_design, cl = cl) # run
stopCluster(cl) # stop the cluster
gc_est_Theta <- laply(gc_cat, function(y) y$thetas) # list of theta
gc_sim_groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(gc_Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(gc_Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(gc_Theta > Theta_groups[i-1] & gc_Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(gc_cat[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(gc_Theta[rows,1], gc_est_Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(gc_Theta[rows,1] - gc_est_Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((gc_Theta[rows,1] - gc_est_Theta[rows])^2))
  nParticipants <- length(rows)
  gc_sim_groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                  "averageItems" = average_items, "correlation" = correlation, 
                                  "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
  }
gc_sim_groups <- bind_rows(gc_sim_groups)
gc_sim_groups
write.table(x = gc_sim_groups, file = paste0("sim_tables/","gc_sim_groups.csv"), sep = ",", row.names = F)

#---                                                                       ---#
################################## Induction  ##################################
#---                                                                       ---#

### --- Get Item List
gf_items <- scan(file = "gfitems.txt", character(), quote = "") # retrieve retained items from ICS
gf_items <- str_remove(gf_items, "score_") # standardise item names consistent with item parameters file

### --- Get ICS Parameters
gf_pars <- read.csv(file = "gf_rasch_parameters.csv") ### Parameters from ICS
rownames(gf_pars) <- gf_pars[, 1] ### put items as row names
gf_pars <- gf_pars[, -1] ### remove item names from df

### --- Clean Data
gf_data_all <- select(filter(data, gfscore > 0), phase, age, gender, nationality, device, all_of(gf_items)) ### raw data for all participants 
gf_data_all <- gf_data_all[c(which(apply(gf_data_all, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gf_data_all <- mutate(gf_data_all, gfscore = rowSums(select(gf_data_all, starts_with("gf")))) ### calculate total score

gf_data_underage <- select(filter(data, gfscore > 0, age < 19), phase, age, gender, nationality, device, all_of(gf_items)) ### raw data for underage participants
gf_data_underage <- gf_data_underage[c(which(apply(gf_data_underage, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gf_data_underage <- mutate(gf_data_underage, gfscore = rowSums(select(gf_data_underage, starts_with("gf")))) ### calculate total score

### --- Explore Data for all data
psych::describe(select(gf_data_all, age, gfscore))

ggplot(data = gf_data_all, mapping = aes(x = age, y = gfscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gf_data_all) + 
  geom_boxplot(mapping = aes(x = phase, y = gfscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gf_data_all$gfscore),5)), limits = c(0,max(gf_data_all$gfscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gf_data_all) + 
  geom_boxplot(mapping = aes(x = nationality, y = gfscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gf_data_all$gfscore),5)), limits = c(0,max(gf_data_all$gfscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gf_data_all) + 
  geom_boxplot(mapping = aes(x = device, y = gfscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gf_data_all$gfscore),5)), limits = c(0,max(gf_data_all$gfscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gf_correct_all <- as.data.frame(sapply(gf_data_all[gf_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gf_correct_all) <- str_remove(row.names(gf_correct_all), "gf")
names(gf_correct_all) <- "Correct"
gf_correct_all <- rownames_to_column(gf_correct_all, var = "Item")
gf_correct_all <- gf_correct_all[order(-gf_correct_all$Correct),]
ggplot(data = gf_correct_all, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Explore data for underage
psych::describe(select(gf_data_underage, age, gfscore))

ggplot(data = gf_data_underage, mapping = aes(x = age, y = gfscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gf_data_underage) + 
  geom_boxplot(mapping = aes(x = phase, y = gfscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gf_data_underage$gfscore),5)), limits = c(0,max(gf_data_underage$gfscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gf_data_underage) + 
  geom_boxplot(mapping = aes(x = nationality, y = gfscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gf_data_underage$gfscore),5)), limits = c(0,max(gf_data_underage$gfscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gf_data_underage) + 
  geom_boxplot(mapping = aes(x = device, y = gfscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gf_data_underage$gfscore),5)), limits = c(0,max(gf_data_underage$gfscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gf_correct_underage <- as.data.frame(sapply(gf_data_underage[gf_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gf_correct_underage) <- str_remove(row.names(gf_correct_underage), "gf")
names(gf_correct_underage) <- "Correct"
gf_correct_underage <- rownames_to_column(gf_correct_underage, var = "Item")
gf_correct_underage <- gf_correct_underage[order(-gf_correct_underage$Correct),]
ggplot(data = gf_correct_underage, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Generate MIRT object using ICS Parameters
gf_mirt_ics <- generate.mirt_object(parameters = gf_pars, itemtype = "Rasch")
plot(gf_mirt_ics, type = "infoSE")
coef(gf_mirt_ics, simplify = TRUE)$item

### -- Simulation for ICS parameters
gf_sim_ics <- list()
for (i in 1:length(sem)) { gf_sim_ics[[i]] <- cat_simulation(chc_mo = gf_mirt_ics, sem = sem[i]) }
gf_sim_ics <- bind_rows(gf_sim_ics)
gf_sim_ics$reliability <- reliability
gf_sim_ics
write.table(x = gf_sim_ics, file = paste0("sim_tables/","gf_sim_ics.csv"), sep = ",", row.names = F)

### --- Generate MIRT object using all age data
gf_mirt_all <- mirt(data = gf_data_all[gf_items], model = 1, itemtype = "Rasch")
plot(gf_mirt_all, type = "infoSE")
coef(gf_mirt_all, simplify = TRUE)$item

### -- Simulation using all age data
gf_sim_all <- list()
for (i in 1:length(sem)) { gf_sim_all[[i]] <- cat_simulation(chc_mo = gf_mirt_all, sem = sem[i]) }
gf_sim_all <- bind_rows(gf_sim_all)
gf_sim_all$reliability <- reliability
gf_sim_all
write.table(x = gf_sim_all, file = paste0("sim_tables/","gf_sim_all.csv"), sep = ",", row.names = F)

### --- Generate MIRT object for underage participants, no multiple imputation required
gf_mirt_underage <- mirt(data = gf_data_underage[gf_items], model = 1, itemtype = "Rasch")
plot(gf_mirt_underage, type = "infoSE")
coef(gf_mirt_underage, simplify = TRUE)$item

### -- Simulation for underage data
gf_sim_underage <- list()
for (i in 1:length(sem)) { gf_sim_underage[[i]] <- cat_simulation(chc_mo = gf_mirt_underage, sem = sem[i]) }
gf_sim_underage <- bind_rows(gf_sim_underage)
gf_sim_underage$reliability <- reliability
gf_sim_underage
write.table(x = gf_sim_underage, file = paste0("sim_tables/","gf_sim_underage.csv"), sep = ",", row.names = F)

### --- Performance of participants by group
gf_Theta <- matrix(rnorm(n = 5000, mean = 0, sd = 1)) # Simulate theta for 5000 people
gf_pattern <- generate_pattern(mo = gf_mirt_underage, Theta = gf_Theta) ### Generate response pattern
gf_cat_design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
gf_cat <- mirtCAT(mo = gf_mirt_underage, local_pattern = gf_pattern, start_item = "MI", method = "EAP", criteria = "MI", design = gf_cat_design, cl = cl) # run
stopCluster(cl) # stop the cluster
gf_est_Theta <- laply(gf_cat, function(y) y$thetas) # list of theta
gf_sim_groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(gf_Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(gf_Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(gf_Theta > Theta_groups[i-1] & gf_Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(gf_cat[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(gf_Theta[rows,1], gf_est_Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(gf_Theta[rows,1] - gf_est_Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((gf_Theta[rows,1] - gf_est_Theta[rows])^2))
  nParticipants <- length(rows)
  gf_sim_groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                  "averageItems" = average_items, "correlation" = correlation, 
                                  "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
gf_sim_groups <- bind_rows(gf_sim_groups)
gf_sim_groups
write.table(x = gf_sim_groups, file = paste0("sim_tables/","gf_sim_groups.csv"), sep = ",", row.names = F)

#---                                                                           ---#
################################## Visualisation  ##################################
#---                                                                           ---#

### --- Get Item List
gv_items <- scan(file = "gvitems.txt", character(), quote = "") # retrieve retained items from ICS
gv_items <- str_remove(gv_items, "score_") # standardise item names consistent with item parameters file

### --- Get ICS Parameters
gv_pars <- read.csv(file = "gv_rasch_parameters.csv") ### Parameters from ICS
rownames(gv_pars) <- gv_pars[, 1] ### put items as row names
gv_pars <- gv_pars[, -1] ### remove item names from df

### --- Clean Data
gv_data_all <- select(filter(data, gvscore > 0), phase, age, gender, nationality, device, all_of(gv_items)) ### raw data for all participants 
gv_data_all <- gv_data_all[c(which(apply(gv_data_all, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gv_data_all <- mutate(gv_data_all, gvscore = rowSums(select(gv_data_all, starts_with("gv")))) ### calculate total score

gv_data_underage <- select(filter(data, gvscore > 0, age < 19), phase, age, gender, nationality, device, all_of(gv_items)) ### raw data for underage participants
gv_data_underage <- gv_data_underage[c(which(apply(gv_data_underage, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gv_data_underage <- mutate(gv_data_underage, gvscore = rowSums(select(gv_data_underage, starts_with("gv")))) ### calculate total score

### --- Explore Data for all data
psych::describe(select(gv_data_all, age, gvscore))

ggplot(data = gv_data_all, mapping = aes(x = age, y = gvscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gv_data_all) + 
  geom_boxplot(mapping = aes(x = phase, y = gvscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gv_data_all$gvscore),5)), limits = c(0,max(gv_data_all$gvscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gv_data_all) + 
  geom_boxplot(mapping = aes(x = nationality, y = gvscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gv_data_all$gvscore),5)), limits = c(0,max(gv_data_all$gvscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gv_data_all) + 
  geom_boxplot(mapping = aes(x = device, y = gvscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gv_data_all$gvscore),5)), limits = c(0,max(gv_data_all$gvscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gv_correct_all <- as.data.frame(sapply(gv_data_all[gv_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gv_correct_all) <- str_remove(row.names(gv_correct_all), "gv")
names(gv_correct_all) <- "Correct"
gv_correct_all <- rownames_to_column(gv_correct_all, var = "Item")
gv_correct_all <- gv_correct_all[order(-gv_correct_all$Correct),]
ggplot(data = gv_correct_all, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Explore data for underage
psych::describe(select(gv_data_underage, age, gvscore))

ggplot(data = gv_data_underage, mapping = aes(x = age, y = gvscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gv_data_underage) + 
  geom_boxplot(mapping = aes(x = phase, y = gvscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gv_data_underage$gvscore),5)), limits = c(0,max(gv_data_underage$gvscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gv_data_underage) + 
  geom_boxplot(mapping = aes(x = nationality, y = gvscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gv_data_underage$gvscore),5)), limits = c(0,max(gv_data_underage$gvscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gv_data_underage) + 
  geom_boxplot(mapping = aes(x = device, y = gvscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gv_data_underage$gvscore),5)), limits = c(0,max(gv_data_underage$gvscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gv_correct_underage <- as.data.frame(sapply(gv_data_underage[gv_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gv_correct_underage) <- str_remove(row.names(gv_correct_underage), "gv")
names(gv_correct_underage) <- "Correct"
gv_correct_underage <- rownames_to_column(gv_correct_underage, var = "Item")
gv_correct_underage <- gv_correct_underage[order(-gv_correct_underage$Correct),]
ggplot(data = gv_correct_underage, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Generate MIRT object using ICS Parameters
gv_mirt_ics <- generate.mirt_object(parameters = gv_pars, itemtype = "Rasch")
plot(gv_mirt_ics, type = "infoSE")
coef(gv_mirt_ics, simplify = TRUE)$item

### -- Simulation for ICS parameters
gv_sim_ics <- list()
for (i in 1:length(sem)) { gv_sim_ics[[i]] <- cat_simulation(chc_mo = gv_mirt_ics, sem = sem[i]) }
gv_sim_ics <- bind_rows(gv_sim_ics)
gv_sim_ics$reliability <- reliability
gv_sim_ics
write.table(x = gv_sim_ics, file = paste0("sim_tables/","gv_sim_ics.csv"), sep = ",", row.names = F)

### --- Generate MIRT object using all age data
gv_mirt_all <- mirt(data = gv_data_all[gv_items], model = 1, itemtype = "Rasch")
plot(gv_mirt_all, type = "infoSE")
coef(gv_mirt_all, simplify = TRUE)$item

### -- Simulation using all age data
gv_sim_all <- list()
for (i in 1:length(sem)) { gv_sim_all[[i]] <- cat_simulation(chc_mo = gv_mirt_all, sem = sem[i]) }
gv_sim_all <- bind_rows(gv_sim_all)
gv_sim_all$reliability <- reliability
gv_sim_all
write.table(x = gv_sim_all, file = paste0("sim_tables/","gv_sim_all.csv"), sep = ",", row.names = F)

### --- Generate MIRT object for underage participants, no multiple imputation required
gv_mirt_underage <- mirt(data = gv_data_underage[gv_items], model = 1, itemtype = "Rasch")
plot(gv_mirt_underage, type = "infoSE")
coef(gv_mirt_underage, simplify = TRUE)$item

### -- Simulation for underage data
gv_sim_underage <- list()
for (i in 1:length(sem)) { gv_sim_underage[[i]] <- cat_simulation(chc_mo = gv_mirt_underage, sem = sem[i]) }
gv_sim_underage <- bind_rows(gv_sim_underage)
gv_sim_underage$reliability <- reliability
gv_sim_underage
write.table(x = gv_sim_underage, file = paste0("sim_tables/","gv_sim_underage.csv"), sep = ",", row.names = F)

### --- Performance of participants by group
gv_Theta <- matrix(rnorm(n = 5000, mean = 0, sd = 1)) # Simulate theta for 5000 people
gv_pattern <- generate_pattern(mo = gv_mirt_underage, Theta = gv_Theta) ### Generate response pattern
gv_cat_design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
gv_cat <- mirtCAT(mo = gv_mirt_underage, local_pattern = gv_pattern, start_item = "MI", method = "EAP", criteria = "MI", design = gv_cat_design, cl = cl) # run
stopCluster(cl) # stop the cluster
gv_est_Theta <- laply(gv_cat, function(y) y$thetas) # list of theta
gv_sim_groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(gv_Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(gv_Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(gv_Theta > Theta_groups[i-1] & gv_Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(gv_cat[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(gv_Theta[rows,1], gv_est_Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(gv_Theta[rows,1] - gv_est_Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((gv_Theta[rows,1] - gv_est_Theta[rows])^2))
  nParticipants <- length(rows)
  gv_sim_groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                  "averageItems" = average_items, "correlation" = correlation, 
                                  "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
gv_sim_groups <- bind_rows(gv_sim_groups)
gv_sim_groups
write.table(x = gv_sim_groups, file = paste0("sim_tables/","gv_sim_groups.csv"), sep = ",", row.names = F)

#---                                                                             ---#
################################## Working Memory  ##################################
#---                                                                             ---#

### --- Get Item List
gwm_items <- scan(file = "gwmitems.txt", character(), quote = "") # retrieve retained items from ICS
gwm_items <- str_remove(gwm_items, "score_") # standardise item names consistent with item parameters file

### --- Get ICS Parameters
gwm_pars <- read.csv(file = "gwm_rasch_parameters.csv") ### Parameters from ICS
rownames(gwm_pars) <- gwm_pars[, 1] ### put items as row names
gwm_pars <- gwm_pars[, -1] ### remove item names from df

### --- Clean Data
gwm_data_all <- select(filter(data, gwmscore > 0), phase, age, gender, nationality, device, all_of(gwm_items)) ### raw data for all participants 
gwm_data_all <- gwm_data_all[c(which(apply(gwm_data_all, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gwm_data_all <- mutate(gwm_data_all, gwmscore = rowSums(select(gwm_data_all, starts_with("gwm")))) ### calculate total score

gwm_data_underage <- select(filter(data, gwmscore > 0, age < 19), phase, age, gender, nationality, device, all_of(gwm_items)) ### raw data for underage participants
gwm_data_underage <- gwm_data_underage[c(which(apply(gwm_data_underage, 1, percentage_missing) == 0)), ] ### keep participants that attempted all retained items
gwm_data_underage <- mutate(gwm_data_underage, gwmscore = rowSums(select(gwm_data_underage, starts_with("gwm")))) ### calculate total score

### --- Explore Data for all data
psych::describe(select(gwm_data_all, age, gwmscore))

ggplot(data = gwm_data_all, mapping = aes(x = age, y = gwmscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gwm_data_all) + 
  geom_boxplot(mapping = aes(x = phase, y = gwmscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gwm_data_all$gwmscore),5)), limits = c(0,max(gwm_data_all$gwmscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gwm_data_all) + 
  geom_boxplot(mapping = aes(x = nationality, y = gwmscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gwm_data_all$gwmscore),5)), limits = c(0,max(gwm_data_all$gwmscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gwm_data_all) + 
  geom_boxplot(mapping = aes(x = device, y = gwmscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gwm_data_all$gwmscore),5)), limits = c(0,max(gwm_data_all$gwmscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gwm_correct_all <- as.data.frame(sapply(gwm_data_all[gwm_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gwm_correct_all) <- str_remove(row.names(gwm_correct_all), "gwm")
names(gwm_correct_all) <- "Correct"
gwm_correct_all <- rownames_to_column(gwm_correct_all, var = "Item")
gwm_correct_all <- gwm_correct_all[order(-gwm_correct_all$Correct),]
ggplot(data = gwm_correct_all, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Explore data for underage
psych::describe(select(gwm_data_underage, age, gwmscore))

ggplot(data = gwm_data_underage, mapping = aes(x = age, y = gwmscore)) +
  geom_point() +
  geom_smooth(method = "loess")

ggplot(data = gwm_data_underage) + 
  geom_boxplot(mapping = aes(x = phase, y = gwmscore, fill = phase), width = .25) +
  scale_x_discrete(name = "Phase") +
  scale_y_continuous(breaks = c(seq(0,max(gwm_data_underage$gwmscore),5)), limits = c(0,max(gwm_data_underage$gwmscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gwm_data_underage) + 
  geom_boxplot(mapping = aes(x = nationality, y = gwmscore, fill = nationality), width = .25) +
  scale_x_discrete(name = "Nationality") +
  scale_y_continuous(breaks = c(seq(0,max(gwm_data_underage$gwmscore),5)), limits = c(0,max(gwm_data_underage$gwmscore)), name = "Total Score") +
  theme(legend.position = "none")

ggplot(data = gwm_data_underage) + 
  geom_boxplot(mapping = aes(x = device, y = gwmscore, fill = device), width = .25) +
  scale_x_discrete(name = "Device") +
  scale_y_continuous(breaks = c(seq(0,max(gwm_data_underage$gwmscore),5)), limits = c(0,max(gwm_data_underage$gwmscore)), name = "Total Score") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

gwm_correct_underage <- as.data.frame(sapply(gwm_data_underage[gwm_items], function(x) { (1 - (table(x)[[1]]/length(x))) * 100 }))
row.names(gwm_correct_underage) <- str_remove(row.names(gwm_correct_underage), "gwm")
names(gwm_correct_underage) <- "Correct"
gwm_correct_underage <- rownames_to_column(gwm_correct_underage, var = "Item")
gwm_correct_underage <- gwm_correct_underage[order(-gwm_correct_underage$Correct),]
ggplot(data = gwm_correct_underage, aes(x = factor(Item, levels = Item), y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  scale_x_discrete(name = "Item Number") + 
  scale_y_continuous(name = "% Correct", breaks = seq(0,100,10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

### --- Generate MIRT object using ICS Parameters
gwm_mirt_ics <- generate.mirt_object(parameters = gwm_pars, itemtype = "Rasch")
plot(gwm_mirt_ics, type = "infoSE")
coef(gwm_mirt_ics, simplify = TRUE)$item

### -- Simulation for ICS parameters
gwm_sim_ics <- list()
for (i in 1:length(sem)) { gwm_sim_ics[[i]] <- cat_simulation(chc_mo = gwm_mirt_ics, sem = sem[i]) }
gwm_sim_ics <- bind_rows(gwm_sim_ics)
gwm_sim_ics$reliability <- reliability
gwm_sim_ics
write.table(x = gwm_sim_ics, file = paste0("sim_tables/","gwm_sim_ics.csv"), sep = ",", row.names = F)

### --- Generate MIRT object using all age data
gwm_mirt_all <- mirt(data = gwm_data_all[gwm_items], model = 1, itemtype = "Rasch")
plot(gwm_mirt_all, type = "infoSE")
coef(gwm_mirt_all, simplify = TRUE)$item

### -- Simulation using all age data
gwm_sim_all <- list()
for (i in 1:length(sem)) { gwm_sim_all[[i]] <- cat_simulation(chc_mo = gwm_mirt_all, sem = sem[i]) }
gwm_sim_all <- bind_rows(gwm_sim_all)
gwm_sim_all$reliability <- reliability
gwm_sim_all
write.table(x = gwm_sim_all, file = paste0("sim_tables/","gwm_sim_all.csv"), sep = ",", row.names = F)

### --- Generate MIRT object for underage participants, no multiple imputation required
gwm_mirt_underage <- mirt(data = gwm_data_underage[gwm_items], model = 1, itemtype = "Rasch")
plot(gwm_mirt_underage, type = "infoSE")
coef(gwm_mirt_underage, simplify = TRUE)$item

### -- Simulation for underage data
gwm_sim_underage <- list()
for (i in 1:length(sem)) { gwm_sim_underage[[i]] <- cat_simulation(chc_mo = gwm_mirt_underage, sem = sem[i]) }
gwm_sim_underage <- bind_rows(gwm_sim_underage)
gwm_sim_underage$reliability <- reliability
gwm_sim_underage
write.table(x = gwm_sim_underage, file = paste0("sim_tables/","gwm_sim_underage.csv"), sep = ",", row.names = F)

### --- Performance of participants by group
gwm_Theta <- matrix(rnorm(n = 5000, mean = 0, sd = 1)) # Simulate theta for 5000 people
gwm_pattern <- generate_pattern(mo = gwm_mirt_underage, Theta = gwm_Theta) ### Generate response pattern
gwm_cat_design <- list(min_SEM = sem[3]) ### the design factors of the simulation
cl <- makeCluster(detectCores()) ### for parallel processing
gwm_cat <- mirtCAT(mo = gwm_mirt_underage, local_pattern = gwm_pattern, start_item = "MI", method = "EAP", criteria = "MI", design = gwm_cat_design, cl = cl) # run
stopCluster(cl) # stop the cluster
gwm_est_Theta <- laply(gwm_cat, function(y) y$thetas) # list of theta
gwm_sim_groups <- list()
for (i in 1:(length(Theta_groups)+1)) { 
  if (i == 1) {
    rows <- which(gwm_Theta < Theta_groups[i])
    group <- paste0("Theta < ", Theta_groups[i])
  } else if (i == 12) {
    rows <- which(gwm_Theta > Theta_groups[i-1])
    group <- paste0("Theta > ", Theta_groups[i-1])
  } else {
    rows <- which(gwm_Theta > Theta_groups[i-1] & gwm_Theta < Theta_groups[i])
    group <- paste0(Theta_groups[i-1], " > Theta > ", Theta_groups[i])
  }
  average_items <- mean(laply(gwm_cat[rows], function(y) length(y$items_answered))) ### average number of items
  correlation <- cor(gwm_Theta[rows,1], gwm_est_Theta[rows]) ### correlation between true theta and estimated theta
  bias <- mean(gwm_Theta[rows,1] - gwm_est_Theta[rows]) ### amount of bias between true theta and estimated beta
  rmsd <- sqrt(mean((gwm_Theta[rows,1] - gwm_est_Theta[rows])^2))
  nParticipants <- length(rows)
  gwm_sim_groups[[i]]<- data.frame("group" = group, "nParticipants" = nParticipants, 
                                  "averageItems" = average_items, "correlation" = correlation, 
                                  "bias" = bias, "rmsd" = rmsd)
  rm(average_items,correlation,bias,rmsd,group,nParticipants)
}
gwm_sim_groups <- bind_rows(gwm_sim_groups)
gwm_sim_groups
write.table(x = gwm_sim_groups, file = paste0("sim_tables/","gwm_sim_groups.csv"), sep = ",", row.names = F)
