# Data Analysis for the Convergent Validity study of the CHC-CAT Test
# Original code developed by Mr Jake Kraksa (Monash University)
# Built on R 3.6.3 with R Studio 1.2.5033
# set the working directory to source file directory manually
# restart R session prior to beginning

#---                                                                      ---#
################################## Libraries ##################################
#---                                                                      ---#

# Data exploration

library(psych) # version 1.9.12

# Missing Data

library(mice) # version 3.8.0
library(VIM) # version 5.1.1
library(BaylorEdPsych) # version 0.5

# Visualise Results

library(ggplot2) # version 3.3.0
library(knitr) # version 1.28
library(latticeExtra)

# Data cleaning

library(plyr) # version 1.8.6
library(dplyr) # version 0.8.5
library(stringr) # version 1.4.0
library(tidyr) # version 1.0.2
library(tibble) # version 2.1.3

#---                                                                      ---#
################################## Set Options ##################################
#---                                                                      ---#

### Set Seed for replication
set.seed(145)

### Set print options
options(max.print = 10000)

#---                                                                                      ---#
################################## Get Items and Load Data ##################################
#---                                                                                      ---#

gc_items <- scan(file = "gcitems.txt", character(), quote = "") # retrieve retained items from ICS
gf_items <- scan(file = "gfitems.txt", character(), quote = "") # retrieve retained items from ICS
gv_items <- scan(file = "gvitems.txt", character(), quote = "") # retrieve retained items from ICS
gwm_items <- scan(file = "gwmitems.txt", character(), quote = "") # retrieve retained items from ICS

wisc_gc_raw <- c("wisc.si.raw","wisc.vc.raw")
wisc_gc_ss <- c("wisc.si.ss","wisc.vc.ss")
wisc_gc_comp <- "wisc.vci.comp"

wisc_gf_raw <- c("wisc.mr.raw","wisc.fw.raw")
wisc_gf_ss <- c("wisc.mr.ss","wisc.fw.ss")
wisc_gf_comp <- "wisc.fri.comp"

wisc_gv_raw <- c("wisc.bd.raw","wisc.vp.raw")
wisc_gv_ss <- c("wisc.bd.ss","wisc.vp.ss")
wisc_gv_comp <- "wisc.vsi.comp"

wisc_gwm_raw <- c("wisc.ds.raw","wisc.ps.raw")
wisc_gwm_ss <- c("wisc.ds.ss","wisc.ps.ss")
wisc_gwm_comp <- "wisc.wmi.comp"

### --- Load main data for first correlation between raw scores
data <- read.csv(file = "simulation_data.csv", stringsAsFactors = FALSE)
data$phase <- recode(data$phase, itos = "ITOS", p2a = "ICS-A", p2u = "ICS-U")
data$phase <- as.factor(data$phase)
data$gender <- as.factor(data$gender)
data$nationality <- as.factor(data$nationality)
data$device <- as.factor(data$device)
data <- select(data, id, phase, age, gender, nationality, device, all_of(gc_items), all_of(gf_items), all_of(gv_items), all_of(gwm_items),
               starts_with("wisc"))
data <- mutate(data,
               gcscore = rowSums(select(data, starts_with("score_gc")), na.rm = TRUE),
               gfscore = rowSums(select(data, starts_with("score_gf")), na.rm = TRUE),
               gvscore = rowSums(select(data, starts_with("score_gv")), na.rm = TRUE),
               gwmscore = rowSums(select(data, starts_with("score_gwm")), na.rm = TRUE))
data <- select(data, -starts_with("score_"))

#---                                                                                     ---#
################################## Participant Demographics ##################################
#---                                                                                      ---#

demographics <- filter(data, wisc.fsiq.comp > 0)
psych::describe(demographics, na.rm = FALSE)

table(demographics$gender)
table(demographics$age)
table(demographics$nationality)
table(demographics$device)

ggplot(demographics, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,30,5), labels = paste(seq(0,30,5)), limits = c(0,30)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

#---                                                                                          ---#
################################## Lexical Knowledge Raw Score ##################################
#---                                                                                          ---#

### --- Load Gc Data for first correlation between raw scores
gc_raw_data <- filter(data, gcscore > 0, wisc.fsiq.comp > 0)

### --- Explore raw data

psych::describe(gc_raw_data)

table(gc_raw_data$gender)
table(gc_raw_data$age)
table(gc_raw_data$nationality)
table(gc_raw_data$device)

ggplot(gc_raw_data, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,20,5), labels = paste(seq(0,20,5)), limits = c(0,20)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

### --- Conduct correlations

si.cor <- paste0("r = ",round(cor(gc_raw_data$wisc.si.raw,gc_raw_data$gcscore),2))

ggplot(data = gc_raw_data, mapping = aes(x = gcscore, y = wisc.si.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Lexical Knowledge Raw Score") +
  scale_y_continuous(name = "WISC-V Similarities Raw Score") +
  annotate(geom = "text", x = 10, y = 40, label = si.cor)

vc.cor <- paste0("r = ",round(cor(gc_raw_data$wisc.vc.raw,gc_raw_data$gcscore),2))

ggplot(data = gc_raw_data, mapping = aes(x = gcscore, y = wisc.vc.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Lexical Knowledge Raw Score") +
  scale_y_continuous(name = "WISC-V Vocabulary Raw Score") +
  annotate(geom = "text", x = 10, y = 40, label = vc.cor)

### --- Conduct raw score analysis

round(cor(select(gc_raw_data, age, gcscore, ends_with("raw"))), 2)[,1:2]

#---                                                                                             ---#
################################## Lexical Knowledge Scaled Score ##################################
#---                                                                                              ---#

### --- Load ICS data for second correlation between fscore and scaled scores

gc_ics_data <- read.csv(file = "gc_fscores.csv", stringsAsFactors = FALSE) # read the data
gc_ics_data <- filter(gc_ics_data, !is.na(fscores)) # remove rows that do not have a fscore
gc_ics_data$gc.fscores <- gc_ics_data$fscores # copy fscores into a gc specific named column

gc_ics_data <- select(gc_ics_data, id, age, all_of(wisc_gc_ss), all_of(wisc_gc_comp), gc.fscores) # select only necessary columns

gc_mean_fscore <- round(aggregate(gc_ics_data$gc.fscores,list(gc_ics_data$age), mean),2) # calculate the mean fscore by age
gc_sd_fscore <- round(aggregate(gc_ics_data$gc.fscores,list(gc_ics_data$age), sd),2) # calculate the sd of fscore by age

gc_ics_data$gc.fscore.age.mean <- NA # create new variable ready for input
gc_ics_data$gc.fscore.age.sd <- NA # create new variable ready for input
for (i in 6:16) gc_ics_data$gc.fscore.age.mean[which(gc_ics_data$age == i)] <- gc_mean_fscore[which(gc_mean_fscore$Group.1 == i),2]
for (i in 6:16) gc_ics_data$gc.fscore.age.sd[which(gc_ics_data$age == i)] <- gc_sd_fscore[which(gc_sd_fscore$Group.1 == i),2]
gc_ics_data <- mutate(gc_ics_data, gc.fscore.z.score = (gc.fscores - gc.fscore.age.mean)/gc.fscore.age.sd)

### --- Explore ICS data

si.ss.cor <- paste0("r = ",round(cor(gc_ics_data$wisc.si.ss,gc_ics_data$gc.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gc_ics_data, mapping = aes(x = gc.fscore.z.score, y = wisc.si.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Lexical Knowledge IRT Theta Z Score") +
  scale_y_continuous(name = "WISC-V Similarities Scaled Score", limits = c(0,20)) +
  annotate(geom = "text", x = -1.5, y = 14, label = si.ss.cor)

vc.ss.cor <- paste0("r = ",round(cor(gc_ics_data$wisc.vc.ss,gc_ics_data$gc.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gc_ics_data, mapping = aes(x = gc.fscore.z.score, y = wisc.vc.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Lexical Knowledge IRT Theta Z Score") +
  scale_y_continuous(name = "WISC-V Vocabulary Scaled Score", limits = c(0,20)) +
  annotate(geom = "text", x = -1.5, y = 15, label = vc.ss.cor)

vci.comp.cor <- paste0("r = ",round(cor(gc_ics_data$wisc.vci.comp,gc_ics_data$gc.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gc_ics_data, mapping = aes(x = gc.fscore.z.score, y = wisc.vci.comp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Lexical Knowledge IRT Theta Z Score") +
  scale_y_continuous(name = "WISC-V Verbal Comprehension Index", limits = c(60,140)) +
  annotate(geom = "text", x = -1.5, y = 120, label = vci.comp.cor)

### --- Conduct scaled score analysis

cor(select(gc_ics_data, wisc.si.ss,wisc.vc.ss,gc.fscore.z.score), use = "complete.obs")

#---                                                                                  ---#
################################## Induction Raw Score ##################################
#---                                                                                  ---#

### --- Load Gc Data for first correlation between raw scores
gf_raw_data <- filter(data, gfscore > 0, wisc.fsiq.comp > 0)

### --- Explore raw data

psych::describe(gf_raw_data)

table(gf_raw_data$gender)
table(gf_raw_data$age)
table(gf_raw_data$nationality)
table(gf_raw_data$device)

ggplot(gf_raw_data, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,25,5), labels = paste(seq(0,25,5)), limits = c(0,25)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

### --- Conduct correlations

mr.cor <- paste0("r = ",round(cor(gf_raw_data$wisc.mr.raw,gf_raw_data$gfscore),2))

ggplot(data = gf_raw_data, mapping = aes(x = gfscore, y = wisc.mr.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Induction Raw Score") +
  scale_y_continuous(name = "WISC-V Matrix Reasoning Raw Score", limits = c(0,30)) +
  annotate(geom = "text", x = 10, y = 10, label = mr.cor)

fw.cor <- paste0("r = ",round(cor(gf_raw_data$wisc.fw.raw,gf_raw_data$gfscore),2))

ggplot(data = gf_raw_data, mapping = aes(x = gfscore, y = wisc.fw.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Induction Raw Score") +
  scale_y_continuous(name = "WISC-V Figure Weights Raw Score", limits = c(0,30)) +
  annotate(geom = "text", x = 10, y = 10, label = fw.cor)

### --- Conduct raw score analysis

round(cor(select(gf_raw_data, age, gfscore, ends_with("raw"))), 2)[,1:2]

#---                                                                                     ---#
################################## Induction Scaled Score ##################################
#---                                                                                      ---#

### --- Load ICS data for second correlation between fscore and scaled scores

gf_ics_data <- read.csv(file = "gf_fscores.csv", stringsAsFactors = FALSE) # read the data
gf_ics_data <- filter(gf_ics_data, !is.na(fscores)) # remove rows that do not have a fscore
gf_ics_data$gf.fscores <- gf_ics_data$fscores # copy fscores into a gf specific named column

gf_ics_data <- select(gf_ics_data, id, age, all_of(wisc_gf_ss), all_of(wisc_gf_comp), gf.fscores) # select only necessary columns

gf_mean_fscore <- round(aggregate(gf_ics_data$gf.fscores,list(gf_ics_data$age), mean),2) # calculate the mean fscore by age
gf_sd_fscore <- round(aggregate(gf_ics_data$gf.fscores,list(gf_ics_data$age), sd),2) # calculate the sd of fscore by age

gf_ics_data$gf.fscore.age.mean <- NA # create new variable ready for input
gf_ics_data$gf.fscore.age.sd <- NA # create new variable ready for input
for (i in 6:16) gf_ics_data$gf.fscore.age.mean[which(gf_ics_data$age == i)] <- gf_mean_fscore[which(gf_mean_fscore$Group.1 == i),2]
for (i in 6:16) gf_ics_data$gf.fscore.age.sd[which(gf_ics_data$age == i)] <- gf_sd_fscore[which(gf_sd_fscore$Group.1 == i),2]
gf_ics_data <- mutate(gf_ics_data, gf.fscore.z.score = (gf.fscores - gf.fscore.age.mean)/gf.fscore.age.sd)

### --- Explore ICS data

mr.ss.cor <- paste0("r = ",round(cor(gf_ics_data$wisc.mr.ss,gf_ics_data$gf.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gf_ics_data, mapping = aes(x = gf.fscore.z.score, y = wisc.mr.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Induction IRT Theta Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Matrix Reasoning Scaled Score", limits = c(0,20)) +
  annotate(geom = "text", x = -1.5, y = 15, label = mr.ss.cor)

fw.ss.cor <- paste0("r = ",round(cor(gf_ics_data$wisc.fw.ss,gf_ics_data$gf.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gf_ics_data, mapping = aes(x = gf.fscore.z.score, y = wisc.fw.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Induction IRT Theta Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Figure Weights Scaled Score", limits = c(0,20)) +
  annotate(geom = "text", x = -1.5, y = 15, label = fw.ss.cor)

fri.comp.cor <- paste0("r = ",round(cor(gf_ics_data$wisc.fri.comp,gf_ics_data$gf.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gf_ics_data, mapping = aes(x = gf.fscore.z.score, y = wisc.fri.comp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Induction IRT Theta Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Fluid Reasoning Index", limits = c(60,140)) +
  annotate(geom = "text", x = -1.5, y = 120, label = fri.comp.cor)


### --- Conduct scaled score analymrs

cor(select(gf_ics_data, wisc.mr.ss,wisc.fw.ss,gf.fscore.z.score), use = "complete.obs")

#---                                                                                      ---#
################################## Visualisation Raw Score ##################################
#---                                                                                      ---#

### --- Load Gc Data for first correlation between raw scores
gv_raw_data <- filter(data, gvscore > 0, wisc.fsiq.comp > 0)

### --- Explore raw data

psych::describe(gv_raw_data)

table(gv_raw_data$gender)
table(gv_raw_data$age)
table(gv_raw_data$nationality)
table(gv_raw_data$device)

ggplot(gv_raw_data, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,20,5), labels = paste(seq(0,20,5)), limits = c(0,20)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

### --- Conduct correlations

bd.cor <- paste0("r = ",round(cor(gv_raw_data$wisc.bd.raw,gv_raw_data$gvscore),2))

ggplot(data = gv_raw_data, mapping = aes(x = gvscore, y = wisc.bd.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Visualisation Raw Score") +
  scale_y_continuous(name = "WISC-V Block Design Raw Score") +
  annotate(geom = "text", x = 22.5, y = 10, label = bd.cor)

vp.cor <- paste0("r = ",round(cor(gv_raw_data$wisc.vp.raw,gv_raw_data$gvscore),2))

ggplot(data = gv_raw_data, mapping = aes(x = gvscore, y = wisc.vp.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Visualisation Raw Score") +
  scale_y_continuous(name = "WISC-V Visual Puzzles Raw Score") +
  annotate(geom = "text", x = 22.5, y = 7.5, label = vp.cor)

### --- Conduct raw score analysis

round(cor(select(gv_raw_data, age, gvscore, ends_with("raw"))), 2)[,1:2]

#---                                                                                          ---#
################################## Visualisation Scaled Score ##################################
#---                                                                                          ---#

### --- Load ICS data for second correlation between fscore and scaled scores

gv_ics_data <- read.csv(file = "gv_fscores.csv", stringsAsFactors = FALSE) # read the data
gv_ics_data <- filter(gv_ics_data, !is.na(fscores)) # remove rows that do not have a fscore
gv_ics_data$gv.fscores <- gv_ics_data$fscores # copy fscores into a gv specific named column

gv_ics_data <- select(gv_ics_data, id, age, all_of(wisc_gv_ss), all_of(wisc_gv_comp), gv.fscores) # select only necessary columns

gv_mean_fscore <- round(aggregate(gv_ics_data$gv.fscores,list(gv_ics_data$age), mean),2) # calculate the mean fscore by age
gv_sd_fscore <- round(aggregate(gv_ics_data$gv.fscores,list(gv_ics_data$age), sd),2) # calculate the sd of fscore by age

gv_ics_data$gv.fscore.age.mean <- NA # create new variable ready for input
gv_ics_data$gv.fscore.age.sd <- NA # create new variable ready for input
for (i in 6:16) gv_ics_data$gv.fscore.age.mean[which(gv_ics_data$age == i)] <- gv_mean_fscore[which(gv_mean_fscore$Group.1 == i),2]
for (i in 6:16) gv_ics_data$gv.fscore.age.sd[which(gv_ics_data$age == i)] <- gv_sd_fscore[which(gv_sd_fscore$Group.1 == i),2]
gv_ics_data <- mutate(gv_ics_data, gv.fscore.z.score = (gv.fscores - gv.fscore.age.mean)/gv.fscore.age.sd)

### --- Explore ICS data

bd.ss.cor <- paste0("r = ",round(cor(gv_ics_data$wisc.bd.ss,gv_ics_data$gv.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gv_ics_data, mapping = aes(x = gv.fscore.z.score, y = wisc.bd.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Visualisation IRT Theta Z Score") +
  scale_y_continuous(name = "WISC-V Block Design Scaled Score", limits = c(0,20)) +
  annotate(geom = "text", x = -1.5, y = 1, label = bd.ss.cor)

vp.ss.cor <- paste0("r = ",round(cor(gv_ics_data$wisc.vp.ss,gv_ics_data$gv.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gv_ics_data, mapping = aes(x = gv.fscore.z.score, y = wisc.vp.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Visualisation IRT Theta Z Score") +
  scale_y_continuous(name = "WISC-V Visual Puzzles Scaled Score", limits = c(0,20)) +
  annotate(geom = "text", x = -1.5, y = 1, label = vp.ss.cor)

vsi.comp.cor <- paste0("r = ",round(cor(gv_ics_data$wisc.vsi.comp,gv_ics_data$gv.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gv_ics_data, mapping = aes(x = gv.fscore.z.score, y = wisc.vsi.comp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Visualisation IRT Theta Z Score") +
  scale_y_continuous(name = "WISC-V Visual Spatial Index", limits = c(60,140)) +
  annotate(geom = "text", x = -1.5, y = 120, label = vsi.comp.cor)

### --- Conduct scaled score analybds

cor(select(gv_ics_data, wisc.bd.ss,wisc.vp.ss,gv.fscore.z.score), use = "complete.obs")

#---                                                                                      ---#
################################## Working Memory Raw Score ##################################
#---                                                                                      ---#

### --- Load Gc Data for first correlation between raw scores
gwm_raw_data <- filter(data, gwmscore > 0, wisc.fsiq.comp > 0)

### --- Explore raw data

psych::describe(gwm_raw_data)

table(gwm_raw_data$gender)
table(gwm_raw_data$age)
table(gwm_raw_data$nationality)
table(gwm_raw_data$device)

ggplot(gwm_raw_data, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,20,5), labels = paste(seq(0,20,5)), limits = c(0,20)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

### --- Conduct correlations

ds.cor <- paste0("r = ",round(cor(gwm_raw_data$wisc.ds.raw,gwm_raw_data$gwmscore),2))

ggplot(data = gwm_raw_data, mapping = aes(x = gwmscore, y = wisc.ds.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Working Memory Raw Score") +
  scale_y_continuous(name = "WISC-V Digit Span Raw Score", limits = c(5,45)) +
  annotate(geom = "text", x = 10, y = 10, label = ds.cor)

ps.cor <- paste0("r = ",round(cor(gwm_raw_data$wisc.ps.raw,gwm_raw_data$gwmscore),2))

ggplot(data = gwm_raw_data, mapping = aes(x = gwmscore, y = wisc.ps.raw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Working Memory Raw Score") +
  scale_y_continuous(name = "WISC-V Picture Span Raw Score", limits = c(5,45)) +
  annotate(geom = "text", x = 10, y = 10, label = ps.cor)

### --- Conduct raw score analysis

round(cor(select(gwm_raw_data, age, gwmscore, ends_with("raw"))), 2)[,1:2]

#---                                                                                        ---#
################################## Working Memory Scaled Score ##################################
#---                                                                                        ---#

### --- Load ICS data for second correlation between fscore and scaled scores

gwm_ics_data <- read.csv(file = "gwm_fscores.csv", stringsAsFactors = FALSE) # read the data
gwm_ics_data <- filter(gwm_ics_data, !is.na(fscores), wisc.ds.ss < 20, wisc.ps.ss < 20) # remove rows that do not have a fscore or have invalid ss
gwm_ics_data$gwm.fscores <- gwm_ics_data$fscores # copy fscores into a gwm specific named column

gwm_ics_data <- select(gwm_ics_data, id, age, all_of(wisc_gwm_ss), all_of(wisc_gwm_comp), gwm.fscores) # select only necessary columns

gwm_mean_fscore <- round(aggregate(gwm_ics_data$gwm.fscores,list(gwm_ics_data$age), mean),2) # calculate the mean fscore by age
gwm_sd_fscore <- round(aggregate(gwm_ics_data$gwm.fscores,list(gwm_ics_data$age), sd),2) # calculate the sd of fscore by age

gwm_ics_data$gwm.fscore.age.mean <- NA # create new variable ready for input
gwm_ics_data$gwm.fscore.age.sd <- NA # create new variable ready for input
for (i in 6:16) gwm_ics_data$gwm.fscore.age.mean[which(gwm_ics_data$age == i)] <- gwm_mean_fscore[which(gwm_mean_fscore$Group.1 == i),2]
for (i in 6:16) gwm_ics_data$gwm.fscore.age.sd[which(gwm_ics_data$age == i)] <- gwm_sd_fscore[which(gwm_sd_fscore$Group.1 == i),2]
gwm_ics_data <- mutate(gwm_ics_data, gwm.fscore.z.score = (gwm.fscores - gwm.fscore.age.mean)/gwm.fscore.age.sd)

### --- Explore ICS data

ds.ss.cor <- paste0("r = ",round(cor(gwm_ics_data$wisc.ds.ss, gwm_ics_data$gwm.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gwm_ics_data, mapping = aes(x = gwm.fscore.z.score, y = wisc.ds.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Working Memory IRT Theta Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Digit Span Scaled Score", limits = c(0, 20)) +
  annotate(geom = "text", x = -1.5, y = 1, label = ds.ss.cor)

ps.ss.cor <- paste0("r = ",round(cor(gwm_ics_data$wisc.ps.ss,gwm_ics_data$gwm.fscore.z.score, use = "complete.obs"),2))

ggplot(data = gwm_ics_data, mapping = aes(x = gwm.fscore.z.score, y = wisc.ps.ss)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Working Memory IRT Theta Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Picture Span Scaled Score", limits = c(0, 20)) +
  annotate(geom = "text", x = -1.5, y = 1, label = ps.ss.cor)

psi.comp.cor <- paste0("r = ", round(cor(gwm_ics_data$wisc.wmi.comp, gwm_ics_data$gwm.fscore.z.score, use = "complete.obs"), 2))

ggplot(data = gwm_ics_data, mapping = aes(x = gwm.fscore.z.score, y = wisc.wmi.comp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT Working Memory Theta Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Working Memory Index", limits = c(60,140)) +
  annotate(geom = "text", x = -1.5, y = 110, label = psi.comp.cor)

### --- Conduct scaled score analydss

cor(select(gwm_ics_data, wisc.ds.ss,wisc.ps.ss,gwm.fscore.z.score), use = "complete.obs")

#---                                                                              ---#
################################## Multidimensional ##################################
#---                                                                              ---#

### --- set up data set
multi.data <- left_join(x = data, y = gc_ics_data)
multi.data <- left_join(x = multi.data, y = gf_ics_data)
multi.data <- left_join(x = multi.data, y = gv_ics_data)
multi.data <- left_join(x = multi.data, y = gwm_ics_data)
multi.data <- filter(multi.data, wisc.fsiq.comp > 0)
multi.data <- filter(multi.data, !is.na(gc.fscores) | !is.na(gf.fscores) | !is.na(gv.fscores) | !is.na(gwm.fscores))

### --- Explore raw data

psych::describe(multi.data)

table(multi.data$gender)
table(multi.data$age)
table(multi.data$nationality)
table(multi.data$device)

ggplot(multi.data, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,25,5), labels = paste(seq(0,25,5)), limits = c(0,25)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

### --- missing data analysis

aggr(multi.data[c("gc.fscores","gf.fscores","gv.fscores","gwm.fscores")], combined = TRUE, labels = c("Gc:VL", "Gf:I", "Gv:Vz", "Gwm"))

### --- reduce down to the three abilities Gc, Gf and Gv
multi.data <- filter(multi.data, !is.na(gc.fscores), !is.na(gf.fscores), !is.na(gv.fscores))
multi.data <- mutate(multi.data, sum.of.z.score = gc.fscore.z.score + gf.fscore.z.score + gv.fscore.z.score)
multi.data$sum.of.z.score.mean <- mean(multi.data$sum.of.z.score)
multi.data$sum.of.z.score.sd <- sd(multi.data$sum.of.z.score)
multi.data$g.z.score <- (multi.data$sum.of.z.score - multi.data$sum.of.z.score.mean)/multi.data$sum.of.z.score.sd

### --- Explore raw data

psych::describe(multi.data)

table(multi.data$gender)
table(multi.data$age)
table(multi.data$nationality)
table(multi.data$device)

ggplot(multi.data, aes(x = age)) +
  geom_bar() +
  scale_x_continuous(name = "Age", breaks = seq(6,16,1)) +
  scale_y_continuous(name = "Count", breaks = seq(0,10,2), labels = paste(seq(0,10,2)), limits = c(0,10)) +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1)

### --- Obtain correlations

round(cor(select(multi.data, g.z.score, contains("comp")), use = "complete.obs", method = "pearson"),2)

fsiq.comp.cor <- paste0("r = ", round(cor(multi.data$wisc.fsiq.comp, multi.data$g.z.score, use = "complete.obs"), 2))

ggplot(data = multi.data, mapping = aes(x = g.z.score, y = wisc.fsiq.comp)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(name = "CHC-CAT g Z Score", limits = c(-2,2)) +
  scale_y_continuous(name = "WISC-V Full Scale Index", limits = c(60,140)) +
  annotate(geom = "text", x = -1.5, y = 110, label = fsiq.comp.cor)
