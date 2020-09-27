# Cleaning Data for the Item Calibration Study (ICS) CHC Test
# Original code developed by Mr Jake Kraksa (Monash University)
# Built on R 3.6.2 with R Studio 1.2.5033
# set the working directory to source file directory manually
# restart R session prior to beginning

### Load libraries

library(tidyr) # version 1.0.2
library(dplyr) # version 0.8.4
library(stringr) # version 1.4.0

######## ITOS Cleanup ######### 

### Clean ITOS Contact
itoscontact <- read.csv("itos_contact.csv", stringsAsFactors = FALSE)
itoscontact <- select(itoscontact, -id, -data_set_id)
itoscontact <- pivot_wider(itoscontact, names_from = name, values_from = value)
itoscontact <- select(itoscontact, session_id)

### Clean ITOS Demographics
itosdemo <- read.csv("itos_demographics.csv", stringsAsFactors = FALSE)
itosdemo <- select(itosdemo, -id, -data_set_id)
itosdemo <- pivot_wider(itosdemo, names_from = name, values_from = value)
itosdemo <- select(itosdemo, session_id, age, gender, nationality)

### clean ITOS Gc Responses
itosgc <- read.csv("itos_gc_results.csv", stringsAsFactors = FALSE)
itosgc$item <- paste0(itosgc$trait,itosgc$item_id)
itosgc <- select(itosgc, -id, -data_set_id, -trait, -item_id, -response)
names(itosgc)[names(itosgc)=="correct"] <- "score"
names(itosgc)[names(itosgc)=="time_taken"] <- "timeTaken"
itosgc <- pivot_wider(itosgc, names_from = item, values_from = c(score,timeTaken))

### clean ITOS Gf Responses
itosgf <- read.csv("itos_gf_results.csv", stringsAsFactors = FALSE)
itosgf$item <- paste0("gf",itosgf$item_id)
itosgf <- select(itosgf, -id, -data_set_id, -trait, -item_id, -response)
names(itosgf)[names(itosgf)=="correct"] <- "score"
names(itosgf)[names(itosgf)=="time_taken"] <- "timeTaken"
itosgf <- pivot_wider(itosgf, names_from = item, values_from = c(score,timeTaken))

### clean ITOS Gv Responses
itosgv <- read.csv("itos_gv_results.csv", stringsAsFactors = FALSE)
itosgv$item <- paste0(itosgv$trait,itosgv$item_id)
itosgv <- select(itosgv, -id, -data_set_id, -trait, -item_id, -response)
names(itosgv)[names(itosgv)=="correct"] <- "score"
names(itosgv)[names(itosgv)=="time_taken"] <- "timeTaken"
itosgv <- pivot_wider(itosgv, names_from = item, values_from = c(score,timeTaken))

### clean ITOS Gwm Responses
itosgwm <- read.csv("itos_gwm_results.csv", stringsAsFactors = FALSE)
itosgwm$item <- paste0(itosgwm$trait,itosgwm$item_id)
itosgwm <- select(itosgwm, -id, -data_set_id, -trait, -item_id, -response)
names(itosgwm)[names(itosgwm)=="correct"] <- "score"
names(itosgwm)[names(itosgwm)=="time_taken"] <- "timeTaken"
itosgwm <- pivot_wider(itosgwm, names_from = item, values_from = c(score,timeTaken))

### Merge ITOS
itos <- full_join(itosdemo,itoscontact, by = "session_id")
itos <- full_join(itos,itosgc, by = "session_id")
itos <- full_join(itos,itosgf, by = "session_id")
itos <- full_join(itos,itosgv, by = "session_id")
itos <- full_join(itos,itosgwm, by = "session_id")
names(itos)[names(itos) == "session_id"] <- "id"

######## Phase 2 Underage Cleanup ######### 

### Clean Phase 2 Underage Sessions
p2u.sessions <- read.csv("stageTwoUnderageSessions.csv", stringsAsFactors = FALSE)
p2u.sessions <- select(p2u.sessions, -id, -test_id, -startedTime, -updateTime, -finished)

### Clean Phase 2 Underage Users
p2u.users <- read.csv("stageTwoUnderageUsers.csv", stringsAsFactors = FALSE)
p2u.users <- select(p2u.users, -password, -enabled)

### Clean Phase 2 Underage Responses
p2u.responses <- read.csv("stageTwoUnderageResponses.csv", stringsAsFactors = FALSE)
p2u.responses$trait <- recode(p2u.responses$trait, gf30 = "gf", gf60 = "gf")
p2u.responses$item <- paste0(p2u.responses$trait,p2u.responses$item_id)
p2u.responses <- select(p2u.responses, -id, -item_id, -response, -trait)
p2u.responses <- pivot_wider(p2u.responses, names_from = item, values_from = c(score,timeTaken))
p2u.responses$session_id <- str_remove(p2u.responses$session_id, "[i]")
p2u.responses$session_id <- as.integer(p2u.responses$session_id)

### Clean Phase 2 Underage WISC
p2u.wisc <- read.csv("stageTwoUnderageWISC.csv", stringsAsFactors = FALSE)
p2u.wisc <- select(p2u.wisc, -Date.of.WISC, -Name, -Practitioner)
names(p2u.wisc) <- str_to_lower(names(p2u.wisc))
p2u.wisc <- rename_all(p2u.wisc, function(x) paste0("wisc.",x))
names(p2u.wisc)[names(p2u.wisc)=="wisc.id"] <- "id"
names(p2u.wisc)[names(p2u.wisc)=="wisc.age.at.testing..months."] <- "age"
names(p2u.wisc)[names(p2u.wisc)=="wisc.gender"] <- "gender"
names(p2u.wisc)[names(p2u.wisc)=="wisc.chc.cat.device"] <- "device"

### Merge Phase 2 Underage
p2u <- full_join(p2u.sessions,p2u.users, by = c("user_id" = "id"))
p2u <- full_join(p2u,p2u.responses, by = c("internal_id" = "session_id"))
p2u <- select(p2u, -user_id, -internal_id)
p2u <- filter(p2u, !login == "test")
p2u <- pivot_longer(p2u, cols = -login, 
                         names_to = "variable", values_to = "result", 
                         values_drop_na = TRUE )
p2u <- distinct(p2u, login, variable, .keep_all = TRUE)
p2u <- pivot_wider(p2u, names_from = variable, values_from = result)
p2u <- inner_join(p2u, p2u.wisc, by = c("login" = "id"))
names(p2u)[names(p2u) == "login"] <- "id"

######## Phase 2 Adult Cleanup ######### 

### Clean Phase 2 Age Sessions
p2a.sessions <- read.csv("stageTwoAdultSessions.csv", stringsAsFactors = FALSE)
p2a.sessions <- select(p2a.sessions, -id, -test_id, -startedTime, -updateTime, -finished, -state)

### Clean Phase 2 Adult Users
p2a.users <- read.csv("stageTwoAdultUsers.csv", stringsAsFactors = FALSE) # 2185 users completed the demographic information
p2a.users <- select(p2a.users, -password, -enabled)

### Clean Phase 2 Adult Demographics
p2a.demo <- read.csv("stageTwoAdultDemographics.csv", stringsAsFactors = FALSE) 
p2a.demo <- select(p2a.demo, -id, -colour, -number, -animal, -lastletter, -firstletter, -session)
p2a.demo <- distinct(p2a.demo)

### Clean Phase 2 Adult Responses
p2a.responses <- read.csv("stageTwoAdultAssessmentResponses.csv", stringsAsFactors = FALSE)
p2a.responses$trait <- recode(p2a.responses$trait, gf30 = "gf", gf60 = "gf")
p2a.responses$item <- paste0(p2a.responses$trait,p2a.responses$item_id)
p2a.responses <- select(p2a.responses, -id, -item_id, -response, -trait, -theta, -sem)
p2a.responses <- pivot_wider(p2a.responses, names_from = item, values_from = c(score,timeTaken))
p2a.responses$session_id <- str_remove(p2a.responses$session_id, "[i]")
p2a.responses$session_id <- as.integer(p2a.responses$session_id)

### Merge Phase 2 Adult
p2a <- full_join(p2a.sessions,p2a.users, by = c("user_id" = "id"))
p2a <- full_join(p2a,p2a.responses, by = c("internal_id" = "session_id"))
p2a <- select(p2a, -user_id, -internal_id)
p2a <- pivot_longer(p2a, cols = -login, 
                         names_to = "variable", values_to = "result", 
                         values_drop_na = TRUE )
p2a <- distinct(p2a, login, variable, .keep_all = TRUE)
p2a <- pivot_wider(p2a, names_from = variable, values_from = result) # 1929 users went past the demographic section
p2a <- left_join(p2a, p2a.demo, by = c("login" = "userid"))
names(p2a)[names(p2a) == "login"] <- "id"

######## Merge ITOS, P2U and P2A ######### 

### give each dataframe a phase column
itos$phase <- "itos"
p2u$phase <- "p2u"
p2a$phase <- "p2a"

### standardise nationality
p2u$nationality <- "au" # all underage participants were from Australian schools
p2a$nationality <- recode(p2a$nationality, af = "nonau", as = "nonau", ca = "nonau", eu = "nonau", me = "nonau", na = "nonau", nz = "nonau", oc = "nonau", sa = "nonau")
p2a$nationality <- recode(na_if(p2a$nationality,""), .missing = "pnts")
itos$nationality <- str_to_lower(itos$nationality)
itos$nationality <- str_replace(itos$nationality, ".*au.*", "au")
temp <- filter(itos, !nationality == "au")
temp2 <- filter(itos, nationality == "au")
temp$nationality <- recode(na_if(temp$nationality,""), .missing = "pnts")
temp3 <- filter(temp, !nationality == "pnts")
temp4 <- filter(temp, nationality == "pnts")
temp3$nationality <- "nonau"
itos <- rbind(temp2,temp3,temp4)
rm(temp,temp2,temp3,temp4)

### standardise gender
itos$gender <- recode(itos$gender, male = "m", female = "f", other = "o")
p2u$gender <- recode(p2u$gender, M = "m", F = "f")
p2u$gender <- recode(na_if(p2u$gender,""), .missing = "pnts")
p2a$gender <- recode(na_if(p2a$gender,""), .missing = "pnts")

### standardise device
p2u$device <- recode(na_if(p2u$device,""), iPad = "ipad", Laptop = "laptop", .missing = "pnts")
p2a$device <- recode(na_if(p2a$device,""), .missing = "pnts")
itos$device <- "unknown"

### standardise age
p2u$age <- floor(p2u$age/12)
p2u$age <- as.integer(p2u$age)
itos$age <- as.integer(itos$age)

### score the subtests and time taken
data <- list(itos = itos, p2a = p2a, p2u = p2u)
data <- lapply(data, function(x) {
  x <- mutate(x, gcscore = rowSums(select(x, starts_with("score_gc")), na.rm = TRUE))
  x <- mutate(x, gfscore = rowSums(select(x, starts_with("score_gf")), na.rm = TRUE))
  x <- mutate(x, gvscore = rowSums(select(x, starts_with("score_gv")), na.rm = TRUE))
  x <- mutate(x, gwmscore = rowSums(select(x, starts_with("score_gwm")), na.rm = TRUE))
  x <- mutate(x, gctime = rowSums(select(x, starts_with("timeTaken_gc")), na.rm = TRUE))
  x <- mutate(x, gftime = rowSums(select(x, starts_with("timeTaken_gf")), na.rm = TRUE))
  x <- mutate(x, gvtime = rowSums(select(x, starts_with("timeTaken_gv")), na.rm = TRUE))
  x <- mutate(x, gwmtime = rowSums(select(x, starts_with("timeTaken_gwm")), na.rm = TRUE))
})               
list2env(data, .GlobalEnv)

### merge the data frames
data <- full_join(itos,p2a)
data <- full_join(data,p2u)

### sort the columns
data <- select(data, id, phase, age, gender, nationality, device, ends_with("score"), ends_with("time"), starts_with("score_gc"),
               starts_with("timeTaken_gc"), starts_with("score_gf"), starts_with("timeTaken_gf"), starts_with("score_gv"), 
               starts_with("timeTaken_gv"), starts_with("score_gwm"), starts_with("timeTaken_gwm"), starts_with("wisc"))
itos <- select(itos, id, phase, age, gender, nationality, device, ends_with("score"), ends_with("time"), starts_with("score_gc"),
               starts_with("timeTaken_gc"), starts_with("score_gf"), starts_with("timeTaken_gf"), starts_with("score_gv"), 
               starts_with("timeTaken_gv"), starts_with("score_gwm"), starts_with("timeTaken_gwm"), starts_with("wisc"))
p2a <- select(p2a, id, phase, age, gender, nationality, device, ends_with("score"), ends_with("time"), starts_with("score_gc"),
              starts_with("timeTaken_gc"), starts_with("score_gf"), starts_with("timeTaken_gf"), starts_with("score_gv"), 
              starts_with("timeTaken_gv"), starts_with("score_gwm"), starts_with("timeTaken_gwm"), starts_with("wisc"))
p2u <- select(p2u, id, phase, age, gender, nationality, device, ends_with("score"), ends_with("time"), starts_with("score_gc"),
              starts_with("timeTaken_gc"), starts_with("score_gf"), starts_with("timeTaken_gf"), starts_with("score_gv"), 
              starts_with("timeTaken_gv"), starts_with("score_gwm"), starts_with("timeTaken_gwm"), starts_with("wisc"))

######## Save Data ######### 

#### save to csv
write.csv(data, "data.csv", row.names = FALSE)
write.csv(itos, "itos_data.csv", row.names = FALSE)
write.csv(p2a, "p2a_data.csv", row.names = FALSE)
write.csv(p2u, "p2u_data.csv", row.names = FALSE)

### Remove all variables from environment
rm(list = ls())

