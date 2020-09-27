# Data Analysis of the Item Tryout Study (ITOS) CHC Test
# Original code developed by Mr Jake Kraksa (Monash University)
# Built on R 3.6.0 with R Studio 1.1.456
# set the working directory to source file directory manually

#### Load required libraries ####

require(ggplot2)
require(car)
require(knitr)
require(lavaan)
require(mirt)
require(mokken)
require(lordif)
require(dplyr)
require(tidyr)
require(stringr)
require(tibble)

#### Set Seed ####

set.seed(145)

#### Set up functions ####

fix.columns <- function(df,x,y) {
  colnames(df) <- sub(x,paste(y,deparse(substitute(df)), sep = ""),colnames(df))
  return(df)
}

cfa.parameters <- function(object, csvname, print = TRUE) {
  y <- parameterEstimates(object, standardized = TRUE) %>%
    filter(op == "=~") %>%
    dplyr::select(Item = rhs, B = est, SE = se, Z = z, "p-value" = pvalue, Beta = std.all)
  
  write.table(y, file = paste0("itos_tables/",csvname,".csv"), sep = ",", quote = FALSE, row.names = F)
    
  y %>%  
    kable(digits = 3, format="pandoc", caption=paste("Factor Loadings for ", deparse(substitute(object)), sep = ""))
}

mokken.outcomes <- function(object, csvname, print = TRUE) {
  x <- coefH(object)
  write.table(x$Hi, file = paste0("itos_tables/",csvname,".csv"), sep = ",", quote = FALSE, row.names = T)
  print(x$Hi)
  print(x$H)
}

aisp.outcomes <- function(object, csvname, print = TRUE) {
  aisp <- aisp(object)
  write.table(aisp, file = paste0("itos_tables/",csvname,".csv"), sep = ",", quote = FALSE, row.names = T)
  print(aisp)
}

rasch.outcomes <- function(object, csvname, print = TRUE) {
  rasch <- mirt(object, model = 1, itemtype = "Rasch")
  rasch.fit <- itemfit(rasch)
  print("Rasch Item Fit")
  print(rasch.fit)
  write.table(rasch.fit, file = paste0("itos_tables/",csvname,".csv"), sep = ",", quote = FALSE, row.names = T)
  print("M2 stats")
  print(M2(rasch))
  print("Marginal Rxx")
  print(marginal_rxx(rasch))
  
  return(list(rasch,rasch.fit))
}


#### Load Raw Data ####

demographics.raw <- read.csv("itos_demographics.csv", header = TRUE) # load in demographics
demographics.raw <- rbind(demographics.raw,read.csv("itos_name.csv", header = TRUE)) # load in participant names
demographics.raw <- rbind(demographics.raw,read.csv("itos_contact.csv", header = TRUE))
gc.raw <- read.csv("itos_gc_results.csv", header = TRUE) # load in Gc
gf.raw <- read.csv("itos_gf_results.csv", header = TRUE) # load in Gf
gv.raw <- read.csv("itos_gv_results.csv", header = TRUE) # load in Gv
gvt.raw <- read.csv("itos_gvt_results.csv", header = TRUE) # load in Gvt
gwm.raw <- read.csv("itos_gwm_results.csv", header = TRUE) # load in Gwm
dfList <- list(demographics.raw,gc.raw,gf.raw,gv.raw,gwm.raw,gvt.raw) # create a list of the above data frames
names(dfList) <- c("demographics.raw","gc.raw","gf.raw","gv.raw","gwm.raw","gvt.raw") # name the data frames
rm("demographics.raw","gc.raw","gf.raw","gv.raw","gwm.raw","gvt.raw") # remove individual data frames

#### Remove Unnecessary Columns ####

dfList  <- lapply(dfList, function(x) {
  x <- x[,!names(x) %in% c("id","data_set_id","trait")]
})

#### Create Data Frames ####

list2env(dfList,envir = .GlobalEnv) # turn the dfList into individual data frames
rm(dfList)

#### Remove Unnecessary Rows ####

demographics.raw <- demographics.raw[!demographics.raw$name %in% c("buttonPressed", "isTimeout", "timeTaken"), ]

#### Reshape Data ####

demographics <- spread(demographics.raw, key = name, value = value)
gc <- reshape(gc.raw, idvar="session_id", timevar = "item_id", direction = "wide")
gf <- reshape(gf.raw, idvar="session_id", timevar = "item_id", direction = "wide")
gv <- reshape(gv.raw, idvar="session_id", timevar = "item_id", direction = "wide")
gv <- gv[,c(1:100,155:157,101:154)]
gwm <- reshape(gwm.raw, idvar="session_id", timevar = "item_id", direction = "wide")
rm("demographics.raw","gc.raw","gf.raw","gv.raw","gwm.raw","gvt.raw")

#### Recode Data Types ####

demographics$age <- as.numeric(as.character(demographics$age))
demographics$email_address <- as.character(demographics$email_address)
demographics$gender <- as.factor(as.character(demographics$gender))
demographics$krongold <- as.factor(as.character(demographics$krongold))
demographics$nationality <- as.character(demographics$nationality)
demographics$postcode <- as.character(demographics$postcode)
demographics$name <- as.character(demographics$name)

#### Rename Columns ####

names(demographics)[names(demographics) == "session_id"] <- "id"
names(demographics)[names(demographics) == "email_address"] <- "email"
names(gc)[names(gc) == "session_id"] <- "id"
names(gf)[names(gf) == "session_id"] <- "id"
names(gv)[names(gv) == "session_id"] <- "id"
names(gwm)[names(gwm) == "session_id"] <- "id"

gc <- fix.columns(gc,"correct.","c.")
gc <- fix.columns(gc,"response.","r.")
gc <- fix.columns(gc,"time_taken.","t.")
gf <- fix.columns(gf,"correct.","c.")
gf <- fix.columns(gf,"response","r.")
gf <- fix.columns(gf,"time_taken.","t.")
gv <- fix.columns(gv,"correct.","c.")
gv <- fix.columns(gv,"response","r.")
gv <- fix.columns(gv,"time_taken.","t.")
gwm <- fix.columns(gwm,"correct.","c.")
gwm <- fix.columns(gwm,"response","r.")
gwm <- fix.columns(gwm,"time_taken.","t.")

#### Gender Outliers ####

filter(demographics, gender == "other")
demographics <- filter(demographics, gender != "other")
demographics$gender <- factor(demographics$gender)

#### Age Outliers ####

demographics <- filter(demographics, age < 91 & age > 17)
ggplot(demographics) + 
  geom_boxplot(mapping = aes(x = factor(0), y = age)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous("Age (years old)", breaks = c(18,30,40,50,60,70,80,90), limits = c(18,90))
age.outliers <- boxplot.stats(demographics$age)$out
age.lower <- boxplot.stats(demographics$age)$stats[1]
age.upper <- boxplot.stats(demographics$age)$stats[5]
age.outliers
print(paste("Number of outliers by age:",length(age.outliers), sep = " "))
demographics <- filter(demographics, age > age.lower - 1 & age < age.upper + 1 )
rm(age.lower,age.upper,age.outliers)

#### Create Age Groups ####

demographics$agegroup <- cut(demographics$age, breaks = c(18,30,40,50,60,70,80,90), right = FALSE)

#### Set Options ####

options(knitr.kable.NA = '')
options(max.print = 3000)

#### Demographics ####

# Plot for Gender
describe(demographics$gender)
ggplot(demographics) +
  geom_bar(mapping = aes(gender)) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
ggsave("itos_plots/Figure 1 - Gender.png")

# Tables and Plot for Age
psych::describe(demographics$age) %>%
  kable(digits=2, format="pandoc", caption="Table 1: Age for Whole Sample")

ggplot(demographics) +
  geom_bar(mapping = aes(agegroup)) +
  xlab("Age Group") +
  scale_x_discrete(labels = c("[18,30)" = "18-29", "[30,40)" = "30-39", "[40,50)" = "40-49",
                              "[50,60)" = "50-59", "[60,70)" = "60-69", "[70,80)" = "70-79",
                              "[80,90)" = "80-90")) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) +
  ylab("Count")
ggsave("itos_plots/Figure 2 - Age.png")

#### Gc Data Preparation ####

gc.items <- merge(select(gc,id,contains("c."),contains("t.")),select(demographics,id,age,gender,agegroup))
gc.items[is.na(gc.items)] <- 0 # Change all NA to 0

#### Gc Column Identifiers ####

item.count <- ncol(select(gc.items,contains("c.")))
c.firstcol <- which(colnames(gc.items) =="c.gc1")
c.lastcol <- which(colnames(gc.items) == paste("c.gc",item.count,sep=""))
t.firstcol <- which(colnames(gc.items) == "t.gc1")
t.lastcol <- which(colnames(gc.items) == paste("t.gc",item.count,sep=""))

#### Gc Outcome Columns ####

gc.items$result <- rowSums(gc.items[,c.firstcol:c.lastcol]) # add a result column
gc.items$time <- rowSums(gc.items[,t.firstcol:t.lastcol]) # add a time taken column

#### Gc Time Outliers ####

gc.time.upper <- quantile(gc.items$time, 0.95)
gc.time.lower <- quantile(gc.items$time, 0.05)
print(paste("95th percentile for time:", gc.time.upper, sep = " "))
print(paste("5th percentile for time:", gc.time.lower, sep = " "))
gc.time.outliers <- nrow(filter(gc.items, time < gc.time.lower | time > gc.time.upper))
print(paste("Number of outliers by time:", gc.time.outliers, sep = " "))
gc.items <- filter(gc.items, time >= gc.time.lower & time <= gc.time.upper )
rm(gc.time.lower,gc.time.upper,gc.time.outliers)

#### Gc Item and Person Removal ####

number.of.items <- ncol(select(gc.items, contains("c.")))
number.of.people <- nrow(gc.items)
gc.items <- filter(gc.items, result !=0 & result != item.count)  # remove persons that achieved a full score or a score of 0
gc.items <- gc.items[vapply(gc.items, function(x) length(unique(x)) > 1, logical(1L))] # remove items that noone got wrong or everyone got right
number.of.items.2 <- ncol(select(gc.items, contains("c.")))
number.of.people.2 <- nrow(gc.items)
print(paste("Number of people removed:",number.of.people - number.of.people.2, sep = " "))
print(paste("Number of items removed:",number.of.items - number.of.items.2, sep = " "))

#### Gc Demographics ####

describe(gc.items$gender)

ggplot(gc.items) +
  geom_bar(mapping = aes(gender)) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
ggsave("itos_plots/Figure 3 - Gc Gender.png")

psych::describe(gc.items$age) %>%
  kable(digits=2, format="pandoc", caption="Table 3: Age for Gc")

ggplot(gc.items) +
  geom_bar(mapping = aes(agegroup)) +
  xlab("Age Group") +
  scale_x_discrete(labels = c("[18,30)" = "18-29", "[30,40)" = "30-39", "[40,50)" = "40-49",
                              "[50,60)" = "50-59", "[60,70)" = "60-69", "[70,80)" = "70-79",
                              "[80,90)" = "80-90")) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) +
  ylab("Count")
ggsave("itos_plots/Figure 4 - Gc Age.png")

#### Gc Graphs ####

ggplot(gc.items, aes(x = time, y = result)) +
  geom_point() +
  xlab("Time Taken (secs)") +
  ylab("Total Score") +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(60,720,30)), limits = c(60, 720)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))
ggsave("itos_plots/Figure 5 - Gc Time Taken.png")

ggplot(gc.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = time)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 6 - Gc Time Boxplot.png")

ggplot(gc.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = result)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 7 - Gc Result Boxplot.png")

ggplot(gc.items) +
  geom_boxplot(mapping = aes(x = gender, y = result)) +
  ylab("Total Score") +
  xlab("Gender")
ggsave("itos_plots/Figure 8 - Gc Gender Result Boxplot.png")

ggplot(gc.items) +
  geom_boxplot(mapping = aes(x = gender, y = time)) +
  ylab("Time Taken") +
  xlab("Gender")
ggsave("itos_plots/Figure 9 - Gc Gender Time Taken Boxplot.png")

#### Gc Gender Differences ####

t.test(gc.items$result ~ gc.items$gender)
t.test(gc.items$time ~ gc.items$gender)

#### Gc Sample Setup ####

gc.eval <- select(gc.items,contains("c."))

#### Gc Reliability ####

psych::alpha(gc.eval)

#### Gc Count Number Correct ####

gc.correct <- as.data.frame(
  sapply(gc.eval, function(x) {
    (1 - (table(x)[1]/length(x))) * 100
  })
)

row.names(gc.correct) <- str_remove(row.names(gc.correct), "c.gc")
row.names(gc.correct) <- str_remove(row.names(gc.correct), "\\.0")

names(gc.correct) <- "Correct"

gc.correct <- rownames_to_column(gc.correct, var = "Item")

gc.correct$Item <- as.numeric(as.vector(gc.correct$Item))

ggplot(data = gc.correct, aes(x = Item, y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  ylab("% Correct") + 
  scale_x_continuous(breaks=c(1:55)) + 
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("itos_plots/Figure 10 - Gc Number of Items Correct.png")

psych::describe(gc.items$result)
psych::describe(gc.items$time)

#### Gc CFA ####

gc.model <- paste("gc =~ ", paste(colnames(gc.eval), collapse = "+"), sep = "")
gc.fit <- cfa(gc.model, data=gc.eval, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gc.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(object = gc.fit, csvname = "Gc CFA 1") # Parameter Estimates for the CFA

#### Gc Item Removal ####

gc.eval.2 <- gc.eval
gc.eval.2$c.gc1 <- NULL # removed because not statistically signficant
gc.eval.2$c.gc3 <- NULL # removed because not statistically signficaint

#### Gc CFA 2 ####

gc.model.2 <- paste("gc =~ ", paste(colnames(gc.eval.2), collapse = "+"), sep = "")
gc.fit.2 <- cfa(gc.model.2, data=gc.eval.2, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gc.fit.2, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gc.fit.2,"Gc CFA 2") # Parameter Estimates for the CFA

#### Gc Mokken ####

mokken.outcomes(gc.eval.2, "Gc Mokken 1")
gc.aisp <- aisp.outcomes(gc.eval.2, "Gc AISP 1")

#### Gc Item Removal 2 ####

gc.eval.3 <- gc.eval.2[, gc.aisp == 1]

#### Gc CFA 3 ####

gc.model.3 <- paste("gc =~ ", paste(colnames(gc.eval.3), collapse = "+"), sep = "")
gc.fit.3 <- cfa(gc.model.3, data=gc.eval.3, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gc.fit.3, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gc.fit.3,"Gc CFA 3") # Parameter Estimates for the CFA

#### Gc Mokken 2 ####

mokken.outcomes(gc.eval.3, "Gc Mokken 2")
gc.aisp.2 <- aisp.outcomes(gc.eval.3, "Gc AISP 2")

#### Gc IRT ####

rasch <- rasch.outcomes(gc.eval.3,"Gc Rasch 1")
gc.rasch <- rasch[[1]]
gc.rasch.fit <- rasch[[2]]

#### Gc Local Dependency ####

residuals(gc.rasch, type = "Q3", digits = 2, suppress = +.2)

#### Gf Item Removal 3 ####

gc.eval.4 <- select(gc.eval.3, -(c.gc4:c.gc33))
gc.eval.4$c.gc35 <- NULL
gc.eval.4$c.gc36 <- NULL
gc.eval.4$c.gc37 <- NULL
gc.eval.4$c.gc38 <- NULL
gc.eval.4$c.gc40 <- NULL
gc.eval.4$c.gc41 <- NULL
gc.eval.4$c.gc43 <- NULL
gc.eval.4$c.gc44 <- NULL
gc.eval.4$c.gc45 <- NULL
gc.eval.4$c.gc47 <- NULL

# Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for 
# dichotomous item response theory models. Applied Psychological Measurement, 24, 50-64.
gc.rasch <- mirt(gc.eval.4, model = 1, itemtype = "Rasch") # conduct Rasch modeling
gc.rasch.fit <- itemfit(gc.rasch)
colnames(select(gc.eval.4, as.vector(filter(gc.rasch.fit, p.S_X2 < 0.01)$item))) # items that are going to get removed
colnames(select(gc.eval.4, as.vector(filter(gc.rasch.fit, p.S_X2 > 0.01)$item))) # items that are going to be kept
gc.eval.4 <- select(gc.eval.4, as.vector(filter(gc.rasch.fit, p.S_X2 > 0.01)$item))

#### Gc CFA 4 ####

gc.model.4 <- paste("gc =~ ", paste(colnames(gc.eval.4), collapse = "+"), sep = "")
gc.fit.4 <- cfa(gc.model.4, data=gc.eval.4, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gc.fit.4, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gc.fit.4,"Gc CFA 4") # Parameter Estimates for the CFA

#### Gc Mokken 3 ####

mokken.outcomes(gc.eval.4, "Gc Mokken 3")
gc.aisp.3 <- aisp.outcomes(gc.eval.4, "Gc AISP 3")

#### Gc IRT 2 ####

rasch <- rasch.outcomes(gc.eval.4, "Gc Rasch 2")
gc.rasch.2 <- rasch[[1]]
gc.rasch.fit.2 <- rasch[[2]]

#### Gc Local Dependency 2 ####

residuals(gc.rasch.2, type = "Q3", digits = 2, suppress = +.2)

#### Gc Reliability 2 ####

psych::alpha(gc.eval.4)

#### Gc Differential Item Functioning ####

sapply(names(gc.eval.4), function(x) {
  table(gc.eval.4[[x]],gc.items[,"gender"])
})

gc.dif <- lordif(gc.eval.4,gc.items[,"gender"])
summary(gc.dif)
plot(gc.dif)

#### Gc Item Characteristic Curves ####

jpeg(file = paste("itos_plots//Figure 11 - ICC of Gc Items.jpeg"))
plot(gc.rasch.2, type="trace", facet_items = FALSE, main ="")
dev.off()

#### Gc Test Information Curve ####

jpeg(file = paste("itos_plots//Figure 12 - Gc Test Information Curve.jpeg"))
plot(gc.rasch.2, type = "info", xlim=c(-5,5), main = "")
dev.off()

#### Gc Item Parameters ####

gc.parameters <- coef(gc.rasch.2, simplify=TRUE, IRTpars = TRUE, printSE = TRUE)
gc.parameters$items

#### Gf Data Preparation ####

gf.items <- merge(select(gf,id,contains("c."),contains("t.")),select(demographics,id,age,gender,agegroup))
gf.items[is.na(gf.items)] <- 0 # Change all NA to 0

#### Gf Column Identifiers ####

item.count <- ncol(select(gf.items,contains("c.")))
c.firstcol <- which(colnames(gf.items) =="c.gf1")
c.lastcol <- which(colnames(gf.items) == paste("c.gf",item.count,sep=""))
t.firstcol <- which(colnames(gf.items) == "t.gf1")
t.lastcol <- which(colnames(gf.items) == paste("t.gf",item.count,sep=""))

#### Gf Outcome Columns ####

gf.items$result <- rowSums(gf.items[,c.firstcol:c.lastcol]) # add a result column
gf.items$time <- rowSums(gf.items[,t.firstcol:t.lastcol]) # add a time taken column

#### Gf Time Outliers ####

gf.time.upper <- quantile(gf.items$time, 0.95)
gf.time.lower <- quantile(gf.items$time, 0.05)
print(paste("95th percentile for time:", gf.time.upper, sep = " "))
print(paste("5th percentile for time:", gf.time.lower, sep = " "))
gf.time.outliers <- nrow(filter(gf.items, time < gf.time.lower | time > gf.time.upper))
print(paste("Number of outliers by time:", gf.time.outliers, sep = " "))
gf.items <- filter(gf.items, time >= gf.time.lower & time <= gf.time.upper )
rm(gf.time.lower,gf.time.upper,gf.time.outliers)

#### Gf Item and Person Removal ####

number.of.items <- ncol(select(gf.items, contains("c.")))
number.of.people <- nrow(gf.items)
gf.items <- filter(gf.items, result !=0 & result != item.count)  # remove persons that achieved a full score or a score of 0
gf.items <- gf.items[vapply(gf.items, function(x) length(unique(x)) > 1, logical(1L))] # remove items that noone got wrong or everyone got right
number.of.items.2 <- ncol(select(gf.items, contains("c.")))
number.of.people.2 <- nrow(gf.items)
print(paste("Number of people removed:",number.of.people - number.of.people.2, sep = " "))
print(paste("Number of items removed:",number.of.items - number.of.items.2, sep = " "))

#### Gf Demographics ####

describe(gf.items$gender)

ggplot(gf.items) +
  geom_bar(mapping = aes(gender)) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
ggsave("itos_plots/Figure 3 - Gf Gender.png")

psych::describe(gf.items$age) %>%
  kable(digits=2, format="pandoc", caption="Table 3: Age for Gf")

ggplot(gf.items) +
  geom_bar(mapping = aes(agegroup)) +
  xlab("Age Group") +
  scale_x_discrete(labels = c("[18,30)" = "18-29", "[30,40)" = "30-39", "[40,50)" = "40-49",
                              "[50,60)" = "50-59", "[60,70)" = "60-69", "[70,80)" = "70-79",
                              "[80,90)" = "80-90")) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) +
  ylab("Count")
ggsave("itos_plots/Figure 4 - Gf Age.png")

#### Gf Graphs ####

ggplot(gf.items, aes(x = time, y = result)) +
  geom_point() +
  xlab("Time Taken (secs)") +
  ylab("Total Score") +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(210,1080,30)), limits = c(210, 1080)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))
ggsave("itos_plots/Figure 5 - Gf Time Taken.png")

ggplot(gf.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = time)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 6 - Gf Time Boxplot.png")

ggplot(gf.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = result)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 7 - Gf Result Boxplot.png")

ggplot(gf.items) +
  geom_boxplot(mapping = aes(x = gender, y = result)) +
  ylab("Total Score") +
  xlab("Gender")
ggsave("itos_plots/Figure 8 - Gf Gender Result Boxplot.png")

ggplot(gf.items) +
  geom_boxplot(mapping = aes(x = gender, y = time)) +
  ylab("Time Taken") +
  xlab("Gender")
ggsave("itos_plots/Figure 9 - Gf Gender Time Taken Boxplot.png")

#### Gf Gender Differences ####

t.test(gf.items$result ~ gf.items$gender)
t.test(gf.items$time ~ gf.items$gender)

#### Gf Sample Setup ####

gf.eval <- select(gf.items,contains("c."))

#### Gf Reliability ####

psych::alpha(gf.eval)

#### Gf Count Number Correct ####

gf.correct <- as.data.frame(
  sapply(gf.eval, function(x) {
    (1 - (table(x)[1]/length(x))) * 100
  })
)

row.names(gf.correct) <- str_remove(row.names(gf.correct), "c.gf")
row.names(gf.correct) <- str_remove(row.names(gf.correct), "\\.0")

names(gf.correct) <- "Correct"

gf.correct <- rownames_to_column(gf.correct, var = "Item")

gf.correct$Item <- as.numeric(as.vector(gf.correct$Item))

ggplot(data = gf.correct, aes(x = Item, y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  ylab("% Correct") + 
  scale_x_continuous(breaks=c(1:33)) + 
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("itos_plots/Figure 10 - Gf Number of Items Correct.png")

psych::describe(gf.items$result)
psych::describe(gf.items$time)

#### Gf Analyse Item Format Differences ####

matrixtype <- c(1,2,1,2,2,3,3,1,2,1,2,2,2,4,3,2,1,1,3,2,4,4,4,5,1,1,4,4,4,5,4,1,5)
responseoptions <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,2,3,2,2,2,2,3,3,2,2,2,2)
itemtype <- c(1,2,1,2,2,3,3,1,2,1,2,2,2,4,3,2,1,1,3,2,4,4,5,6,7,7,4,5,5,6,4,7,6)
timeallowed <- c(rep(1,15),rep(2,18))
items <- names(gf.eval)
lookup <- data.frame(items,matrixtype,responseoptions,itemtype,timeallowed)

lookup$matrixtype <- factor(lookup$matrixtype, 
                            levels = c(1:5),
                            labels = c("5x1","2x2","4x1","3x3","4x2"))

lookup$responseoptions <- factor(lookup$responseoptions,
                                 levels = c(1:3),
                                 labels = c("Four","Five","Six"))
lookup$itemtype <- factor(lookup$itemtype, 
                          levels = c(1:7),
                          labels = c("5x1M 4RO", "2x2M 4RO", "4x1M 4RO", "3x3M 5RO", "3x3M 6RO", "4x2M 5RO", "5x1M 5RO"))

lookup$timeallowed <- factor(lookup$timeallowed,
                             levels = c(1:2),
                             labels = c("30 seconds", "60 seconds"))

gf.wide <- gather(gf.eval, item, score, c.gf1:c.gf33, factor_key=TRUE)
gf.wide$timeallowed <- ifelse(gf.wide$item %in% names(gf.eval)[1:15], yes = "30 seconds", no = "60 seconds")
gf.wide$matrixtype <- sapply(gf.wide$item, function(x) lookup$matrixtype[match(x, lookup$items)])
gf.wide$responseoptions <- sapply(gf.wide$item, function(x) lookup$responseoptions[match(x, lookup$items)])
gf.wide$itemtype <- sapply(gf.wide$item, function(x) lookup$itemtype[match(x, lookup$items)])
gf.wide$score <- factor(gf.wide$score)

# ANOVA

gf.correct$itemtype <- as.factor(itemtype)
gf.correct$timeallowed <- as.factor(timeallowed)
gf.correct$matrixtype <- as.factor(matrixtype)
gf.correct$responseoptions <- as.factor(responseoptions)

table(gf.correct$itemtype, gf.correct$timeallowed)
table(gf.correct$matrixtype, gf.correct$responseoptions)

ggplot(data = gf.correct, aes(x = itemtype, y = Correct, colour = timeallowed)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Correct %") +
  xlab("Item Group") +
  labs(colour = "Time Allowed")

ggplot(data = gf.correct, aes(x = matrixtype, y = Correct, colour = responseoptions)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  ylab("Correct %") +
  xlab("Matrix Type") +
  labs(colour = "Response Options")

item.group.anova.1 <- aov(Correct ~ itemtype * timeallowed, data = gf.correct)
Anova(item.group.anova.1, type = "II")

item.group.anova.2 <- aov(Correct ~ matrixtype * responseoptions * timeallowed, data = gf.correct)
Anova(item.group.anova.2, type = "II")

# Logistic Regression (Matrix Type and RO Separate) 

gf.logistic1 <- glm(score ~ matrixtype + responseoptions + timeallowed, data = gf.wide, family = "binomial")
summary(gf.logistic1)
exp(cbind(OR = coef(gf.logistic1), confint(gf.logistic1)))

# Logistic Regression (Matrix Type and RO Combined)

gf.logistic2 <- glm(score ~ itemtype + timeallowed, data = gf.wide, family = "binomial")
summary(gf.logistic2)
exp(cbind(OR = coef(gf.logistic2), confint(gf.logistic2)))

rm(matrixtype,responseoptions,itemtype,timeallowed,items,lookup,gf.wide,gf.logistic2,gf.correct,gf.logistic1)

#### Gf CFA ####

gf.model <- paste("gf =~ ", paste(colnames(gf.eval), collapse = "+"), sep = "")
gf.fit <- cfa(gf.model, data=gf.eval, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gf.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gf.fit, "Gf CFA 1") # Parameter Estimates for the CFA

#### Gf Item Removal ####

gf.eval.2 <- gf.eval
gf.eval.2$c.gf2 <- NULL # removed because not statistically signficant
gf.eval.2$c.gf9 <- NULL # removed because not statistically signficaint

#### Gf CFA 2 ####

gf.model.2 <- paste("gf =~ ", paste(colnames(gf.eval.2), collapse = "+"), sep = "")
gf.fit.2 <- cfa(gf.model.2, data=gf.eval.2, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gf.fit.2, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gf.fit.2, "Gf CFA 2") # Parameter Estimates for the CFA

#### Gf Mokken ####

mokken.outcomes(gf.eval.2, "Gf Mokken 1")
gf.aisp <- aisp.outcomes(gf.eval.2, "Gf AISP 1")

#### Gf Item Removal ####

gf.eval.3 <- gf.eval.2[, gf.aisp == 1]

#### Gf CFA 3 ####

gf.model.3 <- paste("gf =~ ", paste(colnames(gf.eval.3), collapse = "+"), sep = "")
gf.fit.3 <- cfa(gf.model.3, data=gf.eval.3, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gf.fit.3, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gf.fit.3, "Gf CFA 3") # Parameter Estimates for the CFA

#### Gf Mokken 2 ####

mokken.outcomes(gf.eval.3, "Gf Mokken 2")
gf.aisp.2 <- aisp.outcomes(gf.eval.3, "Gf AISP 2")

#### Gf IRT ####

rasch <- rasch.outcomes(gf.eval.3, "Gf Rasch 1")
gf.rasch <- rasch[[1]]
gf.rasch.fit <- rasch[[2]]

#### Gf Local Dependency ####

residuals(gf.rasch, type = "Q3", digits = 2, suppress = +.2)

#### Gf Item Removal 2 ####

# Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for 
# dichotomous item response theory models. Applied Psychological Measurement, 24, 50-64.
colnames(select(gf.eval.3, as.vector(filter(gf.rasch.fit, p.S_X2 < 0.01)$item))) # items that are going to get removed
colnames(select(gf.eval.3, as.vector(filter(gf.rasch.fit, p.S_X2 > 0.01)$item))) # items that are going to be kept
gf.eval.4 <- select(gf.eval.3, as.vector(filter(gf.rasch.fit, p.S_X2 > 0.01)$item))

#### Gf CFA 4 ####

gf.model.4 <- paste("gf =~ ", paste(colnames(gf.eval.4), collapse = "+"), sep = "")
gf.fit.4 <- cfa(gf.model.4, data=gf.eval.4, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gf.fit.4, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gf.fit.4, "Gf CFA 4") # Parameter Estimates for the CFA

#### Gf Mokken 3 ####

mokken.outcomes(gf.eval.4, "Gf Mokken 3")
aisp.outcomes(gf.eval.4, "Gf AISP 3")

#### Gf IRT 2 ####

rasch <- rasch.outcomes(gf.eval.4, "Gf Rasch 2")
gf.rasch.2 <- rasch[[1]]
gf.rasch.fit.2 <- rasch[[2]]

#### Gf Local Dependency 2 ####

residuals(gf.rasch.2, type = "Q3", digits = 2, suppress = +.2)

#### Gf Reliability 2 ####

psych::alpha(gf.eval.4)

#### Gf Differential Item Functioning ####

sapply(names(gf.eval.4), function(x) {
  table(gf.eval.4[[x]],gf.items[,"gender"])
})

gf.eval.4$c.gf10 <- NULL

gf.dif <- lordif(gf.eval.4,gf.items[,"gender"])
summary(gf.dif)
plot(gf.dif)

#### Gf Item Characteristic Curves ####

jpeg(file = paste("itos_plots//Figure 11 - ICC of Gf Items.jpeg"))
plot(gf.rasch.2, type="trace", facet_items = FALSE, main ="")
dev.off()

#### Gf Test Information Curve ####

jpeg(file = paste("itos_plots//Figure 12 - Gf Test Information Curve.jpeg"))
plot(gf.rasch.2, type = "info", xlim=c(-5,5), main = "")
dev.off()

#### Gf Item Parameters ####

gf.parameters <- coef(gf.rasch.2, simplify=TRUE, IRTpars = TRUE, printSE = TRUE)
gf.parameters$items

#### Gv Data Preparation ####

gv.items <- merge(select(gv,id,contains("c."),contains("t.")),select(demographics,id,age,gender,agegroup))
gv.items[is.na(gv.items)] <- 0 # Change all NA to 0

#### Gv Column Identifiers ####

item.count <- ncol(select(gv.items,contains("c.")))
c.firstcol <- which(colnames(gv.items) =="c.gv1")
c.lastcol <- which(colnames(gv.items) == paste("c.gv",item.count,sep=""))
t.firstcol <- which(colnames(gv.items) == "t.gv1")
t.lastcol <- which(colnames(gv.items) == paste("t.gv",item.count,sep=""))

#### Gv Outcome Columns ####

gv.items$result <- rowSums(gv.items[,c.firstcol:c.lastcol]) # add a result column
gv.items$time <- rowSums(gv.items[,t.firstcol:t.lastcol]) # add a time taken column

#### Gv Time Outliers ####

gv.time.upper <- quantile(gv.items$time, 0.95)
gv.time.lower <- quantile(gv.items$time, 0.05)
print(paste("95th percentile for time:", gv.time.upper, sep = " "))
print(paste("5th percentile for time:", gv.time.lower, sep = " "))
gv.time.outliers <- nrow(filter(gv.items, time < gv.time.lower | time > gv.time.upper))
print(paste("Number of outliers by time:", gv.time.outliers, sep = " "))
gv.items <- filter(gv.items, time >= gv.time.lower & time <= gv.time.upper )
rm(gv.time.lower,gv.time.upper,gv.time.outliers)

#### Gv Item and Person Removal ####

number.of.items <- ncol(select(gv.items, contains("c.")))
number.of.people <- nrow(gv.items)
gv.items <- filter(gv.items, result !=0 & result != item.count)  # remove persons that achieved a full score or a score of 0
gv.items <- gv.items[vapply(gv.items, function(x) length(unique(x)) > 1, logical(1L))] # remove items that noone got wrong or everyone got right
number.of.items.2 <- ncol(select(gv.items, contains("c.")))
number.of.people.2 <- nrow(gv.items)
print(paste("Number of people removed:",number.of.people - number.of.people.2, sep = " "))
print(paste("Number of items removed:",number.of.items - number.of.items.2, sep = " "))

#### Gv Demographics ####

describe(gv.items$gender)

ggplot(gv.items) +
  geom_bar(mapping = aes(gender)) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
ggsave("itos_plots/Figure 3 - Gv Gender.png")

psych::describe(gv.items$age) %>%
  kable(digits=2, format="pandoc", caption="Age for Gv")

ggplot(gv.items) +
  geom_bar(mapping = aes(agegroup)) +
  xlab("Age Group") +
  scale_x_discrete(labels = c("[18,30)" = "18-29", "[30,40)" = "30-39", "[40,50)" = "40-49",
                              "[50,60)" = "50-59", "[60,70)" = "60-69", "[70,80)" = "70-79",
                              "[80,90)" = "80-90")) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) +
  ylab("Count")
ggsave("itos_plots/Figure 4 - Gv Age.png")

#### Gv Graphs ####

ggplot(gv.items, aes(x = time, y = result)) +
  geom_point() +
  xlab("Time Taken (secs)") +
  ylab("Total Score") +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(30,960,30)), limits = c(30, 960)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))
ggsave("itos_plots/Figure 5 - Gv Time Taken.png")

ggplot(gv.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = time)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 6 - Gv Time Boxplot.png")

ggplot(gv.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = result)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 7 - Gv Result Boxplot.png")

ggplot(gv.items) +
  geom_boxplot(mapping = aes(x = gender, y = result)) +
  ylab("Total Score") +
  xlab("Gender")
ggsave("itos_plots/Figure 8 - Gv Gender Result Boxplot.png")

ggplot(gv.items) +
  geom_boxplot(mapping = aes(x = gender, y = time)) +
  ylab("Time Taken") +
  xlab("Gender")
ggsave("itos_plots/Figure 9 - Gv Gender Time Taken Boxplot.png")

#### Gv Gender Differences ####

t.test(gv.items$result ~ gv.items$gender)
t.test(gv.items$time ~ gv.items$gender)

#### Gv Sample Setup ####

gv.eval <- select(gv.items,contains("c."))

#### Gv Reliabiliy ####

psych::alpha(gv.eval)

#### Gv Count Number Correct ####

gv.correct <- as.data.frame(
  sapply(gv.eval, function(x) {
    (1 - (table(x)[1]/length(x))) * 100
  })
)

row.names(gv.correct) <- str_remove(row.names(gv.correct), "c.gv")
row.names(gv.correct) <- str_remove(row.names(gv.correct), "\\.0")

names(gv.correct) <- "Correct"

gv.correct <- rownames_to_column(gv.correct, var = "Item")

gv.correct$Item <- as.numeric(as.vector(gv.correct$Item))

ggplot(data = gv.correct, aes(x = Item, y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  ylab("% Correct") + 
  scale_x_continuous(breaks=c(1:52)) + 
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("itos_plots/Figure 10 - Gv Number of Items Correct.png")

psych::describe(gv.items$result)
psych::describe(gv.items$time)

#### Gv Analyse Item Format Differences ####

pieces <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
border <- c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
items <- names(gv.eval)
lookup <- data.frame(items,pieces,border)

gv.wide <- gather(gv.eval, item, score, c.gv1:c.gv52)
gv.wide$pieces <- sapply(gv.wide$item, function(x) lookup$pieces[match(x, lookup$items)])
gv.wide$border <- sapply(gv.wide$item, function(x) lookup$border[match(x, lookup$items)])
gv.wide$score <- factor(gv.wide$score)

# Logistic Regression 

gv.logistic <- glm(score ~ pieces + border, data = gv.wide, family = "binomial")
summary(gv.logistic)
exp(cbind(OR = coef(gv.logistic), confint(gv.logistic)))

rm(pieces,border,items,lookup,gv.wide,gv.logistic)

#### Gv CFA ####

gv.model <- paste("gv =~ ", paste(colnames(gv.eval), collapse = "+"), sep = "")
gv.fit <- cfa(gv.model, data=gv.eval, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gv.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gv.fit, "Gv CFA 1") # Parameter Estimates for the CFA

#### Gv Mokken ####

mokken.outcomes(gv.eval, "Gv Mokken 1")
gv.aisp <- aisp.outcomes(gv.eval, "Gv AISP 1")

#### Gv Item Removal ####

gv.eval.2 <- gv.eval[, gv.aisp == 1]

#### Gv CFA 2 ####

gv.model.2 <- paste("gv =~ ", paste(colnames(gv.eval.2), collapse = "+"), sep = "")
gv.fit.2 <- cfa(gv.model.2, data=gv.eval.2, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gv.fit.2, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gv.fit.2, "Gv CFA 2") # Parameter Estimates for the CFA

#### Gv Mokken 2 ####

mokken.outcomes(gv.eval.2, "Gv Mokken 2")
gv.aisp.2 <- aisp.outcomes(gv.eval.2, "Gv AISP 2")

#### Gv IRT ####

rasch <- rasch.outcomes(gv.eval.2, "Gv Rasch 1")
gv.rasch <- rasch[[1]]
gv.rasch.fit <- rasch[[2]]

#### Gv Local Dependency ####

residuals(gv.rasch, type = "Q3", digits = 2, suppress = +.2)

#### Gv Item Removal 2 ####

# Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for 
# dichotomous item response theory models. Applied Psychological Measurement, 24, 50-64. 
colnames(select(gv.eval.2, as.vector(filter(gv.rasch.fit, p.S_X2 < 0.01)$item))) # items that are going to get removed
colnames(select(gv.eval.2, as.vector(filter(gv.rasch.fit, p.S_X2 > 0.01)$item))) # items that are going to be kept
gv.eval.3 <- select(gv.eval.2, as.vector(filter(gv.rasch.fit, p.S_X2 > 0.01)$item))

#### Gv CFA 3 ####

gv.model.3 <- paste("gv =~ ", paste(colnames(gv.eval.3), collapse = "+"), sep = "")
gv.fit.3 <- cfa(gv.model.3, data=gv.eval.3, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gv.fit.3, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gv.fit.3, "Gv CFA 3") # Parameter Estimates for the CFA

#### Gv Mokken 3 ####

mokken.outcomes(gv.eval.3, "Gv Mokken 3")
gv.aisp.3 <- aisp.outcomes(gv.eval.3, "Gv AISP 3")

#### Gv IRT 2 ####

rasch <- rasch.outcomes(gv.eval.3, "Gv Rasch 2")
gv.rasch.2 <- rasch[[1]]
gv.rasch.fit.2 <- rasch[[2]]

#### Gv Local Dependency 2 ####

residuals(gv.rasch.2, type = "Q3", digits = 2, suppress = +.2)

#### Gv Item Removal 3 ####

# Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for 
# dichotomous item response theory models. Applied Psychological Measurement, 24, 50-64. 
colnames(select(gv.eval.3, as.vector(filter(gv.rasch.fit.2, p.S_X2 < 0.01)$item))) # items that are going to get removed
colnames(select(gv.eval.3, as.vector(filter(gv.rasch.fit.2, p.S_X2 > 0.01)$item))) # items that are going to be kept
gv.eval.4 <- select(gv.eval.3, as.vector(filter(gv.rasch.fit.2, p.S_X2 > 0.01)$item))

#### Gv CFA 4 ####

gv.model.4 <- paste("gv =~ ", paste(colnames(gv.eval.4), collapse = "+"), sep = "")
gv.fit.4 <- cfa(gv.model.4, data=gv.eval.4, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gv.fit.4, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gv.fit.4, "Gv CFA 4") # Parameter Estimates for the CFA

#### Gv Mokken 4 ####

mokken.outcomes(gv.eval.4, "Gv Mokken 4")
gv.aisp.4 <- aisp.outcomes(gv.eval.4, "Gv AISP 4")

#### Gv IRT 3 ####

rasch <- rasch.outcomes(gv.eval.4, "Gv Rasch 3")
gv.rasch.3 <- rasch[[1]]
gv.rasch.fit.3 <- rasch[[2]]

#### Gv Local Dependency 3 ####

residuals(gv.rasch.3, type = "Q3", digits = 2, suppress = +.2)

#### Gv Reliability 2 ####

psych::alpha(gv.eval.4)

#### Gv Differential Item Functioning ####

sapply(names(gv.eval.4), function(x) {
  table(gv.eval.4[[x]],gv.items[,"gender"])
})

gv.dif <- lordif(gv.eval.4,gv.items[,"gender"])
summary(gv.dif)
plot(gv.dif)

#### Gv Item Characteristic Curves ####

jpeg(file = paste("itos_plots//Figure 11 - ICC of Gv Items.jpeg"))
plot(gv.rasch.3, type="trace", facet_items = FALSE, main ="")
dev.off()

#### Gv Test Information Curve ####

jpeg(file = paste("itos_plots//Figure 12 - Gv Test Information Curve.jpeg"))
plot(gv.rasch.3, type = "info", xlim=c(-5,5), main = "")
dev.off()

#### Gv Item Parameters ####

gv.parameters <- mirt::coef(gv.rasch.3, simplify=TRUE, IRTpars = TRUE, printSE = TRUE)
gv.parameters$items

#### Gwm Data Preparation ####

gwm.items <- merge(select(gwm,id,contains("c."),contains("t.")),select(demographics,id,age,gender,agegroup))
gwm.items[is.na(gwm.items)] <- 0 # Change all NA to 0

#### Gwm Column Identifiers ####

item.count <- ncol(select(gwm.items,contains("c."))) # set a variable that counts how many items are in the data set
c.firstcol <- which(colnames(gwm.items) =="c.gwm1") # set the first item
c.lastcol <- which(colnames(gwm.items) == paste("c.gwm",item.count,sep="")) # set the last item
t.firstcol <- which(colnames(gwm.items) == "t.gwm1") # set the first item for time taken
t.lastcol <- which(colnames(gwm.items) == paste("t.gwm",item.count,sep="")) # set the last item for time taken

#### Gwm Outcome Columns ####

gwm.items$result <- rowSums(gwm.items[,c.firstcol:c.lastcol]) # add a result column
gwm.items$time <- rowSums(gwm.items[,t.firstcol:t.lastcol]) # add a time taken column

#### Gwm Time Outliers ####

gwm.time.upper <- quantile(gwm.items$time, 0.95) # identify the upper time cutoff
gwm.time.lower <- quantile(gwm.items$time, 0.05) # identify the lower time cutoff
print(paste("99th percentile for time:", gwm.time.upper, sep = " ")) # print the upper time cutoff
print(paste("5th percentile for time:", gwm.time.lower, sep = " ")) # print the lower time cutoff
gwm.time.outliers <- nrow(filter(gwm.items, time < gwm.time.lower | time > gwm.time.upper)) # choose the participants outside the time cutoffs
print(paste("Number of outliers by time:", gwm.time.outliers, sep = " ")) # print the number of participants outside the time cutoffs
gwm.items <- filter(gwm.items, time >= gwm.time.lower & time <= gwm.time.upper ) # remove the time outliers from the data set
rm(gwm.time.lower,gwm.time.upper,gwm.time.outliers)

#### Gwm Item and Person Removal ####

number.of.items <- ncol(select(gwm.items, contains("c."))) # set a variable that contains the number of items
number.of.people <- nrow(gwm.items) # set a variable that contains the number of people
gwm.items <- filter(gwm.items, result !=0 & result != item.count)  # remove persons that achieved a full score or a score of 0
gwm.items <- gwm.items[vapply(gwm.items, function(x) length(unique(x)) > 1, logical(1L))] # remove items that noone got wrong or everyone got right
number.of.items.2 <- ncol(select(gwm.items, contains("c."))) # set a new variable that contains the number of items
number.of.people.2 <- nrow(gwm.items) # set a new variable that contains the number of people
print(paste("Number of people removed:",number.of.people - number.of.people.2, sep = " ")) # print the number of people removed
print(paste("Number of items removed:",number.of.items - number.of.items.2, sep = " ")) # print the number of items removed

#### Gwm Demographics ####

describe(gwm.items$gender) # give the descriptives for gender

ggplot(gwm.items) +
  geom_bar(mapping = aes(gender)) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))
ggsave("itos_plots/Figure 3 - Gwm Gender.png")

psych::describe(gwm.items$age) %>%
  kable(digits=2, format="pandoc", caption="Age for Gwm") # give the descriptives for age

ggplot(gwm.items) +
  geom_bar(mapping = aes(agegroup)) +
  xlab("Age Group") +
  scale_x_discrete(labels = c("[18,30)" = "18-29", "[30,40)" = "30-39", "[40,50)" = "40-49",
                              "[50,60)" = "50-59", "[60,70)" = "60-69", "[70,80)" = "70-79",
                              "[80,90)" = "80-90")) +
  theme(axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=1,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain")) +
  ylab("Count")
ggsave("itos_plots/Figure 4 - Gwm Age.png")

#### Gwm Graphs ####

ggplot(gwm.items, aes(x = time, y = result)) +
  geom_point() +
  xlab("Time Taken (secs)") +
  ylab("Total Score") +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = c(seq(90,960,30)), limits = c(90, 960)) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5))
ggsave("itos_plots/Figure 5 - Gwm Time Taken.png")

ggplot(gwm.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = time)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 6 - Gwm Time Boxplot.png")

ggplot(gwm.items) + 
  geom_boxplot(mapping = aes(x = factor(0), y = result)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("itos_plots/Figure 7 - Gwm Result Boxplot.png")

ggplot(gwm.items) +
  geom_boxplot(mapping = aes(x = gender, y = result)) +
  ylab("Total Score") +
  xlab("Gender")
ggsave("itos_plots/Figure 8 - Gwm Gender Result Boxplot.png")

ggplot(gwm.items) +
  geom_boxplot(mapping = aes(x = gender, y = time)) +
  ylab("Time Taken") +
  xlab("Gender")
ggsave("itos_plots/Figure 9 - Gwm Gender Time Taken Boxplot.png")

#### Gwm Gender Differences ####

t.test(gwm.items$result ~ gwm.items$gender)
t.test(gwm.items$time ~ gwm.items$gender)

#### Gwm Sample Setup ####

gwm.eval <- select(gwm.items,contains("c."))

#### Gwm Reliabiliy ####

psych::alpha(gwm.eval)

#### Count Number Correct ####

gwm.correct <- as.data.frame(
  sapply(gwm.eval, function(x) {
    (1 - (table(x)[1]/length(x))) * 100
  })
)

row.names(gwm.correct) <- str_remove(row.names(gwm.correct), "c.gwm")
row.names(gwm.correct) <- str_remove(row.names(gwm.correct), "\\.0")

names(gwm.correct) <- "Correct"

gwm.correct <- rownames_to_column(gwm.correct, var = "Item")

gwm.correct$Item <- as.numeric(as.vector(gwm.correct$Item))

ggplot(data = gwm.correct, aes(x = Item, y = Correct, group = 1)) +
  geom_line(linetype = "dashed") + 
  geom_point() +
  ylab("% Correct") + 
  scale_x_continuous(breaks=c(1:38)) + 
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave("itos_plots/Figure 10 - Gwm Number of Items Correct.png")

psych::describe(gwm.items$result)
psych::describe(gwm.items$time)

#### Gwm CFA ####

gwm.model <- paste("gwm =~ ", paste(colnames(gwm.eval), collapse = "+"), sep = "") # set the factor model
gwm.fit <- cfa(gwm.model, data=gwm.eval, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gwm.fit, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gwm.fit, "Gwm CFA 1") # Parameter Estimates for the CFA

#### Gwm Item Removal ####

gwm.eval.2 <- gwm.eval
gwm.eval.2$c.gwm8 <- NULL

#### Gwm CFA 2 ####

gwm.model.2 <- paste("gwm =~ ", paste(colnames(gwm.eval.2), collapse = "+"), sep = "") # set the factor model
gwm.fit.2 <- cfa(gwm.model.2, data=gwm.eval.2, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gwm.fit.2, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gwm.fit.2, "Gwm CFA 2") # Parameter Estimates for the CFA

#### Gwm Mokken ####

mokken.outcomes(gwm.eval.2, "Gwm Mokken 1")
gwm.aisp <- aisp.outcomes(gwm.eval.2, "Gwm AISP 1")

#### Gwm Item Removal 2 ####

gwm.eval.3 <- gwm.eval.2[, gwm.aisp == 1] # select items within the first Mokken scale

#### Gwm CFA 3 ####

gwm.model.3 <- paste("gwm =~ ", paste(colnames(gwm.eval.3), collapse = "+"), sep = "") # set the factor model
gwm.fit.3 <- cfa(gwm.model.3, data=gwm.eval.3, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gwm.fit.3, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gwm.fit.3, "Gwm CFA 3") # Parameter Estimates for the CFA

#### Gwm Mokken 2 ####

mokken.outcomes(gwm.eval.3, "Gwm Mokken 2")
gwm.aisp.2 <- aisp.outcomes(gwm.eval.3, "Gwm AISP 2")

#### Gwm IRT ####

rasch <- rasch.outcomes(gwm.eval.3, "Gwm Rasch 1")
gwm.rasch <- rasch[[1]]
gwm.rasch.fit <- rasch[[2]]

#### Gwm Local Dependency ####

residuals(gwm.rasch, type = "Q3", digits = 2, suppress = +.2)

#### Gwm Item Removal 3 ####

# Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for 
# dichotomous item response theory models. Applied Psychological Measurement, 24, 50-64. 
colnames(select(gwm.eval.3, as.vector(filter(gwm.rasch.fit, p.S_X2 < 0.01)$item))) # items that are going to get removed
colnames(select(gwm.eval.3, as.vector(filter(gwm.rasch.fit, p.S_X2 > 0.01)$item))) # items that are going to be kept
gwm.eval.4 <- select(gwm.eval.3, as.vector(filter(gwm.rasch.fit, p.S_X2 > 0.01)$item))

gwm.eval.4$c.gwm7 <- NULL # remove due to local dependency
gwm.eval.4$c.gwm12 <- NULL # remove due to local dependency
gwm.eval.4$c.gwm13 <- NULL # remove due to local dependency

#### Gwm CFA 4 ####

gwm.model.4 <- paste("gwm =~ ", paste(colnames(gwm.eval.4), collapse = "+"), sep = "") # set up factor model
gwm.fit.4 <- cfa(gwm.model.4, data=gwm.eval.4, std.lv=TRUE, missing="fiml") # Conduct the CFA
fitMeasures(gwm.fit.4, c("chisq", "df", "pvalue", "cfi","tli","rmsea", "srmr")) # Check the model is unidimensional
cfa.parameters(gwm.fit.4, "Gwm CFA 4") # Parameter Estimates for the CFA

#### Gwm Mokken 3 ####

mokken.outcomes(gwm.eval.4, "Gwm Mokken 3")
gwm.aisp.3 <- aisp.outcomes(gwm.eval.4, "Gwm AISP 3")

#### Gwm IRT 2 ####

rasch <- rasch.outcomes(gwm.eval.4, "Gwm Rasch 2")
gwm.rasch.2 <- rasch[[1]]
gwm.rasch.fit.2 <- rasch[[2]]

#### Gwm Local Dependency 3 ####

residuals(gwm.rasch.2, type = "Q3", digits = 2, suppress = +.2)

#### Gwm Reliability 2 ####

psych::alpha(gwm.eval.4)

#### Gwm Differential Item Functioning ####

sapply(names(gwm.eval.4), function(x) {
  table(gwm.eval.4[[x]],gwm.items[,"gender"])
})

gwm.eval.4$c.gwm2 <- NULL 
gwm.eval.4$c.gwm11 <- NULL

gwm.dif <- lordif(gwm.eval.4,gwm.items[,"gender"])
summary(gwm.dif)
plot(gwm.dif)

#### Gwm Item Characteristic Curves ####

jpeg(file = paste("itos_plots//Figure 11 - ICC of Gwm Items.jpeg"))
plot(gwm.rasch.2, type="trace", facet_items = FALSE, main ="")
dev.off()

#### Gwm Test Information Curve ####

jpeg(file = paste("itos_plots//Figure 12 - Gwm Test Information Curve.jpeg"))
plot(gwm.rasch.2, type = "info", xlim=c(-5,5), main = "")
dev.off()

#### Gwm Item Parameters ####

gwm.parameters <- coef(gwm.rasch.2, simplify=TRUE, IRTpars = TRUE, printSE = TRUE)
gwm.parameters$items
