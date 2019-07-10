### R Script for Analyses in *Supplementary Materials* for Manuscript Titled:
### Does Your Gain Define My Loss?:
### Socially-Defined Counterfactual Loss and Prevention-Focused Decision-Making
###
### Although this R script contains all analyses, the full supplementary materials
### (including the analyses contained here) can be generated using the R Markdown file 
### located on this Open Science Framework page:
### https://osf.io/e4nr3/?view_only=73c8490e19c6449197d2620d29ea8ad2

# Setup
set.seed(1234)

# Load Libraries and Define Functions
# All code used to load libraries and define functions can be found in the .Rmd file used to
# generate this PDF.

# For easy reference, installation code for packages used in this manuscript that are not available on CRAN:
# devtools::install_github("crsh/papaja") # Required to produce APA-formatted manuscript

library(knitr) # Required for knitting
library(papaja) # Required for APA template
library(citr) # Required for easy insertion of citations
library(tidyverse) # Required for data cleaning
library(broman) # Required for myround() function that doesn't truncate digits
library(formatR) # Required for data cleaning (throws error otherwise)
library(psych) # Required to calculate Cronbach's alpha
library(cocron) # Required to statistically compare Cronbach's alphas
library(broom) # Required to extract effects from models
library(rIP) # Required to filter out fraudulent M-Turk participants
library(ggplot2) # Required for plots
library(jtools) # Required for simple slopes and Johnson-Neyman analysis
library(cowplot) # Required for sim_slopes
library(Hmisc) # Required for correlation matrix
library(kableExtra) # Required for table styling
library(boot) # Required for bootstrapping CIs in mediation analysis
library(lavaan) # Required for mediation analyses
library(semPlot) # Required to plot mediation models

# Define functinos
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

pval <- function(x) {
  if (x >= .05) {
    result <- "n.s."
  }
  else if (x < .05 & x >= .01) {
    result <- "p < .05"
  }
  else if (x < .01 & x >= .001) {
    result <- "p < .01"
  }
  else {
    result <- "p < .001"
  }
  return(result)
}

## Import and Clean Data
# All code used to import and clean data can be found in the .Rmd file used to generate this PDF.
# Please note: To reproduce the reported results and ensure all code runs properly, the raw dataset
# should be exported from Qualtrics with the *numeric values* (**not** *choice text*) option selected.

# Loading this RData file so blockornot dataframe is loaded, as it can't be rerun without a user key:
load("data/rfcfloss_clean.RData")
# All other data is reloaded and recleaned below
cfl <- read.csv(file="data/rfcfloss_raw.csv", header = TRUE)

# Deleting two unnecessary additional header rows from Qualtrics export:
cfl <- cfl[-which(cfl$StartDate == "Start Date"), ] # Repeated header row with column names
cfl <- cfl[-1,] # Header row with Qualtrics code in curly brackets

# Renaming columns for ease of use
names(cfl) <- gsub(x = names(cfl), pattern = "\\.", replacement = "")
names(cfl) <- gsub(x = names(cfl), pattern = "rfq", replacement = "rfq_")
names(cfl) <- gsub(x = names(cfl), pattern = "rmq_rmq_", replacement = "rmq_")
colnames(cfl)[colnames(cfl) == "bit_easy"] <- "easy"
colnames(cfl)[colnames(cfl) == "bit_perf_clear"] <- "clear"
colnames(cfl)[colnames(cfl) == "bit_perf_affect"] <- "perfaffected"
colnames(cfl)[colnames(cfl) == "prev_bit"] <- "previnvested"

# Combining Bitcoin allocation (and relevant metadata) for two conditions
cfl$allo <- paste(cfl$bit_ctr_1,cfl$bit_cfl_11)
cfl$allo_timer <- paste(cfl$bit_timer_PageSubmit,cfl$bit_timer_PageSubmit1)
cfl$allo_clicks <- paste(cfl$bit_timer_ClickCount,cfl$bit_timer_ClickCount1)

# Ensuring proper data types for analysis
cfl <- droplevels(cfl) # Drop unneeded factor levels resulting from old Qualtrics headers
cfl.num <- c("Progress", "Durationinseconds", "rmq_1", "rmq_2", "rmq_3", "rmq_4", "rmq_5", "rmq_6", "rmq_7",
             "rmq_8", "rmq_9", "rmq_10", "rmq_11", "rmq_12", "rmq_13", "rmq_14", "rmq_15", "rmq_16", "rmq_17",
             "rmq_18", "rmq_19", "rmq_20", "rmq_21", "rmq_22", "rmq_23", "rmq_24", "rmq_25", "rmq_26", "rmq_27",
             "rmq_28", "rmq_29", "rmq_30", "rfq_1", "rfq_2", "rfq_3", "rfq_4", "rfq_5", "rfq_6",
             "rfq_7", "rfq_8", "rfq_9", "rfq_10", "rfq_11", "easy", "clear", "perfaffected",
             "missedout", "regret", "investearlier", "happier", "relieved", "age", "allo", "allo_timer", "allo_clicks") # Define numeric columns
cfl[cfl.num] <- lapply(cfl[cfl.num], as.character) # Convert numeric columns to character (before numeric) to prevent values from changing
cfl[cfl.num] <- lapply(cfl[cfl.num], as.numeric) # Convert numeric columns to numeric
rm(cfl.num) # Clear cfl.num from workspace as it is no longer needed

# Renaming and ordering factor levels where appropriate
cfl$gender <- as.character(cfl$gender)
cfl$gender[cfl$gender == "1"] <- "male"
cfl$gender[cfl$gender == "2"] <- "female"
cfl$gender[cfl$gender == "3"] <- "other"
cfl$gender <- as.factor(cfl$gender)

cfl$language <- as.character(cfl$language)
cfl$language[cfl$language == "1"] <- "english"
cfl$language[cfl$language == "2"] <- "spanish"
cfl$language[cfl$language == "3"] <- "other"
cfl$language <- as.factor(cfl$language)

cfl$ethnicity <- as.character(cfl$ethnicity)
cfl$ethnicity[!(cfl$ethnicity == "1" | cfl$ethnicity == "2" | cfl$ethnicity == "3" |
                  cfl$ethnicity == "4" | cfl$ethnicity == "5" | cfl$ethnicity == "6" |
                  cfl$ethnicity == "7" | is.na(cfl$ethnicity) == TRUE)] <- "multi"
cfl$ethnicity[cfl$ethnicity == "1"] <- "amerindian"
cfl$ethnicity[cfl$ethnicity == "2"] <- "asian"
cfl$ethnicity[cfl$ethnicity == "3"] <- "black"
cfl$ethnicity[cfl$ethnicity == "4"] <- "hispanic"
cfl$ethnicity[cfl$ethnicity == "5"] <- "pacislander"
cfl$ethnicity[cfl$ethnicity == "6"] <- "white"
cfl$ethnicity[cfl$ethnicity == "7"] <- "other"
cfl$ethnicity <- as.factor(cfl$ethnicity)

cfl$education <- as.character(cfl$education)
cfl$education[cfl$education == "0"] <- "none"
cfl$education[cfl$education == "1"] <- "elementary"
cfl$education[cfl$education == "2"] <- "somehs"
cfl$education[cfl$education == "3"] <- "highschool"
cfl$education[cfl$education == "4"] <- "somecollege"
cfl$education[cfl$education == "5"] <- "associates"
cfl$education[cfl$education == "6"] <- "bachelors"
cfl$education[cfl$education == "7"] <- "graddegree"
cfl$education[cfl$education == "8"] <- "doctorate"
cfl$education <- as.factor(cfl$education)
cfl$education <- factor(cfl$education, ordered = TRUE, levels = c("none", "elementary", "somehs", "highschool", "somecollege", "associates", "bachelors", "graddegree", "doctorate", "noanswer"))

cfl$income <- as.character(cfl$income)
cfl$income[cfl$income == "Prefer not to answer"] <- NA
cfl$income[cfl$income == "0"] <- "$0-$10K"
cfl$income[cfl$income == "1"] <- "$10K-$20K"
cfl$income[cfl$income == "2"] <- "$20K-$40K"
cfl$income[cfl$income == "3"] <- "$40K-$70K"
cfl$income[cfl$income == "4"] <- "$70K-$100K"
cfl$income[cfl$income == "5"] <- "$100K-$250K"
cfl$income[cfl$income == "6"] <- "$250K-$500K"
cfl$income[cfl$income == "7"] <- "$500K+"
cfl$income <- as.factor(cfl$income)
cfl$income <- factor(cfl$income, ordered = TRUE, levels = c("$0-$10K", "$10K-$20K", "$20K-$40K", "$40K-$70K", "$70K-$100K", "$100K-$250K", "$250K-$500K", "$500K+", "noanswer"))

cfl$condition <- as.character(cfl$condition)
cfl$condition[cfl$condition == "Prefer not to answer"] <- NA
cfl$condition[cfl$condition == "control"] <- "Control"
cfl$condition[cfl$condition == "cfloss"] <- "CF Loss"
cfl$condition <- as.factor(cfl$condition)
cfl$condition <- factor(cfl$condition, ordered = TRUE, levels = c("Control", "CF Loss"))

# Calculating prevention, promotion, assessment, and locomotion scores
cfl$prev <- ((6-cfl$rfq_2) + (6-cfl$rfq_4) + cfl$rfq_5 + (6-cfl$rfq_6) + (6-cfl$rfq_8))/5
cfl$prom <- ((6-cfl$rfq_1) + cfl$rfq_3 + cfl$rfq_7 + (6-cfl$rfq_9) + cfl$rfq_10 + (6-cfl$rfq_11))/6
cfl$ass <- ((7-cfl$rmq_2) + cfl$rmq_6 + cfl$rmq_7 + cfl$rmq_9 + (7-cfl$rmq_10) + cfl$rmq_11 + cfl$rmq_15 + cfl$rmq_19 + cfl$rmq_20 + cfl$rmq_22 + (7-cfl$rmq_27) + cfl$rmq_30)/12
cfl$loc <- (cfl$rmq_1 + cfl$rmq_3 + cfl$rmq_4 + cfl$rmq_5 + cfl$rmq_8 + cfl$rmq_16 + cfl$rmq_21 + cfl$rmq_25 + cfl$rmq_28 + cfl$rmq_29 + (7-cfl$rmq_13) + (7-cfl$rmq_24))/12

# Calculating SD within prevention, promotion, assessment, and locomotion scores
cfl$rfq_2r <- 6-cfl$rfq_2
cfl$rfq_4r <- 6-cfl$rfq_4
cfl$rfq_6r <- 6-cfl$rfq_6
cfl$rfq_8r <- 6-cfl$rfq_8
cfl$rfq_1r <- 6-cfl$rfq_1
cfl$rfq_9r <- 6-cfl$rfq_9
cfl$rfq_11r <- 6-cfl$rfq_11
cfl$rmq_2r <- 7-cfl$rmq_2
cfl$rmq_10r <- 7-cfl$rmq_10
cfl$rmq_27r <- 7-cfl$rmq_27
cfl$rmq_13r <- 7-cfl$rmq_13
cfl$rmq_24r <- 7-cfl$rmq_24
cfl$prev.var <- apply(cfl[, c("rfq_2r","rfq_4r","rfq_5","rfq_6r","rfq_8r")], 1, var)
cfl$prom.var <- apply(cfl[, c("rfq_1r","rfq_3","rfq_7","rfq_9r","rfq_10","rfq_11r")], 1, var)
cfl$ass.var <- apply(cfl[,
                         c("rmq_2r","rmq_6","rmq_7","rmq_9","rmq_10r","rmq_11",
                           "rmq_15","rmq_19","rmq_20","rmq_22","rmq_27r","rmq_30")], 1, var)
cfl$loc.var <- apply(cfl[,
                         c("rmq_1","rmq_3","rmq_4","rmq_5","rmq_8","rmq_16",
                           "rmq_21","rmq_25","rmq_28","rmq_29","rmq_13r","rmq_24r")], 1, var)

# Calculating prevention dominance scores and mean-centering them
cfl$prevdom <- cfl$prev - cfl$prom
cfl$prevdomc <- cfl$prevdom
cfl$prevdomc <- scale(cfl$prevdomc, center = TRUE, scale = FALSE)

# Calculating mean-centered prevention, promotion, assessment, and locomotion scores
cfl$premc <- cfl$prev
cfl$premc <- as.numeric(scale(cfl$premc, center = TRUE, scale = FALSE))
cfl$promc <- cfl$prom
cfl$promc <- as.numeric(scale(cfl$promc, center = TRUE, scale = FALSE))
cfl$amc <- cfl$ass
cfl$amc <- as.numeric(scale(cfl$amc, center = TRUE, scale = FALSE))
cfl$lmc <- cfl$loc
cfl$lmc <- as.numeric(scale(cfl$lmc, center = TRUE, scale = FALSE))

# Standardizing scores for potential mediators
cfl$perfaffected.s <- cfl$perfaffected
cfl$perfaffected.s <- scale(cfl$perfaffected.s)
cfl$missedout.s <- cfl$missedout
cfl$missedout.s <- scale(cfl$missedout.s)
cfl$regret.s <- cfl$regret
cfl$regret.s <- scale(cfl$regret.s)
cfl$investearlier.s <- cfl$investearlier
cfl$investearlier.s <- scale(cfl$investearlier.s)
cfl$happier.s <- cfl$happier
cfl$happier.s <- scale(cfl$happier.s)
cfl$relieved.s <- cfl$relieved
cfl$relieved.s <- scale(cfl$relieved.s)

# Creating dummy variables for high/low median splits by prevention, promotion, assesment, and locomotion
cfl$hiloprev <- ifelse(cfl$prev > median(cfl$prev, na.rm=TRUE),"hiprev","loprev")
cfl$hiloprom <- ifelse(cfl$prom > median(cfl$prom, na.rm=TRUE),"hiprom","loprom")
cfl$hiloass <- ifelse(cfl$ass > median(cfl$ass, na.rm=TRUE),"hiass","loass")
cfl$hiloloc <- ifelse(cfl$loc > median(cfl$loc, na.rm=TRUE),"hiloc","loloc")

# Dummy-coding experimental condition variable (for regression):
cfl$condD[cfl$condition == "CF Loss"] <- 1
cfl$condD[cfl$condition == "Control"] <- 0

# Standardizing variables of interest for mediation analyses:
cfl$premc.s <- scale(cfl$premc)
cfl$promc.s <- scale(cfl$promc)
cfl$allo.s <- scale(cfl$allo)
cfl$prevdomc.s <- scale(cfl$prevdomc)

# Effect-coding experimental condition variable (for mediation analyses):
cfl$condE[cfl$condition == "CF Loss"] <- 1
cfl$condE[cfl$condition == "Control"] <- -1

# Saving coded datasets
cflorig <- cfl
save(cflorig, file = "data/rfcfloss_coded.RData")

# Excluding those missing Bitcoin allocation (DV)
cflorig.alloexcl <- cflorig[ which(is.na(cflorig$allo) == FALSE ) , ]

# Excluding those missing prevention (one of two IVs)
cflorig.prevexcl <- cflorig.alloexcl[ which(is.na(cflorig.alloexcl$prev) == FALSE ) , ]
cflorig.ivdvexcl.n <- 500 - nrow(cflorig.prevexcl)

# Identifying participants with fraudulent IP addresses
cflexcl <- cflorig.prevexcl
# blockornot <- getIPinfo(cflexcl1, "IPAddress", "[key]") # Commented out to avoid rerunning
# # each time I run the full chunk, since script will not work without API key
# # Anyone looking to rerun the line of code above will need to secure a free X-key
# # from iphub.info and place it within the quotes

# blockornot <- dplyr::rename(blockornot, IPAddress = ip)
cflexcl <- merge(cflexcl, blockornot, by = "IPAddress")

# Determining participants flagged as block = (1 or 2) vs. 1 only, plus US only
cflexcl$block12 <- ifelse(cflexcl$block == "0", 0, 1)
cflexcl$block <- as.numeric(cflexcl$block)
cflexcl$block.d <- ifelse(cflexcl$block == 1, 1, 0)
cflexcl.blocked <- cflexcl[cflexcl$block.d == 0,] # n = 334
cflexcl.blocked.usonly <- cflexcl.blocked[cflexcl.blocked$countryCode == "US",] # n = 301
cfl.ipexcl <- cflexcl.blocked.usonly
cflexcl.allblocked.n <- nrow(cflorig.prevexcl) - nrow(cfl.ipexcl) # n = 180

# Comparing internal reliability of fraudulent vs. non-fraudulent responses
# Cronbach alpha when blocking 1s only
cflexcl.a <- select(cflexcl, rfq_2r, rfq_4r, rfq_5, rfq_6r, rfq_8r, block.d)
cflexcl.a.noblock <- subset(cflexcl.a, block.d == 0) # n = 334
cflexcl.a.noblock$block.d <- NULL # Removing column so I can calculate Cronbach's alpha below
cflexcl.a.block <- subset(cflexcl.a, block.d == 1) # n = 147
cflexcl.a.block$block.d <- NULL # Removing column so I can calculate Cronbach's alpha below
a.noblock <- psych::alpha(cflexcl.a.noblock)$total$std.alpha # Cronbach's alpha of unblocked group = 0.86
a.block <- psych::alpha(cflexcl.a.block)$total$std.alpha # Cronbach's alpha of blocked group = 0.55
feldttest <- cocron.two.coefficients(alpha = c(a.noblock,a.block), n = c(nrow(cflexcl.a.noblock),nrow(cflexcl.a.block)), dep = FALSE, alternative = "two.sided") # Significantly lower Cronbach's alpha in blocked group (p < .001)

# Comparing variance of fraudulent vs. non-fraudulent responses
# Variance when blocking 1s only
varblocked <- mean(cflexcl$prev.var[cflexcl$block == 1], na.rm = TRUE) # average prev var = 0.951
varunblocked <- mean(cflexcl$prev.var[cflexcl$block != 1], na.rm = TRUE) # average prev var = 0.609
vartest <- t.test(prev.var ~ block.d, data = cflexcl)

# Excluding fraudulent responses
cflexcl.blocked <- cflexcl[cflexcl$block.d == 0,] # Blocking flagged as "1"; n = 334
cflexcl.blocked.usonly <- cflexcl.blocked[cflexcl.blocked$countryCode == "US",] # n = 301
cfl.ipexcl <- cflexcl.blocked.usonly
cfl <- cflexcl.blocked.usonly
cfl.postip <- nrow(cfl)
cfl.excludedbyip <- nrow(cflexcl) - cfl.postip

## Determining participants who did not previously invest in Bitcoin
cfl.noprevinv <- subset(cfl.ipexcl, previnvested != 1)
cfl.noprevinv.pct <- myround((nrow(cfl.ipexcl)-nrow(cfl.noprevinv))/nrow(cfl.ipexcl) * 100, digits = 2)
# Under 10%, so will not exclude

## Determining participants who were not told about study from other M-Turkers
cfl.didntlearn <- subset(cfl.ipexcl, told_test != 1) # n = 295
cfl.didntlearn.n <- nrow(cfl.ipexcl)-nrow(cfl.didntlearn) # n = 6
cfl.didntlearn.pct <- myround((nrow(cfl.ipexcl)-nrow(cfl.didntlearn))/nrow(cfl.ipexcl) * 100, digits = 2)
# Over 10%, so will exclude in primary analyses (and will also rerun including this group)

# Redefining dataset to reflect exclusions (per pre-registration; to be edited
# if either excluded group is >10%)
# cfl <- subset(cfl, previnvested != 1) NOT EXCLUDED BECAUSE >10%
cfl.rerun <- cfl.ipexcl # n = 301
cfl <- subset(cfl.ipexcl, told_test != 1) # n = 295

n <- length(unique(cfl$ResponseId))

## Removing unnecessary columns from cleaned data
cfl$IPAddress <- NULL
cfl$Status <- NULL
cfl$RecipientFirstName <- NULL
cfl$RecipientLastName <- NULL
cfl$RecipientEmail <- NULL
cfl$LocationLatitude <- NULL
cfl$LocationLongitude <- NULL
cfl$ExternalReference <- NULL
cfl$DistributionChannel <- NULL
cfl$welcome <- NULL
cfl$MTurkCode <- NULL
cfl$RecordedDate <- NULL

## Creating subset to facilitate visualization of results within each condition
cfl.cflcond <- subset(cfl, condition == "CF Loss") # n = 142
cfl.ctrlcond <- subset(cfl, condition == "Control") # n = 153

## Saving cleaned dataframes
save(cfl, blockornot, cflorig, cflexcl1, cflexcl, cfl.cflcond, cfl.ctrlcond, file = "data/rfcfloss_clean.RData")
write.csv(cfl, file = "data/rfcfloss_clean.csv")

## Calculating descriptive statistics for demographic variables
meanage <- mean(cfl$age, na.rm = TRUE)
sdage <- sd(cfl$age, na.rm = TRUE)
maxage <- max(cfl$age, na.rm = TRUE)
minage <- min(cfl$age, na.rm = TRUE)
male <- sum(cfl$gender == "male")
female <- sum(cfl$gender == "female")
malepct <- myround(male/nrow(cfl)*100, digits = 2)
femalepct <- myround(female/nrow(cfl)*100, digits = 2)
english <- sum(cfl$language == "english")
engpct <- myround(english/nrow(cfl)*100, digits = 2)
white <- sum(cfl$ethnicity == "white")
whitepct <- myround(white/nrow(cfl)*100, digits = 2)
hisp <- sum(cfl$ethnicity == "hispanic")
hisppct <- myround(hisp/nrow(cfl)*100, digits = 2)
black <- sum(cfl$ethnicity == "black")
blackpct <- myround(black/nrow(cfl)*100, digits = 2)
asian <- sum(cfl$ethnicity == "asian")
asianpct <- myround(asian/nrow(cfl)*100, digits = 2)
cfl$income[is.na(cfl$income)] <- "noanswer"
tentwenty <- sum(cfl$income == "$10K-$20K")
tentwentypct <- myround(tentwenty/nrow(cfl)*100, digits = 2)
twentyforty <- sum(cfl$income == "$20K-$40K")
twentyfortypct <- myround(twentyforty/nrow(cfl)*100, digits = 2)
fortyseventy <- sum(cfl$income == "$40K-$70K")
fortyseventypct <- myround(fortyseventy/nrow(cfl)*100, digits = 2)
seventyhundred <- sum(cfl$income == "$70K-$100K")
seventyhundredpct <- myround(seventyhundred/nrow(cfl)*100, digits = 2)
hundredtwofifty <- sum(cfl$income == "$100K-$250K")
hundredtwofiftypct <- myround(hundredtwofifty/nrow(cfl)*100, digits = 2)
twofiftyfive <- sum(cfl$income == "$250K-$500K")
twofiftyfivepct <- myround(twofiftyfive/nrow(cfl)*100, digits = 2)
fiveplus <- sum(cfl$income == "$500K+")
fivepluspct <- myround(fiveplus/nrow(cfl)*100, digits = 2)
cfl$income[cfl$income == "noanswer"] <- NA
cfl$education[is.na(cfl$education)] <- "noanswer"
elementary <- sum(cfl$education == "elementary")
elementarypct <- myround(elementary/nrow(cfl)*100, digits = 2)
somehs <- sum(cfl$education == "somehs")
somehspct <- myround(somehs/nrow(cfl)*100, digits = 2)
highschool <- sum(cfl$education == "highschool")
highschoolpct <- myround(highschool/nrow(cfl)*100, digits = 2)
somecollege <- sum(cfl$education == "somecollege")
somecollegepct <- myround(somecollege/nrow(cfl)*100, digits = 2)
associates <- sum(cfl$education == "associates")
associatespct <- myround(associates/nrow(cfl)*100, digits = 2)
bachelors <- sum(cfl$education == "bachelors")
bachelorspct <- myround(bachelors/nrow(cfl)*100, digits = 2)
graddegree <- sum(cfl$education == "graddegree")
graddegreepct <- myround(graddegree/nrow(cfl)*100, digits = 2)
doctorate <- sum(cfl$education == "doctorate")
doctoratepct <- myround(doctorate/nrow(cfl)*100, digits = 2)
cfl$education[cfl$education == "noanswer"] <- NA

# Descriptive Plots

## Participant Age
ggplot(data=cfl, aes(x = age)) +
  geom_histogram(bins = 10) +
  labs(title="Histogram of Participant Age", x="Age", y="Participants")

## Participant Ethnicity
cfl$ethnicity.ord <- factor(cfl$ethnicity, ordered = TRUE, levels = c("white", "black", "asian", "multi", "hispanic"))
cfl$ethnicity.ord <- droplevels(cfl$ethnicity.ord)
ggplot(data=subset(cfl, !is.na(ethnicity.ord)), aes(x = ethnicity.ord)) +
  scale_x_discrete(labels = c("White", "Black", "Asian", "Multiracial", "Hispanic")) +
  geom_histogram(stat = "count") +
  labs(title="Histogram of Participant Ethnicity", x="Ethnicity", y="Participants")

## Participant Income
ggplot(subset(cfl, !is.na(income)), aes(x = income)) +
  scale_x_discrete(labels = c("$10K-\n$20K", "$20K-\n$40K", "$40K-\n$70K", "$70K-\n$100K", "$100K-\n$250K", "$250K-\n$500K", "$500K+")) +
  geom_bar(stat = "count") +
  labs(title = "Histogram of Participant Income", x = "Income", y = "Participants")

## Participant Education
ggplot(subset(cfl, !is.na(education)), aes(x = education)) +
  scale_x_discrete(labels = c("High\nSchool", "Some\nCollege", "Associate's\nDegree", "Bachelor's\nDegree", "Graduate\nDegree", "Doctorate")) +
  geom_bar(stat = "count") +
  labs(title = "Histogram of Participant Education", x = "Education", y = "Participants")

## Bitcoin Allocation by Prevention Pride and Condition
# summarySE() function defined at top of document
cflsummary <- summarySE(cfl, measurevar="allo", groupvars=c("hiloprev","condition"))
ggplot(data=cfl, aes(x = hiloprev, y = allo, color = hiloprev)) +
  facet_wrap(. ~ condition) +
  geom_jitter() +
  geom_errorbar(data = cflsummary, aes(ymin = allo - se, ymax = allo + se),
                , col = "black", size = .5) +
  geom_point(data = cflsummary, col = "black", size = 3, alpha = 1) +
  scale_x_discrete(labels = c("High", "Low")) +
  labs(title="Bitcoin Allocation By Counterfactual Loss and Prevention Pride\n(NOTE: Does Not Control for Effects of Promotion Pride)", x="Prevention Pride (Median-Split)", y="Bitcoin Allocation (%)") +
  theme(legend.position = "none")

# Primary linear analysis and table generation
rfcfl.promctr.nonstd <- lm(allo ~ premc * condD + promc * condD, data = cfl)
table1 <- tidy(rfcfl.promctr.nonstd)
table1 <- dplyr::rename(table1, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1$Estimate <- myround(table1$Estimate, digits = 2)
table1$SE <- myround(table1$SE, digits = 2)
table1$t <- myround(table1$t, digits = 2)
table1$p <- myround(table1$p, digits = 3)
table1$p[table1$p < .001] <- "< .001"
table1$Predictor[table1$Predictor == "(Intercept)"] <- "Intercept"
table1$Predictor[table1$Predictor == "premc"] <- "Prevention Pride"
table1$Predictor[table1$Predictor == "condD"] <- "Counterfactual Loss"
table1$Predictor[table1$Predictor == "promc"] <- "Promotion Pride"
table1$Predictor[table1$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1$Predictor[table1$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1[is.na(table1)] <- ""
table1$order <- c(1, 2, 3, 5, 4, 6)
table1 <- arrange(table1, order)
table1$order <- NULL
df <- rfcfl.promctr.nonstd$df.residual
apa_table(table1, caption = "Summary of Linear Regression Analysis", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')

# Plot of primary model
pred.data.ctrl <- data.frame(premc = seq(min(cfl$premc, na.rm=T), 
                                         max(cfl$premc, na.rm=T), .1),
                             promc = 0, condD = 0, condition = "Control")
rfcfl.promctr.pred.ctrl <- cbind(pred.data.ctrl, 
                                 predict(rfcfl.promctr.nonstd, pred.data.ctrl, 
                                         interval = "confidence", 
                                         type = c("response", "terms")))
rfcfl.promctr.pred.ctrl <- dplyr::rename(rfcfl.promctr.pred.ctrl, allo = fit)
pred.data.CFL <- data.frame(premc = seq(min(cfl$premc, na.rm=T), 
                                        max(cfl$premc, na.rm=T), .1),
                            promc = 0, condD = 1, condition = "CF Loss")
rfcfl.promctr.pred.CFL <- cbind(pred.data.CFL, 
                                predict(rfcfl.promctr.nonstd, pred.data.CFL, 
                                        interval = "confidence", 
                                        type = c("response", "terms")))
rfcfl.promctr.pred.CFL <- dplyr::rename(rfcfl.promctr.pred.CFL, allo = fit)
primarymodel <- ggplot(data=cfl, aes(x=premc, y=allo)) +
  scale_color_manual(values=c("gray0", "gray50")) +
  geom_ribbon(data = rfcfl.promctr.pred.ctrl, aes(ymin = lwr, ymax = upr), alpha = .3, fill = "gray") +
  geom_ribbon(data = rfcfl.promctr.pred.CFL, aes(ymin = lwr, ymax = upr), alpha = .3, fill = "gray") +
  geom_line(data = rfcfl.promctr.pred.ctrl, aes(color=condition), size = 1) +
  geom_line(data = rfcfl.promctr.pred.CFL, aes(color=condition), size = 1) +
  labs(title="Prevention Pride and Counterfactual Loss as\nPredictors of Bitcoin Allocation\n(Controlling for Promotion Pride and the Interaction\nBetween Promotion Pride and Counterfactual Loss)", x="Prevention Pride (Mean-Centered)", y="Bitcoin Allocation (%)", color = "Condition")
primarymodel

## Hypothesized Model Without Promotion Pride and With Demographic Covariates

# Not controlling for promotion pride or interaction between promotion pride and CF Loss (mean-centered RF and dummy-coded condition)
rfcfl <- lm(allo ~ premc * condD, data=cfl)
table1.noprom <- tidy(rfcfl)
table1.noprom <- dplyr::rename(table1.noprom, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.noprom$Estimate <- myround(table1.noprom$Estimate, digits = 2)
table1.noprom$p <- myround(table1.noprom$p, digits = 3)
table1.noprom$p[table1.noprom$p < .001] <- "< .001"
table1.noprom$Predictor[table1.noprom$Predictor == "(Intercept)"] <- "Intercept"
table1.noprom$Predictor[table1.noprom$Predictor == "premc"] <- "Prevention Pride"
table1.noprom$Predictor[table1.noprom$Predictor == "condD"] <- "Counterfactual Loss"
table1.noprom$Predictor[table1.noprom$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.noprom[is.na(table1.noprom)] <- ""

# Not controlling for promotion pride (mean-centered RF and dummy-coded condition)
rfcfl.promctrnointxn <- lm(allo ~ premc * condD + promc, data=cfl)
table1.promctrnointxn <- tidy(rfcfl.promctrnointxn)
table1.promctrnointxn <- dplyr::rename(table1.promctrnointxn, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.promctrnointxn$Estimate <- myround(table1.promctrnointxn$Estimate, digits = 2)
table1.promctrnointxn$p <- myround(table1.promctrnointxn$p, digits = 3)
table1.promctrnointxn$p[table1.promctrnointxn$p < .001] <- "< .001"
table1.promctrnointxn$Predictor[table1.promctrnointxn$Predictor == "(Intercept)"] <- "Intercept"
table1.promctrnointxn$Predictor[table1.promctrnointxn$Predictor == "premc"] <- "Prevention Pride"
table1.promctrnointxn$Predictor[table1.promctrnointxn$Predictor == "condD"] <- "Counterfactual Loss"
table1.promctrnointxn$Predictor[table1.promctrnointxn$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.promctrnointxn$Predictor[table1.promctrnointxn$Predictor == "promc"] <- "Promotion Pride"
table1.promctrnointxn[is.na(table1.promctrnointxn)] <- ""

# Controlling for gender (mean-centered RF and dummy-coded condition)
rfcfl.promctr.gender <- lm(allo ~ premc * condD + gender + promc * condD, data=cfl)
table1.gender <- tidy(rfcfl.promctr.gender)
table1.gender <- dplyr::rename(table1.gender, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.gender$Estimate <- myround(table1.gender$Estimate, digits = 2)
table1.gender$p <- myround(table1.gender$p, digits = 3)
table1.gender$p[table1.gender$p < .001] <- "< .001"
table1.gender$Predictor[table1.gender$Predictor == "(Intercept)"] <- "Intercept"
table1.gender$Predictor[table1.gender$Predictor == "premc"] <- "Prevention Pride"
table1.gender$Predictor[table1.gender$Predictor == "condD"] <- "Counterfactual Loss"
table1.gender$Predictor[table1.gender$Predictor == "promc"] <- "Promotion Pride"
table1.gender$Predictor[table1.gender$Predictor == "gendermale"] <- "Gender: Male"
table1.gender$Predictor[table1.gender$Predictor == "genderfemale"] <- "Gender: Female"
table1.gender$Predictor[table1.gender$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.gender$Predictor[table1.gender$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1.gender[is.na(table1.gender)] <- ""

# Controlling for age (mean-centered RF and dummy-coded condition)
rfcfl.promctr.age <- lm(allo ~ premc * condD + age + promc * condD, data=cfl)
table1.age <- tidy(rfcfl.promctr.age)
table1.age <- dplyr::rename(table1.age, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.age$Estimate <- myround(table1.age$Estimate, digits = 2)
table1.age$p <- myround(table1.age$p, digits = 3)
table1.age$p[table1.age$p < .001] <- "< .001"
table1.age$Predictor[table1.age$Predictor == "(Intercept)"] <- "Intercept"
table1.age$Predictor[table1.age$Predictor == "premc"] <- "Prevention Pride"
table1.age$Predictor[table1.age$Predictor == "condD"] <- "Counterfactual Loss"
table1.age$Predictor[table1.age$Predictor == "promc"] <- "Promotion Pride"
table1.age$Predictor[table1.age$Predictor == "age"] <- "Age"
table1.age$Predictor[table1.age$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.age$Predictor[table1.age$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1.age[is.na(table1.age)] <- ""

# Controlling for ethnicity (mean-centered RF and dummy-coded condition)
rfcfl.promctr.ethnicity <- lm(allo ~ premc * condD + ethnicity + promc * condD, data=cfl)
table1.ethnicity <- tidy(rfcfl.promctr.ethnicity)
table1.ethnicity <- dplyr::rename(table1.ethnicity, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.ethnicity$Estimate <- myround(table1.ethnicity$Estimate, digits = 2)
table1.ethnicity$p <- myround(table1.ethnicity$p, digits = 3)
table1.ethnicity$p[table1.ethnicity$p < .001] <- "< .001"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "(Intercept)"] <- "Intercept"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "premc"] <- "Prevention Pride"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "condD"] <- "Counterfactual Loss"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "promc"] <- "Promotion Pride"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicityasian"] <- "Ethnicity: Asian"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicityblack"] <- "Ethnicity: Black"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicityhispanic"] <- "Ethnicity: Hispanic"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicitymulti"] <- "Ethnicity: Multiracial"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicityother"] <- "Ethnicity: Other"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicitypacislander"] <- "Ethnicity: Pacific Islander"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "ethnicitywhite"] <- "Ethnicity: White"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1.ethnicity[is.na(table1.ethnicity)] <- ""

# Controlling for income (mean-centered RF and dummy-coded condition)
cfl$income.unord <- factor(cfl$income, ordered = F)
rfcfl.promctr.income <- lm(allo ~ premc * condD + income.unord + promc * condD, data=cfl)
table1.income <- tidy(rfcfl.promctr.income)
table1.income <- dplyr::rename(table1.income, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.income$Estimate <- myround(table1.income$Estimate, digits = 2)
table1.income$p <- myround(table1.income$p, digits = 3)
table1.income$p[table1.income$p < .001] <- "< .001"
table1.income$Predictor[table1.income$Predictor == "(Intercept)"] <- "Intercept"
table1.income$Predictor[table1.income$Predictor == "premc"] <- "Prevention Pride"
table1.income$Predictor[table1.income$Predictor == "condD"] <- "Counterfactual Loss"
table1.income$Predictor[table1.income$Predictor == "promc"] <- "Promotion Pride"
table1.income$Predictor[table1.income$Predictor == "income.unord$20K-$40K"] <- "Income: $20K-$40K"
table1.income$Predictor[table1.income$Predictor == "income.unord$40K-$70K"] <- "Income: $40K-$70K"
table1.income$Predictor[table1.income$Predictor == "income.unord$70K-$100K"] <- "Income: $70K-$100K"
table1.income$Predictor[table1.income$Predictor == "income.unord$100K-$250K"] <- "Income: $100K-$250K"
table1.income$Predictor[table1.income$Predictor == "income.unord$250K-$500K"] <- "Income: $250K-$500K"
table1.income$Predictor[table1.income$Predictor == "income.unord$500K+"] <- "Income: $500K+"
table1.income$Predictor[table1.income$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.income$Predictor[table1.income$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1.income[is.na(table1.income)] <- ""

# Controlling for education (mean-centered RF and dummy-coded condition)
cfl$edu.unord <- factor(cfl$edu, ordered = F)
rfcfl.promctr.edu <- lm(allo ~ premc * condD + edu.unord + promc * condD, data=cfl)
table1.edu <- tidy(rfcfl.promctr.edu)
table1.edu <- dplyr::rename(table1.edu, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.edu$Estimate <- myround(table1.edu$Estimate, digits = 2)
table1.edu$p <- myround(table1.edu$p, digits = 3)
table1.edu$p[table1.edu$p < .001] <- "< .001"
table1.edu$Predictor[table1.edu$Predictor == "(Intercept)"] <- "Intercept"
table1.edu$Predictor[table1.edu$Predictor == "premc"] <- "Prevention Pride"
table1.edu$Predictor[table1.edu$Predictor == "condD"] <- "Counterfactual Loss"
table1.edu$Predictor[table1.edu$Predictor == "promc"] <- "Promotion Pride"
table1.edu$Predictor[table1.edu$Predictor == "edu.unordsomecollege"] <- "Education: Some College"
table1.edu$Predictor[table1.edu$Predictor == "edu.unordassociates"] <- "Education: Associate's"
table1.edu$Predictor[table1.edu$Predictor == "edu.unordbachelors"] <- "Education: Bachelor's"
table1.edu$Predictor[table1.edu$Predictor == "edu.unordgraddegree"] <- "Education: Grad. Degree"
table1.edu$Predictor[table1.edu$Predictor == "edu.unorddoctorate"] <- "Education: Doctorate"
table1.edu$Predictor[table1.edu$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.edu$Predictor[table1.edu$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1.edu[is.na(table1.edu)] <- ""

# Code to print tables in PDF from R Markdown
apa_table(table1.noprom, caption = "Summary of Linear Regression Analysis Not Controlling for Promotion Pride or the Interaction Between Promotion Pride and Counterfactual Loss", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention pride was mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
apa_table(table1.promctrnointxn, caption = "Summary of Linear Regression Analysis Not Controlling for the Interaction Between Promotion Pride and Counterfactual Loss", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
apa_table(table1.gender, caption = "Summary of Linear Regression Analysis Controlling for Gender", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
apa_table(table1.age, caption = "Summary of Linear Regression Analysis Controlling for Age", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
apa_table(table1.ethnicity, caption = "Summary of Linear Regression Analysis Controlling for Ethnicity", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
apa_table(table1.income, caption = "Summary of Linear Regression Analysis Controlling for Income", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
apa_table(table1.edu, caption = "Summary of Linear Regression Analysis Controlling for Education", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')

# Exploratory Mediation Analyses (Hayes PROCESS Model 8): Counterfactual Thought and Hypothetical Emotion Mediators

# Descriptive statistics for potential mediators
missedoutmean <- round(mean(cfl$missedout, na.rm = TRUE), digits = 2)
missedoutsd <- round(sd(cfl$missedout, na.rm = TRUE), digits = 2)
regretmean <- round(mean(cfl$regret, na.rm = TRUE), digits = 2)
regretsd <- round(sd(cfl$regret, na.rm = TRUE), digits = 2)
happiermean <- round(mean(cfl$happier, na.rm = TRUE), digits = 2)
happiersd <- round(sd(cfl$happier, na.rm = TRUE), digits = 2)
relievedmean <- round(mean(cfl$relieved, na.rm = TRUE), digits = 2)
relievedsd <- round(sd(cfl$relieved, na.rm = TRUE), digits = 2)

## Exploratory Moderated Mediation Analysis: Missing Out

# Specify model
Mod.Med.Lavaan.promctr.Missedout <- '
# Regressions
missedout.s ~ 1 + a1*premc.s + a2*condE + a3*premc.s:condE + a4*promc.s + a5*condE:promc.s
allo.s ~ 1 + cdash1*premc.s + cdash2*condE + b1*missedout.s + cdash3*premc.s:condE + cdash4*promc.s + cdash5*condE:promc.s

#Indirect effects conditional on moderator (a1 + a3*ModValue)*b1
indirect.ctrl := (a1 + a3*(-1))*b1
indirect.CFL := (a1 + a3*(1))*b1

#Direct effects conditional on moderator (cdash1 + cdash3*ModValue)
direct.ctrl := cdash1 + cdash3*(-1)
direct.CFL := cdash1 + cdash3*(1)

#Total effects conditional on moderator
total.ctrl := direct.ctrl + indirect.ctrl
total.CFL := direct.CFL + indirect.CFL

#Proportion mediated conditional on moderator
#To match the output of "mediate" package
prop.mediated.ctrl := indirect.ctrl / total.ctrl
prop.mediated.CFL := indirect.CFL / total.CFL

#Index of moderated mediation
#An alternative way of testing if conditional indirect effects are significantly different from each other
index.mod.med := a3*b1
'

# Fit model - Commented out because this script takes a long time to rerun.
# set.seed(1234)
# Mod.Med.SEM.promctr.missedout <- sem(model = Mod.Med.Lavaan.promctr.Missedout,
#                                      data = cfl,
#                                      se = "bootstrap",
#                                      bootstrap = 5000)
# save(Mod.Med.SEM.promctr.missedout, file = "models/missedoutPromctrFit.RData")

load("models/missedoutPromctrFit.RData")
summary(Mod.Med.SEM.promctr.missedout,
        fit.measures = FALSE,
        standardize = TRUE,
        estimates = TRUE,
        ci = TRUE,
        rsquare = TRUE)

missedoutPromctrCoefs <- parameterEstimates(Mod.Med.SEM.promctr.missedout)
a3.missedout.promctr <- missedoutPromctrCoefs$est[missedoutPromctrCoefs$label == "a3"]
a3.p.missedout.promctr <- missedoutPromctrCoefs$pvalue[missedoutPromctrCoefs$label == "a3"]
b.missedout.promctr <- missedoutPromctrCoefs$est[missedoutPromctrCoefs$label == "b1"]
b.p.missedout.promctr <- missedoutPromctrCoefs$pvalue[missedoutPromctrCoefs$label == "b1"]
condIndEffctrl.missedout.promctr <- missedoutPromctrCoefs$est[missedoutPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.se.missedout.promctr <- missedoutPromctrCoefs$se[missedoutPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.cilower.missedout.promctr <- missedoutPromctrCoefs$ci.lower[missedoutPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.ciupper.missedout.promctr <- missedoutPromctrCoefs$ci.upper[missedoutPromctrCoefs$label == "indirect.ctrl"]
condIndEffCFL.missedout.promctr <- missedoutPromctrCoefs$est[missedoutPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.se.missedout.promctr <- missedoutPromctrCoefs$se[missedoutPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.cilower.missedout.promctr <- missedoutPromctrCoefs$ci.lower[missedoutPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.ciupper.missedout.promctr <- missedoutPromctrCoefs$ci.upper[missedoutPromctrCoefs$label == "indirect.CFL"]

# Generating table of moderated mediation results
table2.missedout <- missedoutPromctrCoefs 
table2.missedout <- dplyr::rename(table2.missedout, Predictor = rhs, Estimate = est, SE = se, p = pvalue, CI.lower = ci.lower, CI.upper = ci.upper)
table2.missedout$label <- NULL
table2.missedout <- table2.missedout[table2.missedout$op != "~~",]
table2.missedout$Predictor[table2.missedout$lhs == "missedout.s" & table2.missedout$op == "~1"] <- "Intercept"
table2.missedout$op[table2.missedout$lhs == "missedout.s" & table2.missedout$op == "~1"] <- "~"
table2.missedout$Predictor[table2.missedout$lhs == "allo.s" & table2.missedout$op == "~1"] <- "Intercept"
table2.missedout$op[table2.missedout$lhs == "allo.s" & table2.missedout$op == "~1"] <- "~"
table2.missedout <- table2.missedout[table2.missedout$op != "~1",]
table2.missedout <- table2.missedout[1:15,]
table2.missedout <- table2.missedout %>%
mutate(section = recode(lhs, missedout.s = 1, allo.s = 2,
indirect.ctrl = 3, indirect.CFL = 4)) %>%
arrange(section,desc(op))
table2.missedout$Predictor[table2.missedout$lhs == "indirect.ctrl"] <- "Control Condition"
table2.missedout$Predictor[table2.missedout$lhs == "indirect.CFL"] <- "Counterfactual Loss Condition"
table2.missedout$Estimate <- myround(table2.missedout$Estimate, digits = 2)
table2.missedout$SE <- myround(table2.missedout$SE, digits = 2)
table2.missedout$z <- myround(table2.missedout$z, digits = 2)
table2.missedout$p <- myround(table2.missedout$p, digits = 3)
table2.missedout$CI.lower <- myround(table2.missedout$CI.lower, digits = 4)
table2.missedout$CI.upper <- myround(table2.missedout$CI.upper, digits = 4)
table2.missedout$p[table2.missedout$p < .001] <- "< .001"
table2.missedout$Predictor[table2.missedout$Predictor == "premc.s"] <- "Prevention Pride"
table2.missedout$Predictor[table2.missedout$Predictor == "condE"] <- "Counterfactual Loss"
table2.missedout$Predictor[table2.missedout$Predictor == "promc.s"] <- "Promotion Pride"
table2.missedout$Predictor[table2.missedout$Predictor == "premc.s:condE"] <- "Prev. Pride x CF Loss"
table2.missedout$Predictor[table2.missedout$Predictor == "condE:promc.s"] <- "Prom. Pride x CF Loss"
table2.missedout$Predictor[table2.missedout$Predictor == "missedout.s"] <- "Missed Out"
table2.missedout[is.na(table2.missedout)] <- ""
table2.missedout$lhs <- NULL
table2.missedout$op <- NULL
table2.missedout$section <- NULL
table2.missedout$order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 10, 11, 12, 14, 15)
table2.missedout <- arrange(table2.missedout, order)
table2.missedout$order <- NULL

medmodel.missedout <- table2.missedout %>%
  mutate(Estimate = cell_spec(Estimate, "latex", align = "r"),
         SE = cell_spec(SE, "latex", align = "r"),
         z = cell_spec(z, "latex", align = "r"),
         p = cell_spec(p, "latex", align = "r"),
         CI.lower = cell_spec(CI.lower, "latex", align = "r"),
         CI.upper = cell_spec(CI.upper, "latex", align = "r")) %>%
  kable("latex", booktabs = TRUE, escape = FALSE, caption = "Summary of Exploratory Mediation Analysis") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  group_rows("Model 1 (DV = Missed Out)", 1, 6) %>%
  group_rows("Model 2 (DV = Bitcoin Allocation)", 7, 13) %>%
  group_rows("Bootstrapped Conditional Indirect Effects\n(Prev. Pride x CF Loss → Missed Out → Bitcoin Allocation)", 14, 15) %>%
  row_spec(0, align = "c") %>%
  footnote(general = "This analysis included an effect-coded variable for the counterfactual loss experience: -1 = control, 1 = counterfactual loss. Additionally, all other variables were standardized (M = 0, SD = 1) within this model. Estimated effect sizes for the models reported here are standardized regression coefficients.", general_title = "Note: ", footnote_as_chunk = T, title_format = c("italic"), threeparttable = TRUE, escape = FALSE)
medmodel.missedout

## Exploratory Moderated Mediation Analysis: Regret

# Specify model
Mod.Med.Lavaan.promctr.regret <- '
# Regressions
regret.s ~ 1 + a1*premc.s + a2*condE + a3*premc.s:condE + a4*promc.s + a5*condE:promc.s
allo.s ~ 1 + cdash1*premc.s + cdash2*condE + b1*regret.s + cdash3*premc.s:condE + cdash4*promc.s + cdash5*condE:promc.s

#Indirect effects conditional on moderator (a1 + a3*ModValue)*b1
indirect.ctrl := (a1 + a3*(-1))*b1
indirect.CFL := (a1 + a3*(1))*b1

#Direct effects conditional on moderator (cdash1 + cdash3*ModValue)
direct.ctrl := cdash1 + cdash3*(-1)
direct.CFL := cdash1 + cdash3*(1)

#Total effects conditional on moderator
total.ctrl := direct.ctrl + indirect.ctrl
total.CFL := direct.CFL + indirect.CFL

#Proportion mediated conditional on moderator
#To match the output of "mediate" package
prop.mediated.ctrl := indirect.ctrl / total.ctrl
prop.mediated.CFL := indirect.CFL / total.CFL

#Index of moderated mediation
#An alternative way of testing if conditional indirect effects are significantly different from each other
index.mod.med := a3*b1
'

# Fit model - Commented out because this script takes a long time to rerun.
# set.seed(1234)
# Mod.Med.SEM.promctr.regret <- sem(model = Mod.Med.Lavaan.promctr.regret,
#                                   data = cfl,
#                                   se = "bootstrap",
#                                   bootstrap = 5000)
# save(Mod.Med.SEM.promctr.regret, file = "models/regretPromctrFit.RData")

load("models/regretPromctrFit.RData")
summary(Mod.Med.SEM.promctr.regret,
        fit.measures = FALSE,
        standardize = TRUE,
        estimates = TRUE,
        ci = TRUE,
        rsquare = TRUE)

regretPromctrCoefs <- parameterEstimates(Mod.Med.SEM.promctr.regret)
a3.regret.promctr <- regretPromctrCoefs$est[regretPromctrCoefs$label == "a3"]
a3.p.regret.promctr <- regretPromctrCoefs$pvalue[regretPromctrCoefs$label == "a3"]
b.regret.promctr <- regretPromctrCoefs$est[regretPromctrCoefs$label == "b1"]
b.p.regret.promctr <- regretPromctrCoefs$pvalue[regretPromctrCoefs$label == "b1"]
condIndEffctrl.regret.promctr <- regretPromctrCoefs$est[regretPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.se.regret.promctr <- regretPromctrCoefs$se[regretPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.cilower.regret.promctr <- regretPromctrCoefs$ci.lower[regretPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.ciupper.regret.promctr <- regretPromctrCoefs$ci.upper[regretPromctrCoefs$label == "indirect.ctrl"]
condIndEffCFL.regret.promctr <- regretPromctrCoefs$est[regretPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.se.regret.promctr <- regretPromctrCoefs$se[regretPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.cilower.regret.promctr <- regretPromctrCoefs$ci.lower[regretPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.ciupper.regret.promctr <- regretPromctrCoefs$ci.upper[regretPromctrCoefs$label == "indirect.CFL"]

# Generating table of moderated mediation results
table2.regret <- regretPromctrCoefs 
table2.regret <- dplyr::rename(table2.regret, Predictor = rhs, Estimate = est, SE = se, p = pvalue, CI.lower = ci.lower, CI.upper = ci.upper)
table2.regret$label <- NULL
table2.regret <- table2.regret[table2.regret$op != "~~",]
table2.regret$Predictor[table2.regret$lhs == "regret.s" & table2.regret$op == "~1"] <- "Intercept"
table2.regret$op[table2.regret$lhs == "regret.s" & table2.regret$op == "~1"] <- "~"
table2.regret$Predictor[table2.regret$lhs == "allo.s" & table2.regret$op == "~1"] <- "Intercept"
table2.regret$op[table2.regret$lhs == "allo.s" & table2.regret$op == "~1"] <- "~"
table2.regret <- table2.regret[table2.regret$op != "~1",]
table2.regret <- table2.regret[1:15,]
table2.regret <- table2.regret %>%
  mutate(section = recode(lhs, regret.s = 1, allo.s = 2,
                          indirect.ctrl = 3, indirect.CFL = 4)) %>%
  arrange(section,desc(op))
table2.regret$Predictor[table2.regret$lhs == "indirect.ctrl"] <- "Control Condition"
table2.regret$Predictor[table2.regret$lhs == "indirect.CFL"] <- "Counterfactual Loss Condition"
table2.regret$Estimate <- myround(table2.regret$Estimate, digits = 2)
table2.regret$SE <- myround(table2.regret$SE, digits = 2)
table2.regret$z <- myround(table2.regret$z, digits = 2)
table2.regret$p <- myround(table2.regret$p, digits = 3)
table2.regret$CI.lower <- myround(table2.regret$CI.lower, digits = 4)
table2.regret$CI.upper <- myround(table2.regret$CI.upper, digits = 4)
table2.regret$p[table2.regret$p < .001] <- "< .001"
table2.regret$Predictor[table2.regret$Predictor == "premc.s"] <- "Prevention Pride"
table2.regret$Predictor[table2.regret$Predictor == "condE"] <- "Counterfactual Loss"
table2.regret$Predictor[table2.regret$Predictor == "promc.s"] <- "Promotion Pride"
table2.regret$Predictor[table2.regret$Predictor == "premc.s:condE"] <- "Prev. Pride x CF Loss"
table2.regret$Predictor[table2.regret$Predictor == "condE:promc.s"] <- "Prom. Pride x CF Loss"
table2.regret$Predictor[table2.regret$Predictor == "regret.s"] <- "Regret"
table2.regret[is.na(table2.regret)] <- ""
table2.regret$lhs <- NULL
table2.regret$op <- NULL
table2.regret$section <- NULL
table2.regret$order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 10, 11, 12, 14, 15)
table2.regret <- arrange(table2.regret, order)
table2.regret$order <- NULL

medmodel.regret <- table2.regret %>%
  mutate(
    Estimate = cell_spec(Estimate, "latex", align = "r"),
    SE = cell_spec(SE, "latex", align = "r"),
    z = cell_spec(z, "latex", align = "r"),
    p = cell_spec(p, "latex", align = "r"),
    CI.lower = cell_spec(CI.lower, "latex", align = "r"),
    CI.upper = cell_spec(CI.upper, "latex", align = "r")) %>%
  kable("latex", booktabs = TRUE, escape = FALSE, caption = "Summary of Exploratory Mediation Analysis") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  group_rows("Model 1 (DV = Regret)", 1, 6) %>%
  group_rows("Model 2 (DV = Bitcoin Allocation)", 7, 13) %>%
  group_rows("Bootstrapped Conditional Indirect Effects\n(Prev. Pride x CF Loss → Regret → Bitcoin Allocation)", 14, 15) %>%
  row_spec(0, align = "c") %>%
  footnote(general = "This analysis included an effect-coded variable for the counterfactual loss experience: -1 = control, 1 = counterfactual loss. Additionally, all other variables were standardized (M = 0, SD = 1) within this model. Estimated effect sizes for the models reported here are standardized regression coefficients.", general_title = "Note: ", footnote_as_chunk = T, title_format = c("italic"), threeparttable = TRUE, escape = FALSE)
medmodel.regret

## Exploratory Moderated Mediation Analysis: Hypothetical Happiness

# Specify model
Mod.Med.Lavaan.promctr.happier <- '
# Regressions
happier.s ~ 1 + a1*premc.s + a2*condE + a3*premc.s:condE + a4*promc.s + a5*condE:promc.s
allo.s ~ 1 + cdash1*premc.s + cdash2*condE + b1*happier.s + cdash3*premc.s:condE + cdash4*promc.s + cdash5*condE:promc.s

#Indirect effects conditional on moderator (a1 + a3*ModValue)*b1
indirect.ctrl := (a1 + a3*(-1))*b1
indirect.CFL := (a1 + a3*(1))*b1

#Direct effects conditional on moderator (cdash1 + cdash3*ModValue)
direct.ctrl := cdash1 + cdash3*(-1)
direct.CFL := cdash1 + cdash3*(1)

#Total effects conditional on moderator
total.ctrl := direct.ctrl + indirect.ctrl
total.CFL := direct.CFL + indirect.CFL

#Proportion mediated conditional on moderator
#To match the output of "mediate" package
prop.mediated.ctrl := indirect.ctrl / total.ctrl
prop.mediated.CFL := indirect.CFL / total.CFL

#Index of moderated mediation
#An alternative way of testing if conditional indirect effects are significantly different from each other
index.mod.med := a3*b1
'

# Fit model  - Commented out because this script takes a long time to rerun.
# set.seed(1234)
# Mod.Med.SEM.promctr.happier <- sem(model = Mod.Med.Lavaan.promctr.happier,
#                                    data = cfl,
#                                    se = "bootstrap",
#                                    bootstrap = 5000)
# save(Mod.Med.SEM.promctr.happier, file = "models/happierPromctrFit.RData")

load("models/happierPromctrFit.RData")
summary(Mod.Med.SEM.promctr.happier,
        fit.measures = FALSE,
        standardize = TRUE,
        estimates = TRUE,
        ci = TRUE,
        rsquare = TRUE)

happierPromctrCoefs <- parameterEstimates(Mod.Med.SEM.promctr.happier)
a3.happier.promctr <- happierPromctrCoefs$est[happierPromctrCoefs$label == "a3"]
a3.p.happier.promctr <- happierPromctrCoefs$pvalue[happierPromctrCoefs$label == "a3"]
b.happier.promctr <- happierPromctrCoefs$est[happierPromctrCoefs$label == "b1"]
b.p.happier.promctr <- happierPromctrCoefs$pvalue[happierPromctrCoefs$label == "b1"]
condIndEffctrl.happier.promctr <- happierPromctrCoefs$est[happierPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.se.happier.promctr <- happierPromctrCoefs$se[happierPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.cilower.happier.promctr <- happierPromctrCoefs$ci.lower[happierPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.ciupper.happier.promctr <- happierPromctrCoefs$ci.upper[happierPromctrCoefs$label == "indirect.ctrl"]
condIndEffCFL.happier.promctr <- happierPromctrCoefs$est[happierPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.se.happier.promctr <- happierPromctrCoefs$se[happierPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.cilower.happier.promctr <- happierPromctrCoefs$ci.lower[happierPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.ciupper.happier.promctr <- happierPromctrCoefs$ci.upper[happierPromctrCoefs$label == "indirect.CFL"]

# Generating table of moderated mediation results for manuscript
table2.happier <- happierPromctrCoefs 
table2.happier <- dplyr::rename(table2.happier, Predictor = rhs, Estimate = est, SE = se, p = pvalue, CI.lower = ci.lower, CI.upper = ci.upper)
table2.happier$label <- NULL
table2.happier <- table2.happier[table2.happier$op != "~~",]
table2.happier$Predictor[table2.happier$lhs == "happier.s" & table2.happier$op == "~1"] <- "Intercept"
table2.happier$op[table2.happier$lhs == "happier.s" & table2.happier$op == "~1"] <- "~"
table2.happier$Predictor[table2.happier$lhs == "allo.s" & table2.happier$op == "~1"] <- "Intercept"
table2.happier$op[table2.happier$lhs == "allo.s" & table2.happier$op == "~1"] <- "~"
table2.happier <- table2.happier[table2.happier$op != "~1",]
table2.happier <- table2.happier[1:15,]
table2.happier <- table2.happier %>%
  mutate(section = recode(lhs, happier.s = 1, allo.s = 2,
                          indirect.ctrl = 3, indirect.CFL = 4)) %>%
  arrange(section,desc(op))
table2.happier$Predictor[table2.happier$lhs == "indirect.ctrl"] <- "Control Condition"
table2.happier$Predictor[table2.happier$lhs == "indirect.CFL"] <- "Counterfactual Loss Condition"
table2.happier$Estimate <- myround(table2.happier$Estimate, digits = 2)
table2.happier$SE <- myround(table2.happier$SE, digits = 2)
table2.happier$z <- myround(table2.happier$z, digits = 2)
table2.happier$p <- myround(table2.happier$p, digits = 3)
table2.happier$CI.lower <- myround(table2.happier$CI.lower, digits = 4)
table2.happier$CI.upper <- myround(table2.happier$CI.upper, digits = 4)
table2.happier$p[table2.happier$p < .001] <- "< .001"
table2.happier$Predictor[table2.happier$Predictor == "premc.s"] <- "Prevention Pride"
table2.happier$Predictor[table2.happier$Predictor == "condE"] <- "Counterfactual Loss"
table2.happier$Predictor[table2.happier$Predictor == "promc.s"] <- "Promotion Pride"
table2.happier$Predictor[table2.happier$Predictor == "premc.s:condE"] <- "Prev. Pride x CF Loss"
table2.happier$Predictor[table2.happier$Predictor == "condE:promc.s"] <- "Prom. Pride x CF Loss"
table2.happier$Predictor[table2.happier$Predictor == "happier.s"] <- "Hypothetical Happiness"
table2.happier[is.na(table2.happier)] <- ""
table2.happier$lhs <- NULL
table2.happier$op <- NULL
table2.happier$section <- NULL
table2.happier$order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 10, 11, 12, 14, 15)
table2.happier <- arrange(table2.happier, order)
table2.happier$order <- NULL

medmodel.happier <- table2.happier %>%
  mutate(Estimate = cell_spec(Estimate, "latex", align = "r"),
         SE = cell_spec(SE, "latex", align = "r"),
         z = cell_spec(z, "latex", align = "r"),
         p = cell_spec(p, "latex", align = "r"),
         CI.lower = cell_spec(CI.lower, "latex", align = "r"),
         CI.upper = cell_spec(CI.upper, "latex", align = "r")) %>%
  kable("latex", booktabs = TRUE, escape = FALSE, caption = "Summary of Exploratory Mediation Analysis") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  group_rows("Model 1 (DV = Hypothetical Happiness)", 1, 6) %>%
  group_rows("Model 2 (DV = Bitcoin Allocation)", 7, 13) %>%
  group_rows("Bootstrapped Conditional Indirect Effects\n(Prev. Pride x CF Loss → Hypothetical Happiness → Bitcoin Allocation)", 14, 15) %>%
  row_spec(0, align = "c") %>%
  footnote(general = "This analysis included an effect-coded variable for the counterfactual loss experience: -1 = control, 1 = counterfactual loss. Additionally, all other variables were standardized (M = 0, SD = 1) within this model. Estimated effect sizes for the models reported here are standardized regression coefficients.", general_title = "Note: ", footnote_as_chunk = T, title_format = c("italic"), threeparttable = TRUE, escape = FALSE)
medmodel.happier

## Exploratory Moderated Mediation Analysis: Hypothetical Relief

# Specify model
Mod.Med.Lavaan.promctr.Relieved <- '
# Regressions
relieved.s ~ 1 + a1*premc.s + a2*condE + a3*premc.s:condE + a4*promc.s + a5*condE:promc.s
allo.s ~ 1 + cdash1*premc.s + cdash2*condE + b1*relieved.s + cdash3*premc.s:condE + cdash4*promc.s + cdash5*condE:promc.s

#Indirect effects conditional on moderator (a1 + a3*ModValue)*b1
indirect.ctrl := (a1 + a3*(-1))*b1
indirect.CFL := (a1 + a3*(1))*b1

#Direct effects conditional on moderator (cdash1 + cdash3*ModValue)
direct.ctrl := cdash1 + cdash3*(-1)
direct.CFL := cdash1 + cdash3*(1)

#Total effects conditional on moderator
total.ctrl := direct.ctrl + indirect.ctrl
total.CFL := direct.CFL + indirect.CFL

#Proportion mediated conditional on moderator
#To match the output of "mediate" package
prop.mediated.ctrl := indirect.ctrl / total.ctrl
prop.mediated.CFL := indirect.CFL / total.CFL

#Index of moderated mediation
#An alternative way of testing if conditional indirect effects are significantly different from each other
index.mod.med := a3*b1
'

# Fit model - Commented out because this script takes a long time to rerun.
# set.seed(1234)
# Mod.Med.SEM.promctr.relieved <- sem(model = Mod.Med.Lavaan.promctr.Relieved,
#                                     data = cfl,
#                                     se = "bootstrap",
#                                     bootstrap = 5000)
# save(Mod.Med.SEM.promctr.relieved, file = "models/relievedPromctrFit.RData")

load("models/relievedPromctrFit.RData")
summary(Mod.Med.SEM.promctr.relieved,
        fit.measures = FALSE,
        standardize = TRUE,
        estimates = TRUE,
        ci = TRUE,
        rsquare = TRUE)

relievedPromctrCoefs <- parameterEstimates(Mod.Med.SEM.promctr.relieved)
a3.relieved.promctr <- relievedPromctrCoefs$est[relievedPromctrCoefs$label == "a3"]
a3.p.relieved.promctr <- relievedPromctrCoefs$pvalue[relievedPromctrCoefs$label == "a3"]
b.relieved.promctr <- relievedPromctrCoefs$est[relievedPromctrCoefs$label == "b1"]
b.p.relieved.promctr <- relievedPromctrCoefs$pvalue[relievedPromctrCoefs$label == "b1"]
condIndEffctrl.relieved.promctr <- relievedPromctrCoefs$est[relievedPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.se.relieved.promctr <- relievedPromctrCoefs$se[relievedPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.cilower.relieved.promctr <- relievedPromctrCoefs$ci.lower[relievedPromctrCoefs$label == "indirect.ctrl"]
condIndEffctrl.ciupper.relieved.promctr <- relievedPromctrCoefs$ci.upper[relievedPromctrCoefs$label == "indirect.ctrl"]
condIndEffCFL.relieved.promctr <- relievedPromctrCoefs$est[relievedPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.se.relieved.promctr <- relievedPromctrCoefs$se[relievedPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.cilower.relieved.promctr <- relievedPromctrCoefs$ci.lower[relievedPromctrCoefs$label == "indirect.CFL"]
condIndEffCFL.ciupper.relieved.promctr <- relievedPromctrCoefs$ci.upper[relievedPromctrCoefs$label == "indirect.CFL"]

# Generating table of moderated mediation results
table2.relieved <- relievedPromctrCoefs 
table2.relieved <- dplyr::rename(table2.relieved, Predictor = rhs, Estimate = est, SE = se, p = pvalue, CI.lower = ci.lower, CI.upper = ci.upper)
table2.relieved$label <- NULL
table2.relieved <- table2.relieved[table2.relieved$op != "~~",]
table2.relieved$Predictor[table2.relieved$lhs == "relieved.s" & table2.relieved$op == "~1"] <- "Intercept"
table2.relieved$op[table2.relieved$lhs == "relieved.s" & table2.relieved$op == "~1"] <- "~"
table2.relieved$Predictor[table2.relieved$lhs == "allo.s" & table2.relieved$op == "~1"] <- "Intercept"
table2.relieved$op[table2.relieved$lhs == "allo.s" & table2.relieved$op == "~1"] <- "~"
table2.relieved <- table2.relieved[table2.relieved$op != "~1",]
table2.relieved <- table2.relieved[1:15,]
table2.relieved <- table2.relieved %>%
  mutate(section = recode(lhs, relieved.s = 1, allo.s = 2,
                          indirect.ctrl = 3, indirect.CFL = 4)) %>%
  arrange(section,desc(op))
table2.relieved$Predictor[table2.relieved$lhs == "indirect.ctrl"] <- "Control Condition"
table2.relieved$Predictor[table2.relieved$lhs == "indirect.CFL"] <- "Counterfactual Loss Condition"
table2.relieved$Estimate <- myround(table2.relieved$Estimate, digits = 2)
table2.relieved$SE <- myround(table2.relieved$SE, digits = 2)
table2.relieved$z <- myround(table2.relieved$z, digits = 2)
table2.relieved$p <- myround(table2.relieved$p, digits = 3)
table2.relieved$CI.lower <- myround(table2.relieved$CI.lower, digits = 4)
table2.relieved$CI.upper <- myround(table2.relieved$CI.upper, digits = 4)
table2.relieved$p[table2.relieved$p < .001] <- "< .001"
table2.relieved$Predictor[table2.relieved$Predictor == "premc.s"] <- "Prevention Pride"
table2.relieved$Predictor[table2.relieved$Predictor == "condE"] <- "Counterfactual Loss"
table2.relieved$Predictor[table2.relieved$Predictor == "promc.s"] <- "Promotion Pride"
table2.relieved$Predictor[table2.relieved$Predictor == "premc.s:condE"] <- "Prev. Pride x CF Loss"
table2.relieved$Predictor[table2.relieved$Predictor == "condE:promc.s"] <- "Prom. Pride x CF Loss"
table2.relieved$Predictor[table2.relieved$Predictor == "relieved.s"] <- "Hypothetical Relief"
table2.relieved[is.na(table2.relieved)] <- ""
table2.relieved$lhs <- NULL
table2.relieved$op <- NULL
table2.relieved$section <- NULL
table2.relieved$order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 10, 11, 12, 14, 15)
table2.relieved <- arrange(table2.relieved, order)
table2.relieved$order <- NULL

medmodel.relieved <- table2.relieved %>%
  mutate(Estimate = cell_spec(Estimate, "latex", align = "r"),
         SE = cell_spec(SE, "latex", align = "r"),
         z = cell_spec(z, "latex", align = "r"),
         p = cell_spec(p, "latex", align = "r"),
         CI.lower = cell_spec(CI.lower, "latex", align = "r"),
         CI.upper = cell_spec(CI.upper, "latex", align = "r")
         ) %>%
  kable("latex", booktabs = TRUE, escape = FALSE, caption = "Summary of Exploratory Mediation Analysis") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  group_rows("Model 1 (DV = Hypothetical Relief)", 1, 6) %>%
  group_rows("Model 2 (DV = Bitcoin Allocation)", 7, 13) %>%
  group_rows("Bootstrapped Conditional Indirect Effects\n(Prev. Pride x CF Loss → Hypothetical Relief → Bitcoin Allocation)", 14, 15) %>%
  row_spec(0, align = "c") %>%
  footnote(general = "This analysis included an effect-coded variable for the counterfactual loss experience: -1 = control, 1 = counterfactual loss. Additionally, all other variables were standardized (M = 0, SD = 1) within this model. Estimated effect sizes for the models reported here are standardized regression coefficients.", general_title = "Note: ", footnote_as_chunk = T, title_format = c("italic"), threeparttable = TRUE, escape = FALSE)
medmodel.relieved

# Rerunning with Exclusion of Participants Who Learned About Study from Other M-Turkers

# Creating subset to facilitate visualization of results within each condition
cfl.rerun.cflcond <- subset(cfl.rerun, condition == "CF Loss") # n = 144
cfl.rerun.ctrlcond <- subset(cfl.rerun, condition == "Control") # n = 157

# SAVING CLEAN DATASETS
save(cfl.rerun, cfl.rerun.cflcond, cfl.rerun.ctrlcond, file = "data/rfcfloss_rerun_clean.RData")
write.csv(cfl.rerun, file = "data/rfcfloss_rerun_clean.csv")

# Descriptive plot
# summarySE() function defined at top of document
cfl.rerun.summary <- summarySE(cfl.rerun, measurevar="allo", groupvars=c("hiloprev","condition"))
rerunplot <- ggplot(data=cfl.rerun, aes(x = hiloprev, y = allo, color = hiloprev)) +
  facet_wrap(. ~ condition) +
  geom_jitter() +
  geom_errorbar(data = cflsummary, aes(ymin = allo - se, ymax = allo + se),
                col = "black", size = .5) +
  geom_point(data = cflsummary, col = "black", size = 3, alpha = 1) +
  scale_x_discrete(labels = c("High", "Low")) +
  labs(title="RERUN: Bitcoin Allocation By Counterfactual Loss and\nPrevention Pride (NOTE: Does Not Control\nfor Effects of Promotion Pride)", x="Prevention Pride (Median-Split)", y="Bitcoin Allocation (%)") +
  theme(legend.position = "none")
rerunplot

## Primary linear analysis and table generation
rfcfl.promctr.rerun <- lm(allo ~ premc * condD + promc * condD, data=cfl.rerun)
summary(rfcfl.promctr.rerun) # premc:condD b = -11.973, p = .0292
table1r <- tidy(rfcfl.promctr.rerun)
table1r <- dplyr::rename(table1r, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1r$Estimate <- myround(table1r$Estimate, digits = 2)
table1r$SE <- myround(table1r$SE, digits = 2)
table1r$t <- myround(table1r$t, digits = 2)
table1r$p <- myround(table1r$p, digits = 3)
table1r$p[table1r$p < .001] <- "< .001"
table1r$Predictor[table1r$Predictor == "(Intercept)"] <- "Intercept"
table1r$Predictor[table1r$Predictor == "premc"] <- "Prevention Pride"
table1r$Predictor[table1r$Predictor == "condD"] <- "Counterfactual Loss"
table1r$Predictor[table1r$Predictor == "promc"] <- "Promotion Pride"
table1r$Predictor[table1r$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1r$Predictor[table1r$Predictor == "condD:promc"] <- "Promotion Pride x Counterfactual Loss"
table1r[is.na(table1r)] <- ""
table1r$order <- c(1, 2, 3, 5, 4, 6)
table1r <- arrange(table1r, order)
table1r$order <- NULL
apa_table(table1r, caption = "RERUN: Summary of Linear Regression Analysis", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Additionally, prevention and promotion pride were mean-centered within this model. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')

# Generating graph with results of rerun primary linear analysis
pred.data.ctrl.rerun <- data.frame(premc = seq(min(cfl.rerun$premc, na.rm=T), 
                                               max(cfl.rerun$premc, na.rm=T), .1),
                                   promc = 0, condD = 0, condition = "Control")
rfcfl.promctr.pred.ctrl.rerun <- cbind(pred.data.ctrl.rerun, 
                                       predict(rfcfl.promctr.rerun, pred.data.ctrl.rerun, 
                                               interval = "confidence", 
                                               type = c("response", "terms")))
rfcfl.promctr.pred.ctrl.rerun <- dplyr::rename(rfcfl.promctr.pred.ctrl.rerun, allo = fit)
pred.data.CFL.rerun <- data.frame(premc = seq(min(cfl.rerun$premc, na.rm=T), 
                                              max(cfl.rerun$premc, na.rm=T), .1),
                                  promc = 0, condD = 1, condition = "CF Loss")
rfcfl.promctr.pred.CFL.rerun <- cbind(pred.data.CFL.rerun, 
                                      predict(rfcfl.promctr.rerun, pred.data.CFL.rerun, 
                                              interval = "confidence", 
                                              type = c("response", "terms")))
rfcfl.promctr.pred.CFL.rerun <- dplyr::rename(rfcfl.promctr.pred.CFL.rerun, allo = fit)

ggplot(data=cfl.rerun, aes(x=premc, y=allo)) +
  scale_color_manual(values=c("gray0", "gray50")) +
  geom_ribbon(data = rfcfl.promctr.pred.ctrl, aes(ymin = lwr, ymax = upr), alpha = .3, fill = "gray") +
  geom_ribbon(data = rfcfl.promctr.pred.CFL, aes(ymin = lwr, ymax = upr), alpha = .3, fill = "gray") +
  geom_line(data = rfcfl.promctr.pred.ctrl, aes(color=condition), size = 1) +
  geom_line(data = rfcfl.promctr.pred.CFL, aes(color=condition), size = 1) +
  labs(title="RERUN: Prevention Pride and Counterfactual Loss as\nPredictors of Bitcoin Allocation\n(Controlling for Promotion Pride and the Interaction\nBetween Promotion Pride and Counterfactual Loss)", x="Prevention Pride (Mean-Centered)", y="Bitcoin Allocation (%)", color = "Condition")

# Simple Slopes and Johnson-Neyman Analyses

jn1.rerun <- sim_slopes(model = rfcfl.promctr.rerun, pred = premc, modx = condD, modx.values = c(0,1),
                        centered = "all", data = cfl, cond.int = FALSE,
                        johnson_neyman = FALSE, jnplot = FALSE, jnalpha = 0.05, robust = FALSE,
                        digits = getOption("jtools-digits", default = 2), pvals = TRUE,
                        confint = TRUE, ci.width = 0.95)
jn1.rerun <- tidy(jn1.rerun)

rfcfl.promctr.nonstd.rerun <- lm(allo ~ prev * condD + promc * condD, data=cfl.rerun)
jn2.rerun <- sim_slopes(model = rfcfl.promctr.nonstd.rerun, pred = condD, modx = prev, modx.values = c(2.50991, 4.572230), 
                        centered = "all", data = cfl, cond.int = FALSE,
                        johnson_neyman = TRUE, jnplot = TRUE, jnalpha = 0.05, robust = FALSE,
                        digits = getOption("jtools-digits", default = 2), pvals = TRUE,
                        confint = TRUE, ci.width = 0.95)
jncutoff_low.rerun <- jn2.rerun$jn[[1]]$bounds[1]
jncutoff_high.rerun <- jn2.rerun$jn[[1]]$bounds[2]
jn2.rerun <- tidy(jn2.rerun)

## Exploratory Mediation Analysis (Hayes PROCESS Model 8): Counterfactual Relief

# Fit model (no changes in model from original specification)
# Commented out because this script takes a long time to rerun.
# set.seed(1234)
# Mod.Med.SEM.promctr.relieved.rerun <- sem(model = Mod.Med.Lavaan.promctr.Relieved,
#                                           data = cfl.rerun,
#                                           se = "bootstrap",
#                                           bootstrap = 5000)
# save(Mod.Med.SEM.promctr.relieved.rerun, file = "models/relievedPromctrFitRerun.RData")

load("models/relievedPromctrFitRerun.RData")
summary(Mod.Med.SEM.promctr.relieved.rerun,
        fit.measures = FALSE,
        standardize = TRUE,
        estimates = TRUE,
        ci = TRUE,
        rsquare = TRUE)

relievedPromctrCoefs.rerun <- parameterEstimates(Mod.Med.SEM.promctr.relieved.rerun)
a3.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$est[relievedPromctrCoefs.rerun$label == "a3"]
a3.p.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$pvalue[relievedPromctrCoefs.rerun$label == "a3"]
b.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$est[relievedPromctrCoefs.rerun$label == "b1"]
b.p.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$pvalue[relievedPromctrCoefs.rerun$label == "b1"]
condIndEffctrl.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$est[relievedPromctrCoefs.rerun$label == "indirect.ctrl"]
condIndEffctrl.se.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$se[relievedPromctrCoefs.rerun$label == "indirect.ctrl"]
condIndEffctrl.cilower.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$ci.lower[relievedPromctrCoefs.rerun$label == "indirect.ctrl"]
condIndEffctrl.ciupper.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$ci.upper[relievedPromctrCoefs.rerun$label == "indirect.ctrl"]
condIndEffCFL.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$est[relievedPromctrCoefs.rerun$label == "indirect.CFL"]
condIndEffCFL.se.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$se[relievedPromctrCoefs.rerun$label == "indirect.CFL"]
condIndEffCFL.cilower.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$ci.lower[relievedPromctrCoefs.rerun$label == "indirect.CFL"]
condIndEffCFL.ciupper.relieved.promctr.rerun <- relievedPromctrCoefs.rerun$ci.upper[relievedPromctrCoefs.rerun$label == "indirect.CFL"]

table2.relieved.rerun <- relievedPromctrCoefs.rerun
table2.relieved.rerun <- dplyr::rename(table2.relieved.rerun, Predictor = rhs, Estimate = est, SE = se, p = pvalue, CI.lower = ci.lower, CI.upper = ci.upper)
table2.relieved.rerun$label <- NULL
table2.relieved.rerun <- table2.relieved.rerun[table2.relieved.rerun$op != "~~",]
table2.relieved.rerun$Predictor[table2.relieved.rerun$lhs == "relieved.s" & table2.relieved.rerun$op == "~1"] <- "Intercept"
table2.relieved.rerun$op[table2.relieved.rerun$lhs == "relieved.s" & table2.relieved.rerun$op == "~1"] <- "~"
table2.relieved.rerun$Predictor[table2.relieved.rerun$lhs == "allo.s" & table2.relieved.rerun$op == "~1"] <- "Intercept"
table2.relieved.rerun$op[table2.relieved.rerun$lhs == "allo.s" & table2.relieved.rerun$op == "~1"] <- "~"
table2.relieved.rerun <- table2.relieved.rerun[table2.relieved.rerun$op != "~1",]
table2.relieved.rerun <- table2.relieved.rerun[1:15,]
table2.relieved.rerun <- table2.relieved.rerun %>%
  mutate(section = recode(lhs, relieved.s = 1, allo.s = 2,
                          indirect.ctrl = 3, indirect.CFL = 4)) %>%
  arrange(section,desc(op))
table2.relieved.rerun$Predictor[table2.relieved.rerun$lhs == "indirect.ctrl"] <- "Control Condition"
table2.relieved.rerun$Predictor[table2.relieved.rerun$lhs == "indirect.CFL"] <- "Counterfactual Loss Condition"
table2.relieved.rerun$Estimate <- myround(table2.relieved.rerun$Estimate, digits = 2)
table2.relieved.rerun$SE <- myround(table2.relieved.rerun$SE, digits = 2)
table2.relieved.rerun$z <- myround(table2.relieved.rerun$z, digits = 2)
table2.relieved.rerun$p <- myround(table2.relieved.rerun$p, digits = 3)
table2.relieved.rerun$CI.lower <- myround(table2.relieved.rerun$CI.lower, digits = 4)
table2.relieved.rerun$CI.upper <- myround(table2.relieved.rerun$CI.upper, digits = 4)
table2.relieved.rerun$p[table2.relieved.rerun$p < .001] <- "< .001"
table2.relieved.rerun$Predictor[table2.relieved.rerun$Predictor == "premc.s"] <- "Prevention Pride"
table2.relieved.rerun$Predictor[table2.relieved.rerun$Predictor == "condE"] <- "Counterfactual Loss"
table2.relieved.rerun$Predictor[table2.relieved.rerun$Predictor == "promc.s"] <- "Promotion Pride"
table2.relieved.rerun$Predictor[table2.relieved.rerun$Predictor == "premc.s:condE"] <- "Prev. Pride x CF Loss"
table2.relieved.rerun$Predictor[table2.relieved.rerun$Predictor == "condE:promc.s"] <- "Prom. Pride x CF Loss"
table2.relieved.rerun$Predictor[table2.relieved.rerun$Predictor == "relieved.s"] <- "Hypothetical Relief"
table2.relieved.rerun[is.na(table2.relieved.rerun)] <- ""
table2.relieved.rerun$lhs <- NULL
table2.relieved.rerun$op <- NULL
table2.relieved.rerun$section <- NULL
table2.relieved.rerun$order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 10, 11, 12, 14, 15)
table2.relieved.rerun <- arrange(table2.relieved.rerun, order)
table2.relieved.rerun$order <- NULL

medmodel.relieved.rerun <- table2.relieved.rerun %>%
  mutate(Estimate = cell_spec(Estimate, "latex", align = "r"),
         SE = cell_spec(SE, "latex", align = "r"),
         z = cell_spec(z, "latex", align = "r"),
         p = cell_spec(p, "latex", align = "r"),
         CI.lower = cell_spec(CI.lower, "latex", align = "r"),
         CI.upper = cell_spec(CI.upper, "latex", align = "r")) %>%
  kable("latex", booktabs = TRUE, escape = FALSE, caption = "Summary of Exploratory Mediation Analysis") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  group_rows("Model 1 (DV = Hypothetical Relief)", 1, 6) %>%
  group_rows("Model 2 (DV = Bitcoin Allocation)", 7, 13) %>%
  group_rows("Bootstrapped Conditional Indirect Effects\n(Prev. Pride x CF Loss → Hypothetical Relief → Bitcoin Allocation)", 14, 15) %>%
  row_spec(0, align = "c") %>%
  footnote(general = "This analysis included an effect-coded variable for the counterfactual loss experience: -1 = control, 1 = counterfactual loss. Additionally, all other variables were standardized (M = 0, SD = 1) within this model. Estimated effect sizes for the models reported here are standardized regression coefficients.", general_title = "Note: ", footnote_as_chunk = T, title_format = c("italic"), threeparttable = TRUE, escape = FALSE)
medmodel.relieved.rerun