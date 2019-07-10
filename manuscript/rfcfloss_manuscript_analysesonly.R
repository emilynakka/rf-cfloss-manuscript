### R Script for Analyses in Manuscript Titled:
### Does Your Gain Define My Loss?:
### Socially-Defined Counterfactual Loss and Prevention-Focused Decision-Making
###
### Although this R script contains all analyses, the full manuscript (including the analyses
### contained here) can be generated using the R Markdown file located on this Open Science
### Framework page: https://osf.io/e4nr3/?view_only=73c8490e19c6449197d2620d29ea8ad2

# Load packages
library(knitr) # Required for knitting
library(papaja) # Required for APA template
library(citr) # Required for easy insertion of citations
library(tidyverse) # Required for data cleaning
library(broman) # Required for myround() function that doesn't truncate digits
library(rIP) # Required to filter out fraudulent M-Turk participants
library(psych) # Required to calculate Cronbach's alpha
library(cocron) # Required to statistically compare Cronbach's alphas
library(broom) # Required to extract effects from models
library(kableExtra) # Required for table styling
library(ggplot2) # Required for plots
library(shiny) # Required to properly cite shiny when referencing our data visualization app
library(jtools) # Required for Johnson-Neyman and simple slopes analyses
library(cowplot) # Required for sim_slopes
library(boot) # Required for bootstrapping CIs in mediation analysis
library(lavaan) # Required for mediation analyses

# Seed for random number generation
set.seed(1234)

# Set knitr options
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, include = FALSE, message = FALSE)
knitr::opts_knit$set(root.dir = normalizePath('../'))

# Define functions
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

# Load cleaned data (full cleaning script can be found in Supplementary Materials)
load(file = "data/rfcfloss_clean.RData")
load(file = "data/rfcfloss_coded.RData")

## Excluding participants who did not submit a Bitcoin allocation percentage (our DV)
cflorig.alloexcl <- cflorig[ which(is.na(cflorig$allo) == FALSE ) , ] # n = 482 remaining

## Excluding participants who did not complete prevention items (our moderator of interest)
cflorig.prevexcl <- cflorig.alloexcl[ which(is.na(cflorig.alloexcl$prev) == FALSE ) , ] # n = 481 remaining
cflorig.ivdvexcl.n <- 500 - nrow(cflorig.prevexcl)

## Identifying participants with fraudulent IP addresses
cflexcl <- cflorig.prevexcl
# blockornot <- getIPinfo(cflexcl1, "IPAddress", "[key]") # Commented out to avoid rerunning
# # each time I run the full chunk, since script will not work without API key
# # Anyone looking to rerun the line of code above will need to secure a free X-key
# # from iphub.info and place it within the quotes

# blockornot <- dplyr::rename(blockornot, IPAddress = ip)
cflexcl <- merge(cflexcl, blockornot, by = "IPAddress")

## Determining participants flagged as block = (1 or 2) vs. 1 only, plus US only
cflexcl$block12 <- ifelse(cflexcl$block == "0", 0, 1)
cflexcl$block <- as.numeric(cflexcl$block)
cflexcl$block.d <- ifelse(cflexcl$block == 1, 1, 0)
cflexcl.blocked <- cflexcl[cflexcl$block.d == 0,] # n = 334
cflexcl.blocked.usonly <- cflexcl.blocked[cflexcl.blocked$countryCode == "US",] # n = 301
cfl.ipexcl <- cflexcl.blocked.usonly
cflexcl.allblocked.n <- nrow(cflorig.prevexcl) - nrow(cfl.ipexcl) # n = 180

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

# Calculating participant information
n <- length(unique(cfl$ResponseId))
meanage <- mean(cfl$age, na.rm = TRUE)
sdage <- sd(cfl$age, na.rm = TRUE)
maxage <- max(cfl$age, na.rm = TRUE)
minage <- min(cfl$age, na.rm = TRUE)
male <- sum(cfl$gender == "male")
female <- sum(cfl$gender == "female")
malepct <- myround(male/nrow(cfl)*100, digits = 2)
femalepct <- myround(female/nrow(cfl)*100, digits = 2)

# Calculating key descriptive statistics
prommean <- myround(mean(cfl$prom), digits = 2)
promsd <- myround(sd(cfl$prom), digits = 2)
prevmean <- myround(mean(cfl$prev), digits = 2)
prevsd <- myround(sd(cfl$prev), digits = 2)
investearliermean <- myround(mean(cfl$investearlier, na.rm = TRUE), digits = 2)
investearliersd <- myround(sd(cfl$investearlier, na.rm = TRUE), digits = 2)
missedoutmean <- myround(mean(cfl$missedout, na.rm = TRUE), digits = 2)
missedoutsd <- myround(sd(cfl$missedout, na.rm = TRUE), digits = 2)
regretmean <- myround(mean(cfl$regret, na.rm = TRUE), digits = 2)
regretsd <- myround(sd(cfl$regret, na.rm = TRUE), digits = 2)
happiermean <- myround(mean(cfl$happier, na.rm = TRUE), digits = 2)
happiersd <- myround(sd(cfl$happier, na.rm = TRUE), digits = 2)
relievedmean <- myround(mean(cfl$relieved, na.rm = TRUE), digits = 2)
relievedsd <- myround(sd(cfl$relieved, na.rm = TRUE), digits = 2)

# Primary linear analysis
rfcfl.promctr.nonstd <- lm(allo ~ premc * condD + promc * condD, data = cfl) # unstandardized
rfcfl.promctr <- lm(scale(allo) ~ scale(premc) * scale(condD) + scale(promc) * scale(condD), data = cfl) # standardized

# Generating table with results of primary linear analysis
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
apa_table(table1, caption = "Summary of Linear Regression Analysis", note = "This analysis included a dummy-coded variable for the counterfactual loss manipulation: 0 = control, 1 = counterfactual loss. Estimated effect sizes reported here are unstandardized regression coefficients.", align = 'lrrrr')
df <- rfcfl.promctr.nonstd$df.residual

# Generating graph with results of primary linear analysis
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
ggsave("img/rfcflplot.png", plot = last_plot(), device = "png", units = "in", width = 6, height = 4.5)
primarymodel

# Original model, rerun to include participants who previously invested in Bitcoin
# (And creation of associated table of model results)
rfcfl.promctr.rerun <- lm(allo ~ premc * condD + promc * condD, data=cfl.rerun)
table1.rerun <- tidy(rfcfl.promctr.rerun)
table1.rerun <- dplyr::rename(table1.rerun, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.rerun$Estimate <- myround(table1.rerun$Estimate, digits = 2)
table1.rerun$p <- myround(table1.rerun$p, digits = 3)
table1.rerun$p[table1.rerun$p < .001] <- "< .001"
table1.rerun$Predictor[table1.rerun$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.rerun[is.na(table1.rerun)] <- ""

# Simple slopes analysis
jn1 <- sim_slopes(model = rfcfl.promctr.nonstd, pred = premc, modx = condD, modx.values = c(0,1),
                  centered = "all", data = cfl, cond.int = FALSE,
                  johnson_neyman = FALSE, jnplot = FALSE, jnalpha = 0.05, robust = FALSE,
                  digits = getOption("jtools-digits", default = 2), pvals = TRUE,
                  confint = TRUE, ci.width = 0.95)
jn1 <- tidy(jn1)

# Johnson-Neyman analysis
rfcfl.promctr.nonctr <- lm(allo ~ prev * condD + promc * condD, data = cfl) # Un-centering prevention...
# ... in this model so the values reflect the actual prevention scale
jn2 <- sim_slopes(model = rfcfl.promctr.nonctr, pred = condD, modx = prev, modx.values = c(2.50991, 3.44, 4.572230),
                  centered = "all", data = cfl, cond.int = FALSE,
                  johnson_neyman = TRUE, jnplot = TRUE, jnalpha = 0.05, robust = FALSE,
                  digits = getOption("jtools-digits", default = 2), pvals = TRUE,
                  confint = TRUE, ci.width = 0.95)
jncutoff_low <- jn2$jn[[1]]$bounds[1]
jncutoff_high <- jn2$jn[[1]]$bounds[2]
jn2 <- tidy(jn2)

# Original model not controlling for promotion pride (mean-centered RF and dummy-coded condition)
rfcfl <- lm(allo ~ premc * condD, data=cfl)
table1.noprom <- tidy(rfcfl)
table1.noprom <- dplyr::rename(table1.noprom, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.noprom$Estimate <- myround(table1.noprom$Estimate, digits = 2)
table1.noprom$p <- myround(table1.noprom$p, digits = 3)
table1.noprom$p[table1.noprom$p < .001] <- "< .001"
table1.noprom$Predictor[table1.noprom$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.noprom[is.na(table1.noprom)] <- ""

# Original model controlling for gender (mean-centered RF and dummy-coded condition)
rfcfl.promctr.gender <- lm(allo ~ premc * condD + gender + promc * condD, data=cfl)
table1.gender <- tidy(rfcfl.promctr.gender)
table1.gender <- dplyr::rename(table1.gender, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.gender$Estimate <- myround(table1.gender$Estimate, digits = 2)
table1.gender$p <- myround(table1.gender$p, digits = 3)
table1.gender$p[table1.gender$p < .001] <- "< .001"
table1.gender$Predictor[table1.gender$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.gender[is.na(table1.gender)] <- ""

# Original model controlling for age (mean-centered RF and dummy-coded condition)
rfcfl.promctr.age <- lm(allo ~ premc * condD + age + promc * condD, data=cfl)
table1.age <- tidy(rfcfl.promctr.age)
table1.age <- dplyr::rename(table1.age, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.age$Estimate <- myround(table1.age$Estimate, digits = 2)
table1.age$p <- myround(table1.age$p, digits = 3)
table1.age$p[table1.age$p < .001] <- "< .001"
table1.age$Predictor[table1.age$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.age[is.na(table1.age)] <- ""

# Original model controlling for ethnicity (mean-centered RF and dummy-coded condition)
rfcfl.promctr.ethnicity <- lm(allo ~ premc * condD + ethnicity + promc * condD, data=cfl)
table1.ethnicity <- tidy(rfcfl.promctr.ethnicity)
table1.ethnicity <- dplyr::rename(table1.ethnicity, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.ethnicity$Estimate <- myround(table1.ethnicity$Estimate, digits = 2)
table1.ethnicity$p <- myround(table1.ethnicity$p, digits = 3)
table1.ethnicity$p[table1.ethnicity$p < .001] <- "< .001"
table1.ethnicity$Predictor[table1.ethnicity$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.ethnicity[is.na(table1.ethnicity)] <- ""

# Original model controlling for income (mean-centered RF and dummy-coded condition)
rfcfl.promctr.income <- lm(allo ~ premc * condD + income + promc * condD, data=cfl)
table1.income <- tidy(rfcfl.promctr.income)
table1.income <- dplyr::rename(table1.income, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.income$Estimate <- myround(table1.income$Estimate, digits = 2)
table1.income$p <- myround(table1.income$p, digits = 3)
table1.income$p[table1.income$p < .001] <- "< .001"
table1.income$Predictor[table1.income$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.income[is.na(table1.income)] <- ""

# Original model controlling for education (mean-centered RF and dummy-coded condition) - MS
rfcfl.promctr.edu <- lm(allo ~ premc * condD + education + promc * condD, data=cfl)
table1.edu <- tidy(rfcfl.promctr.edu)
table1.edu <- dplyr::rename(table1.edu, Predictor = term, Estimate = estimate, SE = std.error, t = statistic, p = p.value)
table1.edu$Estimate <- myround(table1.edu$Estimate, digits = 2)
table1.edu$p <- myround(table1.edu$p, digits = 3)
table1.edu$p[table1.edu$p < .001] <- "< .001"
table1.edu$Predictor[table1.edu$Predictor == "premc:condD"] <- "Prevention Pride x Counterfactual Loss"
table1.edu[is.na(table1.edu)] <- ""

# Exploratory Moderated Mediation Analysis: Specify model
Mod.Med.Lavaan.promctr.Relieved <- '
# Regressions
relieved.s ~ 1 + a1*premc.s + a2*condE + a3*premc.s:condE + a4*promc.s + a5*promc.s:condE
allo.s ~ 1 + cdash1*premc.s + cdash2*condE + b1*relieved.s + cdash3*premc.s:condE + cdash4*promc.s + cdash5*promc.s:condE

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
#                                     meanstructure = TRUE,
#                                     se = "bootstrap",
#                                     bootstrap = 5000)
# save(Mod.Med.SEM.promctr.relieved, file = "models/relievedFitFinal.RData")

load("models/relievedFitFinal.RData")

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

# Generating table of moderated mediation results for manuscript
table2 <- relievedPromctrCoefs 
table2 <- dplyr::rename(table2, Predictor = rhs, Estimate = est, SE = se, p = pvalue, CI.lower = ci.lower, CI.upper = ci.upper)
table2$label <- NULL
table2 <- table2[table2$op != "~~",]
table2$Predictor[table2$lhs == "relieved.s" & table2$op == "~1"] <- "Intercept"
table2$op[table2$lhs == "relieved.s" & table2$op == "~1"] <- "~"
table2$Predictor[table2$lhs == "allo.s" & table2$op == "~1"] <- "Intercept"
table2$op[table2$lhs == "allo.s" & table2$op == "~1"] <- "~"
table2 <- table2[table2$op != "~1",]
table2 <- table2[1:15,]
table2 <- table2 %>%
  mutate(section = recode(lhs, relieved.s = 1, allo.s = 2,
                          indirect.ctrl = 3, indirect.CFL = 4)) %>%
  arrange(section,desc(op))
table2$Predictor[table2$lhs == "indirect.ctrl"] <- "Control Condition"
table2$Predictor[table2$lhs == "indirect.CFL"] <- "Counterfactual Loss Condition"
table2$Estimate <- myround(table2$Estimate, digits = 2)
table2$SE <- myround(table2$SE, digits = 2)
table2$z <- myround(table2$z, digits = 2)
table2$p <- myround(table2$p, digits = 3)
table2$CI.lower <- myround(table2$CI.lower, digits = 4)
table2$CI.upper <- myround(table2$CI.upper, digits = 4)
table2$p[table2$p < .001] <- "< .001"
table2$Predictor[table2$Predictor == "premc.s"] <- "Prevention Pride"
table2$Predictor[table2$Predictor == "condE"] <- "Counterfactual Loss"
table2$Predictor[table2$Predictor == "promc.s"] <- "Promotion Pride"
table2$Predictor[table2$Predictor == "premc.s:condE"] <- "Prev. Pride x CF Loss"
table2$Predictor[table2$Predictor == "promc.s:condE"] <- "Prom. Pride x CF Loss"
table2$Predictor[table2$Predictor == "relieved.s"] <- "Hypothetical Relief"
table2[is.na(table2)] <- ""
table2$lhs <- NULL
table2$op <- NULL
table2$section <- NULL
table2$order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 10, 11, 12, 14, 15)
table2 <- arrange(table2, order)
table2$order <- NULL
medmodel <- table2 %>%
  mutate(Estimate = cell_spec(Estimate, "latex", align = "r"),
         SE = cell_spec(SE, "latex", align = "r"),
         z = cell_spec(z, "latex", align = "r"),
         p = cell_spec(p, "latex", align = "r"),
         CI.lower = cell_spec(CI.lower, "latex", align = "r"),
         CI.upper = cell_spec(CI.upper, "latex", align = "r")) %>%
  kable("latex", booktabs = TRUE, escape = FALSE, caption = "Summary of Mediation Analysis") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  group_rows("Model 1 (DV = Hypothetical Relief)", 1, 6) %>%
  group_rows("Model 2 (DV = Bitcoin Allocation)", 7, 13) %>%
  group_rows("Bootstrapped Conditional Indirect Effects\n(Prev. Pride x CF Loss → Hypothetical Relief → Bitcoin Allocation)", 14, 15) %>%
  row_spec(0, align = "c") %>%
  footnote(general = "This analysis included an effect-coded variable for the counterfactual loss condition: -1 = control, 1 = counterfactual loss. All other variables were standardized (M = 0, SD = 1). Estimated effect sizes for the models reported here are standardized regression coefficients.", general_title = "Note: ", footnote_as_chunk = T, title_format = c("italic"), threeparttable = TRUE, escape = FALSE)
medmodel