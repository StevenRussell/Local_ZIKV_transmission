#-------------------------------------------------------------------------------------------------------------------------------------#
# In this code, we are running a simulation to model the effectiveness of various surveillance strategies in detecting
# local Zika virus transmission. The end result is a data table which can be used for analysis and visualization.
#
# Created by Steven Russell
# Last updated: September 5, 2017
#-------------------------------------------------------------------------------------------------------------------------------------#

# Optional memory management: restart R session (and remove all R objects)

# .rs.restartR()
# rm(list = ls())

# Loading the required packages

require(dplyr)
require(tidyr)
require(data.table)

# Setting the seed and number of samples

set.seed(123)
n.samples <- 10000

# Creating a sequence of incidences and population totals that we are interested in

incidences = 10^c(seq(from = -5, to = -4.09, by = .08), -4, seq(from = -3.88, to = -3, by = .08))
              
pops = c(10000, 100000, 1000000) 

# Creating a variable that lists the different types of surveillance systems: 
#   Pregnant women, blood bank donors, 31 specific symptom combinations, 3 general symptom combinations

types <- factor(c("pregnant", "blood", "arth", "conj", "fever", "head", "rash",
           "arth+conj", "arth+fever", "arth+head", "arth+rash", "conj+fever",
           "conj+head", "conj+rash", "fever+head", "fever+rash", "head+rash",
           "arth+conj+fever", "arth+conj+head", "arth+conj+rash", "arth+fever+head", "arth+fever+rash",
           "arth+head+rash", "conj+fever+head", "conj+fever+rash", "conj+head+rash", "fever+head+rash",
           "arth+conj+fever+head", "arth+fever+head+rash", "arth+conj+head+rash",
           "arth+conj+fever+rash", "conj+fever+head+rash", 
          "arth+conj+fever+head+rash", "2 or more", "3 or more", "rash + 1"))

# Caculating the number of rows in the data table

dt.rows = length(incidences) * n.samples * length(pops) * length(types)
  
# Creating a data table with all the combinations of iteration, type, incidence and population

full <- data.table(expand.grid(iter=1:n.samples, type=types, incidence = incidences, population=pops))

# Adding general variables

testing.vars <- data.table(
  iter = 1:n.samples,
  
  # Pregnancy variables
  sen.ELISA           = runif(n.samples, 0.80, 0.99), # sensitivity of IgM MAC ELISA
  spec.ELISA          = runif(n.samples, 0.80, 0.95), # specificity of IgM MAC ELISA
  detection.days.ab   = runif(n.samples, 56, 112),    # days in which IgM antibodies are detectible 
  preg.rate           = runif(n.samples, .009, .017), # 10-17 per 1,000 in a population per year

  # Blood variables
  sen.NAAT            = runif(n.samples, 0.98, 1), # sensitivity of NAAT test 
  spec.NAAT           = runif(n.samples, .9999, 1),       # specificity of NAAT test
  detection.days.v    = runif(n.samples, 11, 17),   # days in which virus is detectible (in serum) 
  blood.donation.rate = runif(n.samples, .043, .047),   # 43 per 1,000 in a population per year
                                                 # whole blood and apheresis red blood cell units
  p.asymptomatic      = runif(n.samples, .6, .8),
    
  # Symptom variables
  sen.PCR              = runif(n.samples, 0.80, 0.95), # sensitivity of PCR
  spec.PCR             = runif(n.samples, 0.99, 1),    # specificity of PCR
  p.z.symp.seek.er     =                               # What % of people are Zika infected have symptoms and seek care at an ED?
                         runif(n.samples, .20, .40) *  # What % of people who are Zika infected have symptoms?
                         runif(n.samples, .1, .5)   *  # What % of people who are Zika infected w/ symptoms seek care?
                         runif(n.samples, .05, .50),   # What % of people visit the emergency department in a given week?
  p.ed.visit           = runif(n.samples, .007, .010)
  )

# Adding data on surveillance system specific assumptions

full <- merge(full, testing.vars, by="iter")

# Adding data on emergency department use by syndrome

ed.syndromes <- data.table(read.csv("CSV files/ED_Symptoms.csv"), key = "type")

# Adding data on the prevalence of symptoms among Zikv+ individuals who sought care

zika.symptoms <- data.table(read.csv("CSV files/Zika_Symptoms.csv"), key="type")

# Adding emergency department syndrome distributions to the dataset

full <- merge(full, ed.syndromes, by="type", all.x = TRUE) 

full <- full[, p.ed.syndrome := suppressWarnings(runif(dt.rows, l95, u95))]
full <- full[, c("l95", "u95") := NULL ]
dim(filter(full, is.na(p.ed.syndrome) & !(type %in% c("pregnant", "blood"))))[1] # Should be 0

# Adding Zika symptom distributions to the dataset

full <- merge(full, zika.symptoms, by="type", all.x = TRUE)
full <- full[, p.z.symp.seek.er.syndrome := suppressWarnings(runif(dt.rows, l95, u95))]
full <- full[, c("l95", "u95") := NULL ]
dim(filter(full, is.na(p.ed.syndrome) & !(type %in% c("pregnant", "blood"))))[1] # Should be 0

# Creating variables where type == 'pregnant'

full[type == "pregnant", `:=` (
   
      # Prevalence of pregnant women with IgM antibodies at a given time
      detectable.zika.prevalence = incidence * detection.days.ab / 7,
                 
      # Number of new pregnancies in a week
      weekly.pregnancies = preg.rate / 52 * population
      
      )] 

full[type == "pregnant", `:=` (
  
      # Weekly tests on Zikv positive people
      weekly.zikv.ppl.tested = weekly.pregnancies*2*detectable.zika.prevalence,
                 
      # Weekly tests on Zikv negative people
      weekly.notzikv.ppl.tested = weekly.pregnancies*2*(1-detectable.zika.prevalence)

      )] 

full[type == "pregnant", `:=` (
  
      # Probability of detecting Zika (in a given week) if testing all pregnant women twice
      prob.detect.week = 1-(1-sen.ELISA*weekly.zikv.ppl.tested/population)^population, 
      
      num.zikv.ppl.detected = weekly.zikv.ppl.tested * sen.ELISA,
      perc.zikv.ppl.detected = weekly.zikv.ppl.tested * sen.ELISA / (population * incidence), 
      num.zikv.per.week = (population * incidence),
                 
      # Probability of false positive (in a given week) if testing all pregnant women twice
      prob.fp.week = 1 - spec.ELISA^weekly.notzikv.ppl.tested,
                 
      # Number of IgM tests needed per week
      tests.per.week = weekly.pregnancies * 2,
      
      # Expected number of false positives (in a given week)
      # n = expected number of tests on ZIKAV- people , p = probability of false positive on a single test
      expected.false.positives = weekly.notzikv.ppl.tested * (1-spec.ELISA) #based on E(v) of binomial distribution) 

      )]

# Removing variables to conserve memory

full[, c("sen.ELISA", "spec.ELISA", "detection.days.ab", "preg.rate") := NULL]

      #-------------------------------------------------------------------------------------------#           

# Creating variables where type == 'blood'

full[type == "blood", `:=` (   
        
      # Prevalence of people in blood bank with detectable virus at a given time
      detectable.zika.prevalence = incidence * (detection.days.v / 7) * p.asymptomatic ,
                                  
      # Number of blood donations in a week
      weekly.blood.donors = blood.donation.rate * population / 52, 
                                  
      # Number of NAAT tests needed per week
      tests.per.week = blood.donation.rate * population / 52
      
      )]
      
full[type == "blood", `:=` (       
                                 
      # Weekly tests on Zikv positive people
      weekly.zikv.ppl.tested = weekly.blood.donors*detectable.zika.prevalence,
                                  
      # Weekly tests on Zikv negative people
      weekly.notzikv.ppl.tested = weekly.blood.donors*(1-detectable.zika.prevalence)
      
      )]

full[type == "blood", `:=` ( 
                                  
      # Probability of detecting Zika (in a given week) if testing all blood donors
      prob.detect.week = 1-(1-sen.NAAT*weekly.zikv.ppl.tested/population)^population, 
      
      num.zikv.ppl.detected = weekly.zikv.ppl.tested * sen.NAAT,
      perc.zikv.ppl.detected = (weekly.zikv.ppl.tested * sen.NAAT) / (population * incidence), 
      num.zikv.per.week = population * incidence,
      
      # Probability of false positive (in a given week) if testing all blood donors
      prob.fp.week = 1 - spec.NAAT^weekly.notzikv.ppl.tested,
                                  
      # Expected number of false positives (in a given week)
      # n = weekly.notzikv.ppl.tested, #expected number of tests on ZIKAV- people 
      # p = (1-spec.NAAT), #probability of false positive on a single test
      expected.false.positives = weekly.notzikv.ppl.tested * (1-spec.NAAT) # based on E(v) of binomial distribution 
        
      )]

# Removing variables to conserve memory

full[, c("sen.NAAT", "spec.NAAT", "detection.days.v", "blood.donation.rate", "p.asymptomatic") := NULL]

      #-------------------------------------------------------------------------------------------#          

# Creating variables for symptomatic types

full[type != "blood" & type != "pregnant", `:=` ( 
  
      # Number of ZIKV infected people tested per week
      weekly.zikv.ppl.tested = population * incidence * 
      p.z.symp.seek.er * p.z.symp.seek.er.syndrome, 
                                 
      # Number of ZIKV negative people tested per week
      weekly.notzikv.ppl.tested = population * (1-incidence) * p.ed.visit * p.ed.syndrome

      )]

# Removing variables to conserve memory

full[, c("p.z.symp.seek.er", "p.z.symp.seek.er.syndrome", "p.ed.visit", "p.ed.syndrome") := NULL]

# Creating variables for symptomatic types

full[type != "blood" & type != "pregnant", `:=` ( 
  
      # Probability of detecting Zika (in a given week) if testing all symptomatic people (w/ sympx)
      prob.detect.week = 1 - (1 -sen.PCR*weekly.zikv.ppl.tested/population)^population,
      
      num.zikv.ppl.detected = weekly.zikv.ppl.tested * sen.PCR,
      perc.zikv.ppl.detected = (weekly.zikv.ppl.tested * sen.PCR) / (population * incidence), 
      num.zikv.per.week = (population * incidence),
      
      # Probability of false positive (in a given week) if testing all symptomatic people (w/ sympx)
      prob.fp.week = 1 - spec.PCR^weekly.notzikv.ppl.tested,
      
      # Expected number of false positives per week
      expected.false.positives = weekly.notzikv.ppl.tested * (1 - spec.PCR)
      
      )]

# Removing variables to conserve memory

full[, c("sen.PCR", "spec.PCR") := NULL]

# Creating variables for symptomatic types

full[type != "blood" & type != "pregnant", `:=` ( 
  
      tests.per.week = weekly.zikv.ppl.tested + weekly.notzikv.ppl.tested
      
      )]

# Keeping important variables

full <- full[, list(type, population, incidence, prob.detect.week, prob.fp.week, tests.per.week,
                      expected.false.positives, num.zikv.ppl.detected , perc.zikv.ppl.detected, num.zikv.per.week)]

# Creating surveillance variable

full[type == "pregnant", surveillance := "pregnant"]
full[type == "blood", surveillance := "blood"]
full[type != "pregnant" &  type != "blood", surveillance := "symptom"]

# PPV

full <- mutate(full, PPV = num.zikv.ppl.detected / (num.zikv.ppl.detected + expected.false.positives))

# Calculating probability of detection for each incidence, population, and type
summary.stats <- full %>%
  group_by(incidence, population, type, surveillance) %>%
  summarise(p.detect.m = median(prob.detect.week),
            p.detect.l95 = quantile(prob.detect.week, .025),
            p.detect.u95 = quantile(prob.detect.week, .975),
            p.detect.l50 = quantile(prob.detect.week, .25),
            p.detect.u50 = quantile(prob.detect.week, .75),
            tests.m = median(tests.per.week),
            tests.l95 = quantile(tests.per.week, .025),
            tests.u95 = quantile(tests.per.week, .975),
            tests.l50 = quantile(tests.per.week, .25),
            tests.u50 = quantile(tests.per.week, .75),
            fp.m = median(expected.false.positives),
            fp.l95 = quantile(expected.false.positives, .025),
            fp.u95 = quantile(expected.false.positives, .975),
            fp.l50 = quantile(expected.false.positives, .25),
            fp.u50 = quantile(expected.false.positives, .75),
            
            perc.zikv.ppl.detected.m = quantile(perc.zikv.ppl.detected, .5),
            perc.zikv.ppl.detected.l95 = quantile(perc.zikv.ppl.detected, .025),
            perc.zikv.ppl.detected.u95 = quantile(perc.zikv.ppl.detected, .975),
            perc.zikv.ppl.detected.l50 = quantile(perc.zikv.ppl.detected, .25),
            perc.zikv.ppl.detected.u50 = quantile(perc.zikv.ppl.detected, .75),
            
            num.zikv.ppl.detected.m = quantile(num.zikv.ppl.detected, .5),
            num.zikv.ppl.detected.l95 = quantile(num.zikv.ppl.detected, .025),
            num.zikv.ppl.detected.u95 = quantile(num.zikv.ppl.detected, .975),
            num.zikv.ppl.detected.l50 = quantile(num.zikv.ppl.detected, .25),
            num.zikv.ppl.detected.u50 = quantile(num.zikv.ppl.detected, .75),
            
            
            num.zikv.per.week.m = quantile(num.zikv.per.week, .5),
            
            PPV.m = quantile(PPV, .5),
            PPV.l95 = quantile(PPV, .025),
            PPV.u95 = quantile(PPV, .975),
            PPV.l50 = quantile(PPV, .25),
            PPV.u50 = quantile(PPV, .75)
            
  ) %>%
  ungroup()

summary.stats <- mutate(summary.stats,
                        log10.incidence = log10(incidence))

# Calculating probability of detection for each incidence, population, and syndrome

summary.stats2 <- full %>%
  filter(type != "blood" & type != "pregnant") %>%
  group_by(incidence, population, type) %>%
  summarise(p.detect.m = median(prob.detect.week),
            p.detect.l95 = quantile(prob.detect.week, .025),
            p.detect.u95 = quantile(prob.detect.week, .975),
            p.detect.l50 = quantile(prob.detect.week, .25),
            p.detect.u50 = quantile(prob.detect.week, .75),
            tests.m = median(tests.per.week),
            tests.l95 = quantile(tests.per.week, .025),
            tests.u95 = quantile(tests.per.week, .975),
            tests.l50 = quantile(tests.per.week, .25),
            tests.u50 = quantile(tests.per.week, .75),
            fp.m = median(expected.false.positives),
            fp.l95 = quantile(expected.false.positives, .025),
            fp.u95 = quantile(expected.false.positives, .975),
            fp.l50 = quantile(expected.false.positives, .25),
            fp.u50 = quantile(expected.false.positives, .75),
            
            perc.zikv.ppl.detected.m = quantile(perc.zikv.ppl.detected, .5),
            num.zikv.ppl.detected.m = quantile(num.zikv.ppl.detected, .5)
  ) %>%
  ungroup()

summary.stats2 <- mutate(summary.stats2,
                         log10.incidence = log10(incidence))

#------------------------------------- Optional: Save resulting datasets ------------------------------------------------------#

#save(full, file='')
#save(summary.stats, file='')
#save(summary.stats2, file='')
