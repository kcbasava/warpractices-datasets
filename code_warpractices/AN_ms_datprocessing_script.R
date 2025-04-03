library(ape)
library(caper)
library(bayestraitr)
library(stringr)
library(dplyr)
setwd("~/appendix/ch4_AN")

#read in tree
ch4_ANphy <- read.nexus("an.treeS_hh_new.nex") #set of 1000 phylogenies from Gray et al. 
gray2009_summary.trees <- read.tree("summary.trees") #consensus tree for phylogenetic signal analysis

#read in data
ch4dat <- read.csv("ch4_data.csv")

length(ch4dat[ch4dat$headhunting == 0,]$headhunting) #74
length(ch4dat[ch4dat$headhunting == 1,]$headhunting) #46

#signal analysis----
hh.sig <- phylo.d(ch4dat, gray2009_summary.trees, names.col=taxon, binvar=headhunting, permut=1000)
summary(hh.sig)

setwd("/Appendix_2/appendix_chapter4")
#multistate analysis for ancestral state
ms_logs <- bt_read.log("hh.txt.log.txt")
View(ms_logs)

asr_columns = str_detect(colnames(ms_logs), "Root") 
asr <- colMeans(ms_logs[,asr_columns])
asr

ms_columns = str_detect(colnames(ms_logs), "q") 
hhtransitions <- colMeans(ms_logs[,ms_columns])
View(hhtransitions)

#analyzing logs----
#independent, agriculture----
agr_ind_logs = bt_read.log("IND_agr.Log.txt")

agr_ind_rates_columns = str_detect(colnames(agr_ind_logs), "alpha|beta") 

colSums(is.na(agr_ind_logs[,agr_ind_rates_columns]))
agr_ind_rates <- colMeans(agr_ind_logs[,agr_ind_rates_columns])
View(agr_ind_rates)

agr_abs <- c("gain hh", "lose hh", "gain agriculture", "lose agriculture")
View(agr_abs)

#dependent, agriculture----
agr_dep_logs = bt_read.log("DEP_agr.Log.txt")
agr_dep_rates = str_detect(colnames(agr_dep_logs), "q")
colSums(is.na(agr_dep_logs[,agr_dep_rates]))
agr_dep_rates <- colMeans(agr_dep_logs[,agr_dep_rates])
agr_dep_rates <- as.vector(agr_dep_rates)
View(agr_dep_rates)
class(agr_dep_rates)

qs_agr <- c("no hh, gain agriculture",
            "gain hh, no agriculture", 
            "no hh, lose agriculture",
            "gain hh, keep agriculture", 
            "lose hh, no agriculture", 
            "keep hh, gain agriculture",
            "lose hh, keep agriculture",
            "keep hh, lose agriculture")
qs_agr <- as.data.frame(qs)
View(qs_agr)

qs_agr <- cbind(qs_agr, agr_dep_rates)

#independent model stones lml -143.828470
#dependent model stones lml -142.198810 

2*(-142.198810 - -143.828470)
#BF = 3.26

#for alternative models (stricter hh def)
#independent model stones lml -138.738979
#dependent model stones lml -135.321978
2*(-135.321978 - -138.738979)
#BF = 6.83 

agr_long_df = agr_dep_logs %>% 
  data.frame(.) %>%
  dplyr::select(q12, q13, q21, q24, q31, q34, q42, q43) %>% 
  reshape2::melt()

ggplot(agr_long_df, aes(y = value, x = variable)) + 
  geom_boxplot()

plot(agr_dep_logs$Lh, type = "l")

#independent, pol.bin2----
pol2_ind_logs = bt_read.log("IND_pol2.Log.txt")
pol2_ind_rates_columns = str_detect(colnames(pol2_ind_logs), "alpha|beta") 

colSums(is.na(pol2_ind_logs[,pol2_ind_rates_columns]))
pol2_ind_rates <- colMeans(pol2_ind_logs[,pol2_ind_rates_columns])
View(pol2_ind_rates)

pol2_abs <- c("gain hh", "lose hh", "gain complexity", "lose complexity")
View(pol2_abs)

#dependent, pol.bin2
pol2_dep_logs = bt_read.log("DEP_pol2.Log.txt")

pol2_dep_rates = str_detect(colnames(pol2_dep_logs), "q")

colSums(is.na(pol2_dep_logs[,pol2_dep_rates])) 
pol2_dep_rates <- colMeans(pol2_dep_logs[,pol2_dep_rates])
View(pol2_dep_rates)

qs_pol2 <- c("no hh, gain complexity", "gain hh, no complexity", "no hh, lose complexity", "gain hh, retain complexity", "lose hh, no complexity", "keep hh, retain complexity", "lose hh, retain complexity", "keep hh, lose complexity")
View(qs_pol2)

#independent model stones lml -137.479502
#dependent model stones lml -139.055402
2*(-137.479502 - -139.055402) 
#BF 3.1518 in favor of independent model

#for alternative models (stricter hh def)
#independent model stones lml -132.292999
#dependent model stones lml -133.819796
2*(-132.292999 - -133.819796)
#BF 3.053594 in favor of independent model

#independent, pol.bin3----
pol3_ind_logs = bt_read.log("IND_pol3.Log.txt")
pol3_ind_rates_columns = str_detect(colnames(pol3_ind_logs), "alpha|beta") 

colSums(is.na(pol3_ind_logs[,pol3_ind_rates_columns]))
pol3_ind_rates <- colMeans(pol3_ind_logs[,pol3_ind_rates_columns])
View(pol3_ind_rates)

pol3_abs <- c("gain hh", "lose hh", "gain complexity", "lose complexity")
View(pol3_abs)

#dependent, pol.bin3
pol3_dep_logs = bt_read.log("DEP_pol3.Log.txt")

pol3_dep_rates = str_detect(colnames(pol3_dep_logs), "q")

colSums(is.na(pol3_dep_logs[,pol3_dep_rates])) 
pol3_dep_rates <- colMeans(pol3_dep_logs[,pol3_dep_rates])
View(pol3_dep_rates)
pol3_dep_rates <- c(pol3_dep_rates)

qs_pol3 <- c("no hh, gain complexity",
             "gain hh, no complexity",
             "no hh, lose complexity",
             "gain hh, keep complexity",
             "lose hh, no complexity",
             "keep hh, gain complexity",
             "lose hh, keep complexity",
             "keep hh, lose complexity")
View(qs_pol3)

qs_pol3 <- cbind(qs_pol3, pol3_dep_rates)
View(qs_pol3)

pol3_long_df = pol3_dep_logs %>% 
  data.frame(.) %>%
  dplyr::select(q12, q13, q21, q24, q31, q34, q42, q43) %>% 
  reshape2::melt()

ggplot(pol3_long_df, aes(y = value, x = variable)) + 
  geom_boxplot()

#independent model stones lml -123.575012
#dependent model stones lml -122.937065
2*(-122.937065 - -123.575012)
#BF 1.275894

#for alternative models (stricter hh def)
#independent model stones lml -117.079791
#dependent model stones lml -117.231802
2*(-117.079791 - -117.231802)
#BF -0.304022