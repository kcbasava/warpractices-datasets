options(scipen=999)

library(ggdag)
library(dagitty)
library(brms)
library(cmdstanr)
library(sjPlot)
library(rethinking)

#DAGs----
#dag w/ agr -> expansion + agr -> tt
ttdag_at <- dagify(trophies ~ complexity + expansion + agriculture,
                   complexity ~ agriculture,
                   expansion ~ complexity + agriculture,
                   labels=c("troph"="trophy\n taking", "complx"="political\n complexity", "exp"="political\n expansion", "agr"="agriculture"
                   ))
ggdag(ttdag_at, text = FALSE, use_labels = "label", edge_type = "diagonal", text_size = 6,
      label_size = text_size, label_col = "blue", node=F, text_col="black") + theme_dag_blank()
drawdag(ttdag_at)
plot(ttdag_at)
dagitty

adjustmentSets(ttdag_at, exposure="complexity", outcome="trophies") 
#when predicting trophy-taking from complexity, agriculture is fork backdoor path 

adjustmentSets(ttdag_at, exposure="expansion", outcome="trophies") 
#complexity and agriculture are backdoor paths

#dag w/ agr !-> tt
ttdag_cet <- dagify(trophies ~ complexity + expansion,
                    complexity ~ agriculture,
                    expansion ~ complexity + agriculture,
                    labels=c("troph"="trophy\n taking", "complx"="political\n complexity", "exp"="political\n expansion", "agr"="agriculture"
                    ))
plot(ttdag_cet)
adjustmentSets(ttdag_cet, exposure="complexity", outcome="trophies") 
#when predicting trophy-taking from complexity, agriculture is fork backdoor path 

adjustmentSets(ttdag_cet, exposure="expansion", outcome="trophies") 
#complexity only is backdoor paths

#dag w/ exp !-> tt
ttdag_ct <- dagify(trophies ~ complexity,
                   complexity ~ agriculture,
                   expansion ~ complexity + agriculture,
                   labels=c("troph"="trophy\n taking", "complx"="political\n complexity", "exp"="political\n expansion", "agr"="agriculture"
                   ))
drawdag(ttdag_ct)
adjustmentSets(ttdag_ct, exposure="complexity", outcome="trophies") 
#none

adjustmentSets(ttdag_ct, exposure="expansion", outcome="trophies") 
#complexity only is backdoor path

#read in data----
trophies_clean <- read.csv(here(file="trophies_clean.csv"))

#match to SPC measures from turchin 2018
# turchin.trim <- turchin2018
# colnames(turchin.trim)[2] <- "society"
# turchin.trim <- merge(trophies_clean, turchin.trim[c(2,3,4,13)], by="society", all.x=T, all.y=F)
# turchin.trim <- unique(turchin.trim)
# turchin.trim <- turchin.trim[-c(16)]
# turchin.trim <- unique(turchin.trim)
#wrote to csv and manually averaged the multiple ones/match dates
#write.csv(turchin.trim, "~/dphil_evanth/datasets/turchinmatch.csv")
trophies231022 <- read.csv(here("turchinmatch.csv"))

#split into societies w/ turchin2018 SPC and those without----
troph_sesh <- subset(trophies231022, trophies231022$SPC1 != 'NA')
troph_nosesh <- subset(trophies231022, is.na(trophies231022$SPC1)==TRUE)
troph_nosesh <- troph_nosesh[-c(15,16)]


#PCA for non-seshat----
tt.pca1 <- prcomp(data.frame(C=standardize(troph_nosesh$pol.cent), P=standardize(troph_nosesh$pop.res1), H=standardize(troph_nosesh$ext_hier), row.names=troph_nosesh$society))
screeplot(tt.pca1, type="lines")
summary(tt.pca1)

#proportion of variance explained by each component
pve_tt1 <- (tt.pca1$sdev^2)/sum(tt.pca1$sdev^2) 
pve_tt1 
#1st component = 0.85

#and plot
plot(pve_tt1, xlab = "PC",
     ylab = "Proportion of variance",
     ylim = c(0, 1), type = "b")

ttpca1df <- as.data.frame(tt.pca1$x)

troph_nosesh$PC1 <- ttpca1df$PC1 #add to dataset
troph_nosesh$PC1 <- -1*(troph_nosesh$PC1)
#check correlation w/ population
cor(troph_nosesh$pop.res1, troph_nosesh$PC1) 
#0.87

troph_nosesh$PC1.st <- standardize(troph_nosesh$PC1)

#this is complete dataset, 06.11.22----
# trophies61122 <- merge(trophies231022, troph_nosesh[c(1,16)], by="society", all.x=TRUE, all.y=TRUE)
# trophies61122$pop.st <- standardize(trophies61122$pop.res1)
# trophies61122$SPC1.st <- standardize(trophies61122$SPC1)
# write.csv(trophies61122, "~/dphil_evanth/datasets/trophies61122.csv")

trophies <- read.csv(file="/Users/kiranbasava/dphil_evanth/Appendix_2/publishing/trophies.csv")

#sanity checks----
summary(glm(ter.res~agr.st, data=trophies_clean, family="binomial")) #agriculture increases likelihood of expansion
summary(glm(ter.res~pop.st, data=trophies_clean, family="binomial")) #population increases likelihood of expansion
summary(lm(pop.st~agr.st, data=trophies_clean)) #agriculture and population positively correlated
summary(lm(PC1.st~agr.st, data=troph_nosesh1)) #agriculture +~ complexity w/in non-seshat societies
summary(lm(SPC1~agr.st, data=troph_sesh)) #agriculture +~ complexity w/in seshat societies
summary(glm(ter.res~PC1, data=troph_nosesh1)) #complexity +~ expansion w/in nonseshat societies
summary(glm(ter.res~SPC1, data=troph_sesh)) #complexity +~ expansion w/in seshat societies

#non-Seshat societies analysis----
#complexity total effect, quadratic, non-seshat data----
m.ttethn_p2a <- brm(data=troph_nosesh, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + PC1.st + I(PC1.st^2) + agr.st,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttethn_p2a) 
report(m.ttethn_p2a)

#expansion, adjusting for complexity, non-seshat data----
m.ttethn_p2e <- brm(data=troph_nosesh, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + PC1.st + I(PC1.st^2) + ter.res,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttethn_p2e) 
#ter.res estimate centered around 0

#complexity total effect, linear, non-seshat data----
m.ttethn_pa <- brm(data=troph_nosesh, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + PC1.st + agr.st,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttethn_pa) 

#expansion effect, linear, non-seshat data----
m.ttethn_pe <- brm(data=troph_nosesh, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + PC1.st + ter.res,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttethn_pe) 

#complexity only, quadratic, non-seshat----
m.ttethn_p2 <- brm(data=troph_nosesh, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + PC1.st + I(PC1.st^2),
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttethn_p2) 

#complexity only, linear, non-seshat----
m.ttethn_p <- brm(data=troph_nosesh, family=bernoulli,
                  trophies.bin ~ 1 + (1|Region) + PC1.st,
                  prior = c(prior(normal(0, 1.5), class = Intercept),
                            prior(normal(0, 1.5), class = b)),
                  iter = 4000, warmup = 2000, chains = 3, cores = 3,
                  backend="cmdstanr")
summary(m.ttethn_p) 

#all vars, non-seshat----
m.ttethn_pae <- brm(data=troph_nosesh, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + agr.st + PC1.st + ter.res,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttethn_pae) 

#non-sesh model comparison----
model_weights(m.ttethn_p2a, m.ttethn_p2, m.ttethn_p2e, m.ttethn_pa, m.ttethn_pe, m.ttethn_p, m.ttethn_pae, weights="loo")
#m.ttethn_p 0.56; m.ttethn_pa 0.14; m.ttethn_pe 0.13
tab_model(m.ttethn_p, m.ttethn_pa, m.ttethn_pe, transform=NULL, dv.labels=c("m_PC1", "m_PC1.agr", "m_PC1.exp"))

#Seshat societies analysis----
troph_sesh$SPC1 <- rethinking::standardize(troph_sesh$SPC1)


#complexity total effect, quadratic, seshat data----
m.ttsesh_p2a <- brm(data=troph_sesh, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + SPC1 + I(SPC1^2) + agr.st,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttsesh_p2a) 

#complexity total effect, linear, seshat data----
m.ttsesh_pa <- brm(data=troph_sesh, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + SPC1 + agr.st,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttsesh_pa) 

#expansion effect, linear, seshat data----
m.ttsesh_pe <- brm(data=troph_sesh, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + SPC1 + ter.res,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttsesh_pe)

#complexity only, linear, seshat data----
m.ttsesh_p <- brm(data=troph_sesh, family=bernoulli,
                  trophies.bin ~ 1 + (1|Region) + SPC1,
                  prior = c(prior(normal(0, 1.5), class = Intercept),
                            prior(normal(0, 1.5), class = b)),
                  iter = 4000, warmup = 1000, chains = 3, cores = 3,
                  backend="cmdstanr")
summary(m.ttsesh_p)

#complexity^2 only, seshat data---
m.ttsesh_p2 <- brm(data=troph_sesh, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + SPC1 + I(SPC1^2),
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 1000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttsesh_p2)

#all vars, seshat data---
m.ttsesh_pae <- brm(data=troph_sesh, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + SPC1 + agr.st + ter.res,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 4000, warmup = 1000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttsesh_pae)

#seshat data model comparison----
model_weights(m.ttsesh_pa, m.ttsesh_p2a, m.ttsesh_p2, m.ttsesh_p, m.ttsesh_pe, m.ttsesh_pae, weights="loo")
#m.ttsesh_p2 0.30; m.ttsesh_p 0.27; m.ttsesh_pe 0.16
tab_model(m.ttsesh_p2, m.ttsesh_p, m.ttsesh_pe, transform=NULL, dv.labels=c("m_PC1^2", "m_PC1", "m_PC1.exp"))

#whole dataset analyses----
#population^2 whole dataset----
m.ttwhole_p2a <- brm(data=trophies, family=bernoulli,
                     trophies.bin ~ 1 + (1|Region) + pop.st + I(pop.st^2) + agr.st,
                     prior = c(prior(normal(0, 1.5), class = Intercept),
                               prior(normal(0, 1.5), class = b)),
                     iter = 4000, warmup = 1000, chains = 3, cores = 3,
                     backend="cmdstanr")
summary(m.ttwhole_p2a) 

#population linear whole dataset----
m.ttwhole_pa <- brm(data=trophies, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + pop.st + agr.st,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 2000, warmup = 1000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttwhole_pa)

#population^2 alone whole dataset----
m.ttwhole_p2 <- brm(data=trophies, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + pop.st + I(pop.st^2) ,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 2000, warmup = 1000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttwhole_p2)

#population alone whole dataset----
m.ttwhole_p <- brm(data=trophies, family=bernoulli,
                   trophies.bin ~ 1 + (1|Region) + pop.st ,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 1.5), class = b)),
                   iter = 4000, warmup = 1000, chains = 3, cores = 3,
                   backend="cmdstanr")
summary(m.ttwhole_p)

#expansion whole dataset----
m.ttwhole_pe <- brm(data=trophies, family=bernoulli,
                    trophies.bin ~ 1 + (1|Region) + pop.st + ter.res,
                    prior = c(prior(normal(0, 1.5), class = Intercept),
                              prior(normal(0, 1.5), class = b)),
                    iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    backend="cmdstanr")
summary(m.ttwhole_pe)

#all vars whole dataset----
m.ttwhole_pae <- brm(data=trophies, family=bernoulli,
                     trophies.bin ~ 1 + (1|Region) + pop.st + agr.st + ter.res,
                     prior = c(prior(normal(0, 1.5), class = Intercept),
                               prior(normal(0, 1.5), class = b)),
                     iter = 4000, warmup = 2000, chains = 3, cores = 3,
                     backend="cmdstanr")
summary(m.ttwhole_pae)

#whole dataset model comparison----
model_weights(m.ttwhole_pa, m.ttwhole_p2a, m.ttwhole_p, m.ttwhole_p2, m.ttwhole_pe, m.ttwhole_pae, weights="loo")
#m.ttwhole_p 0.38; m.ttwhole_pe 0.24; m.ttwhole_pa 0.13; m.ttwhole_p2 0.11
tab_model(m.ttwhole_p, m.ttwhole_pe, m.ttwhole_pa, transform=NULL, dv.labels=c("m_pop", "m_pop.exp", "m_pop.agr"))

#phyloglms----
#distance matrix -> tree----
library(ape)
library(phangorn)
library(geodist)
library(caper)
library(phylolm)

trophies_ll3 <- trophies_clean[c(1,10,11)]
ttmat <- geodist(trophies_ll3, measure="geodesic")
colnames(ttmat) <- rownames(ttmat) <- trophies_ll3[['society']] 
ttmat <- ttmat/1000000
tttree <- upgma(ttmat)
plot(tttree, cex=0.7)

ttcomp <- comparative.data(tttree, crania[c(1,6,9,13,18,20)], names.col=society, warn.dropped=T)
phylo.d(ttcomp, binvar=trophies.bin)
#D = 0.64, significantly different from both Brownian and random distribution

ttphy_p <- phyloglm(formula= trophies.bin ~ pop.st, data=ttcomp$data, phy=ttcomp$phy, boot=100)
summary(ttphy_p)
summary(m.ttwhole_p)

ttphy_pa <- phyloglm(formula= trophies.bin ~ pop.st + agr.st, data=ttcomp$data, phy=ttcomp$phy, boot=100)
summary(ttphy_pa)
summary(m.ttwhole_pa)

ttphy_pae <- phyloglm(formula= trophies.bin ~ pop.st + ter.res + agr.st, data=ttcomp$data, phy=ttcomp$phy, boot=100)
summary(ttphy_pae)
summary(m.ttwhole_pae)

ttphy_pe <- phyloglm(formula= trophies.bin ~ pop.st + ter.res, data=ttcomp$data, phy=ttcomp$phy, boot=100)
summary(ttphy_pe)
summary(m.ttwhole_pe)

#crania----
crphy_a <- phyloglm(formula= crania ~ agr.st, data=ttcomp$data, phy=ttcomp$phy, boot=100)
summary(crphy_a)
summary(m.ttwhole_a)

#violin plot trophy-taking and political centralization
ggplot(ttcomp$data, aes(x=as.factor(crania), y=agr.st)) + 
  geom_violin(scale="count") + labs(x="head-taking", y="agriculture")


#violin plot trophy-taking and political centralization
ggplot(ttcomp$data, aes(x=as.factor(trophies.bin), y=pop.st)) + 
  geom_violin(scale="count") + labs(x="trophy-taking", y="population")
