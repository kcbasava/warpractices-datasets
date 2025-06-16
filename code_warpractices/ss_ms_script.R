getwd()

library(brms)
library(cmdstanr)
library(dagitty)
library(ggdag)
library(sjPlot)
library(mice)
library(here)

setwd(here("code_warpractices"))

#DAGs
ss_dag <- dagify(sacrifice ~ complexity + military + subsistence,
                 complexity ~ subsistence,
                 military ~ complexity,
                 labels=c("sacrifice"="self-sacrifice", "complexity"="social\n complexity", "subsistence"="agricultural dependence", "military"="military formalization"
                 ))
drawdag(ss_dag)

adjustmentSets(ss_dag, exposure="complexity", outcome="sacrifice") 
#to estimate effect of complexity on sacrifice, control for subsistence

adjustmentSets(ss_dag, exposure="subsistence", outcome="sacrifice") 
#none

adjustmentSets(ss_dag, exposure="military", outcome="sacrifice") 
#to estimate effect of military on sacrifice, control for complexity

#read in data----
ss_orig <- read.csv("dataset1.1_preimputation.csv")
ss_orig$agr.st <- standardize(ss_orig$agr_ins)

#imputation----
ss_imp <- complete(mice(ss_orig, method="pmm"))

#PCA for political complexity----
ss.pca <- prcomp(data.frame(C=standardize(ss_imp$pol_cent), P=standardize(ss_imp$pop_res1), H=standardize(ss_imp$ext_hier), S=standardize(ss_imp$commsize_scale), row.names=ss_imp$society))
screeplot(ss.pca, type="lines")
summary(ss.pca)

#proportion of variance explained by each component
pve_ss <- (ss.pca$sdev^2)/sum(ss.pca$sdev^2) 
#1st component = 0.81

#and plot
plot(pve_ss, xlab = "PC",
     ylab = "Proportion of variance",
     ylim = c(0, 1), type = "b")

df.pcs_ss <- data.frame(ss.pca$x) 
ss_imp$PC1 <- df.pcs_ss$PC1 #add to imputated dataset
#check correlation w/ population
cor(ss_imp$pop_res1, ss_imp$PC1) 

#sanity checks
summary(lm(PC1 ~ agr.st, data=ss_imp))
summary(lm(PC1 ~ as.factor(mil_org), data=ss_imp))

#read in imputed and standardized data----
ss_imp <- read.csv("dataset1.2_imputed.csv")

ss_imp$sac_scale <- as.ordered(ss_imp$sac_scale)
ss_imp$mil_org <- as.ordered(ss_imp$mil_org)

#analyses----
#with all predictors
ss_PC1am <- brm(data=ss_imp, family=cumulative,
               sac_scale ~ 1 + (1|Region) + PC1 + agr + mil_org,
               prior = c(prior(normal(0, 1.5), class = Intercept),
                         prior(normal(0, 0.5), class = b)),
               iter = 4000, warmup = 2000, chains = 4, cores = 4,
               backend="cmdstanr")
print(ss_PC1am)



#complexity----
ss_PC1 <- brm(data=ss_imp, family=cumulative,
                  sac_scale ~ 1 + (1|Region) + PC1,
                  prior = c(prior(normal(0, 1.5), class = Intercept),
                            prior(normal(0, 0.5), class = b)),
                  iter = 4000, warmup = 2000, chains = 4, cores = 4,
                  backend="cmdstanr")
print(ss_PC1)

#controlling for agriculture
ss_PC1a <- brm(data=ss_imp, family=cumulative,
                   sac_scale ~ 1 + (1|Region) + PC1 + agr,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 0.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 4, cores = 4,
                   backend="cmdstanr")
print(ss_PC1a)

#agriculture
ss_a <- brm(data=ss_imp, family=cumulative,
                sac_scale ~ 1 + (1|Region) + agr,
                prior = c(prior(normal(0, 1.5), class = Intercept),
                          prior(normal(0, 0.5), class = b)),
                iter = 4000, warmup = 2000, chains = 4, cores = 4,
                backend="cmdstanr")
print(ss_a)

#military----
ss_PC1m <- brm(data=ss_imp, family=cumulative,
                   sac_scale ~ 1 + (1|Region) + PC1 + mil_org,
                   prior = c(prior(normal(0, 1.5), class = Intercept),
                             prior(normal(0, 0.5), class = b)),
                   iter = 4000, warmup = 2000, chains = 4, cores = 4,
                   backend="cmdstanr")
print(ss_PC1m)

ss_m <- brm(data=ss_imp, family=cumulative,
               sac_scale ~ 1 + (1|Region) + mil_org,
               prior = c(prior(normal(0, 1.5), class = Intercept),
                         prior(normal(0, 0.5), class = b)),
               iter = 4000, warmup = 2000, chains = 4, cores = 4,
               backend="cmdstanr")
print(ss_m)

#compare models----
model_weights(ss_PC1am, ss_PC1, ss_a, ss_PC1a, ss_PC1m, ss_m, weights="loo")
model_weights(ss_PC1am, ss_PC1, ss_PC1a, ss_PC1m, weights="loo")


tab_model(ss_PC1a, ss_PC1am, ss_PC1m, transform=NULL, show.intercept=FALSE, show.re.var=FALSE, dv.labels=c("m_PC1a", "m_PC1am", "m_PC1m"))
