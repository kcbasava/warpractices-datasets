library(brms)
library(cmdstanr)
library(dagitty)
library(ggdag)
library(sjPlot)
library(mice)
library(ICC)
library(here)
options(scipen=999) #turn off scientific notation

here()

#DAGs----
enemy_dag <- dagify(kill_ext ~ mil.org + exp + pol.complx,
                    exp ~ pol.complx,
                    mil.org ~ pol.complx,
                    labels=c("kill_ext"="indiscriminate\n killing", "mil.org"="military\n formalization", "pol.complx"="social\n complexity", "exp"="political\n expansion"
                    ))

plot(enemy_dag)

#backdoor paths
adjustmentSets(enemy_dag, exposure="pol.complx", outcome="kill_ext") 
#empty set, so kill ~ PC

adjustmentSets(enemy_dag, exposure="mil.org", outcome="kill_ext") 
#political complexity, so kill ~ military + PC

adjustmentSets(enemy_dag, exposure="exp", outcome="kill_ext") 
#political complexity, so kill ~ expansion + PC

impliedConditionalIndependencies(enemy_dag) 
#expansion and military organization independent of each other conditional on political complexity

#read in data----
ch2_preimp <- read.csv(here("inkill_preimputation.csv"))
ch2_imp <- ch2_preimp[c(1,4,9,10,13,14,15,18,19)]

#impute (default should be predictive mean matching)----
ch2_mice <- mice(ch2_imp)
ch2_imp <- complete(ch2_mice)

#sanity checks----
summary(lm(ch2_imp$pop_res1~as.factor(ch2_imp$ter_res))) #societies with higher populations more likely to have political expansion
summary(lm(pop_res1~ as.factor(mil.res), data=ch2_imp)) #and more likely to have formal militaries 
summary(lm(pop_res1 ~ pol_cent, data=ch2_imp)) #population and centralization correlated
summary(lm(pop_res1 ~ ext_hier, data=ch2_imp)) #population and hierarchy correlated
summary(lm(pol_cent ~ ext_hier, data=ch2_imp)) #hierarchy and centralization correlated

#PCA for political complexity----
ch2pc.list <- list(cent=rethinking::standardize(ch2_imp$pol_cent), pop=rethinking::standardize(ch2_imp$pop_res1), hier=rethinking::standardize(ch2_imp$ext_hier))
ch2pc_df <- data.frame(ch2pc.list$cent, ch2pc.list$pop, ch2pc.list$hier, row.names=ch2_imp$society)
ch2pc <- prcomp(ch2pc_df)
biplot(ch2pc)
screeplot(ch2pc, type="lines")
ch2pc$sdev^2 #lambda for PC1 > 1
pve.ch2 <- (ch2pc$sdev^2)/sum(ch2pc$sdev^2) #proportion of variance explained by each PC
summary(ch2pc) #PC1 explains 0.88 of variance

ch2df.pcs <- data.frame(ch2pc$x) #save PCAs in a dataframe
ch2_imp$PC1 <- -1*(ch2df.pcs$PC1)
ch2_imp$PC1.st <- rethinking::standardize(ch2_imp$PC1)

#data list for analyses----
#can't have sd = 0 for the brm, so adding making se for societies w/ no error se = 0.00000001
se_01 <- ifelse (ch2_imp$kill_se == 0.0, ch2_imp$kill_se + 0.00000000000001, ch2_imp$kill_se)
ch2_imp$k_sd <- se_01 / sd(ch2_imp$kill_mean)

ch2_prioronly <-  brm(data = ch2_imp, family = gaussian,
                      kill_mean | se(k_sd) ~ PC1.st + ter_res,
                      prior = c(prior(normal(4, 1.3), class = Intercept),
                                prior(normal(0, 0.5), class=b)),
                      sample_prior = "only", 
                      iter = 2000, warmup = 1000, cores = 2, chains = 2,
                      backend="cmdstanr")

summary(ch2_prioronly)

#complexity only model
m2_p <- 
  brm(data = ch2_imp, family = gaussian,
      kill_mean | se(k_sd, sigma=TRUE) ~ PC1.st,
      prior = c(prior(normal(4, 1.3), class = Intercept),
                prior(normal(0, 0.5), class=b),
                prior(exponential(1), class=sigma)), 
      iter = 4000, warmup = 2000, cores = 3, chains = 3,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 13),
      backend="cmdstanr")
summary(m2_p) #CI centered around 0, likely no effect
mcmc_plot(m2_p, variable="b_PC1.st", type="intervals")

#effect of military after accounting for complexity and complexity -> exp + effect of comp including effects through exp, independent of effects thorugh military
m2_pm <- brm(data = ch2_imp, family = gaussian,
             kill_mean | se(k_sd, sigma=TRUE) ~ PC1.st + mil_res,
             prior = c(prior(normal(4,1.3), class = Intercept),
                       prior(normal(0,0.7), class=b),
                       prior(exponential(1), class=sigma)),
             iter = 4000, warmup = 2000, cores = 3, chains = 3,
             control = list( adapt_delta = 0.99,
                             max_treedepth = 13),
             backend="cmdstanr")
summary(m2_pm)
mcmc_plot(m2_pm, variable=c("b_PC1.st", "b_mil_res"), type="intervals")

#expansion, controlling for complexity
m2_ep <- brm(data = ch2_imp, family = gaussian,
             kill_mean | se(k_sd, sigma=TRUE) ~ ter_res + PC1.st,
             prior = c(prior(normal(4,1.3), class = Intercept),
                       prior(normal(0,0.7), class=b),
                       prior(exponential(1), class=sigma)),
             iter = 4000, warmup = 2000, cores = 3, chains = 3,
             control = list( adapt_delta = 0.99,
                             max_treedepth = 13),
             backend="cmdstanr")
summary(m2_ep) #expansion positive but CI crosses 0
mcmc_plot(m2_ep, variable=c("b_PC1.st", "b_ter_res"), type="intervals")

#all vars
m2_all <- brm(data = ch2_imp, family = gaussian,
              kill_mean | se(k_sd, sigma=TRUE) ~ PC1.st + ter_res + mil_res,
              prior = c(prior(normal(4,1.3), class = Intercept),
                        prior(normal(0,0.7), class=b),
                        prior(cauchy(0,5), class = sigma)),
              iter = 4000, warmup = 2000, cores = 3, chains = 3,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 13),
              backend="cmdstanr") 
summary(m2_all)
mcmc_plot(m2_all, variable=c("b_PC1.st", "b_ter_res", "b_mil_res"), type="intervals")

#model comparison----
model_weights(m2_p, m2_ep, m2_pm, m2_all, weights="loo")
model_weights(m2_all, m2_ep, m2_p, weights="loo")

tab_model(m2_ep, m2_pm, m2_p,  transform=NULL, dv.labels=c("m_complexity_expansion", "m_complexity_military", "m_complexity_only"))

tab_model(m2_all, m2_ep, m2_p,  transform=NULL, dv.labels=c("m_all", "m_expansions", "m_complexity_only"))


#regional autocorrelation----
icc_kill <- ch2_imp[c(1,2,8)]
ICCbare(Region, kill_mean, data=ch2_imp) #0.01
ICCbare(Region, kill_mean+kill_se, data = ch2_imp) #0.05
ICCbare(Region, kill_mean-kill_se, data = ch2_imp) #-0.02

summary(brm(data = ch2_imp, family = gaussian,
            kill_mean ~ PC1.st + (1|Region),
            prior = c(prior(normal(4,1.3), class = Intercept),
                      prior(normal(0,0.5), class=b)),
            iter = 2000, warmup = 1000, cores = 3, chains = 3,
            backend="cmdstanr"))
#PC1.st = 0.04, 95% CI = -0.42, 0.48

summary(brm(data = ch2_imp, family = gaussian,
            kill_mean ~ PC1.st,
            prior = c(prior(normal(4,1.3), class = Intercept),
                      prior(normal(0,0.5), class=b)),
            iter = 2000, warmup = 1000, cores = 3, chains = 3,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 13),
            backend="cmdstanr"))
#PC1.st = 0.03, 95% CI = -0.41, 0.46

#plots----
hist(ch2_imp$PC1, breaks="FD", border="lightgrey", main="", xlab="PC1")
#case studies PC1----
ch2_imp[32, 10] #kamakura
ch2_imp[50, 10] #nuer
ch2_imp[34, 10] #kapauku
ch2_imp[62, 10] #tlingit
ch2_imp[57, 10] #qarmat
ch2_imp[13, 10] #rustamid

killranges_plot <- 
  ggplot() +
  geom_errorbarh(data=ch2_imp, mapping=aes(y=PC1, xmin=(kill_mean-kill_se), xmax=(kill_mean+kill_se)), height=0.06, size=0.6, color="blue") +
  labs(x='kill_ext ranges by society', y = 'PC1') +
  theme_dark(base_size=17)
killranges_plot

#ranges by society, ordered top to bottom by complexity
ordsoc <- ch2_imp[order(ch2_imp$PC1),]
ggplot(data=ordsoc, aes(y=1:74, x=kill_mean, xmin=(kill_mean-kill_se), xmax=(kill_mean+kill_se))) +
  geom_errorbarh(height=.5) +
  scale_y_continuous(breaks=1:nrow(ordsoc), labels=ordsoc$society) +
  scale_x_continuous(breaks = scales::breaks_width(1) ) +
  labs(x='ranges of enemies killed', y = 'societies ordered by ascending political complexity') +
  theme_minimal(base_size=10)
