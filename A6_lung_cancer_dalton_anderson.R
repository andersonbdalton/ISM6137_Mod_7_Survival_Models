#dalton anderson

rm(list=ls())

library(dplyr)
library(ggplot2)
library(corrplot)
library(stargazer)
library(PerformanceAnalytics)
library(stargazer)
library(lme4)



#import data
df_master <- read.table("LungCancer.txt", skip = 15)
#create columns names 
df_master <- df_master %>%
  rename(treatment = V1,
         cell_type = V2,
         survival = V3,
         status = V4,
         karnofsky = V5,
         diagnosis_period =V6,
         age = V7,
         prior_chemo = V8)

#reorder data
df <- df_master %>%
  select(survival,
         treatment,
         prior_chemo,
         cell_type,
         age,
         status,
         diagnosis_period,
         karnofsky)

df$treatment <- relevel(df$treatment, ref ="1")
#kaplan-meier non-parametric
#devtools::install_github("kassambara/survminer", build_vignettes = FALSE)
#library(survminer)
library(survival)

km <- survfit(Surv(survival,status) ~ treatment, data = df) 
summary(km)
#doesn't look like the 'new' medicine works very well at the 6 and 12 month mark.

#survival summary for six months
summary(km, times = 182)

#survival summary for one year
summary(km, times = 365)

#can not install survminer using plot for now
#plot of survival
plot(km, xlab="Days", ylab="Survival Probability")

#feature engineering 

#cell type? After doing some reading there are two types of lung cancer are classified as small cell and non-small cell 
#https://www.cancer.org/cancer/lung-cancer/about/what-is.html
#small cell lung cancer spreads faster than non-small cell lung cancer

#Small cell lung cancer (SCLC) and Non-small cell lung cancer (NSCLC)
tempdf=df %>%
  mutate(
    cancer_type = case_when(
      cell_type == 1 ~ "nsclc",
      cell_type == 2 ~ "sclc",
      cell_type == 3 ~ "nsclc",
      cell_type == 4 ~ "nsclc"
    )
  )
df=tempdf
#treatment
tempdf=df %>%
  mutate(
    treatment_type = case_when(
      treatment == 1 ~ "standard",
      treatment == 2 ~ "test"
    )
  )
df=tempdf

df=tempdf
#prior chemotherapy
tempdf=df %>%
  mutate(
    prior_chemo = case_when(
      prior_chemo == 0 ~ "No",
      prior_chemo == 10 ~ "Yes"
    )
  )

df=tempdf

as.factor(df$cancer_type)
as.factor(df$treatment)
as.factor(df$status)
as.factor(df$prior_chemo)


ggplot(df, aes(x=treatment_type, y=survival)) +
  geom_boxplot() +
  labs(title="Survival by Treatment",x="Treatment", y = "Survival (in days)") +
  geom_boxplot(fill= 'darkgreen', outlier.shape=15,
               outlier.size=2)
#treatment data
df_treat <- df %>%
  select(survival,
         treatment_type) %>%
  where(treatment == 2)

#create df for median tests
df_test<- filter(df_treat, treatment_type == "test")
df_standard<- filter(df_treat, treatment_type == "standard")

#test vs model
median(df_treat$survival)
median(df_standard$survival)

#parametric models

#does treatment and age have effect?
exp <-  survreg(Surv(survival,status) ~ treatment + age + cancer_type + prior_chemo + diagnosis_period +treatment*age
                , data = df, dist ="exponential")
summary(exp)
#no treatment and age are not statistically significantly 
#it might be good to keep in the model to show the drs

#does the company's treatment perform better on a centain type of cancer?
exp <-  survreg(Surv(survival,status) ~ treatment + age + cancer_type + prior_chemo + diagnosis_period
                + treatment*age
                , data = df, dist ="exponential")
summary(exp)

#the treatment performs 66% better on NSCLC than on SCLC
weibull <-  survreg(Surv(survival,status) ~ treatment + age + cancer_type + prior_chemo + diagnosis_period
                    + treatment*age
                    ,data = df, dist ="weibull")
summary(weibull)

loglogistic <-  survreg(Surv(survival,status) ~ treatment + age, data = df, dist ="loglogistic")
summary(loglogistic)

#semi-parametric 

#cox hazard model
cox <- coxph(Surv(survival,status) ~ treatment + age + cancer_type + prior_chemo + diagnosis_period
             + treatment*age
             , data = df, method ="breslow")
summary(cox)

#compare models
library(stargazer)
stargazer(weibull, exp, cox, type="text", 
          title="A6. Lung Cancer")


