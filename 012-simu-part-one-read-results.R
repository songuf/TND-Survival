
{
library(tidyverse)
library(survival)

rm(list = ls())

source("~/rwd2023/TND-Cox/R/000-func-definition.R")

## Read rds files
{
  list.fit.cohort.cox.frailty.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.1.dose.2/")
  df.fit.cohort.cox.frailty.order.1.dose.2 <- create_VE_sd_df(list.fit.cohort.cox.frailty.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.cohort.cox.frailty.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.1.dose.3/")
  df.fit.cohort.cox.frailty.order.1.dose.3 <- create_VE_sd_df(list.fit.cohort.cox.frailty.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.cohort.cox.frailty.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.2.dose.2/")
  df.fit.cohort.cox.frailty.order.2.dose.2 <- create_VE_sd_df(list.fit.cohort.cox.frailty.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.cohort.cox.frailty.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.2.dose.3/")
  df.fit.cohort.cox.frailty.order.2.dose.3 <- create_VE_sd_df(list.fit.cohort.cox.frailty.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  ### 
  
  list.fit.cohort.cox.robust.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.1.dose.2/")
  df.fit.cohort.cox.robust.order.1.dose.2 <- create_VE_sd_df(list.fit.cohort.cox.robust.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.cohort.cox.robust.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.1.dose.3/")
  df.fit.cohort.cox.robust.order.1.dose.3 <- create_VE_sd_df(list.fit.cohort.cox.robust.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.cohort.cox.robust.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.2.dose.2/")
  df.fit.cohort.cox.robust.order.2.dose.2 <- create_VE_sd_df(list.fit.cohort.cox.robust.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.cohort.cox.robust.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.2.dose.3/")
  df.fit.cohort.cox.robust.order.2.dose.3 <- create_VE_sd_df(list.fit.cohort.cox.robust.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  ### 
  
  list.fit.tnd.cond.logis.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.1.dose.2/")
  df.fit.tnd.cond.logis.order.1.dose.2 <- create_VE_sd_df(list.fit.tnd.cond.logis.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cond.logis.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.1.dose.3/")
  df.fit.tnd.cond.logis.order.1.dose.3 <- create_VE_sd_df(list.fit.tnd.cond.logis.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cond.logis.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.2.dose.2/")
  df.fit.tnd.cond.logis.order.2.dose.2 <- create_VE_sd_df(list.fit.tnd.cond.logis.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cond.logis.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.2.dose.3/")
  df.fit.tnd.cond.logis.order.2.dose.3 <- create_VE_sd_df(list.fit.tnd.cond.logis.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  ### 
  
  list.fit.tnd.cox.full.frailty.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.1.dose.2/")
  df.fit.tnd.cox.full.frailty.order.1.dose.2 <- create_VE_sd_df(list.fit.tnd.cox.full.frailty.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cox.full.frailty.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.1.dose.3/")
  df.fit.tnd.cox.full.frailty.order.1.dose.3 <- create_VE_sd_df(list.fit.tnd.cox.full.frailty.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cox.full.frailty.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.2.dose.2/")
  df.fit.tnd.cox.full.frailty.order.2.dose.2 <- create_VE_sd_df(list.fit.tnd.cox.full.frailty.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cox.full.frailty.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.2.dose.3/")
  df.fit.tnd.cox.full.frailty.order.2.dose.3 <- create_VE_sd_df(list.fit.tnd.cox.full.frailty.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])

  ### 
  
  list.fit.tnd.cox.full.robust.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.1.dose.2/")
  df.fit.tnd.cox.full.robust.order.1.dose.2 <- create_VE_sd_df(list.fit.tnd.cox.full.robust.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cox.full.robust.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.1.dose.3/")
  df.fit.tnd.cox.full.robust.order.1.dose.3 <- create_VE_sd_df(list.fit.tnd.cox.full.robust.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cox.full.robust.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.2.dose.2/")
  df.fit.tnd.cox.full.robust.order.2.dose.2 <- create_VE_sd_df(list.fit.tnd.cox.full.robust.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.cox.full.robust.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.2.dose.3/")
  df.fit.tnd.cox.full.robust.order.2.dose.3 <- create_VE_sd_df(list.fit.tnd.cox.full.robust.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  ### 
  
  list.fit.tnd.uncond.logis.has.time.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.1.dose.2/")
  df.fit.tnd.uncond.logis.has.time.order.1.dose.2 <- create_VE_sd_df(list.fit.tnd.uncond.logis.has.time.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.uncond.logis.has.time.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.1.dose.3/")
  df.fit.tnd.uncond.logis.has.time.order.1.dose.3 <- create_VE_sd_df(list.fit.tnd.uncond.logis.has.time.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.uncond.logis.has.time.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.2.dose.2/")
  df.fit.tnd.uncond.logis.has.time.order.2.dose.2 <- create_VE_sd_df(list.fit.tnd.uncond.logis.has.time.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.uncond.logis.has.time.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.2.dose.3/")
  df.fit.tnd.uncond.logis.has.time.order.2.dose.3 <- create_VE_sd_df(list.fit.tnd.uncond.logis.has.time.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  
  ### 
  
  list.fit.tnd.uncond.logis.no.time.order.1.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.1.dose.2/")
  df.fit.tnd.uncond.logis.no.time.order.1.dose.2 <- create_VE_sd_df(list.fit.tnd.uncond.logis.no.time.order.1.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.uncond.logis.no.time.order.1.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.1.dose.3/")
  df.fit.tnd.uncond.logis.no.time.order.1.dose.3 <- create_VE_sd_df(list.fit.tnd.uncond.logis.no.time.order.1.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.uncond.logis.no.time.order.2.dose.2 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.2.dose.2/")
  df.fit.tnd.uncond.logis.no.time.order.2.dose.2 <- create_VE_sd_df(list.fit.tnd.uncond.logis.no.time.order.2.dose.2)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.tnd.uncond.logis.no.time.order.2.dose.3 <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.2.dose.3/")
  df.fit.tnd.uncond.logis.no.time.order.2.dose.3 <- create_VE_sd_df(list.fit.tnd.uncond.logis.no.time.order.2.dose.3)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
}





df.all <- list(df.fit.cohort.cox.frailty.order.1.dose.2, 
                df.fit.cohort.cox.frailty.order.1.dose.3, 
                df.fit.cohort.cox.frailty.order.2.dose.2,    
                df.fit.cohort.cox.frailty.order.2.dose.3,
                
                df.fit.cohort.cox.robust.order.1.dose.2,
                df.fit.cohort.cox.robust.order.1.dose.3,
                df.fit.cohort.cox.robust.order.2.dose.2,
                df.fit.cohort.cox.robust.order.2.dose.3,
                
                df.fit.tnd.cond.logis.order.1.dose.2,
                df.fit.tnd.cond.logis.order.1.dose.3,
                df.fit.tnd.cond.logis.order.2.dose.2,
                df.fit.tnd.cond.logis.order.2.dose.3,
                
                df.fit.tnd.cox.full.frailty.order.1.dose.2,
                df.fit.tnd.cox.full.frailty.order.1.dose.3,
                df.fit.tnd.cox.full.frailty.order.2.dose.2,
                df.fit.tnd.cox.full.frailty.order.2.dose.3,
                
                df.fit.tnd.cox.full.robust.order.1.dose.2,
                df.fit.tnd.cox.full.robust.order.1.dose.3,
                df.fit.tnd.cox.full.robust.order.2.dose.2,
                df.fit.tnd.cox.full.robust.order.2.dose.3,
                
                df.fit.tnd.uncond.logis.has.time.order.1.dose.2,
                df.fit.tnd.uncond.logis.has.time.order.1.dose.3,
                df.fit.tnd.uncond.logis.has.time.order.2.dose.2,
                df.fit.tnd.uncond.logis.has.time.order.2.dose.3,
                
                df.fit.tnd.uncond.logis.no.time.order.1.dose.2,
                df.fit.tnd.uncond.logis.no.time.order.1.dose.3,
                df.fit.tnd.uncond.logis.no.time.order.2.dose.2,
                df.fit.tnd.uncond.logis.no.time.order.2.dose.3)

df.all.name <- Hmisc::Cs(cohort.cox.frailty.order.1.dose.2, 
                             cohort.cox.frailty.order.1.dose.3, 
                             cohort.cox.frailty.order.2.dose.2,    
                             cohort.cox.frailty.order.2.dose.3,
                             
                             cohort.cox.robust.order.1.dose.2,
                             cohort.cox.robust.order.1.dose.3,
                             cohort.cox.robust.order.2.dose.2,
                             cohort.cox.robust.order.2.dose.3,
                             
                             tnd.cond.logis.order.1.dose.2,
                             tnd.cond.logis.order.1.dose.3,
                             tnd.cond.logis.order.2.dose.2,
                             tnd.cond.logis.order.2.dose.3,
                             
                             tnd.cox.full.frailty.order.1.dose.2,
                             tnd.cox.full.frailty.order.1.dose.3,
                             tnd.cox.full.frailty.order.2.dose.2,
                             tnd.cox.full.frailty.order.2.dose.3,

                             tnd.cox.full.robust.order.1.dose.2,
                             tnd.cox.full.robust.order.1.dose.3,
                             tnd.cox.full.robust.order.2.dose.2,
                             tnd.cox.full.robust.order.2.dose.3,
                             
                             tnd.uncond.logis.has.time.order.1.dose.2,
                             tnd.uncond.logis.has.time.order.1.dose.3,
                             tnd.uncond.logis.has.time.order.2.dose.2,
                             tnd.uncond.logis.has.time.order.2.dose.3,
                             
                             tnd.uncond.logis.no.time.order.1.dose.2,
                             tnd.uncond.logis.no.time.order.1.dose.3,
                             tnd.uncond.logis.no.time.order.2.dose.2,
                             tnd.uncond.logis.no.time.order.2.dose.3)

}


model.order <- c('PGF-TND', 
                 'PGR-TND', 
                 'CLR', 
                 'ULR-T', 
                 'ULR', 
                 'PGF-Cohort', 
                 'PGR-Cohort')

non.summ.df.all <- data.frame()

for (i in 1:length(df.all)) {
  
  df.all.i <- df.all[[i]]
  df.all.i$model <- str_extract(df.all.name[i], ".*(?=\\.order)") 
  df.all.i$order <- str_extract(df.all.name[i], "order\\.\\d+") 
  df.all.i$dose <- str_extract(df.all.name[i], "dose\\.\\d+") 
  df.all.i$order.dose <- str_extract(df.all.name[i], "order\\.\\d+.dose\\.\\d+") 
  
  df.all.i$study <- unlist(strsplit(df.all.i$model, "\\."))[1]
  df.all.i$study <- ifelse(df.all.i$study == "tnd", "TND", df.all.i$study)
  df.all.i$study <- ifelse(df.all.i$study == "cohort", "Cohort", df.all.i$study)
  
  non.summ.df.all <- rbind(non.summ.df.all, df.all.i)
  
}

non.summ.df.all <- non.summ.df.all %>% 
  mutate(
    model = case_when(
      model == "cohort.cox.frailty"         ~ "PGF-Cohort",
      model == "cohort.cox.robust"          ~ "PGR-Cohort",
      model == "tnd.cond.logis"             ~ "CLR",
      model == "tnd.cox.full.frailty"       ~ "PGF-TND",
      model == "tnd.cox.full.robust"        ~ "PGR-TND",
      model == "tnd.uncond.logis.has.time"  ~ "ULR-T",
      model == "tnd.uncond.logis.no.time"   ~ "ULR",
      TRUE ~ NA_character_
    ) 
  ) %>% 
  mutate(
    order = case_when(
      order == "order.1"          ~ "1st Inf",
      order == "order.2"          ~ "2nd Inf",
      TRUE ~ NA_character_
    ),
    dose = case_when(
      dose == "dose.2"          ~ "2 Doses",
      dose == "dose.3"          ~ "3 Doses",
      TRUE ~ NA_character_
    )
  ) 


# Calculate mean and quantiles
quantile.df.all <- non.summ.df.all %>%
  group_by(model, order, dose, order.dose, study) %>%
  summarise(mean_VE = mean(VE),
            ci_upp = quantile(VE, 0.975),
            ci_low = quantile(VE, 0.025)) 

## summ.df.all
nloop       <- nrow(df.all[[1]])
nmethod     <- length(df.all)

summ.df.all <- matrix(NA, nrow = nmethod, ncol = 6, 
                      dimnames = list(df.all.name,
                                      c("VE.avg", "VE.bias.emp", "VE.sd.avg", "VE.sd.emp", "CI.cov.95", "MSE")))


for (i in 1:length(df.all)) {
  imethod <- df.all.name[i]
  iod     <- substr(imethod, nchar(imethod) - 13, nchar(imethod))
  
  VE.true      <- case_when(
    iod == "order.1.dose.2"    ~ 0.5,
    iod == "order.1.dose.3"    ~ 0.6,
    iod == "order.2.dose.2"    ~ 0.7,
    iod == "order.2.dose.3"    ~ 0.8,
    TRUE ~ NA_real_
  )
  
  VE.list    <- df.all[[i]]$VE
  VE.sd.list <- df.all[[i]]$VE.sd
  
  summ.df.all[i, 1] <- mean(VE.list)
  summ.df.all[i, 2] <- (mean(VE.list) - VE.true)
  summ.df.all[i, 3] <- mean(VE.sd.list)
  summ.df.all[i, 4] <- sd(VE.list)  
  
  VE.true
  
  VE.ci.low <- VE.list - 1.96*VE.sd.list
  VE.ci.upp <- VE.list + 1.96*VE.sd.list
  
  VE.cover   <- rep(NA, len = nloop)
  for (jloop in 1:nloop) {
    VE.cover[jloop] <- 
      VE.true >= VE.ci.low[jloop] & VE.true <= VE.ci.upp[jloop]
  }
  
  summ.df.all[i, 5] <- mean(VE.cover) 
  
  summ.df.all[i, 6] <- mean((VE.list-VE.true)^2)
  
}

summ.df.all <- summ.df.all %>% 
  as.data.frame() %>% 
  rownames_to_column("rowname") %>% 
  mutate(
    model = str_extract(rowname, ".*(?=\\.order)"), 
    order = str_extract(rowname, "order\\.\\d+"), 
    dose  = str_extract(rowname, "dose\\.\\d+"), 
    order.dose = str_extract(rowname, "order\\.\\d+.dose\\.\\d+"),
    study = str_extract(rowname, "\\b\\w+\\b(?=\\.)"),
    study = ifelse(study == "tnd",    "TND", study),
    study = ifelse(study == "cohort", "Cohort", study)
  ) %>% 
  filter(!model %in% c("tnd.cox.part.frailty", "tnd.cox.part.nonrobust", 
                       "tnd.cox.part.robust", "tnd.cox.full.nonrobust",
                       "cohort.cox.nonrobust")) %>% 
  mutate(
    model = case_when(
      model == "cohort.cox.frailty"         ~ "PGF-Cohort",
      model == "cohort.cox.robust"          ~ "PGR-Cohort",
      model == "tnd.cond.logis"             ~ "CLR",
      model == "tnd.cox.full.frailty"       ~ "PGF-TND",
      model == "tnd.cox.full.robust"        ~ "PGR-TND",
      model == "tnd.uncond.logis.has.time"  ~ "ULR-T",
      model == "tnd.uncond.logis.no.time"   ~ "ULR",
      TRUE ~ NA_character_
    ) 
  ) %>% 
  mutate(
    order = case_when(
      order == "order.1"          ~ "1st Inf",
      order == "order.2"          ~ "2nd Inf",
      TRUE ~ NA_character_
    ),
    dose = case_when(
      dose == "dose.2"          ~ "2 Doses",
      dose == "dose.3"          ~ "3 Doses",
      TRUE ~ NA_character_
    )
  ) %>% 
  arrange(factor(model, levels = model.order)) %>% 
  select(-rowname)
  

summ.df.all.order.1.dose.2 <- summ.df.all %>% 
  filter(
    order == "1st Inf" & dose == "2 Doses"
  ) %>% 
  select(
    c("model", "VE.avg", "VE.bias.emp", "VE.sd.avg", "VE.sd.emp", 
      "CI.cov.95", "MSE")
  ) %>% 
  mutate(
    VE.avg      = paste0(sprintf("%.2f", round(VE.avg, 4)*100), "%"),
    VE.bias.emp = paste0(sprintf("%.2f", round(VE.bias.emp, 4)*100), "%"),
    VE.sd.avg   = paste0(sprintf("%.2f", round(VE.sd.avg, 4)*100), "%"),
    VE.sd.emp   = paste0(sprintf("%.2f", round(VE.sd.emp, 4)*100), "%"),
    CI.cov.95   = paste0(sprintf("%.2f", round(CI.cov.95, 4)*100), "%"),
    MSE         = paste0(sprintf("%.2f", round(MSE, 4)*100), "%")
  ) 

summ.df.all.order.1.dose.3 <- summ.df.all %>% 
  filter(
    order == "1st Inf" & dose == "3 Doses"
  ) %>% 
  select(
    c("model", "VE.avg", "VE.bias.emp", "VE.sd.avg", "VE.sd.emp", 
      "CI.cov.95", "MSE")
  ) %>% 
  mutate(
    VE.avg      = paste0(sprintf("%.2f", round(VE.avg, 4)*100), "%"),
    VE.bias.emp = paste0(sprintf("%.2f", round(VE.bias.emp, 4)*100), "%"),
    VE.sd.avg   = paste0(sprintf("%.2f", round(VE.sd.avg, 4)*100), "%"),
    VE.sd.emp   = paste0(sprintf("%.2f", round(VE.sd.emp, 4)*100), "%"),
    CI.cov.95   = paste0(sprintf("%.2f", round(CI.cov.95, 4)*100), "%"),
    MSE         = paste0(sprintf("%.2f", round(MSE, 4)*100), "%")
  ) 


summ.df.all.order.2.dose.2 <- summ.df.all %>% 
  filter(
    order == "2nd Inf" & dose == "2 Doses"
  ) %>% 
  select(
    c("model", "VE.avg", "VE.bias.emp", "VE.sd.avg", "VE.sd.emp", 
      "CI.cov.95", "MSE")
  ) %>% 
  mutate(
    VE.avg      = paste0(sprintf("%.2f", round(VE.avg, 4)*100), "%"),
    VE.bias.emp = paste0(sprintf("%.2f", round(VE.bias.emp, 4)*100), "%"),
    VE.sd.avg   = paste0(sprintf("%.2f", round(VE.sd.avg, 4)*100), "%"),
    VE.sd.emp   = paste0(sprintf("%.2f", round(VE.sd.emp, 4)*100), "%"),
    CI.cov.95   = paste0(sprintf("%.2f", round(CI.cov.95, 4)*100), "%"),
    MSE         = paste0(sprintf("%.2f", round(MSE, 4)*100), "%")
  ) 


summ.df.all.order.2.dose.3 <- summ.df.all %>% 
  filter(
    order == "2nd Inf" & dose == "3 Doses"
  ) %>% 
  select(
    c("model", "VE.avg", "VE.bias.emp", "VE.sd.avg", "VE.sd.emp", 
      "CI.cov.95", "MSE")
  ) %>% 
  mutate(
    VE.avg      = paste0(sprintf("%.2f", round(VE.avg, 4)*100), "%"),
    VE.bias.emp = paste0(sprintf("%.2f", round(VE.bias.emp, 4)*100), "%"),
    VE.sd.avg   = paste0(sprintf("%.2f", round(VE.sd.avg, 4)*100), "%"),
    VE.sd.emp   = paste0(sprintf("%.2f", round(VE.sd.emp, 4)*100), "%"),
    CI.cov.95   = paste0(sprintf("%.2f", round(CI.cov.95, 4)*100), "%"),
    MSE         = paste0(sprintf("%.2f", round(MSE, 4)*100), "%")
  ) 


summary_stats <- left_join(quantile.df.all, summ.df.all, 
                           by = c("model", "order", "dose", "study", "order.dose"))

# Plotting
ggplot() +
  geom_errorbar(data = summary_stats,
                aes(x = factor(model, level = model.order),
                    ymin = ci_low,
                    ymax = ci_upp), #color = study),
                width = 0.2) +  # Add error bars
  geom_point(data = summary_stats, 
             aes(x = factor(model, level = model.order), 
                 y = mean_VE), #color = study), 
             shape = 8, 
             size = 2) + 
  geom_text(data = summary_stats,
            aes(x = factor(model, level = model.order),
                y = VE.avg,
                label = paste0(
                  paste(format(round(mean_VE*100,   1), nsmall = 1), "%"),
                  "\n(",
                  paste(format(round(ci_low*100,   1), nsmall = 1), "%"),
                  ", ",
                  paste(format(round(ci_upp*100,   1), nsmall = 1), "%"),
                  ")\n",
                  paste(format(round(CI.cov.95*100, 1), nsmall = 1), "%")
                )
                ),
            vjust = 1.5, hjust = 0.0, size = 3, color = "black") +  # Add labels for 95% CI coverage
  scale_x_discrete(labels= c('PGF-TND', 
                             'PGR-TND', 
                             'CL-TND', 
                             'ULA-TND', 
                             'UL-TND', 
                             'PGF-Cohort ', 
                             'PGR-Cohort') ) + 
  scale_y_continuous(labels = scales::percent) +  
  facet_wrap(~order.dose, #scales = "free", 
             labeller = as_labeller(
               c(order.1.dose.2 = "Initial Infection - Full Doses",
                 order.1.dose.3 = "Initial Infection - Booster Dose",
                 order.2.dose.2 = "Reinfection - Full Doses",
                 order.2.dose.3 = "Reinfection - Booster Dose")
             )) +
  geom_hline(data = subset(non.summ.df.all, order.dose == "order.1.dose.2"), aes(yintercept = 0.5), linetype = "dotdash", color = "black") +
  geom_hline(data = subset(non.summ.df.all, order.dose == "order.1.dose.3"), aes(yintercept = 0.6), linetype = "dashed", color = "black") +
  geom_hline(data = subset(non.summ.df.all, order.dose == "order.2.dose.2"), aes(yintercept = 0.7), linetype = "dotted", color = "black") +
  geom_hline(data = subset(non.summ.df.all, order.dose == "order.2.dose.3"), aes(yintercept = 0.8), linetype = "twodash", color = "black") +
  theme_bw() + 
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.6, vjust =0),
        axis.title.x = element_text(angle = 0, hjust = 0.5, vjust =-2)) + 
  xlab("Estimation Scenarios") + 
  ylab("VE") 


ggsave("~/rwd2023/TND-Cox/plot/TND-Simu-Results.pdf", 
       width = 30,
       height = 20,
       units = "cm")






