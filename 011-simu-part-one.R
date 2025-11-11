library(survival)
library(tidyverse)
library(coxme)

# rm(list = ls())

source("~/rwd2023/TND-Cox/R/000-func-definition.R")


( array.id.index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) )

## Parallel Computing among array.ids
for (array.id in (6*array.id.index+1-6):(6*array.id.index)) {
  
  set.seed(array.id+1234) 
  
  N         =  30000
  
  # Set up confounders 'C' in the DAGs
  AgeGroup  =  rep(c("Children", "Adult", "Elderly"), 
                   times = N*c(0.2, 0.6, 0.2)) %>% sample()
  Sex       =  rep(c("Male", "Female"), 
                   times = N*c(0.5, 0.5)) %>% sample()
  RiskGroup =  rep(c("HighRisk", "LowRisk"), 
                   times = N*c(0.1, 0.9)) %>% sample()
  
  dat <- data.frame(
    id       = 1:N,
    AgeGroup = AgeGroup, 
    Sex      = Sex, 
    RiskGroup = RiskGroup
  ) 
  
  dat.dummy <- model.matrix(
    ~ AgeGroup + Sex + RiskGroup
  ) %>% data.frame() %>% 
    select(-X.Intercept.)
  
  
  # Set up health-care seeking behaviors 'H' in the DAGs
  HSB.comp <- rnorm(N) # HSB independent components
  HSB.lp <- with(
    dat.dummy,
    2 + 0.5*AgeGroupElderly + -0.5*SexMale + -2*RiskGroupLowRisk + 0.5*HSB.comp
  )
  HSB.p     = 1/(1+exp(-HSB.lp)) 
  HSB       = rbinom(N, 1, HSB.p)
  dat.dummy$HSB.p <- dat$HSB.p <- HSB.p
  dat.dummy$HSB   <- dat$HSB   <- HSB
  
  
  vac.develop.time <- 0
  boost.gap.time   <- 20
  vac.lambda.0  <- sample( c(1/25, 1/50, 1/100, 1/200, 1/300, 1/500), size = N, replace = TRUE )
  vac.hazrd.lp  <-  with(
    dat.dummy,
    0.2*AgeGroupChildren + 0.5*AgeGroupElderly + -0.5*SexMale + -1*RiskGroupLowRisk + 2*HSB 
  )

  vac.prima.time <- vac.develop.time + gener_exp_t(vac.lambda.0, vac.hazrd.lp) %>% ceiling()
  vac.boost.time <- vac.prima.time + boost.gap.time + gener_exp_t(vac.lambda.0, vac.hazrd.lp) %>% ceiling()
  
  dat.dummy$vac.prima.time <- dat$vac.prima.time <- vac.prima.time
  dat.dummy$vac.boost.time <- dat$vac.boost.time <- vac.boost.time
  dat.dummy <- mutate(dat.dummy, id = 1:N, .before = everything())
  
  
  {
    #VE of the primary and booster against the 1st and 2nd infection
    VE.inf1.prima <- 0.5 
    VE.inf1.boost <- 0.6 
    VE.inf2.prima <- 0.7 
    VE.inf2.boost <- 0.8 
    
    #Event constant parameters
    beta.AgeGroupChildren        =  0.3
    beta.AgeGroupElderly         =  0.5
    beta.SexMale                 = -0.2
    beta.RiskGroupLowRisk        = -1.0
    beta.HSB                     = -0.5
    
    ## event 1 parameters
    lambda.1                     =  1/30
    beta.inf1.prima = log(1-VE.inf1.prima)
    beta.inf1.boost = log(1-VE.inf1.boost) - log(1-VE.inf1.prima) ## Relative HR
    
    ## event 2 parameters
    lambda.2       = 1/60
    beta.inf2.prima = log(1-VE.inf2.prima)
    beta.inf2.boost = log(1-VE.inf2.boost) - log(1-VE.inf2.prima) ## Relative HR
  }
  
  for (i in 1:nrow(dat.dummy)) {
    
    AgeGroupChildren.i <- dat.dummy$AgeGroupChildren[i]
    AgeGroupElderly.i  <- dat.dummy$AgeGroupElderly[i]
    SexMale.i          <- dat.dummy$SexMale[i]
    RiskGroupLowRisk.i <- dat.dummy$RiskGroupLowRisk[i]
    HSB.i              <- dat.dummy$HSB[i]
    d2t.i              <- dat.dummy$vac.prima.time[i]
    d3t.i              <- dat.dummy$vac.boost.time[i]
    
    betax.i            <- 
      beta.AgeGroupChildren*AgeGroupChildren.i + 
      beta.AgeGroupElderly*AgeGroupElderly.i + 
      beta.SexMale*SexMale.i + 
      beta.RiskGroupLowRisk*RiskGroupLowRisk.i + 
      beta.HSB*HSB.i
    
    dat[i, c("et1", "et2", "gene.case")] <- dat.dummy[i, c("et1", "et2", "gene.case")] <- 
      gener_rec_tdc_data(lambda.1     = lambda.1, 
                         lambda.2     = lambda.2, 
                         betax        = betax.i, 
                         alpha.2.inf1 = beta.inf1.prima, 
                         alpha.2.inf2 = beta.inf2.prima, 
                         alpha.3.inf1 = beta.inf1.boost, 
                         alpha.3.inf2 = beta.inf2.boost, 
                         d2           = d2t.i, 
                         d3           = d3t.i)
    
    dat[i, c("coet1")] <- dat.dummy[i, c("coet1")] <- 
      gener_exp_t(1/100, 0) %>% ceiling()
    
    dat[i, c("coet2")] <- dat.dummy[i, c("coet2")] <- 
      gener_exp_t(1/200, 0) %>% ceiling()
    
    dat[i, c("coet3")] <- dat.dummy[i, c("coet3")] <- 
      gener_exp_t(1/300, 0) %>% ceiling()
    
    dat[i, c("coet4")] <- dat.dummy[i, c("coet4")] <- 
      gener_exp_t(1/400, 0) %>% ceiling()
    
    dat[i, c("coet5")] <- dat.dummy[i, c("coet5")] <- 
      gener_exp_t(1/500, 0) %>% ceiling()
    
    dat[i, c("coet6")] <- dat.dummy[i, c("coet6")] <- 
      gener_exp_t(1/500, 0) %>% ceiling()
    
    dat[i, c("coet7")] <- dat.dummy[i, c("coet7")] <- 
      gener_exp_t(1/600, 0) %>% ceiling()
    
    dat[i, c("coet8")] <- dat.dummy[i, c("coet8")] <- 
      gener_exp_t(1/600, 0) %>% ceiling()
    
    dat[i, c("coet9")] <- dat.dummy[i, c("coet9")] <- 
      gener_exp_t(1/800, 0) %>% ceiling()
    
    dat[i, c("coet10")] <- dat.dummy[i, c("coet10")] <- 
      gener_exp_t(1/800, 0) %>% ceiling()
    
  }
  
  ## Set all infections are symptomatic
  Symp.p    = 1
  Symp      = rbinom(N, 1, Symp.p)
  
  dat.dummy$Symp.p  <- dat$Symp.p  <- Symp.p
  dat.dummy$Symp    <- dat$Symp    <- Symp
  
  ## Testing
  Test.lp <- with(
    dat.dummy,
    -100000 + 1*Symp + 100000*HSB
  )
  Test.p     = 1/(1+exp(-Test.lp)) 
  Test      = rbinom(N, 1, Test.p)
  
  dat.dummy$Test.p  <- dat$Test.p  <- Test.p
  dat.dummy$Test    <- dat$Test    <- Test
  
  ## Censoring
  
  max.time <- 365*2
  
  dat.dummy$censoring.time <- max.time
  dat.dummy$stop.time      <- pmin(dat.dummy$censoring.time, dat.dummy$et2)
  
  ## Remove Tied et1 and et2
  dat.dummy <- dat.dummy %>% filter(!et1 == et2)
  

  ### TND datasets  ###
  
  dat.dummy.TND.logis.dup.Y1 <- data.frame() ## store 1st infection data
  dat.dummy.TND.logis.dup.Y2 <- data.frame() ## store 2nd infection data
  
  for (i in 1:max.time ) {
    
    ## Data for 1st infection data
    
    dat.case.i.Y1          <- filter(dat.dummy,  Test==1 & et1 == i)
    dat.control.i.Y1       <- filter(dat.dummy,  Test==1 & et1 > i & (coet1 == i | coet2 == i | coet3 == i | coet4 == i | coet5 == i | coet6 == i | coet7 == i | coet8 == i | coet9 == i | coet10 == i ))
    
    if ( nrow(dat.case.i.Y1)>0 & nrow(dat.control.i.Y1)>0 ) {
      
      dat.case.i.Y1$ctype     <- glue::glue("Y1 case at {i}")
      dat.control.i.Y1$ctype  <- glue::glue("Y1 control at {i}")
      dat.case.i.Y1$test.time <- dat.control.i.Y1$test.time <- i
      
      print(glue::glue("{i}/{max.time}; case:control {nrow(dat.case.i.Y1)}:{nrow(dat.control.i.Y1)}"))
      
      
      dat.dummy.TND.logis.dup.Y1 <- 
        rbind(dat.dummy.TND.logis.dup.Y1,
              dat.case.i.Y1,
              dat.control.i.Y1)
    }
    
    ## Data for 2nd infection data
    
    dat.case.i.Y2          <- filter(dat.dummy,  Test==1 & et1 < i & et2 == i)
    dat.control.i.Y2       <- filter(dat.dummy,  Test==1 & et1 < i & et2 > i & (coet1 == i | coet2 == i | coet3 == i | coet4 == i | coet5 == i | coet6 == i | coet7 == i | coet8 == i | coet9 == i | coet10 == i ))
    
    if ( nrow(dat.case.i.Y2)>0 & nrow(dat.control.i.Y2)>0 ) {
      
      dat.case.i.Y2$ctype     <- glue::glue("Y2 case at {i}")
      dat.control.i.Y2$ctype  <- glue::glue("Y2 control at {i}")
      dat.case.i.Y2$test.time <- dat.control.i.Y2$test.time <- i
      
      print(glue::glue("{i}/{max.time}; case:control {nrow(dat.case.i.Y2)}:{nrow(dat.control.i.Y2)}"))
      
      
      dat.dummy.TND.logis.dup.Y2 <- 
        rbind(dat.dummy.TND.logis.dup.Y2,
              dat.case.i.Y2,
              dat.control.i.Y2)
    }
  }
  
  dat.dummy.TND.logis.dup.Y1$Y1            <- as.numeric(dat.dummy.TND.logis.dup.Y1$test.time >= dat.dummy.TND.logis.dup.Y1$et1)
  dat.dummy.TND.logis.dup.Y1$Y2            <- as.numeric(dat.dummy.TND.logis.dup.Y1$test.time >= dat.dummy.TND.logis.dup.Y1$et2)
  dat.dummy.TND.logis.dup.Y1$dose2         <- as.numeric(dat.dummy.TND.logis.dup.Y1$test.time >= dat.dummy.TND.logis.dup.Y1$vac.prima.time)
  dat.dummy.TND.logis.dup.Y1$dose3         <- as.numeric(dat.dummy.TND.logis.dup.Y1$test.time >= dat.dummy.TND.logis.dup.Y1$vac.boost.time)
  dat.dummy.TND.logis.dup.Y1.sub           <- dat.dummy.TND.logis.dup.Y1
  
  dat.dummy.TND.logis.dup.Y2$Y1            <- as.numeric(dat.dummy.TND.logis.dup.Y2$test.time >= dat.dummy.TND.logis.dup.Y2$et1)
  dat.dummy.TND.logis.dup.Y2$Y2            <- as.numeric(dat.dummy.TND.logis.dup.Y2$test.time >= dat.dummy.TND.logis.dup.Y2$et2)
  dat.dummy.TND.logis.dup.Y2$dose2         <- as.numeric(dat.dummy.TND.logis.dup.Y2$test.time >= dat.dummy.TND.logis.dup.Y2$vac.prima.time)
  dat.dummy.TND.logis.dup.Y2$dose3         <- as.numeric(dat.dummy.TND.logis.dup.Y2$test.time >= dat.dummy.TND.logis.dup.Y2$vac.boost.time)
  dat.dummy.TND.logis.dup.Y2.sub           <- dat.dummy.TND.logis.dup.Y2
  
  
  #### Cohort Study - Count-Process Format Data for Survival Models ####
  
  dat.dummy.cp <- tmerge(dat.dummy,    dat.dummy,    
                         id = id, 
                         tstop = stop.time, 
                         status = event(et1), 
                         status = event(et2)) 
  
  dat.dummy.cp <- tmerge(dat.dummy.cp, dat.dummy.cp, id = id,
                         dose2 = tdc(vac.prima.time), dose3 = tdc(vac.boost.time))
  
  dat.dummy.cp <- dat.dummy.cp %>% 
    group_by(id) %>% 
    mutate(order = ifelse(tstop <= et1, 1, 2)) %>% 
    select(id, vac.prima.time, vac.boost.time, et1, et2, 
           tstart, tstop, status, order,
           AgeGroupChildren, AgeGroupElderly, SexMale, RiskGroupLowRisk, 
           HSB.p, HSB, Symp.p, Symp, Test.p, Test, dose2, dose3) %>% 
    ungroup() %>%
    arrange(id, order, tstart) %>%
    group_by(id, order) %>%
    mutate(
      tstart_gap = tstart - first(tstart),
      tstop_gap = tstop - first(tstart)
    )  %>%
    ungroup()
  
  
  #### Model Fitting ####
  
  ## Cohort Study - PGR model 
  
  fit.cohort.cox.robust <- coxph(
    Surv(tstart_gap, tstop_gap, status) ~ AgeGroupChildren + 
      AgeGroupElderly + SexMale + RiskGroupLowRisk + cluster(id) + 
      strata(order) / (dose2 + dose3), robust = TRUE,
    data = dat.dummy.cp,
    control = coxph.control(timefix = FALSE,
                            iter.max = 1000, 
                            eps = 1e-15, 
                            toler.chol = 2e-16)
  )
  

  save.fit.esti.cox(fit.cohort.cox.robust, 
                    order.1.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.1.dose.2/fit.cohort.cox.robust.order.1.dose.2.seed.{array.id}.rds"), 
                    order.1.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.1.dose.3/fit.cohort.cox.robust.order.1.dose.3.seed.{array.id}.rds"), 
                    order.2.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.2.dose.2/fit.cohort.cox.robust.order.2.dose.2.seed.{array.id}.rds"), 
                    order.2.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-robust/order.2.dose.3/fit.cohort.cox.robust.order.2.dose.3.seed.{array.id}.rds"))
  
  ## Cohort Study - PGF model 
  
  fit.cohort.cox.frailty <- coxme(Surv(tstart_gap, tstop_gap, status) ~ AgeGroupChildren + 
                                    AgeGroupElderly + SexMale + RiskGroupLowRisk+ (1|id) + 
                                    strata(order) / (dose2 + dose3), 
                                  dat.dummy.cp)
  
  save.fit.esti.cox(fit.cohort.cox.frailty, 
                    order.1.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.1.dose.2/fit.cohort.cox.frailty.order.1.dose.2.seed.{array.id}.rds"), 
                    order.1.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.1.dose.3/fit.cohort.cox.frailty.order.1.dose.3.seed.{array.id}.rds"), 
                    order.2.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.2.dose.2/fit.cohort.cox.frailty.order.2.dose.2.seed.{array.id}.rds"), 
                    order.2.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/cohort-cox-frailty/order.2.dose.3/fit.cohort.cox.frailty.order.2.dose.3.seed.{array.id}.rds"))
  

  ## TND - UnConditional Logistic Regression w/ calendar time
  
  fit.lgr.dup.has.time.Y1 <- glm(Y1 ~ AgeGroupChildren + AgeGroupElderly + SexMale + RiskGroupLowRisk + 
                                   dose2 + dose3 + test.time, dat.dummy.TND.logis.dup.Y1.sub, 
                                 family = binomial(), model = FALSE) ## model = FALSE is used to save storage.
  
  save.fit.esti.lgr(fit.lgr.dup.has.time.Y1, 
                    dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.1.dose.2/fit.lgr.dup.has.time.order.1.dose.2.seed.{array.id}.rds"), 
                    dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.1.dose.3/fit.lgr.dup.has.time.order.1.dose.3.seed.{array.id}.rds"))
  
  
  fit.lgr.dup.has.time.Y2 <- glm(Y2 ~ AgeGroupChildren + AgeGroupElderly + SexMale + RiskGroupLowRisk + 
                                   dose2 + dose3 + test.time, dat.dummy.TND.logis.dup.Y2.sub, 
                                 family = binomial(), model = FALSE)
  
  save.fit.esti.lgr(fit.lgr.dup.has.time.Y2, 
                    dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.2.dose.2/fit.lgr.dup.has.time.order.2.dose.2.seed.{array.id}.rds"), 
                    dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-has-time/order.2.dose.3/fit.lgr.dup.has.time.order.2.dose.3.seed.{array.id}.rds"))
  
  
  ## TND - UnConditional Logistic Regression w/o calendar time
  
  fit.lgr.dup.no.time.Y1 <- glm(Y1 ~ AgeGroupChildren + AgeGroupElderly + SexMale + RiskGroupLowRisk + 
                                   dose2 + dose3, dat.dummy.TND.logis.dup.Y1.sub, 
                                 family = binomial(), model = FALSE)
  
  save.fit.esti.lgr(fit.lgr.dup.no.time.Y1, 
                    dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.1.dose.2/fit.lgr.dup.no.time.order.1.dose.2.seed.{array.id}.rds"), 
                    dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.1.dose.3/fit.lgr.dup.no.time.order.1.dose.3.seed.{array.id}.rds"))
  
  
  fit.lgr.dup.no.time.Y2 <- glm(Y2 ~ AgeGroupChildren + AgeGroupElderly + SexMale + RiskGroupLowRisk + 
                                   dose2 + dose3, dat.dummy.TND.logis.dup.Y2.sub, 
                                 family = binomial(), model = FALSE)
  
  save.fit.esti.lgr(fit.lgr.dup.no.time.Y2, 
                    dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.2.dose.2/fit.lgr.dup.no.time.order.2.dose.2.seed.{array.id}.rds"), 
                    dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-uncond-logis-no-time/order.2.dose.3/fit.lgr.dup.no.time.order.2.dose.3.seed.{array.id}.rds"))
  
  
  ## TND - Conditional Logistic Regression
  
  fit.clgr.dup.Y1 <- clogit(Y1 ~ AgeGroupChildren + AgeGroupElderly + SexMale + RiskGroupLowRisk + 
                              dose2 + dose3 + strata(test.time), method = "exact", 
                            dat.dummy.TND.logis.dup.Y1.sub, model = FALSE)
  
  save.fit.esti.lgr(fit.clgr.dup.Y1, 
                    dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.1.dose.2/fit.clgr.dup.order.1.dose.2.seed.{array.id}.rds"), 
                    dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.1.dose.3/fit.clgr.dup.order.1.dose.3.seed.{array.id}.rds"))
  
  
  fit.clgr.dup.Y2 <- clogit(Y2 ~ AgeGroupChildren + AgeGroupElderly + SexMale + RiskGroupLowRisk + 
                              dose2 + dose3 + strata(test.time), method = "exact", 
                            dat.dummy.TND.logis.dup.Y2.sub, model = FALSE)
  
  save.fit.esti.lgr(fit.clgr.dup.Y2, 
                    dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.2.dose.2/fit.clgr.dup.order.2.dose.2.seed.{array.id}.rds"), 
                    dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cond-logis/order.2.dose.3/fit.clgr.dup.order.2.dose.3.seed.{array.id}.rds"))
  
  
  #### TND - Count-Process Format Data for Survival Models ####
  
  dat.dummy.TND.dup.full <- dat.dummy %>% 
    filter(id %in% unique(c(dat.dummy.TND.logis.dup.Y1$id, dat.dummy.TND.logis.dup.Y2$id)))
  
  dat.dummy.TND.dup.full.cp <- tmerge(dat.dummy.TND.dup.full,    dat.dummy.TND.dup.full,    
                                      id = id, 
                                      tstop = stop.time, 
                                      status = event(et1), 
                                      status = event(et2)) 
  
  dat.dummy.TND.dup.full.cp <- tmerge(dat.dummy.TND.dup.full.cp, dat.dummy.TND.dup.full.cp, id = id,
                                      dose2 = tdc(vac.prima.time), dose3 = tdc(vac.boost.time))
  
  dat.dummy.TND.dup.full.cp <- dat.dummy.TND.dup.full.cp %>% 
    group_by(id) %>% 
    mutate(order = ifelse(tstop <= et1, 1, 2)) %>% 
    select(id, vac.prima.time, vac.boost.time, et1, et2, 
           tstart, tstop, status, order,
           AgeGroupChildren, AgeGroupElderly, SexMale, RiskGroupLowRisk, 
           HSB.p, HSB, Symp.p, Symp, Test.p, Test, dose2, dose3) %>% 
    ungroup() %>%
    arrange(id, order, tstart) %>%
    group_by(id, order) %>%
    mutate(
      tstart_gap = tstart - first(tstart),
      tstop_gap = tstop - first(tstart)
    )  %>%
    ungroup()
  
  #### Model Fitting ####
  
  ## TND - PGR model 
  
  fit.TND.cox.full.robust <- coxph(
    Surv(tstart_gap, tstop_gap, status) ~ AgeGroupChildren + 
      AgeGroupElderly + SexMale + RiskGroupLowRisk + cluster(id) + 
      strata(order) / (dose2 + dose3), robust = T,
    data = dat.dummy.TND.dup.full.cp,
    control = coxph.control(timefix = FALSE,
                            iter.max = 1000, 
                            eps = 1e-15, 
                            toler.chol = 2e-16)
  )

  save.fit.esti.cox(fit.TND.cox.full.robust, 
                    order.1.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.1.dose.2/fit.TND.cox.full.robust.order.1.dose.2.seed.{array.id}.rds"), 
                    order.1.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.1.dose.3/fit.TND.cox.full.robust.order.1.dose.3.seed.{array.id}.rds"), 
                    order.2.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.2.dose.2/fit.TND.cox.full.robust.order.2.dose.2.seed.{array.id}.rds"), 
                    order.2.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-robust/order.2.dose.3/fit.TND.cox.full.robust.order.2.dose.3.seed.{array.id}.rds"))
  
  ## TND - PGF model 
  
  fit.TND.cox.full.frailty <- coxme(Surv(tstart_gap, tstop_gap, status) ~ AgeGroupChildren + 
                                      AgeGroupElderly + SexMale + RiskGroupLowRisk+ (1|id) + 
                                      strata(order) / (dose2 + dose3), 
                                    dat.dummy.TND.dup.full.cp)
  
  save.fit.esti.cox(fit.TND.cox.full.frailty, 
                    order.1.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.1.dose.2/fit.TND.cox.full.frailty.order.1.dose.2.seed.{array.id}.rds"), 
                    order.1.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.1.dose.3/fit.TND.cox.full.frailty.order.1.dose.3.seed.{array.id}.rds"), 
                    order.2.dose.2.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.2.dose.2/fit.TND.cox.full.frailty.order.2.dose.2.seed.{array.id}.rds"), 
                    order.2.dose.3.path = glue::glue("~/rwd2023/TND-Cox/rds/comp-TND-cohort/tnd-cox-full-frailty/order.2.dose.3/fit.TND.cox.full.frailty.order.2.dose.3.seed.{array.id}.rds"))
  
  
}





