library(survival)
library(tidyverse)

# rm(list = ls())

source("~/rwd2023/TND-Cox/R/000-func-definition.R")


( array.id.index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) )

for (array.id in (3*array.id.index+1-3):(3*array.id.index)) {

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
  vac.lambda.0  <- sample( c(1/25, 1/50, 1/100, 1/200, 1/300, 1/500), size = N, replace = TRUE )
  vac.hazrd.lp  <-  with(
    dat.dummy,
    0.2*AgeGroupChildren + 0.5*AgeGroupElderly + -0.5*SexMale + -1*RiskGroupLowRisk + 2*HSB 
  )
  vac.prima.time <- vac.develop.time + gener_exp_t(vac.lambda.0, vac.hazrd.lp) %>% ceiling()
  dat.dummy$vac.prima.time <- dat$vac.prima.time <- vac.prima.time
  dat.dummy <- mutate(dat.dummy, id = 1:N, .before = everything())
  
  
  {
    #VE of the primary against the 1st infection
    VE.inf1.prima.short  <- 0.5 
    VE.inf1.prima.medium <- 0.4 
    VE.inf1.prima.long   <- 0.3
    
    VE.inf2.prima.short  <- 0.8 
    VE.inf2.prima.medium <- 0.7 
    VE.inf2.prima.long   <- 0.6 
    
    
    #Event constant parameters
    beta.AgeGroupChildren        =  0.3
    beta.AgeGroupElderly         =  0.5
    beta.SexMale                 = -0.2
    beta.RiskGroupLowRisk        = -1.0
    beta.HSB                     = -0.5
    
    ## event 1 parameters
    lambda.1                     =  1/30
    beta.inf1.prima.short  = log(1-VE.inf1.prima.short)
    beta.inf1.prima.medium = log(1-VE.inf1.prima.medium) - beta.inf1.prima.short
    beta.inf1.prima.long   = log(1-VE.inf1.prima.long) - beta.inf1.prima.short - beta.inf1.prima.medium
    
    ## event 2 parameters
    lambda.2       = 1/60
    beta.inf2.prima.short  = log(1-VE.inf2.prima.short)
    beta.inf2.prima.medium = log(1-VE.inf2.prima.medium) - beta.inf2.prima.short
    beta.inf2.prima.long   = log(1-VE.inf2.prima.long) - beta.inf2.prima.short - beta.inf2.prima.medium
    
  }
  
  for (i in 1:nrow(dat.dummy)) {
    
    AgeGroupChildren.i <- dat.dummy$AgeGroupChildren[i]
    AgeGroupElderly.i  <- dat.dummy$AgeGroupElderly[i]
    SexMale.i          <- dat.dummy$SexMale[i]
    RiskGroupLowRisk.i <- dat.dummy$RiskGroupLowRisk[i]
    HSB.i              <- dat.dummy$HSB[i]
    d2.short.i         <- dat.dummy$vac.prima.time[i]
    d2.medium.i        <- dat.dummy$vac.prima.time[i]+90
    d2.long.i          <- dat.dummy$vac.prima.time[i]+180
    

    betax.i            <- 
      beta.AgeGroupChildren*AgeGroupChildren.i + 
      beta.AgeGroupElderly*AgeGroupElderly.i + 
      beta.SexMale*SexMale.i + 
      beta.RiskGroupLowRisk*RiskGroupLowRisk.i + 
      beta.HSB*HSB.i
    
    dat[i, c("et1", "et2", "gene.case")] <- dat.dummy[i, c("et1", "et2", "gene.case")] <- 
      gener_step_waned_VE_rec_tdc_data(lambda.1     = lambda.1, 
                                       lambda.2     = lambda.2,
                                       betax        = betax.i, 
                                       alpha.short.inf1  = beta.inf1.prima.short,
                                       alpha.short.inf2  = beta.inf2.prima.short,
                                       alpha.medium.inf1 = beta.inf1.prima.medium,
                                       alpha.medium.inf2 = beta.inf2.prima.medium,
                                       alpha.long.inf1   = beta.inf1.prima.long,
                                       alpha.long.inf2   = beta.inf2.prima.long,
                                       d2.short     = d2.short.i, 
                                       d2.medium    = d2.medium.i, 
                                       d2.long      = d2.long.i) %>% 
      ceiling()
    
    
    
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
  
  ## Symptoms
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
  
  mean(dat.dummy$stop.time < dat.dummy$et1)
  mean(dat.dummy$stop.time < dat.dummy$et2)
  
  ## Remove Tied et1 and et2
  dat.dummy <- dat.dummy %>% filter(!et1 == et2)
  
  
  # TND Datasets
  
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
  
  
  ## TND - Count-Process Format Data for Survival Models ####
  
  dat.dummy.TND.dup.full <- dat.dummy %>% 
    filter(id %in% unique(c(dat.dummy.TND.logis.dup.Y1$id, dat.dummy.TND.logis.dup.Y2$id)))
  
  dat.dummy.TND.dup.full.cp <- tmerge(dat.dummy.TND.dup.full,    dat.dummy.TND.dup.full,    
                                      id = id, 
                                      tstop = stop.time, 
                                      status = event(et1), 
                                      status = event(et2)) 
  
  dat.dummy.TND.dup.full.cp <- tmerge(dat.dummy.TND.dup.full.cp, dat.dummy.TND.dup.full.cp, id = id,
                                      dose2.short  = tdc(vac.prima.time),
                                      dose2.medium = tdc(vac.prima.time+90),
                                      dose2.long   = tdc(vac.prima.time+180))
  
  dat.dummy.TND.dup.full.cp <- dat.dummy.TND.dup.full.cp %>% 
    group_by(id) %>% 
    mutate(order = ifelse(tstop <= et1, 1, 2)) %>% 
    select(id, vac.prima.time, et1, et2, 
           tstart, tstop, status, order,
           AgeGroupChildren, AgeGroupElderly, SexMale, RiskGroupLowRisk, 
           HSB.p, HSB, Symp.p, Symp, Test.p, Test, dose2.short, dose2.medium, dose2.long) %>% 
    ungroup() %>%
    arrange(id, order, tstart) %>%
    group_by(id, order) %>%
    mutate(
      tstart_gap = tstart - first(tstart),
      tstop_gap = tstop - first(tstart)
    )  %>%
    ungroup()
  
  
  library(coxme)
  
  fit.TND.cox.full.frailty <- coxme(Surv(tstart_gap, tstop_gap, status) ~ AgeGroupChildren + 
                                      AgeGroupElderly + SexMale + RiskGroupLowRisk+ (1|id) + 
                                      strata(order) / (dose2.short + dose2.medium + dose2.long), 
                                    dat.dummy.TND.dup.full.cp)

  
  save.fit.esti.cox.step(fit.TND.cox.full.frailty, 
  order.1.dose.2.short.path  = glue::glue("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.1.dose.2.short/fit.order.1.dose.2.short.seed.{array.id}.rds"),
  order.1.dose.2.medium.path = glue::glue("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.1.dose.2.medium/fit.order.1.dose.2.medium.seed.{array.id}.rds"),
  order.1.dose.2.long.path   = glue::glue("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.1.dose.2.long/fit.order.1.dose.2.long.seed.{array.id}.rds"),
  order.2.dose.2.short.path  = glue::glue("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.2.dose.2.short/fit.order.2.dose.2.short.seed.{array.id}.rds"),
  order.2.dose.2.medium.path = glue::glue("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.2.dose.2.medium/fit.order.2.dose.2.medium.seed.{array.id}.rds"),
  order.2.dose.2.long.path   = glue::glue("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.2.dose.2.long/fit.order.2.dose.2.long.seed.{array.id}.rds")
  )
  
}













