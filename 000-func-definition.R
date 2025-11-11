## function definitions

### Save the important info in logistic regression results 
save.fit.esti.lgr <- function(fit.obj, 
                              dose.2.name="dose2", 
                              dose.3.name="dose3", 
                              dose.2.path,
                              dose.3.path){
  
  ind.dose.2         <- which(names(coef(fit.obj)) == dose.2.name)
  ind.dose.3         <- which(names(coef(fit.obj)) == dose.3.name)
  
  beta.list  <- coef(fit.obj)
  VE.dose.2     <- 1-exp(beta.list[ind.dose.2])
  VE.dose.3     <- 1-exp(beta.list[ind.dose.2]+beta.list[ind.dose.3])
  
  SIGMA      <- vcov(fit.obj)
  
  VE.sd.dose.2  <- sqrt( diag(SIGMA)[ind.dose.2] * exp(beta.list[ind.dose.2])^2 )
  
  VE.sd.dose.3 <- 
    sqrt( c(exp(beta.list[ind.dose.2]+beta.list[ind.dose.3]), 
            exp(beta.list[ind.dose.2]+beta.list[ind.dose.3])) %*% 
            SIGMA[c(ind.dose.2,ind.dose.3), c(ind.dose.2,ind.dose.3)] %*%
            c(exp(beta.list[ind.dose.2]+beta.list[ind.dose.3]), 
              exp(beta.list[ind.dose.2]+beta.list[ind.dose.3])) ) %>% as.vector()
  
  
  o.dose.2          <- list(VE.dose.2    = VE.dose.2, 
                            VE.sd.dose.2 = VE.sd.dose.2)
  saveRDS(o.dose.2, file = dose.2.path)
  
  o.dose.3          <- list(VE.dose.3 = VE.dose.3, 
                            VE.sd.dose.3 = VE.sd.dose.3)
  
  saveRDS(o.dose.3, file = dose.3.path)
  
  
  list(o.dose.2, o.dose.3)
}

### Save the important info in Cox regression results (Non-Step)
save.fit.esti.cox <- function(fit.obj, 
                              order.1.dose.2.name="strata(order)order=1:dose2", 
                              order.1.dose.3.name="strata(order)order=1:dose3", 
                              order.2.dose.2.name="strata(order)order=2:dose2", 
                              order.2.dose.3.name="strata(order)order=2:dose3", 
                              order.1.dose.2.path,
                              order.1.dose.3.path,
                              order.2.dose.2.path,
                              order.2.dose.3.path){

  ind.order.1.dose.2         <- which(names(coef(fit.obj)) == order.1.dose.2.name)
  ind.order.1.dose.3         <- which(names(coef(fit.obj)) == order.1.dose.3.name)
  ind.order.2.dose.2         <- which(names(coef(fit.obj)) == order.2.dose.2.name)
  ind.order.2.dose.3         <- which(names(coef(fit.obj)) == order.2.dose.3.name)
  
  beta.list  <- coef(fit.obj)
  VE.order.1.dose.2     <- 1-exp(beta.list[ind.order.1.dose.2])
  VE.order.1.dose.3     <- 1-exp(beta.list[ind.order.1.dose.2]+beta.list[ind.order.1.dose.3])
  VE.order.2.dose.2     <- 1-exp(beta.list[ind.order.2.dose.2])
  VE.order.2.dose.3     <- 1-exp(beta.list[ind.order.2.dose.2]+beta.list[ind.order.2.dose.3])
  
  SIGMA      <- vcov(fit.obj)
  
  VE.sd.order.1.dose.2  <- sqrt( diag(SIGMA)[ind.order.1.dose.2] * exp(beta.list[ind.order.1.dose.2])^2 )
  
  VE.sd.order.1.dose.3 <- 
    sqrt( c(exp(beta.list[ind.order.1.dose.2]+beta.list[ind.order.1.dose.3]), 
            exp(beta.list[ind.order.1.dose.2]+beta.list[ind.order.1.dose.3])) %*% 
            SIGMA[c(ind.order.1.dose.2,ind.order.1.dose.3), c(ind.order.1.dose.2,ind.order.1.dose.3)] %*%
            c(exp(beta.list[ind.order.1.dose.2]+beta.list[ind.order.1.dose.3]), 
              exp(beta.list[ind.order.1.dose.2]+beta.list[ind.order.1.dose.3])) ) %>% as.vector()
  
  VE.sd.order.2.dose.2  <- sqrt( diag(SIGMA)[ind.order.2.dose.2] * exp(beta.list[ind.order.2.dose.2])^2 )
  
  VE.sd.order.2.dose.3 <- 
    sqrt( c(exp(beta.list[ind.order.2.dose.2]+beta.list[ind.order.2.dose.3]), 
            exp(beta.list[ind.order.2.dose.2]+beta.list[ind.order.2.dose.3])) %*% 
            SIGMA[c(ind.order.2.dose.2,ind.order.2.dose.3), c(ind.order.2.dose.2,ind.order.2.dose.3)] %*%
            c(exp(beta.list[ind.order.2.dose.2]+beta.list[ind.order.2.dose.3]), 
              exp(beta.list[ind.order.2.dose.2]+beta.list[ind.order.2.dose.3])) ) %>% as.vector()

  
  o.order.1.dose.2          <- list(VE.order.1.dose.2    = VE.order.1.dose.2, 
                                    VE.sd.order.1.dose.2 = VE.sd.order.1.dose.2)
  saveRDS(o.order.1.dose.2, file = order.1.dose.2.path)
  
  o.order.1.dose.3          <- list(VE.order.1.dose.3 = VE.order.1.dose.3, 
                                    VE.sd.order.1.dose.3 = VE.sd.order.1.dose.3)
  saveRDS(o.order.1.dose.3, file = order.1.dose.3.path)  
  
  o.order.2.dose.2          <- list(VE.order.2.dose.2 = VE.order.2.dose.2, 
                                    VE.sd.order.2.dose.2 = VE.sd.order.2.dose.2)
  saveRDS(o.order.2.dose.2, file = order.2.dose.2.path)
  
  o.order.2.dose.3          <- list(VE.order.2.dose.3 = VE.order.2.dose.3, 
                                    VE.sd.order.2.dose.3 = VE.sd.order.2.dose.3)
  saveRDS(o.order.2.dose.3, file = order.2.dose.3.path)
  
  list(o.order.1.dose.2, o.order.1.dose.3, o.order.2.dose.2, o.order.2.dose.3)
}


### Save the important info in Cox regression results (Step)

save.fit.esti.cox.step <- function(fit.obj, 
                              order.1.dose.2.short.name ="strata(order)order=1:dose2.short", 
                              order.1.dose.2.medium.name="strata(order)order=1:dose2.medium", 
                              order.1.dose.2.long.name  ="strata(order)order=1:dose2.long", 
                              order.2.dose.2.short.name ="strata(order)order=2:dose2.short", 
                              order.2.dose.2.medium.name="strata(order)order=2:dose2.medium", 
                              order.2.dose.2.long.name  ="strata(order)order=2:dose2.long", 
                              order.1.dose.2.short.path,
                              order.1.dose.2.medium.path,
                              order.1.dose.2.long.path,
                              order.2.dose.2.short.path,
                              order.2.dose.2.medium.path,
                              order.2.dose.2.long.path){
  
  ind.order.1.dose.2.short          <- which(names(coef(fit.obj)) == order.1.dose.2.short.name)
  ind.order.1.dose.2.medium         <- which(names(coef(fit.obj)) == order.1.dose.2.medium.name)
  ind.order.1.dose.2.long           <- which(names(coef(fit.obj)) == order.1.dose.2.long.name)
  ind.order.2.dose.2.short          <- which(names(coef(fit.obj)) == order.2.dose.2.short.name)
  ind.order.2.dose.2.medium         <- which(names(coef(fit.obj)) == order.2.dose.2.medium.name)
  ind.order.2.dose.2.long           <- which(names(coef(fit.obj)) == order.2.dose.2.long.name)
  
  
  beta.list  <- coef(fit.obj)
  VE.order.1.dose.2.short      <- 1-exp(beta.list[ind.order.1.dose.2.short])
  VE.order.1.dose.2.medium     <- 1-exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium])
  VE.order.1.dose.2.long       <- 1-exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long])
    
  VE.order.2.dose.2.short      <- 1-exp(beta.list[ind.order.2.dose.2.short])
  VE.order.2.dose.2.medium     <- 1-exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium])
  VE.order.2.dose.2.long       <- 1-exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long])
  
  SIGMA      <- vcov(fit.obj)
  
  VE.sd.order.1.dose.2.short  <- sqrt( diag(SIGMA)[ind.order.1.dose.2.short] * exp(beta.list[ind.order.1.dose.2.short])^2 )
  
  VE.sd.order.1.dose.2.medium <- 
    sqrt( c(exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]), 
            exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium])) %*% 
            SIGMA[c(ind.order.1.dose.2.short,ind.order.1.dose.2.medium), c(ind.order.1.dose.2.short,ind.order.1.dose.2.medium)] %*%
            c(exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]), 
              exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium])) ) %>% as.vector()
  
  VE.sd.order.1.dose.2.long <- 
    sqrt( c(exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long]), 
            exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long]),
            exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long])) %*% 
            SIGMA[c(ind.order.1.dose.2.short,ind.order.1.dose.2.medium, ind.order.1.dose.2.long), 
                  c(ind.order.1.dose.2.short,ind.order.1.dose.2.medium, ind.order.1.dose.2.long)] %*%
            c(exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long]), 
              exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long]),
              exp(beta.list[ind.order.1.dose.2.short]+beta.list[ind.order.1.dose.2.medium]+beta.list[ind.order.1.dose.2.long])) ) %>% as.vector()
  
  VE.sd.order.2.dose.2.short  <- sqrt( diag(SIGMA)[ind.order.2.dose.2.short] * exp(beta.list[ind.order.2.dose.2.short])^2 )
  
  VE.sd.order.2.dose.2.medium <- 
    sqrt( c(exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]), 
            exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium])) %*% 
            SIGMA[c(ind.order.2.dose.2.short,ind.order.2.dose.2.medium), c(ind.order.2.dose.2.short,ind.order.2.dose.2.medium)] %*%
            c(exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]), 
              exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium])) ) %>% as.vector()
  
  VE.sd.order.2.dose.2.long <- 
    sqrt( c(exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long]), 
            exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long]),
            exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long])) %*% 
            SIGMA[c(ind.order.2.dose.2.short,ind.order.2.dose.2.medium, ind.order.2.dose.2.long), 
                  c(ind.order.2.dose.2.short,ind.order.2.dose.2.medium, ind.order.2.dose.2.long)] %*%
            c(exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long]), 
              exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long]),
              exp(beta.list[ind.order.2.dose.2.short]+beta.list[ind.order.2.dose.2.medium]+beta.list[ind.order.2.dose.2.long])) ) %>% as.vector()
  
  
  
  o.order.1.dose.2.short          <- list(VE.order.1.dose.2.short    = VE.order.1.dose.2.short, 
                                          VE.sd.order.1.dose.2.short = VE.sd.order.1.dose.2.short)
  saveRDS(o.order.1.dose.2.short, file = order.1.dose.2.short.path)
  
  o.order.1.dose.2.medium          <- list(VE.order.1.dose.2.medium = VE.order.1.dose.2.medium, 
                                    VE.sd.order.1.dose.2.medium = VE.sd.order.1.dose.2.medium)
  saveRDS(o.order.1.dose.2.medium, file = order.1.dose.2.medium.path)  
  
  o.order.1.dose.2.long          <- list(VE.order.1.dose.2.long = VE.order.1.dose.2.long, 
                                           VE.sd.order.1.dose.2.long = VE.sd.order.1.dose.2.long)
  saveRDS(o.order.1.dose.2.long, file = order.1.dose.2.long.path)  
  
  o.order.2.dose.2.short          <- list(VE.order.2.dose.2.short = VE.order.2.dose.2.short, 
                                    VE.sd.order.2.dose.2.short = VE.sd.order.2.dose.2.short)
  saveRDS(o.order.2.dose.2.short, file = order.2.dose.2.short.path)
  
  o.order.2.dose.2.medium          <- list(VE.order.2.dose.2.medium = VE.order.2.dose.2.medium, 
                                    VE.sd.order.2.dose.2.medium = VE.sd.order.2.dose.2.medium)
  saveRDS(o.order.2.dose.2.medium, file = order.2.dose.2.medium.path)
  
  o.order.2.dose.2.long          <- list(VE.order.2.dose.2.long = VE.order.2.dose.2.long, 
                                           VE.sd.order.2.dose.2.long = VE.sd.order.2.dose.2.long)
  saveRDS(o.order.2.dose.2.long, file = order.2.dose.2.long.path)
  
  list(o.order.1.dose.2.short, o.order.1.dose.2.medium, o.order.1.dose.2.long, 
       o.order.2.dose.2.short, o.order.2.dose.2.medium, o.order.2.dose.2.long)
}


### Generating infection time etc by exponential disctributions
gener_exp_t <- function(lambda.0, lp){
  u <- runif(length(lp))
  (-log(u)) / (lambda.0*exp(lp))
}

gener_rec_tdc_data <-
  function(lambda.1, lambda.2,
           betax,
           alpha.2.inf1, alpha.2.inf2, alpha.3.inf1, alpha.3.inf2,
           d2, d3){
    theta <- rgamma(1, shape = 5, scale = 0.2) #5, 0.2) # Frailty term
    et <- rep(NA, 3)
    t1 <-  ( -log(runif(1)) )/( lambda.1*exp(betax)*theta )
    if (t1 <= d2){
      et[1] <- t1
      t2    <-  ( -log(runif(1)) )/( lambda.2*exp(betax)*theta )
      if (t1 + t2 <= d2){
        et[2] <- t1+t2; et[3] <- 1  ## inf1, inf2, || dose2, dose3
      } else {
        t3 <-  ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.2.inf2)*theta )
        if (d2+t3<=d3){
          et[2] <- d2+t3; et[3] <- 2 ## inf1, dose2, inf2, || dose3
        } else {
          t4 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.2.inf2+alpha.3.inf2)*theta )
          et[2] <- d3+t4; et[3] <- 3 ## inf1, dose2, dose3, inf2 ||
        }
      }
    } else {
      t5 <- ( -log(runif(1)) )/( lambda.1*exp(betax+alpha.2.inf1)*theta )
      if (d2+t5 <= d3){
        et[1] <- d2+t5
        t6 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.2.inf2)*theta )
        if (d2+t5+t6<=d3){
          et[2] <- d2+t5+t6; et[3] <- 4 ## dose2, inf1, inf2, || dose3
        } else {
          t7 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.2.inf2+alpha.3.inf2)*theta )
          et[2] <- d3+t7; et[3] <- 5 ## dose2, inf1, dose3, inf2 ||
        }
      } else {
        t8 <- ( -log(runif(1)) )/( lambda.1*exp(betax+alpha.2.inf1+alpha.3.inf1)*theta )
        et[1] <- d3+t8
        t9 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.2.inf2+alpha.3.inf2)*theta )
        et[2] <- d3+t8+t9; et[3] <- 6 ## dose2, dose3, inf1, inf2 ||
      }
    }
    ceiling(et) # Lead to tied infection times
  }


gener_step_waned_VE_one_inf_data <-
  function(lambda,
           betax,
           alpha.short, alpha.medium, alpha.long, ## additive effect, NOT stand-alone
           d2.short,    d2.medium,    d2.long){ 
    
    ## event times - single
    et <- rep(NA, 1)
    
    t1 <-  ( -log(runif(1)) )/( lambda*exp(betax) )
    
    if (t1 <= d2.short){ 
      et[1] <- t1 
    } else {
      t2 <- ( -log(runif(1)) )/( lambda*exp(betax+alpha.short) )
      if (d2.short+t2 <= d2.medium) {
        et[1] <- d2.short+t2 
      } else {
        t3 <- ( -log(runif(1)) )/( lambda*exp(betax+alpha.short + alpha.medium) )
        if (d2.medium+t3 <= d2.long) {
          et[1] <- d2.medium+t3
        } else {
          t4 <- ( -log(runif(1)) )/( lambda*exp(betax+alpha.short + alpha.medium + alpha.long) )
          et[1] <- d2.long+t4
        }
      }
    }
    et 
  }

gener_step_waned_VE_rec_tdc_data <-
  function(lambda.1, lambda.2,
           betax,
           alpha.short.inf1,  alpha.short.inf2, 
           alpha.medium.inf1, alpha.medium.inf2,
           alpha.long.inf1,   alpha.long.inf2,
           d2.short,    d2.medium,    d2.long){
    theta <- rgamma(1, shape = 5, scale = 0.2) #5, 0.2) # Frailty term
    et <- rep(NA, 3)
    t1 <-  ( -log(runif(1)) )/( lambda.1*exp(betax)*theta )
    if (t1 <= d2.short){
      et[1] <- t1
      t2    <-  ( -log(runif(1)) )/( lambda.2*exp(betax)*theta )
      if (t1 + t2 <= d2.short){
        et[2] <- t1+t2; et[3] <- 1  ## inf1, inf2, || dose2, dose3, dose4
      } else {
        t3 <-  ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2)*theta )
        if (d2.short+t3<=d2.medium){
          et[2] <- d2.short+t3; et[3] <- 2 ## inf1, dose2, inf2, || dose3, dose4
        } else {
          t4 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2)*theta )
          if (d2.medium+t4<=d2.long) {
            et[2] <- d2.medium+t4; et[3] <- 3 ## inf1, dose2, dose3, inf2 || dose 4
          } else {
            t5 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2+alpha.long.inf2)*theta )
            et[2] <- d2.long+t5; et[3] <- 4 ## inf1, dose2, dose3, dose 4, inf2 || 
          }
        }
      }
    } else {
      t6 <- ( -log(runif(1)) )/( lambda.1*exp(betax+alpha.short.inf1)*theta )
      if (d2.short+t6 <= d2.medium){
        et[1] <- d2.short+t6
        t7 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2)*theta )
        if (d2.short+t6+t7<=d2.medium){
          et[2] <- d2.short+t6+t7; et[3] <- 5 ## dose2, inf1, inf2, || dose3, dose4
        } else {
          t8 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2)*theta )
          if (d2.medium+t8 <=d2.long){
            et[2] <- d2.medium+t8; et[3] <- 6 ## dose2, inf1, dose3, inf2 || dose4
          } else {
            t9 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2+alpha.long.inf2)*theta )
            et[2] <- d2.long+t9; et[3] <- 7 ## dose2, inf1, dose3, dose4, inf2 || 
          }
        }
      } else {
        t10 <- ( -log(runif(1)) )/( lambda.1*exp(betax+alpha.short.inf1+alpha.medium.inf1)*theta )
        if (d2.medium+t10 <= d2.long) {
          et[1] <- d2.medium+t10
          t11 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2)*theta )
          if (d2.medium+t10+t11 <= d2.long) {
            et[2] <- d2.medium+t10+t11; et[3] <- 8 ## dose2, dose3, inf1, inf2 || dose4
          } else {
            t12 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2+alpha.long.inf2)*theta )
            et[2] <- d2.long+t12; et[3] <- 9 ## dose2, dose3, inf1, dose4, inf2 ||
          }
        } else {
          t13 <- ( -log(runif(1)) )/( lambda.1*exp(betax+alpha.short.inf1+alpha.medium.inf1+alpha.long.inf1)*theta )
          et[1] <- d2.long+t13
          t14 <- ( -log(runif(1)) )/( lambda.2*exp(betax+alpha.short.inf2+alpha.medium.inf2+alpha.long.inf2)*theta )
          et[2] <- d2.long+t13+t14; et[3] <- 10 ## dose2, dose3, dose4, inf1, inf2 ||
        }
      }
    }
    ceiling(et) # Lead to tied infection times
  }


## Misc functions for results summary
create_VE_sd_df <- function(x){
  o <- data.frame(
    seed  = as.numeric(regmatches(x, regexpr("(?<=seed\\.)\\d+", x, perl=TRUE))),
    VE    = sapply(x, function(x) get(x)[[1]]),
    VE.sd = sapply(x, function(x) get(x)[[2]]), row.names = NULL
  ) %>% arrange(seed)
  
  o
}

mse_matrix <- function(actual_vector, predicted_matrix) {
  if (length(actual_vector) != ncol(predicted_matrix)) {
    stop("Length of the actual vector must match the number of columns in the predicted matrix.")
  }
  
  # Compute MSE for each parameter across all instances
  return(apply((t(predicted_matrix) - actual_vector)^2, 1, mean))
}


