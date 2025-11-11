
{
library(tidyverse)
library(survival)

rm(list = ls())

source("~/rwd2023/TND-Cox/R/000-func-definition.R")

## Read rds files
{
  
  ### 
  
  list.fit.order.1.dose.2.short <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.1.dose.2.short/")
  df.fit.order.1.dose.2.short <- create_VE_sd_df(list.fit.order.1.dose.2.short)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.order.1.dose.2.medium <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.1.dose.2.medium/")
  df.fit.order.1.dose.2.medium <- create_VE_sd_df(list.fit.order.1.dose.2.medium)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.order.1.dose.2.long <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.1.dose.2.long/")
  df.fit.order.1.dose.2.long <- create_VE_sd_df(list.fit.order.1.dose.2.long)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  
  list.fit.order.2.dose.2.short <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.2.dose.2.short/")
  df.fit.order.2.dose.2.short <- create_VE_sd_df(list.fit.order.2.dose.2.short)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.order.2.dose.2.medium <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.2.dose.2.medium/")
  df.fit.order.2.dose.2.medium <- create_VE_sd_df(list.fit.order.2.dose.2.medium)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
  list.fit.order.2.dose.2.long <-
    ssfunc::readRDS_dir("~/rwd2023/TND-Cox/rds/simu-step-wane/tnd-cox-full-frailty/order.2.dose.2.long/")
  df.fit.order.2.dose.2.long <- create_VE_sd_df(list.fit.order.2.dose.2.long)
  rm(list = ls()[ls() %>% str_starts("fit.")])
  
}

{
df.fit.order.1.dose.2.short$trueVE  <- 0.5
df.fit.order.1.dose.2.medium$trueVE <- 0.4
df.fit.order.1.dose.2.long$trueVE   <- 0.3

df.fit.order.2.dose.2.short$trueVE  <- 0.8
df.fit.order.2.dose.2.medium$trueVE <- 0.7
df.fit.order.2.dose.2.long$trueVE   <- 0.6

df.fit.order.1.dose.2.short$order <- 
  df.fit.order.1.dose.2.medium$order <- 
  df.fit.order.1.dose.2.long$order <- "Initial Infection"

df.fit.order.2.dose.2.short$order  <- 
  df.fit.order.2.dose.2.medium$order  <- 
  df.fit.order.2.dose.2.long$order  <- "Reinfection"

df.fit.order.1.dose.2.short$term <- 
  df.fit.order.2.dose.2.short$term <- "Short Term"

df.fit.order.1.dose.2.medium$term <- 
  df.fit.order.2.dose.2.medium$term <- "Medium Term"

df.fit.order.1.dose.2.long$term <- 
  df.fit.order.2.dose.2.long$term <- "Long Term"

df.fit.order.1.dose.2.short$xdays <- 
  df.fit.order.2.dose.2.short$xdays <- 45

df.fit.order.1.dose.2.medium$xdays <- 
  df.fit.order.2.dose.2.medium$xdays <- 135

df.fit.order.1.dose.2.long$xdays <- 
  df.fit.order.2.dose.2.long$xdays <- 225


df.all <- rbind(df.fit.order.1.dose.2.short,
                df.fit.order.1.dose.2.medium,
                df.fit.order.1.dose.2.long,
                df.fit.order.2.dose.2.short,
                df.fit.order.2.dose.2.medium,
                df.fit.order.2.dose.2.long)

df.all$VE.ci.low <- df.all$VE - 1.96*df.all$VE.sd 
df.all$VE.ci.upp <- df.all$VE + 1.96*df.all$VE.sd 
df.all$VE.cover  <- df.all$trueVE >= df.all$VE.ci.low & df.all$trueVE <= df.all$VE.ci.upp



# Calculate mean and quantiles
quantile.df.all <- df.all %>%
  group_by(order, term) %>%
  summarise(xdays = mean(xdays),
            mean_VE = mean(VE),
            ci_upp = quantile(VE, 0.975),
            ci_low = quantile(VE, 0.025)) 

# Calculate Summ
summ.df.all <- df.all %>%
  group_by(order, term) %>%
  summarise(xdays = mean(xdays),
            trueVE = mean(trueVE),
            VE.avg = mean(VE),
            VE.bias.emp = (mean(VE) - trueVE),
            VE.se.avg = mean(VE.sd),
            VE.sd.emp = sd(VE),  
            CI.cov.95 = mean(VE.cover),
            MSE = mean((VE-trueVE)^2)) 


}

quantile.df.all$CI.cov.95 <- summ.df.all$CI.cov.95
# Define the data for the piecewise function
trueVE.dat <- data.frame(
  x = c(0,  90,  180,  300,    0,  90,  180,  300),
  y = c(0.5, 0.4,  0.3,  0.3,  0.8, 0.7,  0.6,  0.6),
  order = c(rep("Initial Infection", 4),
            rep("Reinfection", 4))
)

# Create the plot
ggplot() +
  geom_step(data = trueVE.dat, 
            aes(x = x, y = y, linetype = order)) + 
  scale_linetype_manual(values = c("dotdash", "dashed")) +  
  geom_errorbar(data = quantile.df.all,
                aes(x = xdays,
                    ymin = ci_low,
                    ymax = ci_upp),
                width = 10) +  # Add error bars
  geom_point(data = quantile.df.all, 
             aes(x = xdays, 
                 y = mean_VE), 
             shape = 8, 
             size = 2) + 
  geom_text(data = quantile.df.all,
            aes(x = xdays,
                y = mean_VE,
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
            vjust = 2, hjust = 0.0, size = 3, color = "black") +  # Add labels for 95% CI coverage
  scale_x_continuous(breaks = c(0, 90, 180, 270), limits = c(0, 300)) +
  scale_y_continuous(breaks = c(0, 0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9), limits = c(0.15, 0.84), labels = scales::percent) +  
  # facet_wrap(~order, scales = "free", 
  #            labeller = as_labeller(
  #              c(`1st Inf` = "Initial Infection",
  #                `2nd Inf` = "Reinfection")
  #            )) +
  theme_bw() +
  labs(x = "Days Since Vacciation",
       y = "VE") + 
  theme(
         legend.position = c(.95, .95),
         legend.justification = c("right", "top"),
         legend.box.background = element_rect(fill="white", color="black"),
         legend.box.just = "right",
         legend.margin = margin(6, 6, 6, 6),
         legend.title=element_blank()
       )


ggsave("~/rwd2023/TND-Cox/plot/TND-Simu-Stepwane-Results.pdf", 
       width = 20,
       height = 10,
       units = "cm")  
  
}
  
summ.df.all %>%
  mutate(
    VE.avg      = paste0(sprintf("%.2f", round(VE.avg, 4)*100), "%"),
    VE.bias.emp = paste0(sprintf("%.2f", round(VE.bias.emp, 4)*100), "%"),
    VE.se.avg   = paste0(sprintf("%.2f", round(VE.se.avg, 4)*100), "%"),
    VE.sd.emp   = paste0(sprintf("%.2f", round(VE.sd.emp, 4)*100), "%"),
    CI.cov.95   = paste0(sprintf("%.2f", round(CI.cov.95, 4)*100), "%"),
    MSE         = paste0(sprintf("%.2f", round(MSE, 4)*100), "%")
  )
  


