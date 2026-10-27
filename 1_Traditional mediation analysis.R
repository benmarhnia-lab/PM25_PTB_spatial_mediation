##### Spatially Varying Contribution of Fine Particulate Matter to Racial Inequities #####
##### in Preterm Birth Across California Neighborhoods ###################################
##### Authors: Mathilde Le Vu, Sharan Kalkat, Anais Teyton, Tarik Benmarhnia, Yiqun Ma ###

library(dplyr)
library(sf)
library(sp)
library(spmoran)
library(spgwr)
library("GWmodel")

### Traditional mediation analysis

## 1. mediator variables
# full-pregnancy period and trimester-specific PM2.5 concentration)

mediator.vars <- c("pm25_lmp_full", "pm25_lmp_t1", "pm25_lmp_t2", "pm25_lmp_t3")


## 2. function to extract estimates
estimate.effects <- function(data) {
  
  model1 <- lm(outcome ~ exposure + noedu.pct + house.income, data = data)
  
  model2 <- lm(mediator ~ exposure, data = data)
  
  model3 <- lm(
    outcome ~ exposure + mediator + noedu.pct + house.income, data = data
  )
  
  a <- coef(model2)[["exposure"]]
  
  c(
    total = coef(model1)[["exposure"]],
    exposure.mediator = a,
    direct = coef(model3)[["exposure"]],
    mediator.outcome = coef(model3)[["mediator"]],
    indirect = a * coef(model3)[["mediator"]]
  )
}

## 3. Bootstrap settings
set.seed(2026)
n.boot <- 1000L

result.list <- list()
i <- 0L

## 4. run models for all mediators
for (mediator.var in mediator.vars) {
  
  data <- data.all %>%
    transmute(
      outcome = outcome,
      exposure = exposure,
      mediator = .data[[mediator.var]],
      noedu.pct = noedu.pct,
      house.income = house.income
    )
  

  # point estimates from the original sample
  estimates <- estimate.effects(data)
  
  # Bootstrap
  boot.estimates <- replicate(n.boot, {
    index <- sample.int(n, size = n, replace = TRUE)
    estimate.effects(data[index, , drop = FALSE])
  })
  
  # 95% confidence intervals
  ci <- t(apply(
    boot.estimates,
    MARGIN = 1,
    FUN = quantile,
    probs = c(0.025, 0.975),
    names = FALSE
  ))
  
  # summarize results
  result.row <- tibble(
    mediator = mediator.var,
    n = n,
    n.boot = n.boot
  )
  
  for (effect in names(estimates)) {
    result.row[[effect]] <- unname(estimates[[effect]])
    result.row[[paste0(effect, ".lower")]] <- ci[effect, 1]
    result.row[[paste0(effect, ".upper")]] <- ci[effect, 2]
  }
  
  i <- i + 1L
  result.list[[i]] <- result.row
}

result.table <- bind_rows(result.list)

print(result.table, n = Inf, width = Inf)
