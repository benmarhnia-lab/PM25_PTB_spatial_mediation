##### Spatially Varying Contribution of Fine Particulate Matter to Racial Inequities #####
##### in Preterm Birth Across California Neighborhoods ###################################
##### Authors: Mathilde Le Vu, Sharan Kalkat, Anais Teyton, Tarik Benmarhnia, Yiqun Ma ###

library(dplyr)
library(sf)
library(sp)
library(spmoran)
library(spgwr)
library("GWmodel")

### Spatial mediation analysis

## 1. mediator variables
# full-pregnancy period and trimester-specific PM2.5 concentration)

mediator.vars <- c("pm25_lmp_full", "pm25_lmp_t1", "pm25_lmp_t2", "pm25_lmp_t3")


## 2. function to extract estimates
get.coef <- function(model, column, n, expected.cols) {
  
  b <- as.matrix(model$b_vc)
  
  if (nrow(b) != n || ncol(b) != expected.cols) {
    stop(
      "Unexpected b_vc dimensions: ",
      paste(dim(b), collapse = " x "),
      "; expected ", n, " x ", expected.cols,
      ". Column names: ", paste(colnames(b), collapse = ", ")
    )
  }
  
  as.numeric(b[, column])
}

rmse <- function(observed, predicted) {
  sqrt(mean((as.numeric(predicted) - observed)^2))
}

## 3. Run for each mediator
summary.list <- list()
i <- 0

for (mediator.var in mediator.vars) {
  
  i <- i + 1
  
  data <- data.all %>%
    transmute(
      ZCTA,
      outcome,
      exposure,
      mediator = .data[[mediator.var]],
      noedu.pct,
      house.income,
      Longitude,
      Latitude
    )

  points <- st_as_sf(
    data,
    coords = c("Longitude", "Latitude"),
    crs = 4326,
    remove = FALSE
  ) %>%
    st_transform(3310)
  
  coords <- st_coordinates(points)
  
  x.exposure1 <- as.matrix(
    data[, c(
      "exposure", "noedu.pct", "house.income"
    ), drop = FALSE]
  )
  
  x.exposure2 <- as.matrix(
    data[, "exposure", drop = FALSE]
  )
  
  x.outcome <- as.matrix(
    data[, c(
      "exposure", "mediator", "noedu.pct", "house.income"
    ), drop = FALSE]
  )
  
  ## Fit models
  t_mod <- besf_vc(
    y = data$outcome,
    x = x.exposure1,
    coords = coords,
    x_sel = FALSE,
    ncores = 10
  )
  
  m_mod <- besf_vc(
    y = data$mediator,
    x = x.exposure2,
    coords = coords,
    x_sel = FALSE,
    ncores = 10
  )
  
  y_mod <- besf_vc(
    y = data$outcome,
    x = x.outcome,
    coords = coords,
    x_sel = FALSE,
    ncores = 10
  )
  
  # get ZCTA-specific coefficients and error
  n <- nrow(data)
  
  data$TE <- get.coef(t_mod, column = 2, n = n, expected.cols = 4)
  data$a  <- get.coef(m_mod, column = 2, n = n, expected.cols = 2)

  data$DE <- get.coef(y_mod, column = 2, n = n, expected.cols = 5)
  data$b  <- get.coef(y_mod, column = 3, n = n, expected.cols = 5)
  
  data$IE <- data$a * data$b
  
  data$Y_pred <- as.numeric(y_mod$pred)
  data$error <- (data$outcome - data$Y_pred) / sd(data$outcome)
  
  data$exposure.var <- exposure.var
  data$mediator.var <- mediator.var
  
  
  ## model-fit summary
  summary.list[[i]] <- tibble(
    mediator = mediator.var,
    n = n,
    rmse.total = rmse(data$outcome, t_mod$pred),
    rmse.mediator = rmse(data$mediator, m_mod$pred),
    rmse.outcome = rmse(data$outcome, y_mod$pred)
  )
  
  # Update after every combination to retain progress
  result.summary <- bind_rows(summary.list)
  
}

print(result.summary, n = Inf)