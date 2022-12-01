library(tidyverse)

source("ci_ld_functions.R")
load("simulated_population.RData")
load("true_L50.RData")

AICc <- function (obj) {
  aic <- AIC(obj)
  n <- obj$df.null + 1
  p <- n - obj$df.residual
  aicc <- aic + (2 * p^2 + 2 * p)/(n - p - 1)
  return(aicc)
}

BICc <- function (obj) {
  bic <- BIC(obj)
  n <- obj$df.null + 1
  p <- n - obj$df.residual
  bicc <- bic + (log(n) * (p + 1) * p)/(n - p - 1)
  return(bicc)
}

draw_sample <- function(data = walleye_population,
                        effect_size = c("5_perc","10_perc","25_perc"),
                        n) {
  effect_size <- match.arg(effect_size)
  x_range <- switch(effect_size,
                    "5_perc" = 10001:20000,
                    "10_perc" = 20001:30000,
                    "25_perc" = 30001:40000)
  sample_id1 <- sample(x = 1:10000, size = n, replace = FALSE)
  sample_id2 <- sample(x = x_range, size = n, replace = FALSE)
  sample_id <- c(sample_id1,sample_id2)
  sample_data <- data[sample_id,]
  sample_data$pop <- droplevels(sample_data$pop)
  return(sample_data)
}

fit_model <- function(data) {
  fit <- glm(maturity ~ 0 + pop + TL, family = binomial(probit), data = data)
  return(fit)
}

get_p_value <- function(model) {
  data <- model$data
  fit1 <- glm(maturity ~ 0 + TL, family = binomial(probit), data = data)
  comparison <- anova(fit1, model, test = "Chisq")
  return(list(p_value = comparison$`Pr(>Chi)`[2],
              AICc_ADD = AICc(model),
              AICc_TLO = AICc(fit1),
              BICc_ADD = BICc(model),
              BICc_TLO = BICc(fit1)))
}

get_L50 <- function(obj, p = 0.5, level = 0.95, nboot = 10000,
                    method, interval_type) {
  L50_1 <- try(confint_L(obj, p = p, cf = c(1,3), level = level, nboot = nboot,
                          method = method, interval_type = interval_type), silent = TRUE)
  L50_2 <- try(confint_L(obj, p = p, cf = c(2,3), level = level, nboot = nboot,
                          method = method, interval_type = interval_type), silent = TRUE)
  if(class(L50_1)[1] == "try-error") {
    L50_1 <- c(NA,NA,NA)
  }
  if(class(L50_2)[1] == "try-error") {
    L50_2 <- c(NA,NA,NA)
  }
  return(list("original" = L50_1,
              "effect" = L50_2))
}

get_coverage <- function(interval,
                         effect = c("5_perc","10_perc","25_perc"),
                         truth = true_L50) {
  effect <- match.arg(effect)
  orig_lower <- interval$original[1]
  orig_upper <- interval$original[3]
  eff_lower <- interval$effect[1]
  eff_upper <- interval$effect[3]
  orig_true <- true_L50[1]
  eff_true <- true_L50[effect]
  orig_in <- orig_true >= orig_lower & orig_true <= orig_upper
  eff_in <- eff_true >= eff_lower & eff_true <= eff_upper
  return(c("original" = as.logical(orig_in),
           "effect" = as.logical(eff_in)))
}

get_overlap <- function(interval) {
  return(as.logical(interval$effect[3] >= interval$original[1]))
}

get_perc_overlap <- function(interval) {
  if(as.logical(interval$effect[3] >= interval$original[1])) {
    total_width <- interval$original[3] - interval$effect[1]
    overlap_width <- interval$effect[3] - interval$original[1]
    return(overlap_width/total_width)
  } else {
    return(NA)
  }
}

get_width <- function(interval) {
  return(as.numeric(c(interval$original[3] - interval$original[1],
                      interval$effect[3] - interval$effect[1])))
}

get_bias <- function(interval,
                     effect = c("5_perc","10_perc","25_perc"),
                     truth = true_L50) {
  effect <- match.arg(effect)
  bias_original <- as.numeric(interval$original[2] - true_L50["original"])
  bias_effect <- as.numeric(interval$effect[2] - true_L50[effect])
  return(c(bias_original, bias_effect))
}