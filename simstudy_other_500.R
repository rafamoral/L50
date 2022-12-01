sample_size <- 500

source("simstudy_functions.R")
load("simulated_datasets500.RData")

n_sim <- 1000
effect_sizes <- c("5_perc","10_perc","25_perc")
written_levels <- c("ninety_five","eighty_four")
ci_levels <- c("ninety_five" = .95, "eighty_four" = .84)
method <- c("delta","fieller","proflik","bayesian")

sublist <- list("5_perc" = list("ninety_five" = matrix(NA, ncol = 2, nrow = n_sim),
                                "eighty_four" = matrix(NA, ncol = 2, nrow = n_sim)),
                "10_perc" = list("ninety_five" = matrix(NA, ncol = 2, nrow = n_sim),
                                 "eighty_four" = matrix(NA, ncol = 2, nrow = n_sim)),
                "25_perc" = list("ninety_five" = matrix(NA, ncol = 2, nrow = n_sim),
                                 "eighty_four" = matrix(NA, ncol = 2, nrow = n_sim)))

interval_width <- bias_L50 <- perc_overlap <-
  coverage <- overlap <- list("delta" = sublist,
                              "fieller" = sublist,
                              "proflik" = sublist,
                              "bayesian" = sublist)

set.seed(2022)

for(i in 1:n_sim) {
  for(e in effect_sizes) {
    
    the_sample <- all_datasets[[e]][[i]]
    model <- fit_model(the_sample)
    
    for(m in method) {
      for(l in written_levels) {
        
        L50 <- get_L50(model, level = ci_levels[l], method = m)
        
        coverage[[m]][[e]][[l]][i,] <- get_coverage(interval = L50, effect = e)
        overlap[[m]][[e]][[l]][i,1] <- get_overlap(interval = L50)
        perc_overlap[[m]][[e]][[l]][i,1] <- get_perc_overlap(interval = L50)
        interval_width[[m]][[e]][[l]][i,] <- get_width(interval = L50)
        bias_L50[[m]][[e]][[l]][i,] <- get_bias(interval = L50, effect = e)
        
        cat(l, e, m, i, "\n")
      }
    }
  }
}

left_side <- expand.grid(what = c("original","effect"),
                         level = written_levels,
                         effect_size = effect_sizes)[,3:1]

## number of times the Fieller method didn't find a root
n_problems <- cbind(left_side, sapply(coverage, function(x)
  sapply(x, function(y)
    sapply(y, function(z)
      apply(z, 2, function(w)
        sum(is.na(w)))))))

## proportion of times the L50 CIs for the two populations overlapped
overlap_result <- cbind(left_side, sapply(overlap, function(x)
  sapply(x, function(y)
    sapply(y, function(z)
      apply(z, 2, function(w)
        sum(na.omit(w))/sum(!is.na(w)))))))

## percentage overlap when the L50 CIs for the two populations overlapped relative to the width of both CIs
perc_overlap_result <- cbind(left_side, sapply(perc_overlap, function(x)
  sapply(x, function(y)
    sapply(y, function(z)
      apply(z, 2, function(w)
        sum(na.omit(w))/sum(!is.na(w)))))))

## CI coverage
coverage_result <- cbind(left_side, sapply(coverage, function(x)
  sapply(x, function(y)
    sapply(y, function(z)
      apply(z, 2, function(w)
        sum(na.omit(w))/sum(!is.na(w)))))))

## CI width
width_result <- cbind(left_side, sapply(interval_width, function(x)
  sapply(x, function(y)
    sapply(y, function(z)
      apply(z, 2, function(w)
        sum(na.omit(w))/sum(!is.na(w)))))))

## L50 average bias
bias_result <- cbind(left_side, sapply(bias_L50, function(x)
  sapply(x, function(y)
    sapply(y, function(z)
      apply(z, 2, function(w)
        sum(na.omit(w))/sum(!is.na(w))))))) %>%
  as.data.frame %>%
  filter(level == "ninety_five") %>%
  dplyr::select(effect_size, what, delta, bayesian)
names(bias_result)[3:4] <- c("bias_others", "bias_bayesian")

save.image("500other.RData")