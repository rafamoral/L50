sample_size <- 100

source("simstudy_functions.R")
load("simulated_datasets100.RData")

n_sim <- 1000
effect_sizes <- c("5_perc","10_perc","25_perc")
written_levels <- c("ninety_five","eighty_four")
ci_levels <- c("ninety_five" = .95, "eighty_four" = .84)
method <- c("parboot","nonparboot","montecarlo")

sublist <- list("5_perc" = list("ninety_five" = matrix(NA, ncol = 2, nrow = n_sim),
                                "eighty_four" = matrix(NA, ncol = 2, nrow = n_sim)),
                "10_perc" = list("ninety_five" = matrix(NA, ncol = 2, nrow = n_sim),
                                 "eighty_four" = matrix(NA, ncol = 2, nrow = n_sim)),
                "25_perc" = list("ninety_five" = matrix(NA, ncol = 2, nrow = n_sim),
                                 "eighty_four" = matrix(NA, ncol = 2, nrow = n_sim)))

interval_width <- bias_L50 <- perc_overlap <-
  coverage <- overlap <- list("bootstrap_eti" = sublist,
                              "bootstrap_hdi" = sublist,
                              "bootstrap_bca" = sublist,
                              "nonpar_bootstrap_eti" = sublist,
                              "nonpar_bootstrap_hdi" = sublist,
                              "nonpar_bootstrap_bca" = sublist,
                              "montecarlo_eti" = sublist,
                              "montecarlo_hdi" = sublist,
                              "montecarlo_bca" = sublist)

set.seed(2022)

for(i in 1:n_sim) {
  for(e in effect_sizes) {
    
    the_sample <- all_datasets[[e]][[i]]
    model <- fit_model(the_sample)
    
    for(m in method) {
      for(l in written_levels) {
        
        L50 <- get_L50(model, level = ci_levels[l], method = m, interval_type = "all")
        
        m_eti <- paste0(m, "_eti")
        m_hdi <- paste0(m, "_hdi")
        m_bca <- paste0(m, "_bca")
        
        L50_eti <- list(original = L50$original$eti,
                               effect = L50$effect$eti)
        L50_hdi <- list(original = L50$original$hdi,
                        effect = L50$effect$hdi)
        L50_bca <- list(original = L50$original$bca,
                        effect = L50$effect$bca)
        
        coverage[[m_eti]][[e]][[l]][i,] <- get_coverage(interval = L50_eti, effect = e)
        coverage[[m_hdi]][[e]][[l]][i,] <- get_coverage(interval = L50_hdi, effect = e)
        coverage[[m_bca]][[e]][[l]][i,] <- get_coverage(interval = L50_bca, effect = e)
        
        overlap[[m_eti]][[e]][[l]][i,1] <- get_overlap(interval = L50_eti)
        overlap[[m_hdi]][[e]][[l]][i,1] <- get_overlap(interval = L50_hdi)
        overlap[[m_bca]][[e]][[l]][i,1] <- get_overlap(interval = L50_bca)
        
        perc_overlap[[m_eti]][[e]][[l]][i,1] <- get_perc_overlap(interval = L50_eti)
        perc_overlap[[m_hdi]][[e]][[l]][i,1] <- get_perc_overlap(interval = L50_hdi)
        perc_overlap[[m_bca]][[e]][[l]][i,1] <- get_perc_overlap(interval = L50_bca)
        
        interval_width[[m_eti]][[e]][[l]][i,] <- get_width(interval = L50_eti)
        interval_width[[m_hdi]][[e]][[l]][i,] <- get_width(interval = L50_hdi)
        interval_width[[m_bca]][[e]][[l]][i,] <- get_width(interval = L50_bca)
        
        bias_L50[[m_eti]][[e]][[l]][i,] <- get_bias(interval = L50_eti, effect = e)
        bias_L50[[m_hdi]][[e]][[l]][i,] <- get_bias(interval = L50_hdi, effect = e)
        bias_L50[[m_bca]][[e]][[l]][i,] <- get_bias(interval = L50_bca, effect = e)
        
        cat(l, e, m, i, "\n")
      }
    }
  }
}

left_side <- expand.grid(what = c("original","effect"),
                         level = written_levels,
                         effect_size = effect_sizes)[,3:1]

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
  dplyr::select(effect_size, what, bootstrap_eti)
names(bias_result)[3] <- "bias"

save.image("100boot.RData")