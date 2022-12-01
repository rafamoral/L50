sample_size <- 100

source("simstudy_functions.R")

n_sim <- 1000
effect_sizes <- c("5_perc","10_perc","25_perc")

all_datasets <- all_models <- convergence_flag <-
  p_values <- prop_one <- TL_overlap <- 
  model_AICc <- model_BICc <- list("5_perc" = list(),
                                   "10_perc" = list(),
                                   "25_perc" = list())

set.seed(2022)

for(i in 1:n_sim) {
  for(e in effect_sizes) {
    
    the_sample <- draw_sample(n = sample_size, effect_size = e)
    model <- fit_model(the_sample)
    convergence_flag <- model$converged
    
    while(!convergence_flag) {
      the_sample <- draw_sample(n = sample_size, effect_size = e)
      model <- fit_model(the_sample)
      convergence_flag <- model$converged
    }
    
    all_datasets[[e]][[i]] <- the_sample
    
    inf_criteria <- get_p_value(model)
    p_values[[e]][[i]] <- inf_criteria$p_value
    
    model_AICc[[e]][[i]] <- c("ADD" = inf_criteria$AICc_ADD,
                              "TLO" = inf_criteria$AICc_TLO)
    model_BICc[[e]][[i]] <- c("ADD" = inf_criteria$BICc_ADD,
                              "TLO" = inf_criteria$BICc_TLO)
    
    prop_one[[e]][[i]] <- sum(the_sample$maturity[the_sample$pop != "original"])/sample_size
    
    TL_overlap[[e]][[i]] <- (max(the_sample$TL[the_sample$maturity == 0 & the_sample$pop != "original"]) -
                             min(the_sample$TL[the_sample$maturity == 1 & the_sample$pop != "original"])) /
                            (diff(range(the_sample$TL[the_sample$pop != "original"])))
    
  }
  cat(i,"\n")
}

## mean p-values for the LRT between TLO and ADD and percentage it was < 0.05
p_values_result <- data.frame("mean" = sapply(lapply(p_values, unlist), mean),
                              "perc_less_0.05" = sapply(lapply(p_values, unlist),
                                                        function(x) sum(x < 0.05)/length(x)*100))

## mean AICc and percentage the ADD model had a lower AICc than TLO
AICc_mean_results <- sapply(lapply(model_AICc, function(x) do.call(rbind, x)), colMeans)
AICc_perc_results <- sapply(lapply(model_AICc, function(x) do.call(rbind, x)),
                            function(y) sum(y[,1] < y[,2])/nrow(y)*100)

## mean BICc and percentage the ADD model had a lower BICc than TLO
BICc_mean_results <- sapply(lapply(model_BICc, function(x) do.call(rbind, x)), colMeans)
BICc_perc_results <- sapply(lapply(model_BICc, function(x) do.call(rbind, x)),
                            function(y) sum(y[,1] < y[,2])/nrow(y)*100)

## mean proportion of 1s in the sample
prop_one_results <- lapply(prop_one, unlist)
prop_one_means <- sapply(prop_one_results, mean)

## ratio overlap between 0s and 1s wrt TL range
TL_overlap_results <- lapply(TL_overlap, unlist)
TL_overlap_means <- sapply(TL_overlap_results, mean)

save.image("simulated_datasets100.RData")