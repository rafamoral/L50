## Assessing the adequacy of binomial models and their predicted uncertainty
## for the monitoring of the L50 from sampled fish
## by Mainguy et al.

## this script reproduces the simulated "true" population of fish used in the
## simulation study within the paper

library(hnp)
library(lme4)
library(readxl)
library(tidyverse)

pendj <- read_xlsx("PENDJ_2022-02-23.xlsx")
pendj$YEAR <- as.factor(pendj$YEAR)
pendj$NAME <- as.factor(pendj$NAME)

fit0 <- glm(MATURITY ~ BIN5,
            family = binomial(probit),
            data = pendj)

- coef(fit0)[1] / coef(fit0)[2]

fit <- glmer(MATURITY ~ BIN5 + (0 + YEAR || NAME),
             family = binomial(probit),
             data = pendj,
             nAGQ = 0)

- fixef(fit)[1] / fixef(fit)[2]

## Mesh size distribution for resampling
pendj2 <- pendj[pendj$MESH != "NA",]
pendj2$MESH <- as.factor(pendj2$MESH)

j <- 1
the_density <- dens_fun <- dens_fun2 <- list()

for(i in levels(pendj2$MESH)) {
  the_density[[j]] <- pendj2 %>%
    filter(MESH == i) %>%
    pull(TL) %>%
    density(x = .,
            kernel = "gaussian",
            bw = "SJ",
            n = 2^15)
  the_density[[j]]$x <- c(the_density[[j]]$x,
                          min(the_density[[j]]$x) - 1e-10,
                          max(the_density[[j]]$x) + 1e-10)
  the_density[[j]]$y <- c(the_density[[j]]$y, 0, 0)
  dens_fun[[j]] <- approxfun(the_density[[j]])
  j <- j + 1
}

dens_fun2[[1]] <- Vectorize(function(x) {
  y <- dens_fun[[1]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[2]] <- Vectorize(function(x) {
  y <- dens_fun[[2]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[3]] <- Vectorize(function(x) {
  y <- dens_fun[[3]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[4]] <- Vectorize(function(x) {
  y <- dens_fun[[4]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[5]] <- Vectorize(function(x) {
  y <- dens_fun[[5]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[6]] <- Vectorize(function(x) {
  y <- dens_fun[[6]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[7]] <- Vectorize(function(x) {
  y <- dens_fun[[7]](x)
  if(is.na(y)) return(0) else return(y)
})
dens_fun2[[8]] <- Vectorize(function(x) {
  y <- dens_fun[[8]](x)
  if(is.na(y)) return(0) else return(y)
})

fitted_dist <- tibble(TL = rep(95:999, 8),
                      MESH = factor(rep(levels(pendj2$MESH), each = length(95:999)),
                                    levels = levels(pendj2$MESH)),
                      y = c(dens_fun2[[1]](TL[MESH == "25"]),
                            dens_fun2[[2]](TL[MESH == "38"]),
                            dens_fun2[[3]](TL[MESH == "51"]),
                            dens_fun2[[4]](TL[MESH == "64"]),
                            dens_fun2[[5]](TL[MESH == "76"]),
                            dens_fun2[[6]](TL[MESH == "102"]),
                            dens_fun2[[7]](TL[MESH == "127"]),
                            dens_fun2[[8]](TL[MESH == "152"])))

pendj2 %>%
  ggplot(aes(x = TL, fill = MESH, colour = MESH)) +
  theme_bw() +
  geom_histogram(aes(y = ..density..), alpha = .5) +
  facet_wrap(~ MESH) +
  geom_line(data = fitted_dist, aes(y = y)) +
  ylab("Estimated density") +
  xlab("Fish length (mm)") +
  ggtitle("Individual densities")

fitted_mix <- tibble(TL = 95:999,
                     y = 1/8 * dens_fun2[[1]](TL) +
                         1/8 * dens_fun2[[2]](TL) +
                         1/8 * dens_fun2[[3]](TL) +
                         1/8 * dens_fun2[[4]](TL) +
                         1/8 * dens_fun2[[5]](TL) +
                         1/8 * dens_fun2[[6]](TL) +
                         1/8 * dens_fun2[[7]](TL) +
                         1/8 * dens_fun2[[8]](TL))

fitted_dist$MESH <- as.numeric(as.character(fitted_dist$MESH))
fitted_dist$MESH <- factor(fitted_dist$MESH, levels = sort(unique(fitted_dist$MESH)))

png("fig1b.png", res = 800, units = "in", w = 8, h = 6)
fitted_mix %>%
  ggplot(aes(x = TL, y = y)) +
  theme_bw() +
  geom_ribbon(aes(ymin = 0, ymax = y),
              alpha = .2, col = 1) +
  geom_ribbon(data = fitted_dist,
              aes(ymin = 0, ymax = y * 1/8, fill = MESH),
              alpha = .4) +
  ylab("Estimated density") +
  xlab("Fish length (mm)") +
  ggtitle("Gaussian kernel density estimate") +
  labs(fill = "Mesh size (mm)")
dev.off()

draw_nonpar_mix <- function(n, weights = rep(1/8, 8), data = pendj2) {
  y <- numeric(n)
  for(i in 1:n) {
    which_dens <- sample(1:length(weights), 1, prob = weights)
    the_mesh <- levels(data$MESH)[which_dens]
    the_limits <- range(data$TL[data$MESH == the_mesh])
    the_grid <- seq(the_limits[1], the_limits[2], length.out = 1e3)
    y[i] <- sample(x = the_grid, size = 1, prob = dens_fun2[[which_dens]](the_grid))
    cat("sample", i, "of", n, "\n")
  }
  return(y)
}

y_sim <- tibble(TL = draw_nonpar_mix(n = 1000))

y_sim %>%
  ggplot(aes(x = TL)) +
  theme_bw() +
  geom_histogram(aes(y = ..density..), alpha = .4, bins = 50, col = 1, cex = .1) +
  geom_line(data = fitted_mix, aes(y = y)) +
  ylab("Density") +
  xlab("Fish length (mm)") +
  ggtitle("Simulated from mixture with equal weights") +
  geom_ribbon(data = fitted_dist,
              aes(ymin = 0, ymax = y * 1/8, fill = MESH),
              alpha = .4) +
  labs(fill = "Mesh size (mm)")

## creating fictive population

set.seed(2022)
walleye_population <- tibble(TL = draw_nonpar_mix(n = 40000))

true_beta1 <- true_beta2 <- true_beta3 <- true_beta4 <- fixef(fit)
true_beta2[1] <- true_beta2[1] * .95
true_beta3[1] <- true_beta3[1] * .9
true_beta4[1] <- true_beta4[1] * .75

true_L50_1 <- - true_beta1[1] / true_beta1[2]
true_L50_2 <- - true_beta2[1] / true_beta2[2]
true_L50_3 <- - true_beta3[1] / true_beta3[2]
true_L50_4 <- - true_beta4[1] / true_beta4[2]

true_L50 <- c("original" = as.numeric(true_L50_1),
              "5_perc" = as.numeric(true_L50_2),
              "10_perc" = as.numeric(true_L50_3),
              "25_perc" = as.numeric(true_L50_4))

X1 <- model.matrix(~ TL, data = walleye_population[1:10000,])
X2 <- model.matrix(~ TL, data = walleye_population[10001:20000,])
X3 <- model.matrix(~ TL, data = walleye_population[20001:30000,])
X4 <- model.matrix(~ TL, data = walleye_population[30001:40000,])

p1 <- pnorm(as.numeric(X1 %*% true_beta1))
p2 <- pnorm(as.numeric(X2 %*% true_beta2))
p3 <- pnorm(as.numeric(X3 %*% true_beta3))
p4 <- pnorm(as.numeric(X4 %*% true_beta4))

walleye_population$pop <- gl(4, 10000, labels = c("original","5_perc","10_perc","25_perc"))

set.seed(2022)
walleye_population$maturity <- c(rbinom(10000, 1, p1),
                                 rbinom(10000, 1, p2),
                                 rbinom(10000, 1, p3),
                                 rbinom(10000, 1, p4))

save(true_L50, file = "true_L50.RData")
save(walleye_population, file = "simulated_population.RData")
save.image("original_fit_all_objects.RData")