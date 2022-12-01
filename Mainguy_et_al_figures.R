## Assessing the adequacy of binomial models and their predicted uncertainty
## for the monitoring of the L50 from sampled fish
## by Mainguy et al.

## this script reproduces all figures in the paper except

library(tidyverse)
library(readxl)
library(ggpubr)

## Figure 1
fig1a <- read_excel("Figure_1A.xlsx")

pendj <- read_xlsx("PENDJ_2022-02-23.xlsx")
pendj$YEAR <- as.factor(pendj$YEAR)
pendj$NAME <- as.factor(pendj$NAME)
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

fig1a_plot <- fig1a %>%
  ggplot(aes(x = BIN5, y = Percent)) +
  theme_bw() +
  geom_bar(stat = "identity", fill = "gray90", col = 1, lwd = .3) +
  geom_bar(data = fig1a %>%
             filter(BIN5 == 465),
           stat = "identity", fill = 1, col = 1, lwd = .3, width = 4) +
  xlab("Total Length (TL) by 5-mm class (BIN5)") +
  ylab("Frequency (%)") +
  xlim(95, 860) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("A")

fig1b_plot <- fitted_mix %>%
  ggplot(aes(x = TL, y = y)) +
  theme_bw() +
  geom_ribbon(aes(ymin = 0, ymax = y),
              alpha = .2, col = 1) +
  geom_ribbon(data = fitted_dist,
              aes(ymin = 0, ymax = y * 1/8, fill = MESH),
              alpha = .4) +
  ylab("Estimated density") +
  xlab("Fish length (mm)") +
  xlim(95, 860) +
  ggtitle("B") +
  labs(fill = "Mesh size (mm)") +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 16))

png("figure1.png", res = 800, units = "in", w = 6, h = 8)
ggarrange(fig1a_plot, fig1b_plot, ncol = 1)
dev.off()

## Figure 2
fig2 <- read_excel("Figure_2.xlsx")

fig2$method <- factor(fig2$method, levels = c("Monte Carlo [BCa]",
                                              "Monte Carlo [ETI]",
                                              "Fieller",
                                              "Profile-likelihood",
                                              "Non-parametric bootstrapping [BCa]",
                                              "Parametric bootstrapping [BCa]",
                                              "Parametric bootstrapping [ETI]",
                                              "Monte Carlo [HDI]",
                                              "Delta",
                                              "Parametric bootstrapping [HDI]",
                                              "Non-parametric bootstrapping [HDI]",
                                              "Non-parametric bootstrapping [ETI]",
                                              "Bayesian") %>% rev)

fig2_plot <- fig2 %>%
  ggplot(aes(x = coverage, y = method)) +
  theme_bw() +
  geom_boxplot() +
  geom_point(data = fig2 %>%
               group_by(method) %>%
               summarise(mean = mean(coverage)),
             aes(x = mean),
             pch = "+", size = 5) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12, face = 2),
        axis.text.x = element_text(size = 12)) +
  xlab("Coverage probability") +
  ylab("") +
  geom_vline(xintercept = .84, lty = 2) +
  scale_x_continuous(breaks = seq(.77, .87, .01))

png("figure2.png", res = 800, units = "in", w = 8, h = 6)
fig2_plot
dev.off()

## Figure 3
surv1_pred <- read_excel("Figure_3.xlsx", sheet = 1) 
surv1_obs <- read_excel("Figure_3.xlsx", sheet = 2)
surv2_pred <- read_excel("Figure_3.xlsx", sheet = 3)
surv2_obs <- read_excel("Figure_3.xlsx", sheet = 4)

names(surv1_pred)[2:5] <- names(surv2_pred)[2:5] <- c("logit","probit","cloglog","cauchit")

surv1_pred <- surv1_pred %>%
  pivot_longer(2:5,
               names_to = "link",
               values_to = "pred")
surv1_pred$link <- factor(surv1_pred$link, levels = c("probit","logit","cloglog","cauchit"))

surv2_pred <- surv2_pred %>%
  pivot_longer(2:5,
               names_to = "link",
               values_to = "pred")
surv2_pred$link <- factor(surv2_pred$link, levels = c("probit","logit","cloglog","cauchit"))

fig3_plot <- surv1_pred  %>%
  ggplot() +
  theme_bw() +
  geom_line(aes(x = TL, y = pred, lty = link),
            col = "gray60") +
  geom_line(data = surv2_pred,
            aes(x = TL, y = pred, lty = link)) +
  geom_point(data = surv2_obs,
             aes(x = TL, y = PROP, size = TOT),
             shape = 1) +
  geom_point(data = surv1_obs,
             aes(x = TL, y = PROP, size = TOT),
             col = "gray40", alpha = .5) +
  geom_segment(aes(x = 181, y = .5, xend = 555.9428, yend = .5),
               lty = 3, size = .4) +
  geom_segment(aes(x = 458.1764, y = 0, xend = 458.1764, yend = .5),
               lty = 3, size = .4) +
  geom_segment(aes(x = 555.9428, y = 0, xend = 555.9428, yend = .5),
               lty = 3, size = .4) +
  geom_segment(aes(x = 444.1271, y = .25, xend = 477.6581, yend = .25),
               lty = 1, size = .8) +
  geom_segment(aes(x = 541.552, y = .25, xend = 572.4319, yend = .25),
               lty = 1, size = .8) +
  xlab("Total Length (mm)") +
  ylab("Probability to observe developed gonads") +
  scale_size(name = "Sample size") +
  scale_linetype_manual(name = "Link function",
                        values = c("solid","longdash","dotted","dotdash")) +
  theme(panel.grid = element_blank(),
        legend.position = "top")

png("figure3.png", res = 800, units = "in", w = 8, h = 6)
fig3_plot
dev.off()

## Figure 4
over_obs <- read_excel("Figure_4.xlsx", sheet = 1)
over_pred <- read_excel("Figure_4.xlsx", sheet = 2)
pone_obs <- read_excel("Figure_4.xlsx", sheet = 3)
pone_pred <- read_excel("Figure_4.xlsx", sheet = 4)

fig4a_plot <- ggplot() +
  theme_bw() +
  geom_point(data = over_obs,
             aes(x = OVERLAP, y = WIDTH_84CI)) +
  geom_line(data = over_pred,
            aes(x = OVERLAP, y = pred_O),
            lwd = 1) +
  geom_ribbon(data = over_pred,
              aes(x = OVERLAP, ymin = LL_pred_O, ymax = UL_pred_O),
              fill = "gray60", alpha = .4) +
  xlab(expression(Percentage~of~italic(overlap)~between~lenghts~with~0~and~1)) +
  ylab("84% CI width (mm)") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("A")

fig4b_plot <- ggplot() +
  theme_bw() +
  geom_point(data = pone_obs,
             aes(x = PERCONE, y = WIDTH_84CI)) +
  geom_line(data = pone_pred,
            aes(x = PERCONE, y = pred_P),
            lwd = 1) +
  geom_ribbon(data = pone_pred,
              aes(x = PERCONE, ymin = LL_pred_P, ymax = UL_pred_P),
              fill = "gray60", alpha = .4) +
  xlab(expression(Percentage~of~occurrences~"("~italic(percOne)~")")) +
  ylab("84% CI width (mm)") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("B")

png("figure4.png", res = 800, units = "in", w = 6, h = 8)
ggarrange(fig4a_plot, fig4b_plot, ncol = 1)
dev.off()