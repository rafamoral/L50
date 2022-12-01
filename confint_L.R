library(MASS)
library(R2jags)
library(HDInterval)
library(mvtnorm)

confint_L <- function(object, p = 0.5, cf = 1:2, level = 0.95, nboot = 10000,
                      method = c("delta","fieller","proflik",
                                 "parboot","nonparboot",
                                 "bayesian","montecarlo"),
                      interval_type = c("eti","hdi","bca","all"), ...) {
  method <- match.arg(method)
  if(method %in% c("parboot","nonparboot","montecarlo")) {
    interval_type <- match.arg(interval_type)
  }
  switch(method,
         "delta" = ld_delta(object = object, p = p, cf = cf, level = level),
         "fieller" = ld_fieller(object = object, cf = cf, p = p, level = level),
         "proflik" = ld_proflik(object = object, cf = cf, p = p, level = level, ...),
         "parboot" = ld_boot(object = object, cf = cf, p = p, level = level, nboot = nboot, interval_type = interval_type),
         "nonparboot" = ld_boot_nonpar(object = object, cf = cf, p = p, level = level, nboot = nboot, interval_type = interval_type),
         "bayesian" = ld_bayesian(object = object, cf = cf, p = p, level = level, ...),
         "montecarlo" = ld_montecarlo(object = object, cf = cf, p = p, level = level, nboot = nboot, interval_type = interval_type))
}

## Delta Method
ld_delta <- function(object, p, cf, level) {
  dose_object <- dose.p(object, p = p, cf = cf)
  parm <- seq_along(dose_object)
  nam <- names(dose_object)[parm]
  se <- attr(dose_object, "SE")[parm]
  p <- attr(dose_object, "p")[parm]
  dose_object <- as.vector(dose_object[parm])
  z <- sqrt(qchisq(level, 1))
  res <- cbind(lower = dose_object - z*se,
               estimate = dose_object,
               upper = dose_object + z*se)
  row.names(res) <- paste("p = ", format(p), ":", sep = "")
  return(res)
}

## Fieller's Method
fieller_ofun <- function(xt, b, xi, v, chi2) {
  (b[1] + b[2]*xt - xi)^2/(v[1,1] + 2*xt*v[1,2] + xt^2*v[2,2]) - chi2
}

fieller_ci <- function(b, xi, v, chi2) {
  xhat <- (xi - b[1])/b[2]
  xl <- xhat - 1
  maxit <- 10000
  iter <- 1
  while(fieller_ofun(xl, b, xi, v, chi2) < 0) {
    xl <- xl - 1
    iter <- iter + 1
    if(iter > maxit) stop("maximum number of iterations exceeded")
  }
  low <- uniroot(fieller_ofun, interval = c(xl, min(xhat, xl + 1)),
                 b = b, xi = xi, v = v, chi2 = chi2)$root
  xu <- xhat + 1
  iter <- 1
  while(fieller_ofun(xu, b, xi, v, chi2) < 0) {
    xu <- xu + 1
    iter <- iter + 1
    if(iter > maxit) stop("maximum number of iterations exceeded")
  }
  upp <- uniroot(fieller_ofun, interval = c(max(xhat, xu - 1), xu),
                 b = b, xi = xi, v = v, chi2 = chi2)$root
  return(c(lower = as.vector(low), estimate = as.vector(xhat), upper = as.vector(upp)))
}

ld_fieller <- function(object, cf = 1:2, p = 0.5, level = 0.95) {
  b <- coef(object)[cf]
  V <- vcov(object)[cf, cf]
  xiv <- family(object)$linkfun(p)
  chi2 <- qchisq(level, df = 1)
  
  R <- NULL
  for(xi in xiv) {
    R <- rbind(R, fieller_ci(b, xi, V, chi2))
  }
  row.names(R) <- paste("p = ", format(p), ":", sep = "")
  structure(R, p = p, class = "Fieller")
}

print.Fieller <- function(x, ...) {
  attr(x, "p") <- class(x) <- NULL
  NextMethod("print", x, ...)
}

## Profile likelihood

prof_ofun <- Vectorize(function(theta_p, fam, Y, X0, x, etastart, wts, control, chi2, D0, off) {
  glm.fit(x = cbind(X0, x-theta_p), y = Y, weights = wts,
          etastart = etastart, offset = off, family = fam,
          control = control, intercept = TRUE)$deviance - D0 - chi2
}, "theta_p")

prof_ci <- function(b, xi, fam, Y, X0, x, etastart, wts, control, chi2, D0, off) {
  theta_p <- (xi - b[1])/b[2]
  theta_p_l <- theta_p-1
  iter <- 1
  maxit <- 10000
  while(prof_ofun(theta_p_l, fam, Y, X0, x, etastart, wts, control, chi2, D0, off) < 0) {
    theta_p_l <- theta_p_l - 1
    iter <- iter + 1
    if(iter > maxit) stop("maximum number of iterations exceeded")
  }
  low <- uniroot(prof_ofun, interval = c(theta_p_l, min(theta_p, theta_p_l+1)),
                 fam = fam, Y = Y, X0 = X0, x = x, etastart = etastart,
                 wts = wts, control = control, chi2 = chi2, D0 = D0, off = off)$root
  theta_p_u <- theta_p+1
  iter <- 1
  while(prof_ofun(theta_p_u, fam, Y, X0, x, etastart, wts, control, chi2, D0, off) < 0) {
    theta_p_u <- theta_p_u + 1
    iter <- iter + 1
    if(iter > maxit) stop("maximum number of iterations exceeded")
  }
  upp <- uniroot(prof_ofun, interval = c(max(theta_p, theta_p_u-1), theta_p_u),
                 fam = fam, Y = Y, X0 = X0, x = x, etastart = etastart,
                 wts = wts, control = control, chi2 = chi2, D0 = D0, off = off)$root
  c(lower = as.vector(low), estimate = as.vector(theta_p), upper = as.vector(upp))
}

prof <- function(b, xi, R, fam, Y, X0, x, etastart, wts, control, chi2, D0, off, ...) {
  R <- R[1,]
  theta_p <- (xi - b[1])/b[2]
  inc <- diff(range(R))
  dc <- D0 + chi2
  dev.x <- seq(theta_p - inc, theta_p + inc, length=50)
  dev.y <- prof_ofun(dev.x, fam, Y, X0, x, etastart, wts, control, chi2, D0, off) + dc
  plot(dev.x, dev.y, type="l", ylab="Deviance", las=1, ...)
  abline(h=dc, lty=2)
  arrows(R, c(0,0,0), R, c(dc,min(dev.y),dc), 
         code = 3, length = 0, lty = 2)
  points(R[c(1,3)], rep(par("usr")[3], 2), xpd = TRUE, pch = 16, cex = .8)
  points(R[2], par("usr")[3], xpd = TRUE, pch = "*", cex = 1.8)
  #axis(1, at = R[c(1,3)], labels = FALSE)
  text(x = R[c(1,3)], par("usr")[3], cex = .7,
       labels = c("lower limit","upper limit"), 
       xpd = TRUE, srt = 45, pos = 1)
}

ld_proflik <- function(object, cf = 1:2, p = 0.5, level = 0.95, profile = FALSE, ...) {
    fam <- family(object)
    Y <- object$y
    X <- model.matrix(object)
    X0 <- X[, -cf]
    x <- X[, cf[2]]
    b <- as.vector(coef(object)[cf])
    etastart <- object$linear.predictors
    wts <- weights(object)
    originalOffset <- if(is.null(o <- object$offset)) {
      0
    } else {
      o
    }
    control <- object$control
    xiv <- fam$linkfun(p)
    chi2 <- qchisq(level, 1)
    D0 <- deviance(object)
    R <- NULL
    for(xi in xiv) {
      off <- originalOffset + xi
      R <- rbind(R, prof_ci(b, xi, fam, Y, X0, x, etastart, wts, control, chi2, D0, off))
    }
    if(profile) {
      if(length(p) > 1) cat("Profile produced only for p = ", p[1], sep="", "\n")
      xi <- xiv[1]
      off <- originalOffset + xi
      prof(b, xi, R, fam, Y, X0, x, etastart, wts, control, chi2, D0, off, ...)
    }
    row.names(R) <- paste("p = ", format(p), ":", sep = "")
    structure(R, p = p, class = "LR_glm_dose")
  }

print.LR_glm_dose <- function(x, ...) {
  attr(x, "p") <- class(x) <- NULL
  NextMethod("print", x, ...)
}

## Parametric Bootstrap

ld_boot <- function(object, p = 0.5, cf = 1:2, level = 0.95, nboot = 1000,
                    interval_type = c("eti","hdi")) {
  data <- object$data
  fmla <- formula(delete.response(terms(formula(object))))
  X <- model.matrix(fmla, data = data)
  original_d_hat <- dose.p(object, p = p, cf = cf)
  
  d_hat <- matrix(NA, ncol = length(p), nrow = nboot)
  
  for(i in 1:nboot) {
    new_y <- as.matrix(simulate(object))
    new_fit <- glm(new_y ~ 0 + X, family = object$family)
    d_hat[i,] <- as.numeric(dose.p(new_fit, p = p, cf = cf))
  }
  
  res <- switch(interval_type,
                "eti" = cbind(lower = apply(d_hat, 2, quantile, (1-level)/2),
                                     estimate = original_d_hat,
                                     upper = apply(d_hat, 2, quantile, (1+level)/2)),
                "hdi" = {
                  if(length(level) == 1) {
                    hdi_interval <- apply(d_hat, 2, hdi, level)
                    return(
                      cbind(lower = hdi_interval[1,],
                            estimate = original_d_hat,
                            upper = hdi_interval[2,])
                    )
                  } else {
                    hdi_interval <- list()
                    for(l in 1:length(level)) {
                      hdi_interval[[l]] <- as.numeric(apply(d_hat, 2, hdi, level[l]))
                    }
                    ret <- do.call(rbind, hdi_interval)
                    ret <- cbind(ret, rep(original_d_hat, length(level)))
                    ret <- ret[,c(1,3,2)]
                    colnames(ret) <- c("lower","estimate","upper")
                    return(ret)}
                },
                "bca" = {
                  bca_interval <- get_bca(object = object, p = p, cf = cf,
                                          d_hat = d_hat, original_d_hat = original_d_hat,
                                          nboot = nboot, level = level)
                  bca_interval <- matrix(bca_interval, ncol = 2, nrow = length(level), byrow = FALSE)
                  return(cbind(lower = bca_interval[,1],
                               estimate = original_d_hat,
                               upper = bca_interval[,2]))
                },
                "all" = {
                  perc <- cbind(lower = apply(d_hat, 2, quantile, (1-level)/2),
                                estimate = original_d_hat,
                                upper = apply(d_hat, 2, quantile, (1+level)/2))
                  if(length(level) == 1) {
                    hdi_interval <- apply(d_hat, 2, hdi, level)
                    hdi_int <- cbind(lower = hdi_interval[1,],
                                     estimate = original_d_hat,
                                     upper = hdi_interval[2,])
                  } else {
                    hdi_interval <- list()
                    for(l in 1:length(level)) {
                      hdi_interval[[l]] <- as.numeric(apply(d_hat, 2, hdi, level[l]))
                    }
                    ret <- do.call(rbind, hdi_interval)
                    ret <- cbind(ret, rep(original_d_hat, length(level)))
                    ret <- ret[,c(1,3,2)]
                    colnames(ret) <- c("lower","estimate","upper")
                    hdi_int <- ret
                  }
                  bca_interval <- get_bca(object = object, p = p, cf = cf,
                                          d_hat = d_hat, original_d_hat = original_d_hat,
                                          nboot = nboot, level = level)
                  bca_interval <- matrix(bca_interval, ncol = 2, nrow = length(level), byrow = FALSE)
                  bca_int <- cbind(lower = bca_interval[,1],
                                   estimate = original_d_hat,
                                   upper = bca_interval[,2])
                  return(list("eti" = perc, "hdi" = hdi_int, "bca" = bca_int))
                })
  return(res)
}

## Non-parametric Bootstrap

ld_boot_nonpar <- function(object, p = 0.5, cf = 1:2, level = 0.95, nboot = 1000,
                           interval_type = c("eti","hdi","bca","all")) {
  data <- object$data
  fmla <- formula(delete.response(terms(formula(object))))
  X <- model.matrix(fmla, data = data)
  original_d_hat <- dose.p(object, p = p, cf = cf)
  
  d_hat <- matrix(NA, ncol = length(p), nrow = nboot)
  
  for(i in 1:nboot) {
    sampled_rows <- sample(1:nrow(data), nrow(data), replace = TRUE)
    new_data <- data[sampled_rows,]
    new_fit <- update(object, data = new_data)
    d_hat[i,] <- as.numeric(dose.p(new_fit, p = p, cf = cf))
  }
  
  res <- switch(interval_type,
                "eti" = cbind(lower = apply(d_hat, 2, quantile, (1-level)/2),
                                     estimate = original_d_hat,
                                     upper = apply(d_hat, 2, quantile, (1+level)/2)),
                "hdi" = {
                  if(length(level) == 1) {
                    hdi_interval <- apply(d_hat, 2, hdi, level)
                    return(
                      cbind(lower = hdi_interval[1,],
                            estimate = original_d_hat,
                            upper = hdi_interval[2,])
                    )
                  } else {
                    hdi_interval <- list()
                    for(l in 1:length(level)) {
                      hdi_interval[[l]] <- as.numeric(apply(d_hat, 2, hdi, level[l]))
                    }
                    ret <- do.call(rbind, hdi_interval)
                    ret <- cbind(ret, rep(original_d_hat, length(level)))
                    ret <- ret[,c(1,3,2)]
                    colnames(ret) <- c("lower","estimate","upper")
                    return(ret)}
                },
                "bca" = {
                  bca_interval <- get_bca(object = object, p = p, cf = cf,
                                          d_hat = d_hat, original_d_hat = original_d_hat,
                                          nboot = nboot, level = level)
                  bca_interval <- matrix(bca_interval, ncol = 2, nrow = length(level), byrow = FALSE)
                  return(cbind(lower = bca_interval[,1],
                               estimate = original_d_hat,
                               upper = bca_interval[,2]))
                },
                "all" = {
                  perc <- cbind(lower = apply(d_hat, 2, quantile, (1-level)/2),
                                estimate = original_d_hat,
                                upper = apply(d_hat, 2, quantile, (1+level)/2))
                  if(length(level) == 1) {
                    hdi_interval <- apply(d_hat, 2, hdi, level)
                    hdi_int <- cbind(lower = hdi_interval[1,],
                                     estimate = original_d_hat,
                                     upper = hdi_interval[2,])
                  } else {
                    hdi_interval <- list()
                    for(l in 1:length(level)) {
                      hdi_interval[[l]] <- as.numeric(apply(d_hat, 2, hdi, level[l]))
                    }
                    ret <- do.call(rbind, hdi_interval)
                    ret <- cbind(ret, rep(original_d_hat, length(level)))
                    ret <- ret[,c(1,3,2)]
                    colnames(ret) <- c("lower","estimate","upper")
                    hdi_int <- ret
                  }
                  bca_interval <- get_bca(object = object, p = p, cf = cf,
                                          d_hat = d_hat, original_d_hat = original_d_hat,
                                          nboot = nboot, level = level)
                  bca_interval <- matrix(bca_interval, ncol = 2, nrow = length(level), byrow = FALSE)
                  bca_int <- cbind(lower = bca_interval[,1],
                                   estimate = original_d_hat,
                                   upper = bca_interval[,2])
                  return(list("eti" = perc, "hdi" = hdi_int, "bca" = bca_int))
                })
  return(res)
}

## Bayesian credible intervals

ld_jags <- function(object, cf, p, level, progress.bar = "none",
                    n.chains = 3, n.burnin = 500, n.iter = 1000, n.thin = 5) {
  
  the_link <- object$family$link
  
  if(the_link == "logit") {
    model_code <- "
model
{
  # likelihood
  for (i in 1:N) {
    y[i] ~ dbinom(p[i], m[i])
    logit(p[i]) <- inprod(X[i,], beta)
  }
  
  # priors
  for (i in 1:N_betas) {
    beta[i] ~ dnorm(0, 0.01)
  }
  
  ld <- (p_const - beta[cf[1]]) / beta[cf[2]]
}
"
  } else if(the_link == "probit") {
    model_code <- "
model
{
  # likelihood
  for (i in 1:N) {
    y[i] ~ dbinom(p[i], m[i])
    probit(p[i]) <- inprod(X[i,], beta)
  }
  
  # priors
  for (i in 1:N_betas) {
    beta[i] ~ dnorm(0, 0.01)
  }
  
  ld <- (p_const - beta[cf[1]]) / beta[cf[2]]
}
"
  } else if(the_link == "cloglog") {
    model_code <- "
model
{
  # likelihood
  for (i in 1:N) {
    y[i] ~ dbinom(p[i], m[i])
    p[i] <- 1 - exp(-exp(inprod(X[i,], beta)))
  }
  
  # priors
  for (i in 1:N_betas) {
    beta[i] ~ dnorm(0, 0.01)
  }
  
  ld <- (p_const - beta[cf[1]]) / beta[cf[2]]
}
"
  } else {
    stop("only implemented for logit, probit or cloglog") 
  }
  
  linkfun <- object$family$linkfun
  p_const <- linkfun(p)
  
  model_data <- list(y = as.numeric(object$y * object$prior.weights),
                     m = as.numeric(object$prior.weights),
                     X = model.matrix(object),
                     N = nrow(model.matrix(object)),
                     N_betas = ncol(model.matrix(object)),
                     cf = cf, p_const = p_const)
  
  model_parameters <- "ld"
  
  model_run <- jags(data = model_data,
                    parameters.to.save = model_parameters,
                    model.file = textConnection(model_code),
                    progress.bar = progress.bar,
                    n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter, n.thin = n.thin,
                    inits = list(list(beta = coef(object) + rnorm(model_data$N_betas, 0, .001)),
                                 list(beta = coef(object) + rnorm(model_data$N_betas, 0, .001)),
                                 list(beta = coef(object) + rnorm(model_data$N_betas, 0, .001))))
  
  ld_post <- model_run$BUGSoutput$sims.list$ld
  ld_quantile <- as.numeric(apply(ld_post, 2, quantile, prob = c((1-level)/2, (1+level)/2)))
  
  return(c(ld_quantile[1], mean(ld_post), ld_quantile[2]))
}

ld_bayesian <- function(object, cf, p, level, ...) {
  R <- matrix(NA, ncol = 3, nrow = length(p))
  colnames(R) <- c("lower","estimate","upper")
  for(i in 1:length(p)) {
    R[i,] <- ld_jags(object = object, cf = cf, p = p[i], level = level, ...)
  }
  row.names(R) <- paste("p = ", format(p), ":", sep = "")
  structure(R, p = p, class = "ld_bayesian")
}

print.ld_bayesian <- function(x, ...) {
  attr(x, "p") <- class(x) <- NULL
  NextMethod("print", x, ...)
}

## Monte Carlo

ld_montecarlo <- function(object, p = 0.5, cf = 1:2, level = 0.95, nboot = 1000,
                          interval_type = c("eti","hdi","bca","all")) {
  original_d_hat <- dose.p(object, p = p, cf = cf)
  varcovar_mat <- vcov(object)[cf,cf]
  mean_vec <- coef(object)[cf]
  
  beta_sim <- rmvnorm(n = nboot, mean = mean_vec, sigma = varcovar_mat)
  
  linkfun <- object$family$linkfun
  p_const <- linkfun(p)
  
  d_hat <- matrix(NA, ncol = length(p), nrow = nboot)
  
  for(i in 1:length(p)) {
    d_hat[,i] <- (p_const[i] - beta_sim[,1]) / beta_sim[,2]
  }
  
  res <- switch(interval_type,
                "eti" = cbind(lower = apply(d_hat, 2, quantile, (1-level)/2),
                                     estimate = original_d_hat,
                                     upper = apply(d_hat, 2, quantile, (1+level)/2)),
                "hdi" = {
                  if(length(level) == 1) {
                    hdi_interval <- apply(d_hat, 2, hdi, level)
                    return(
                      cbind(lower = hdi_interval[1,],
                            estimate = original_d_hat,
                            upper = hdi_interval[2,])
                    )
                  } else {
                    hdi_interval <- list()
                    for(l in 1:length(level)) {
                      hdi_interval[[l]] <- as.numeric(apply(d_hat, 2, hdi, level[l]))
                    }
                    ret <- do.call(rbind, hdi_interval)
                    ret <- cbind(ret, rep(original_d_hat, length(level)))
                    ret <- ret[,c(1,3,2)]
                    colnames(ret) <- c("lower","estimate","upper")
                    return(ret)}
                },
                "bca" = {
                  bca_interval <- get_bca(object = object, p = p, cf = cf,
                                          d_hat = d_hat, original_d_hat = original_d_hat,
                                          nboot = nboot, level = level)
                  bca_interval <- matrix(bca_interval, ncol = 2, nrow = length(level), byrow = FALSE)
                  return(cbind(lower = bca_interval[,1],
                               estimate = original_d_hat,
                               upper = bca_interval[,2]))
                },
                "all" = {
                  perc <- cbind(lower = apply(d_hat, 2, quantile, (1-level)/2),
                                estimate = original_d_hat,
                                upper = apply(d_hat, 2, quantile, (1+level)/2))
                  if(length(level) == 1) {
                    hdi_interval <- apply(d_hat, 2, hdi, level)
                    hdi_int <- cbind(lower = hdi_interval[1,],
                                     estimate = original_d_hat,
                                     upper = hdi_interval[2,])
                  } else {
                    hdi_interval <- list()
                    for(l in 1:length(level)) {
                      hdi_interval[[l]] <- as.numeric(apply(d_hat, 2, hdi, level[l]))
                    }
                    ret <- do.call(rbind, hdi_interval)
                    ret <- cbind(ret, rep(original_d_hat, length(level)))
                    ret <- ret[,c(1,3,2)]
                    colnames(ret) <- c("lower","estimate","upper")
                    hdi_int <- ret
                  }
                  bca_interval <- get_bca(object = object, p = p, cf = cf,
                                          d_hat = d_hat, original_d_hat = original_d_hat,
                                          nboot = nboot, level = level)
                  bca_interval <- matrix(bca_interval, ncol = 2, nrow = length(level), byrow = FALSE)
                  bca_int <- cbind(lower = bca_interval[,1],
                                   estimate = original_d_hat,
                                   upper = bca_interval[,2])
                  return(list("eti" = perc, "hdi" = hdi_int, "bca" = bca_int))
                })
  return(res)
}

get_bca <- function(object, p, cf, d_hat, original_d_hat, nboot, level) {
  alpha <- c((1-level)/2, level + (1-level)/2)
  n <- nrow(object$data)
  u <- rep(0, n)
  for(i in 1:n) {
    u[i] <- as.numeric(dose.p(update(object, data = object$data[-i,]), p = p, cf = cf))
  }
  z0 <- qnorm(sum(d_hat < original_d_hat)/nboot)
  uu <- mean(u) - u
  acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
  zalpha <- qnorm(alpha)
  tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
  confpoints <- as.numeric(quantile(x = d_hat, probs = tt, type = 1))
  return(confpoints)
}