##### Main Functions

#' Create a Poisson .stan model
#'
#' This is an internal function used to create a Stan model.
#' @param P The highest power of WBS signal in regression
#' @return The .stan model
#' @examples
#' cat(createStanPoisson(2))

createStanPoisson = function(P = 2) {
  scode = "data {
    int<lower=0> N;
    int<lower=0> L;
    int<lower=0> Q;
    vector[L] evalues;
    vector[N] WWmean;
    matrix[N, L] WW;
    matrix[N, Q] covariate;
    vector[N] offset;
    int<lower=0> y[N];
  }

  parameters {
    real beta0;
    vector[L] betaL;
    vector[2] betaP;
    vector[Q] betaQ;
    vector[N] phi;
    real <lower=0> sigma_beta;
    real <lower=0> sigma_betaL;
    real alpha1;
    real <lower=0> sigma_phi;
  }

  transformed parameters  {
    vector[N] mu;
    vector[N] wbe;
    wbe = WWmean + WW * betaL;
    mu = offset + beta0 + phi + sqrt(covariate) * betaQ;
    mu = mu + betaP[1] * wbe + betaP[2] * wbe .* wbe;
  }

  model {
    beta0 ~ normal(0, 10000);
    for (l in 1:L) {
      betaL[l] ~ normal(0, evalues[l]);
    }
    for (p in 1:2) {
      betaP[p] ~ normal(0, sigma_beta);
    }
    for (q in 1:Q) {
      betaQ[q] ~ normal(0, sigma_beta);
    }
    sigma_beta ~ inv_gamma(0.1, 0.1);
    phi[1] ~ normal(0, sigma_phi);
    for (i in 2:N) {
      phi[i] ~ normal(alpha1*phi[i-1], sigma_phi);
    }
    alpha1 ~ uniform(-1, 1);
    sigma_phi ~ inv_gamma(0.1, 0.1);
    y ~ poisson(exp(mu));
  }"

  return(scode)
}


#' Create a Negative Binomial .stan model
#'
#' This is an internal function used to create a Stan model.
#' @param P The highest power of WBS signal in regression
#' @return The .stan model
#' @examples
#' cat(createStanNB(2))

createStanNB = function(P = 2) {
  scode = "data {
    int<lower=0> N;
    int<lower=0> L;
    int<lower=0> Q;
    vector[L] evalues;
    vector[N] WWmean;
    matrix[N, L] WW;
    matrix[N, Q] covariate;
    vector[N] offset;
    int<lower=0> y[N];
  }

  parameters {
    real beta0;
    vector[L] betaL;
    vector[2] betaP;
    vector[Q] betaQ;
    vector[N] phi;
    real <lower=0> sigma_beta;
    real <lower=0> sigma_betaL;
    real alpha1;
    real <lower=0> sigma_phi;
  }

  transformed parameters  {
    vector[N] mu;
    vector[N] wbe;
    wbe = WWmean + WW * betaL;
    mu = offset + beta0 + phi + sqrt(covariate) * betaQ;
    mu = mu + betaP[1] * wbe + betaP[2] * wbe .* wbe;
  }

  model {
    beta0 ~ normal(0, 10000);
    for (l in 1:L) {
      betaL[l] ~ normal(0, evalues[l]);
    }
    for (p in 1:2) {
      betaP[p] ~ normal(0, sigma_beta);
    }
    for (q in 1:Q) {
      betaQ[q] ~ normal(0, sigma_beta);
    }
    sigma_beta ~ inv_gamma(0.1, 0.1);
    phi[1] ~ normal(0, sigma_phi);
    for (i in 2:N) {
      phi[i] ~ normal(alpha1*phi[i-1], sigma_phi);
    }
    alpha1 ~ uniform(-1, 1);
    sigma_phi ~ inv_gamma(0.1, 0.1);
    y ~ neg_binomial_2(exp(mu), 100);
  }"

  return(scode)
}


#' Main Poisson regression function
#'
#' This is the main function for fitting Bayesian Poisson regression model.
#' @details
#' This is the main function used to fit a Bayesian regression model with Poisson-distributed response
#' to predict clinical cases.
#' See ?predictPoisson for how to make predictions of clinical cases from Bayesian regression model.
#' @param modeldata The long-format data frame/table containing date, WBS observation, and case count
#' @param date Name of date column
#' @param p.cases Name of case count column
#' @param p.rate Name of positivity rate (# of cases / # of people tested) column
#' @param ww.signal Names of WBS signal column
#' @param lag The time lag between WBS signal and case count
#' @param start_date Starting date of training data used to build model
#' @param end_date End date of training data used to build model
#' @param StanModel Created Stan model
#' @param iteration A positive integer specifying the number of iterations for each chain (including burnin).
#' @param distribution Distribution of outcome variable: "Poisson" or "NB".
#' @return fit The fitted MCMC object
#' @examples
#' modeldata = as.data.table(readRDS("modeldata.rds"))
#' modeldata = modeldata[modeldata$date <= as.Date("2021-12-20"), ]
#' modelres = PoissonReg(modeldata, date, p.cases, p.rate, ww.signal, lag,
#'                       start_date, end_date, iteration)

PoissonReg = function(modeldata, date, p.cases, p.rate, ww.signal, lag = 1,
                      start_date, end_date, StanModel, iteration, distribution = "Poisson") {
  ldata = as.data.frame(modeldata[, names(modeldata) %in% c(date, ww.signal), with = FALSE])
  cdata = as.data.frame(modeldata[, names(modeldata) %in% c(date, p.cases, p.rate), with = FALSE])
  cdata$offset = round(cdata$clinical.cases / cdata$positivity.rate * 100)

  vc = as.data.frame(t(modeldata[, names(modeldata) %in% ww.signal, with = FALSE]))
  names(vc) = modeldata$date
  # Fit FPCA to impute NAs
  fpca.fit = fpca.sc(Y = as.matrix(vc), npc = 2)
  L = fpca.fit$npc
  Ypred = as.data.frame(cbind(fpca.fit$mu, fpca.fit$efunctions))
  names(Ypred) = c("ww_mean", paste("ww_signal", 1:L, sep = ""))
  ldata[, "collect_date"] = ldata$date + lag
  ldata = cbind(ldata, Ypred)
  finaldata = cdata
  finaldata = merge(finaldata, ldata[, names(ldata) %in% c("collect_date", "ww_mean",
                                                           paste("ww_signal", 1:L, sep = ""))],
                    by.x = "date", by.y = "collect_date", all = TRUE)
  names(finaldata) = c("date", "clinical.cases", "positivity.rate", "offset", "ww_mean",
                       paste("ww_signal", 1:L, sep = ""))
  finaldata = finaldata[finaldata$date >= start_date & finaldata$date <= end_date, ]
  finaldata = na.omit(finaldata)

  data = list(P = 2,
              N = dim(finaldata)[1],
              L = L,
              evalues = fpca.fit$evalues,
              WWmean = finaldata$ww_mean,
              WW = as.data.frame(finaldata[, grep("ww_signal", names(finaldata))]),
              offset = log(finaldata$offset),
              y = finaldata[, which(names(finaldata) %in% p.cases)])
  if (distribution = "Poisson") {
    regression_model = stan_model(createStanPoisson(P = data$P))
  } else if (distribution = "NB") {
    regression_model = stan_model(createStanNB(P = data$P))
  } else {
    stop("Outcome distribution not supported")
  }
  fit = rstan::sampling(regression_model, data = data, chains = 2, iter = iteration, refresh = 0)

  return(list(lag = lag,
              fpca.fit = fpca.fit,
              fit = fit,
              fitdata = finaldata))
}


#' Main prediction function
#'
#' This is the main prediction function for a fitted Poisson regression model.
#' @details
#' This function is used to make predictions using a fitted Poisson regression model.
#' @param model_res The model result from PoissonReg()
#' @param modeldata The long-format data frame/table containing date, WBS observation, and case count
#' @param pred_start Starting date of testing data used to make predictions
#' @param pred_end End date of testing data used to make predictions
#' @param iteration A positive integer specifying the number of iterations for each chain (including burnin).
#' @return Ypred The forecast of observed virus concentration
#' @examples
#' predres = predictPoisson(modelres, modeldata, pred_start, pred_end, iteration)

predictPoisson = function(modelres, modeldata, pred_start, pred_end, iteration) {
  lag = modelres$lag
  fpca.fit = modelres$fpca.fit
  fit = modelres$fit
  draws = extract(fit, permuted = FALSE, inc_warmup = FALSE)

  ldata = as.data.frame(modeldata[, names(modeldata) %in% c(date, ww.signal), with = FALSE])
  cdata = as.data.frame(modeldata[, names(modeldata) %in% c(date, p.cases, p.rate), with = FALSE])
  cdata$offset = round(cdata$clinical.cases / cdata$positivity.rate * 100)

  vc = as.data.frame(t(modeldata[, names(modeldata) %in% ww.signal, with = FALSE]))
  names(vc) = modeldata$date
  # Fit FPCA to impute NAs
  fpca.fit = fpca.sc(Y = as.matrix(vc), npc = 2)
  L = fpca.fit$npc
  Ypred = as.data.frame(cbind(fpca.fit$mu, fpca.fit$efunctions))
  names(Ypred) = c("ww_mean", paste("ww_signal", 1:L, sep = ""))
  ldata[, "collect_date"] = ldata$date + lag
  ldata = cbind(ldata, Ypred)
  finaldata = cdata
  finaldata = merge(finaldata, ldata[, names(ldata) %in% c("collect_date", "ww_mean",
                                                           paste("ww_signal", 1:L, sep = ""))],
                    by.x = "date", by.y = "collect_date", all = TRUE)
  names(finaldata) = c("date", "clinical.cases", "positivity.rate", "offset", "ww_mean",
                       paste("ww_signal", 1:L, sep = ""))
  finaldata = finaldata[finaldata$date >= pred_start & finaldata$date <= pred_end, ]
  finaldata = na.omit(finaldata)
  finalX = finaldata[, grep("ww_signal", names(finaldata))]
  ww_mean = finaldata$ww_mean

  chain = as.data.frame(draws[(dim(draws)[1] - 2500):(dim(draws)[1]), 1, ])
  predmat = matrix(NA, dim(finaldata)[1], dim(chain)[1])
  for (iter in 1:dim(predmat)[2]) {
    beta0 = chain[iter, "beta0"]
    betaL = as.numeric(chain[iter, grep("betaL[[].[]]", names(chain))])
    betaP = as.numeric(chain[iter, grep("betaP[[].[]]", names(chain))])
    phi = as.numeric(chain[iter, grep("phi", names(chain))])
    predvec = as.numeric(forecast(auto.arima(phi), h = dim(finaldata)[1])$mean)
    predvec = predvec + betaP[1] * (ww_mean + as.matrix(finalX) %*% betaL) +
      betaP[2] * (ww_mean + as.matrix(finalX) %*% betaL)^2 +
      beta0 + log(finaldata$offset)
    predmat[, iter] = exp(predvec)
  }

  return(list(actual = finaldata,
              prediction = predmat))
}

