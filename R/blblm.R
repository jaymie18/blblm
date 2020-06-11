#' @import purrr
#' @import stats
#' @import parallel
#' @import bench
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Linear Regression with Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("fit"))

#' Blb lm
#'
#' linear regression that creates subsample and bootstrap
#'
#'@param formula the interest user want from data
#'@param data the dataset
#'@param m numbers of subsamples created
#'@param B the number of bootstraps user want
#' @export
blblm <- function(formula, data, m = 10, B = 5000) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' blb parallel
#'
#' It allows user to do bag of little bootstrap linear regression applying cluster to speed the process
#'
#' @param formula the interest user want from data
#' @param data the dataset
#' @param m numbers of subsamples created
#' @param B the number of bootstraps user want
#' @param cl the cluster user create
#' @export
blblm_par <- function(formula, data, m, B, cl) {
  data_list <- split_data(data, m)
  n <- nrow(data)
  estimates <- parLapply(cl, data_list, function(data, formula,n, B){
    blblm::lm_each_subsample(
      formula = formula, data = data , n = n, B = B)
  }, formula = formula, n = n , B = B)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  res
}

#' split data
#'
#' @param m numbers of subsamples created
#' @param data the dataset
#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' create subsample
#'
#' compute the estimates
#'
#' @param formula the interest user want from data
#' @param data the dataset
#' @param n number of subsample
#' @param B the number of bootstraps user want
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' applied regression
#'
#' compute the regression estimates for a blb dataset
#' @param formula the interest user want from data
#' @param data the dataset
#' @param n numbers of subsamples created
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' linear regression with rcpp
#'
#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula the interest user want from data
#' @param data the dataset
#' @param freqs the times user want the regression weight on
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  x <- model.matrix(formula , data) #transform independent variable
  y <- model.response(model.frame(formula, data)) #transform dependent variable
  fit <- fastlm(x, y, wi = freqs)
  fit$coefficients <- unlist(as.list(t(fit$coefficients))) #rename the coefficients
  names(fit$coefficients) <- c(colnames(x))
  list(coef = blbcoef(fit), sigma = fit$s) #list out the information
}


#' coefficient of blblm
#'
#' compute the coefficients from fit
#'
#' @param fit the linear regression function
#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}

#' sigma of blblm
#'
#' compute sigma from fit
#'
#' @param fit the linear regression function
blbsigma <- function(fit) {
  sigma(fit)
}

#' print blblm model
#'
#' @param x blblm model
#' @param ... extra arguments other than this function
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' compute the condifence interval for sigma
#'
#' @param object the variable or parameter
#' @param level level of confidence interval
#' @param confidence if user wants confidence interval
#' @param ... extra arguments other than this function
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' blblm coef calculation
#'
#' @param object the variable or parameter
#' @param ... extra arguments other than this function
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' confint interval
#'
#' a function that computes the bootstrap confidence interval
#'
#' @param object the regression model
#' @param level level of confidence
#' @param parm the input of parameter
#' @param ... extra arguments other than this function
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' prediction using blblm
#'
#' it predicts new data with the use of blblm
#'
#' @param object the variable or parameter
#' @param new_data data of interest
#' @param level the percentage of confidence level
#' @param confidence pick if they want the confidence level
#' @param ... extra arguments that is not related to this function
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}



mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}


map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}



map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}


map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
