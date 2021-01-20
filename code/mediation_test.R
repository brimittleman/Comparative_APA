#' @title Mediation analysis
#'
#' @description Test for the effect of a single mediator, determining whether the reduction in the
#'     indepenent variable is significant, after including the mediating variable in the model. 
#'     In other words, we would like to evaluate the indirect effect of the independent variable 
#'     on the dependent variable, if the indirect effect is sufficiently large, we can say 
#'     that a large amount of the effect of indepenent variable on the depenetn variable is 
#'     transmitted through the mediating variable, and thus a significant mediating effect. 
#'     This function computes the indirect effect.  Some important details about this analysis.
#'     \itemize{
#'       \item We are interested in the effect sizes (regression coefficient) and hence no need to 
#'            concern about common issues in sequencing data anaysis, such as overdispersion.
#'       \item For significance testing, we recommend bootstrapping. In contrast to the popular
#'             Sobel's test approach, boostrapping doesn't require normality assumption.
#'     }
#'     
#' @param exprs gene by sample expression matrix (G by N)
#' @param fixed_covariate length-N vector of sample condition labels (numeric vector,
#'   such as 0, 1 or 1,2).
#' @param varing_covariate gene by sample covariate matrix (G by N)
#'
#' @return a list of estimates
#'   d: estimated indirect effect
#' @export
#'
test_mediation <- function(exprs, fixed_covariates=list(), varying_covariate) {
  
  library(limma)
  
  # Fitting model 1
  #  design_1 <- model.matrix(~fixed_covariates[[1]])
  design_1 <- model.matrix(~fixed_covariates[[1]]+fixed_covariates[[2]])
  colnames(design_1)[-1] <- c(names(fixed_covariates))
  model_1 <- lmFit(exprs, design_1)
  
  # Fitting model 2
  design_2 <- model.matrix(~fixed_covariates[[1]])
  colnames(design_2)[-1] <- names(fixed_covariates)[1]
  model_2 <- lmFit(varying_covariate, design_2)
  
  # Fitting model 3
  model_3 <- lm_varying_covariate(exprs,
                                  fixed_covariates,
                                  varying_covariate)
  
  tau <- coef(model_1)[,2]
  tau_prime <- model_3$coefs[,2]
  alpha <- coef(model_2)[,2]
  beta <- model_3$coefs[,3]
  sigma2_alpha <- (model_2$stdev.unscaled[,2]*model_2$sigma)^2
  sigma2_beta <- (model_3$stdev.unscaled[,3]*model_3$sigma)^2
  
  d <- tau-tau_prime
  se <- sqrt((alpha^2)*sigma2_beta + (beta^2)*sigma2_alpha)
  return(list(tau=tau,
              tau_prime=tau_prime,
              alpha=alpha,
              beta=beta,
              sigma2_alpha=sigma2_alpha,
              sigma2_beta=sigma2_beta,
              d=d,
              se=se))
}
#################
#' @title linear model for gene expression datasets with gene-specific covariates
#'
#' @param exprs sample by gene expressin matrix.
#' @param fixed_covariate length-N covariate vector constant across genes (such as phenotype).
#' @param varying_covariate sample by gene covariate measurement matrix.
#'
#' @export
lm_varying_covariate <- function(exprs, fixed_covariates=list(), varying_covariate) {
  
  assertthat::assert_that(all.equal(dim(exprs), dim(varying_covariate)))
  
  G <- dim(exprs)[1]
  
  # assume the number of fixed covariate is 1
  est <- lapply(1:G, function(g) {
    cov <- unlist(varying_covariate[g,])
    y <- unlist(exprs[g,])
    #    design <- model.matrix(~fixed_covariates[[1]]+fixed_covariates[[2]] + cov)
    #    colnames(design)[c(2,3)] <- names(fixed_covariates)
    design <- model.matrix(~fixed_covariates[[1]] + fixed_covariates[[2]] + cov)
    colnames(design)[2:3] <- c(names(fixed_covariates))
    
    fit <- lm.fit(y=y, x=design)
    return(fit)
  })
  
  coefs <- do.call(rbind, lapply(est, "coef"))
  rownames(coefs) <- rownames(exprs)
  
  cov.coefficient <- lapply(est, function(x) {
    chol2inv(x$qr$qr, size = x$qr$rank) })
  
  stdev.unscaled <- do.call(rbind, lapply(cov.coefficient, function(x) {
    sqrt(diag(x)) }))
  colnames(stdev.unscaled) <- colnames(coefs)
  
  sigma <- sapply(est, function(x) {
    sqrt(sum(x$residuals^2)/x$df.residual)
  })
  return(list(coefs=coefs,
              stdev.unscaled=stdev.unscaled,
              sigma=sigma))
}
