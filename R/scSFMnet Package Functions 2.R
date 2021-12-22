##############################################
## HMC functions for SFM-SHS and SFM-DHS
##############################################

#' HMC parameter estimation wrapper function for differential network analysis
#'
#' Wrapper function used to implement the Hamiltonian Monte Carlo (HMC) version of the hierarchical
#' Bayesian factor model defined in the "Single-cell Differential Network Analysis with Sparse Bayesian
#' Factor Models" manuscript.
#'
#'
#' @param Y data.frame or matrix of gene expression counts where the rows correspond to genes and
#'          columns correspond to cells; Y must contain integers and have row names
#' @param treatGroup vector or factor of length ncol(Y) indicating the
#'        treatment assignments (control or treatment) of the cells.
#' @param Fac number of factors to consider in the model; only a single number is accepted
#' @param mod.type type of model to implement; choose one of either "SHS" for Single HorseShoe or
#'                 "DHS" for Double HorseShoe
#' @param chains a positive integer specifying the number of Markov chains. The rstan default is 4.
#' @param n.cores number of cores to use when executing the stan chains in parallel
#' @param seed seed for random number generation in rstan
#' @param ... additional parameters to be passed to sampling() function in rstan
#'
#' @return sfm.stan-class object containing:
#'
#' \itemize{
#'   \item{Y}{: data.frame or matrix of gene expression counts}
#'   \item{Fac}{: number of factors considered in the model}
#'   \item{stan.fit}{: an object of S4 class stanfit representing the fitted results}
#'   \item{corr.d}{: symmetric estimated differential correlation matrix generated from optimization}
#'   \item{CI.low.d}{: symmetric matrix representing the lower bound of the 95% credible intervals (CrI) generated from HMC}
#'   \item{CI.upp.d}{: symmetric matrix representing the upper bound of the 95% credible intervals (CrI) generated from HMC}
#'   \item{CI.eval.d}{: symmetric logical matrix generated from 95% credible intervals (CrI) from HMC;
#'                    TRUE means 95% CrI does not include 0 (i.e., significant correlation)}
#'   \item{p.val.d}{: symmetric matrix consisting of approximate "p-values" generated from HMC for differential network analysis}
#'   \item{corr.0}{: symmetric estimated correlation matrix generated from optimization for control group}
#'   \item{CI.low.0}{: symmetric matrix representing the lower bound of the 95% credible intervals (CrI) generated from HMC}
#'   \item{CI.upp.0}{: symmetric matrix representing the upper bound of the 95% credible intervals (CrI) generated from HMC}
#'   \item{CI.eval.0}{: symmetric logical matrix generated from 95% credible intervals (CrI) from HMC;
#'                      TRUE means 95% CrI does not include 0 (i.e., significant correlation)}
#'   \item{p.val.0}{: symmetric matrix consisting of approximate "p-values" generated from HMC for control group}
#'   \item{corr.1}{: symmetric estimated correlation matrix generated from optimization for treatment group}
#'   \item{CI.low.1}{: symmetric matrix representing the lower bound of the 95% credible intervals (CrI) generated from HMC}
#'   \item{CI.upp.1}{: symmetric matrix representing the upper bound of the 95% credible intervals (CrI) generated from HMC}
#'   \item{CI.eval.1}{: symmetric logical matrix generated from 95% credible intervals (CrI) from HMC;
#'                      TRUE means 95% CrI does not include 0 (i.e., significant correlation)}
#'   \item{p.val.1}{: symmetric matrix consisting of approximate "p-values" generated from HMC for treatment group}
#' }
#'
#' @examples
#' \dontrun{
#' ## Load dataset
#' data("diff.gene.mat")
#'
#' ## First 100 columns are the control group, the next 100 columns are the treatment group
#' t_i = rep(0:1,each=100)
#'
#' ## Run HMC via stan
#' res1 = sfm.hmc(Y = diff.gene.mat, treatGroup = t_i, Fac = 7, mod.type = "SHS", seed = 123)
#'
#' ## Output differential network results
#' summary(res1)
#'
#' ## Output control network correlation
#' controlCorr(res1)
#'
#' ## Output treatment network correlation
#' treatCorr(res1)
#'
#' }
#'
#' @export
#'

sfm.hmc = function(Y, treatGroup, Fac, mod.type = "SHS", chains = 4, n.cores = 1, seed = 123, ...){

  ###### Last updated 12/20/2021

  ## Ensure mod.type is in correct format
  if(length(mod.type)>1){
    stop("mod.type must either be 'SHS' or 'DHS'")
  }else if(!mod.type %in% c("SHS","DHS")){
    stop("mod.type must either be 'SHS' or 'DHS'")
  }

  ## Get Stan Model
  if(mod.type == "SHS"){
    mod <- stanmodels$sfm_shs
  }else if(mod.type == "DHS"){
    mod <- stanmodels$sfm_dhs
  }

  ## Check for count data
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  isCountData <- all(is.wholenumber(Y))
  if(identical(isCountData, FALSE)){
    stop("Y must consist of only integers")
  }

  ## Check for number of groups
  if(length(unique(treatGroup)) != 2){
    stop("treatGroup must consist of only two groups")
  }

  ## Obtain information from data
  gene_names = rownames(Y)
  G = nrow(Y)
  N = ncol(Y)
  treatGroup = as.factor(treatGroup)
  t_i = unname(model.matrix(~ treatGroup)[,2])
  treat0 = which(t_i == 0)
  treat1 = which(t_i == 1)

  ## Begin HMC
  addArgs = list(...)

  gene_data <- list(
    G = nrow(Y),
    N = ncol(Y),
    Fac = Fac,
    Y = Y,
    t_i = t_i)

  ## Run optimizing function from rstan
  fit1 = do.call(rstan::sampling,c(list(
    object = mod,
    data = gene_data,    # named list of data
    chains = chains,
    seed = seed,
    cores = n.cores), addArgs)
  )

  fit_dat <- rstan::extract(fit1)

  corr_mat_0_samps <- fit_dat$corr_mu_mu0
  corr_mat_1_samps <- fit_dat$corr_mu_mu1
  corr_mat_d_samps <- corr_mat_0_samps - corr_mat_1_samps

  corr_mat_0_est <- apply(corr_mat_0_samps,2:3,mean)
  CI_0_low_all <- apply(corr_mat_0_samps,2:3,quantile, probs = c(0.025), na.rm=TRUE)
  CI_0_up_all <- apply(corr_mat_0_samps,2:3,quantile, probs = c(0.975), na.rm=TRUE)
  CI_0_eval_all <- (CI_0_low_all * CI_0_up_all) > 0
  p_val_0 <- apply(corr_mat_0_samps,2:3, p.approx)
  diag(p_val_0) <- diag(CI_0_eval_all) <-  NA

  corr_mat_1_est <- apply(corr_mat_1_samps,2:3,mean)
  CI_1_low_all <- apply(corr_mat_1_samps,2:3,quantile, probs = c(0.025), na.rm=TRUE)
  CI_1_up_all <- apply(corr_mat_1_samps,2:3,quantile, probs = c(0.975), na.rm=TRUE)
  CI_1_eval_all <- (CI_1_low_all * CI_0_up_all) > 0
  p_val_1 <- apply(corr_mat_1_samps,2:3, p.approx)
  diag(p_val_1) <- diag(CI_1_eval_all) <-  NA

  corr_mat_d_est <- apply(corr_mat_d_samps,2:3,mean)
  CI_d_low_all <- apply(corr_mat_d_samps,2:3,quantile, probs = c(0.025), na.rm=TRUE)
  CI_d_up_all <- apply(corr_mat_d_samps,2:3,quantile, probs = c(0.975), na.rm=TRUE)
  CI_d_eval_all <- (CI_d_low_all * CI_d_up_all) > 0
  p_val_d <- apply(corr_mat_d_samps,2:3, p.approx)
  diag(p_val_d) <- diag(CI_d_eval_all) <-  NA

  rownames(corr_mat_0_est) <- colnames(corr_mat_0_est) <- gene_names
  rownames(corr_mat_1_est) <- colnames(corr_mat_1_est) <- gene_names
  rownames(corr_mat_d_est) <- colnames(corr_mat_d_est) <- gene_names
  rownames(CI_0_eval_all) <- colnames(CI_0_eval_all) <- gene_names
  rownames(CI_0_low_all) <- colnames(CI_0_low_all) <- gene_names
  rownames(CI_0_up_all) <- colnames(CI_0_up_all) <- gene_names
  rownames(p_val_0) <- colnames(p_val_0) <- gene_names
  rownames(CI_1_eval_all) <- colnames(CI_1_eval_all) <- gene_names
  rownames(CI_1_low_all) <- colnames(CI_1_low_all) <- gene_names
  rownames(CI_1_up_all) <- colnames(CI_1_up_all) <- gene_names
  rownames(p_val_1) <- colnames(p_val_1) <- gene_names
  rownames(CI_d_eval_all) <- colnames(CI_d_eval_all) <- gene_names
  rownames(CI_d_low_all) <- colnames(CI_d_low_all) <- gene_names
  rownames(CI_d_up_all) <- colnames(CI_d_up_all) <- gene_names
  rownames(p_val_d) <- colnames(p_val_d) <- gene_names

  res <- list(Y=Y, Fac=Fac, stan.fit = fit1,
              corr.0 = corr_mat_0_est, CI.low.0 = CI_0_low_all, CI.upp.0 = CI_0_up_all, CI.eval.0 = CI_0_eval_all, p.val.0 = p_val_0,
              corr.1 = corr_mat_1_est, CI.low.1 = CI_1_low_all, CI.upp.1 = CI_1_up_all, CI.eval.1 = CI_1_eval_all, p.val.1 = p_val_1,
              corr.d = corr_mat_d_est, CI.low.d = CI_d_low_all, CI.upp.d = CI_d_up_all, CI.eval.d = CI_d_eval_all, p.val.d = p_val_d)
  class(res) <- "sfm.stan"
  return(res)
}






#####################################
## S3 methods - sfm.stan
#####################################

#' Prints sfm.corr object
#' @export
#' @noRd
`print.sfm.stan` <-
  function(x, ...){
    x.frame <- data.frame(gene.name1 = character(),
                          gene.name2 = character(),
                          corrD = double(),
                          CI.lowD = double(),
                          CI.uppD = double(),
                          CI.sigD = logical(),
                          p.valD = character(),
                          stringsAsFactors = FALSE)
    k<-1
    for(i in 1:(nrow(x$corr.d)-1)){
      for(j in (i+1):ncol(x$corr.d)){
        x.frame[k,1] <- rownames(x$corr.d)[i]
        x.frame[k,2] <- colnames(x$corr.d)[j]
        x.frame[k,3] <- round(x$corr.d[i,j],4)
        x.frame[k,4] <- round(x$CI.low.d[i,j],4)
        x.frame[k,5] <- round(x$CI.upp.d[i,j],4)
        x.frame[k,6] <- ifelse(x$CI.eval.d[i,j],"  *  ","")
        x.frame[k,7] <- formatC(x$p.val.d[i,j],format = "e", digits = 2)
        k<-k+1
      }
    }
    print(x.frame)
  }

#' Summary for sfm.corr object
#' @export
#' @noRd
`summary.sfm.stan` <-
  function(object, ...){
    print(object)
  }

#' Prints control network results from sfm.stan object
#' @export
#' @noRd
`controlCorr` <-
  function(sfm.stan){
    x.frame <- data.frame(gene.name1 = character(),
                          gene.name2 = character(),
                          corr0 = double(),
                          CI.low0 = double(),
                          CI.upp0 = double(),
                          CI.sig0 = logical(),
                          p.val0 = character(),
                          stringsAsFactors = FALSE)
    k<-1
    for(i in 1:(nrow(sfm.stan$corr.0)-1)){
      for(j in (i+1):ncol(sfm.stan$corr.0)){
        x.frame[k,1] <- rownames(sfm.stan$corr.0)[i]
        x.frame[k,2] <- colnames(sfm.stan$corr.0)[j]
        x.frame[k,3] <- round(sfm.stan$corr.0[i,j],4)
        x.frame[k,4] <- round(sfm.stan$CI.low.0[i,j],4)
        x.frame[k,5] <- round(sfm.stan$CI.upp.0[i,j],4)
        x.frame[k,6] <- ifelse(sfm.stan$CI.eval.0[i,j],"  *  ","")
        x.frame[k,7] <- formatC(sfm.stan$p.val.0[i,j],format = "e", digits = 2)
        k<-k+1
      }
    }
    print(x.frame)
  }


#' Prints treatment network results from sfm.stan object
#' @export
#' @noRd
`treatCorr` <-
  function(sfm.stan){
    x.frame <- data.frame(gene.name1 = character(),
                          gene.name2 = character(),
                          corr1 = double(),
                          CI.low1 = double(),
                          CI.upp1 = double(),
                          CI.sig1 = logical(),
                          p.val1 = character(),
                          stringsAsFactors = FALSE)
    k<-1
    for(i in 1:(nrow(sfm.stan$corr.1)-1)){
      for(j in (i+1):ncol(sfm.stan$corr.1)){
        x.frame[k,1] <- rownames(sfm.stan$corr.1)[i]
        x.frame[k,2] <- colnames(sfm.stan$corr.1)[j]
        x.frame[k,3] <- round(sfm.stan$corr.1[i,j],4)
        x.frame[k,4] <- round(sfm.stan$CI.low.1[i,j],4)
        x.frame[k,5] <- round(sfm.stan$CI.upp.1[i,j],4)
        x.frame[k,6] <- ifelse(sfm.stan$CI.eval.1[i,j],"  *  ","")
        x.frame[k,7] <- formatC(sfm.stan$p.val.1[i,j],format = "e", digits = 2)
        k<-k+1
      }
    }
    print(x.frame)
  }
