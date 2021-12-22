##############################################
## R package functions for SFM-SHS and SFM-DHS
##############################################

#' Optimization parameter estimation for differential network analysis
#'
#' Function used to implement the hierarchical Bayesian factor model defined in the
#' "Single-cell Differential Network Analysis with Sparse Bayesian Factor Models"
#' manuscript.
#'
#'
#' @param Y data.frame or matrix of gene expression counts where the rows correspond to genes and
#'          columns correspond to cells; Y must contain integers and have row names
#' @param treatGroup vector or factor of length ncol(Y) indicating the
#'        treatment assignments (control or treatment) of the cells.
#' @param Fac number of factors to consider in the model; only a single number is accepted
#' @param mod.type type of model to implement; choose one of either "SHS" for Single HorseShoe or
#'                 "DHS" for Double HorseShoe
#' @param boot.iter total number of bootstrap iterations
#' @param ori.data if TRUE, the first iteration will consist of the original (non-sampled) data; to be used
#'                 when interested in non-bootstrapped results or in the first optimization when bootstrapping
#' @param seed seed for random number generation
#' @param stan.iter maximum number of optimization iterations
#' @param stan.history.size number of update vectors to use in Hessian approximations; for more details refer to
#'                          the Stan manual
#' @param ... additional parameters to be passed to optimize() function in rstan
#'
#' @return sfm.fit-class object containing:
#'
#' \itemize{
#'   \item{Y}{: data.frame or matrix of gene expression counts}
#'   \item{Fac}{: number of factors considered in the model}
#'   \item{corr.est.0}{: array of estimated correlation matrices for the control group}
#'   \item{corr.est.1}{: array of estimated correlation matrices for the treatment group}
#'   \item{corr.est.d}{: array of estimated correlation difference (control - treatment) matrices}
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
#' ## For first set of bootstrap optimizations set ori.data = TRUE to optimize on original data in first iteration
#' res1 = sfm(Y = diff.gene.mat, treatGroup = t_i, Fac = 7, mod.type = "SHS", boot.iter = 5, ori.data = TRUE, seed = 123,
#' stan.iter = 100000, stan.history.size = 10)
#'
#' }
#'
#' @export
#'

sfm = function(Y, treatGroup, Fac, mod.type = "SHS", boot.iter = 1, ori.data = FALSE, seed = 123,
               stan.iter = 100000, stan.history.size = 10, ...){

  ###### Last updated 11/01/2021

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
  G = nrow(Y)
  N = ncol(Y)
  treatGroup = as.factor(treatGroup)
  t_i = unname(model.matrix(~ treatGroup)[,2])
  treat0 = which(t_i == 0)
  treat1 = which(t_i == 1)

  ## Begin bootstrap
  addArgs = list(...)
  corr_est0_mat = array(NA, dim = c(G,G,boot.iter))
  corr_est1_mat = array(NA, dim = c(G,G,boot.iter))
  for(i in 1:boot.iter){
    set.seed(seed)
    if(ori.data & i == 1){
      ## Run optimize on original data
      Y_boot = Y
      t_i_boot = t_i

    }else{
      ## Randomly sample from within each treatment
      boot_cols = c(sample(treat0, length(treat0), replace=TRUE),sample(treat1, length(treat1), replace=TRUE))
      Y_boot = Y[,boot_cols]
      t_i_boot = t_i[boot_cols]
    }

    gene_data <- list(
      G = nrow(Y_boot),
      N = ncol(Y_boot),
      Fac = Fac,
      Y = Y_boot,
      t_i = t_i_boot)

    ## Run optimizing function from rstan
    fit1 = do.call(rstan::optimizing,c(list(
      object = mod,
      data = gene_data,    # named list of data
      seed = seed+i,
      iter = stan.iter,
      history_size = stan.history.size),addArgs)
    )

    ind_0 = grep("corr_mu_mu0",names(fit1$par))
    ind_1 = grep("corr_mu_mu1",names(fit1$par))

    corr_samps0 = matrix(fit1$par[ind_0], nrow = G, ncol = G)
    corr_samps0[which(is.infinite(corr_samps0))] = 0
    corr_est0_mat[,,i] = corr_samps0

    corr_samps1 = matrix(fit1$par[ind_1], nrow = G, ncol = G)
    corr_samps1[which(is.infinite(corr_samps1))] = 0
    corr_est1_mat[,,i] = corr_samps1

    message(paste("Optimization",i,"completed."))

  }# End bootstrap loop

  corr_diff_mat = corr_est0_mat - corr_est1_mat

  res = list(Y=Y, Fac=Fac, corr.est.0 = corr_est0_mat, corr.est.1 = corr_est1_mat, corr.est.d = corr_diff_mat)
  class(res) <- "sfm.fit"
  return(res)
}





#' Estimate the gene-gene correlation matrices from bootstrapped samples
#'
#' Function used to combine and analyze bootstrapped results from one or more sets of
#' correlation samples generated by the \code{sfm} function.
#'
#' @param sfm.list list where each element contains an sfm.fit-class object; each element of the list
#'                  contains an object from a different parallel bootstrap chain
#'@param ori.samp location of list element that contains optimization on non-sampled data; default is 1
#' @return sfm.corr-class object containing:
#'
#'  \itemize{
#'   \item{corr.d}{: symmetric estimated differential correlation matrix generated from optimization}
#'   \item{CI.low.d}{: symmetric matrix representing the lower bound of the 95% confidence intervals (CI) generated from the bootstrap optimizations for differential network analysis}
#'   \item{CI.upp.d}{: symmetric matrix representing the upper bound of the 95% confidence intervals (CI) generated from the bootstrap optimizations for differential network analysis}
#'   \item{CI.eval.d}{: symmetric logical matrix generated from 95% confidence intervals (CI) from bootstrap optimizations for differential network analysis;
#'                    TRUE means 95% CI does not include 0 (i.e., significant correlation)}
#'   \item{p.val.d}{: symmetric matrix consisting of approximate "p-values" generated from bootstrap optimizations for differential network analysis}
#'   \item{corr.0}{: symmetric estimated correlation matrix generated from optimization for control group}
#'   \item{CI.low.0}{: symmetric matrix representing the lower bound of the 95% confidence intervals (CI) generated from the bootstrap optimizations for control group}
#'   \item{CI.upp.0}{: symmetric matrix representing the upper bound of the 95% confidence intervals (CI) generated from the bootstrap optimizations for control group}
#'   \item{CI.eval.0}{: symmetric logical matrix generated from 95% confidence intervals (CI) from bootstrap optimizations for control group;
#'                      TRUE means 95% CI does not include 0 (i.e., significant correlation)}
#'   \item{p.val.0}{: symmetric matrix consisting of approximate "p-values" generated from bootstrap optimizations for control group}
#'   \item{corr.1}{: symmetric estimated correlation matrix generated from optimization for treatment group}
#'   \item{CI.low.1}{: symmetric matrix representing the lower bound of the 95% confidence intervals (CI) generated from the bootstrap optimizations for treatment group}
#'   \item{CI.upp.1}{: symmetric matrix representing the upper bound of the 95% confidence intervals (CI) generated from the bootstrap optimizations for treatment group}
#'   \item{CI.eval.1}{: symmetric logical matrix generated from 95% confidence intervals (CI) from bootstrap optimizations for treatment group;
#'                      TRUE means 95% CI does not include 0 (i.e., significant correlation)}
#'   \item{p.val.1}{: symmetric matrix consisting of approximate "p-values" generated from bootstrap optimizations for treatment group}
#' }
#'
#' @details
#' The correlation estimates are determined from original optimization on original (non-sampled) data for each gene-gene pair.
#'
#' To determine whether the correlation is significant, a 95% confidence interval (CI) is determined from the bootstrap samples. If the CI
#' includes 0, the correlation is deemed to be non-significant and the "CI.eval" element is FALSE. If the CI does not include 0, the correlation is
#' considered significant and the "CI.eval" element is TRUE.
#'
#' An alternative measurement of significance is the approximate "p-value" calculation, which is output in the "p.val" matrix. This "p-value" is determined
#' by finding the smallest "a" value such that the 100(1-a)% CI contains 0. The corresponding "a" value represents the proportion of posterior distribution
#' that is outside the smallest confidence interval that contains 0.
#'
#' In most cases a significant 95% CI from "CI.eval" will correspond with an approximate "p-value" < 0.05 from "p.val".
#'
#' @examples
#' \dontrun{
#' ## Load dataset
#' data("diff.gene.mat")
#'
#' ## First 100 columns are the control group, the next 100 columns are the treatment group
#' t_i = rep(0:1,each=100)
#'
#' ## For first set of bootstrap optimizations set ori.data = TRUE to optimize on original data in first iteration
#' res1 = sfm(Y = diff.gene.mat, treatGroup = t_i, Fac = 7, mod.type = "SHS", boot.iter = 5, ori.data = TRUE, seed = 123,
#' stan.iter = 100000, stan.history.size = 10)
#'
#' ## For other sets of bootstrap optimizations set ori.data = FALSE to optimize on resampled data
#' res2 = sfm(Y = diff.gene.mat, treatGroup = t_i, Fac = 7, mod.type = "SHS", boot.iter = 5, ori.data = FALSE, seed = 1234,
#' stan.iter = 100000, stan.history.size = 10)
#'
#' ## Combine bootstrap results from multiple sfm.fit-class objects
#' ## Note: The first object has the optimization estimates from original data hence ori.samp = 1
#' res.list = list(res1,res2)
#' boot.res = sfm.results(res.list, ori.samp=1)
#'
#' ## Output differential network results
#' summary(boot.res)
#'
#' ## Output control network correlation
#' controlCorr(boot.res)
#'
#' ## Output treatment network correlation
#' treatCorr(boot.res)
#'
#' }
#'
#' @export

sfm.results <- function(sfm.list,ori.samp=1){
  G <- nrow(sfm.list[[1]]$Y)
  gene_names <- rownames(sfm.list[[1]]$Y)

  ## Identify the total number of samples
  M.samps <- numeric()
  for(l in 1:length(sfm.list)){
    M.samps <- c(M.samps, dim(sfm.list[[l]]$corr.est.d)[3])
  }
  corr_mat_0_samps <- array(NA, dim = c(G,G,sum(M.samps)))
  corr_mat_1_samps <- array(NA, dim = c(G,G,sum(M.samps)))
  corr_mat_d_samps <- array(NA, dim = c(G,G,sum(M.samps)))

  M.start <- 1
  for(l in 1:length(sfm.list)){
    if(l == ori.samp){
      corr_mat_0_est <- sfm.list[[l]]$corr.est.0[,,1]
      corr_mat_1_est <- sfm.list[[l]]$corr.est.1[,,1]
      corr_mat_d_est <- sfm.list[[l]]$corr.est.d[,,1]
    }
    corr_mat_0_samps[,,M.start:(M.start - 1 + M.samps[l])] <- sfm.list[[l]]$corr.est.0
    corr_mat_1_samps[,,M.start:(M.start - 1 + M.samps[l])] <- sfm.list[[l]]$corr.est.1
    corr_mat_d_samps[,,M.start:(M.start - 1 + M.samps[l])] <- sfm.list[[l]]$corr.est.d
    M.start <- M.start + M.samps[l]
  }

  CI_0_low_all <- apply(corr_mat_0_samps,1:2,quantile, probs = c(0.025), na.rm=TRUE)
  CI_0_up_all <- apply(corr_mat_0_samps,1:2,quantile, probs = c(0.975), na.rm=TRUE)
  CI_0_eval_all <- (CI_0_low_all * CI_0_up_all) > 0
  p_val_0 <- apply(corr_mat_0_samps,1:2, p.approx)
  diag(p_val_0) <- diag(CI_0_eval_all) <-  NA

  CI_1_low_all <- apply(corr_mat_1_samps,1:2,quantile, probs = c(0.025), na.rm=TRUE)
  CI_1_up_all <- apply(corr_mat_1_samps,1:2,quantile, probs = c(0.975), na.rm=TRUE)
  CI_1_eval_all <- (CI_1_low_all * CI_0_up_all) > 0
  p_val_1 <- apply(corr_mat_1_samps,1:2, p.approx)
  diag(p_val_1) <- diag(CI_1_eval_all) <-  NA

  CI_d_low_all <- apply(corr_mat_d_samps,1:2,quantile, probs = c(0.025), na.rm=TRUE)
  CI_d_up_all <- apply(corr_mat_d_samps,1:2,quantile, probs = c(0.975), na.rm=TRUE)
  CI_d_eval_all <- (CI_d_low_all * CI_d_up_all) > 0
  p_val_d <- apply(corr_mat_d_samps,1:2, p.approx)
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

  res <- list(corr.0 = corr_mat_0_est, CI.low.0 = CI_0_low_all, CI.upp.0 = CI_0_up_all, CI.eval.0 = CI_0_eval_all, p.val.0 = p_val_0,
              corr.1 = corr_mat_1_est, CI.low.1 = CI_1_low_all, CI.upp.1 = CI_1_up_all, CI.eval.1 = CI_1_eval_all, p.val.1 = p_val_1,
              corr.d = corr_mat_d_est, CI.low.d = CI_d_low_all, CI.upp.d = CI_d_up_all, CI.eval.d = CI_d_eval_all, p.val.d = p_val_d)
  class(res) <- "sfm.corr"
  return(res)
}





#####################################
## S3 methods
#####################################

#' Prints sfm.corr object
#' @export
#' @noRd
`print.sfm.corr` <-
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
`summary.sfm.corr` <-
  function(object, ...){
    print(object)
  }

#' Prints control network results from sfm.corr object
#' @export
#' @noRd
`controlCorr` <-
  function(sfm.corr){
    x.frame <- data.frame(gene.name1 = character(),
                          gene.name2 = character(),
                          corr0 = double(),
                          CI.low0 = double(),
                          CI.upp0 = double(),
                          CI.sig0 = logical(),
                          p.val0 = character(),
                          stringsAsFactors = FALSE)
    k<-1
    for(i in 1:(nrow(sfm.corr$corr.0)-1)){
      for(j in (i+1):ncol(sfm.corr$corr.0)){
        x.frame[k,1] <- rownames(sfm.corr$corr.0)[i]
        x.frame[k,2] <- colnames(sfm.corr$corr.0)[j]
        x.frame[k,3] <- round(sfm.corr$corr.0[i,j],4)
        x.frame[k,4] <- round(sfm.corr$CI.low.0[i,j],4)
        x.frame[k,5] <- round(sfm.corr$CI.upp.0[i,j],4)
        x.frame[k,6] <- ifelse(sfm.corr$CI.eval.0[i,j],"  *  ","")
        x.frame[k,7] <- formatC(sfm.corr$p.val.0[i,j],format = "e", digits = 2)
        k<-k+1
      }
    }
    print(x.frame)
  }


#' Prints treatment network results from sfm.corr object
#' @export
#' @noRd
`treatCorr` <-
  function(sfm.corr){
    x.frame <- data.frame(gene.name1 = character(),
                          gene.name2 = character(),
                          corr1 = double(),
                          CI.low1 = double(),
                          CI.upp1 = double(),
                          CI.sig1 = logical(),
                          p.val1 = character(),
                          stringsAsFactors = FALSE)
    k<-1
    for(i in 1:(nrow(sfm.corr$corr.1)-1)){
      for(j in (i+1):ncol(sfm.corr$corr.1)){
        x.frame[k,1] <- rownames(sfm.corr$corr.1)[i]
        x.frame[k,2] <- colnames(sfm.corr$corr.1)[j]
        x.frame[k,3] <- round(sfm.corr$corr.1[i,j],4)
        x.frame[k,4] <- round(sfm.corr$CI.low.1[i,j],4)
        x.frame[k,5] <- round(sfm.corr$CI.upp.1[i,j],4)
        x.frame[k,6] <- ifelse(sfm.corr$CI.eval.1[i,j],"  *  ","")
        x.frame[k,7] <- formatC(sfm.corr$p.val.1[i,j],format = "e", digits = 2)
        k<-k+1
      }
    }
    print(x.frame)
  }





#####################################
## Internal calculation functions
#####################################
#' Approximate p-values from correlation posterior
#'
#' Find the smallest "a" such that the 100(1-a)% confidence interval
#' contains 0.
#'
#' @param x samples from posterior
#'
#' @return
#' approximate p-value
#' @export
#' @noRd


p.approx <- function(x){
  if(median(x) < 0){
    n <- length(which(x > 0))
  } else if(median(x) > 0){
    n <- length(which(x < 0))
  } else{
    # median is 0
    n <- sum(abs(x) > 0)/2
  }

  # Check for 0 on boundaries
  if(n == 0 & min(abs(x)) == 0){
    n <- sum(abs(x) == 0)
  }

  if(n == 0){
    p.val1 <- .01/(length(x)+2)
  }else{
    p.val1 <- n/(length(x))
  }
  p.val <- 2*p.val1
  return(p.val)
}





#####################################
## Datasets
#####################################

#' Example gene expression matrix
#'
#' An example simulated dataset consisting of gene expression counts for G = 50 genes (rows)
#' across N = 200 cells (columns). The first 100 columns correspond to the control group
#' and the remaining 100 columns correspond to the treatment group
#'
#' @format A matrix with 50 rows and 200 columns
#'
#' @examples
#' data(diff.gene.mat)
#' rownames(diff.gene.mat)
#' colnames(diff.gene.mat)
#'
"diff.gene.mat"

#' Mouse microglia cell (MMC) scRNA-seq data
#'
#' Subset of data obtained from the GEO database under accession number GSE90975 and contains the gene expressions from single-cell
#' analysis of neurodegeneration in microglia cells of mice (Tay et al., 2018). This data includes N = 944 cells G = 98 differentially
#' expressed genes taken from Figure S1 of the Tay et al. manuscript.
#'
#' @format A matrix with 98 rows (genes) and 944 columns (cells)
#'
#' @references
#' Tay, T. L., Dautzenberg, J., Grün, D., & Prinz, M. (2018). Unique microglia recovery population revealed by single-cell RNA-seq following
#' neurodegeneration. Acta neuropathologica communications, 6(1), 1-11.
#'
#'
#' @examples
#' data(mmc.gene.mat)
#' rownames(mmc.gene.mat)
#' colnames(mmc.gene.mat)
#'
"mmc.gene.mat"
