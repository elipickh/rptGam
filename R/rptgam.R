#' @title Repeatability estimates for random effect GAM models
#'
#' @description Estimate repeatabilties for random effect terms in mgcv GAM Gaussian models
#'
#' @param formula A GAM formula. If NULL (default), then a fitted gam object should be
#' passed to gamObj. Must contain at least one random term, and none of the random terms
#' can be specified with interactions with other terms.
#' @param data A data frame or list containing the model response variable and covariates
#' required by the formula. If NULL (default), then a fitted gam object should be passed
#' to gamObj.
#' @param gamObj A fitted gam object as produced by mgcv::gam or mgcv::bam. Must contain at
#' least one random term, and none of the random terms can be specified with interactions
#' with other terms. Only Gaussian models are currently supported. If NULL, then formula
#' and data should both be non-NULL.
#' @param rterms Character string or a vector of character strings of random term(s) to
#' include in repeatability estimations. If NULL (default), include all random terms.
#' @param gam_pars List of parameters for mgcv::gam model specification, which will be
#' used when formula and data are both non-NULL. If TRUE, then the default parameter
#' values in mgcv::gam will be used. Should not contain the control(nthreads) parameter,
#' as this will be set by rptgam. Only Gaussian models are currently supported.
#' @param bam_pars List of parameters for mgcv::bam model specification, which will be
#' used when formula and data are both non-NULL. If TRUE, then the default parameter values
#' in mgcv::bam will be used. Should not contain the nthreads (which is used when
#' DISCRETE = TRUE) or the cluster parameters, as they will be set by rptgam. Only Gaussian
#' models are currently supported.
#' @param nboot Number of bootstrap replications per bootstrap type. Default is 0.
#' @param boot_type Type of bootstraps to perform for the repeatability estimates. Can be
#' one or more of the following: 'case', 'param' (default), 'resid', 'cgr', 'reb0', 'reb1',
#' 'reb2'. If 'all', then all methods will be used. Note that types 'reb0', 'reb1', and
#' 'reb2' are available when the GAM model contains exactly one random term.
#' See \strong{Details} below for more information on the various types.
#' @param ci Confidence level for the bootstrap interval(s). Default to 0.95.
#' @param ci_type Type of ci for the boostraps. Can be 'perc', which uses R's quantile
#' method or 'bca' which uses the bias-corrected and accelerated method. 'all' (default)
#' returns both of these types.
#' @param case_resample Vector of booleans. Required when the 'case' bootstrap method is
#' used, and specifies whether each level of the model should be resampled. "The levels
#' should be specified from the highest level (largest cluster) of the hierarchy to the
#' lowest (observation-level); for example for students within a school, specify the school
#' level first, then the student level" (lmeresampler::case_bootstrap). Length of vector
#' should be one more than the number of random terms in the model (the extra, and last,
#' one being the row unit). See \strong{Details} below for more information.
#' @param case_tries A numeric indicating the maximum number of resampling tries for the
#' 'case' bootstrap to get a usable sample (see case_minrows). Default is 100.
#' @param case_minrows A numeric indicating the minimum number of rows allowable in a
#' 'case' type bootstrap. This number will be ignored if it is below the number of rows
#' allowable to run the gam model. If NULL (default), then will use no less than number
#' of rows allowable to run the gam model.
#' @param nperm Number of permutations for calcualting asymptotic p-values for the
#' repeatability estimates. Default is 0.
#' @param aic Boolean (default is TRUE) indicating whether to calculate the AIC(s) of t
#' he models with and without the random effect(s).
#' @param select Boolean (default is TRUE) indicating whether to calcualte coefficients
#' for random effects with and without selection penalties.
#' @param saveboot Boolean (default is TRUE) indicating whether to save the bootstraped
#' repeatability estimates in the returned object.
#' @param savepermute Boolean (default is TRUE) indicating whether to save the permutated
#' repeatability estimates in the returned object.
#' @param seed Numeric (which is converted to integer) to be used in set.seed to allow
#' reproducible reults of bootstraps and permutations. Default is NULL.
#' @param verbose Boolean (default is TRUE) indicating if messages should be printed.
#' @param parallel Boolean (default is TRUE) indicating if the models, bootstraps, and
#' permutations should be run in parallel
#' @param ncores Integer indicating how many cores to use for parallel processing.
#' Positive integers specify the number of cores. -1 means using all processors, -2 means
#' using 1 less than all processores, etc. Default is -1.
#' @param logical Boolean indicating if virtual CPU cores should be counted when ncores
#' is negative (the same as running parallel::detectCores(logical = TRUE)). FALSE (default)
#' only counts physical cores.
#' @param cl_type One of 'PSOCK' (default) or 'FORK', indicating the type of cluster to
#' use for parallel processing. 'FORK' can be faster than 'PSOCK' but can be unstable
#' and isn't available in Windows.
#' @param blas_threads Integer indicating the number of BLAS threads to use. For
#' multi-threaded BLAS libraries, such as MKL, OpenBLAS and Apple BLAS, and when parallel
#' is set to TRUE, is can be faster to limit BLAS threads to 1. This is accomplished
#' via RhpcBLASctl::blas_set_num_threads(1), which will be installed if it is not already.
#' R's default BLAS library is single threaded, and so this parameter isn't neccesary.
#' NULL (default) skips this process.
#' @return Returns an object of class \code{rptgam}.
#'
#' @examples
#' library(mgcv)
#' dat <- gamSim(1,n=100,scale=2)
#' fac <- sample(1:5,100,replace=TRUE)
#' b <- rnorm(20)*.5
#' dat$y <- dat$y + b[fac]
#' dat$fac <- as.factor(fac)
#'
#' # GAM model with one random terms
#' rm1 <- gam(y ~ s(fac,bs="re")+s(x0)+s(x1)+s(x2)+s(x3),
#' data=dat,method="REML")
#'
#' # Pass the fitted GAM object into rptgam
#' # nboot and nperm of 100 is for illustration purposes
#' # and would typically be set higher.
#' out = rptgam(gamObj=rm1, parallel=TRUE, nboot = 100,
#' nperm = 100, aic=TRUE, select=TRUE, verbose=TRUE, seed=1,
#'             boot_type = 'all', ci_type = 'all',
#'             case_resample = c(TRUE,FALSE))
#'
#' # Alternatively, run the GAM model in rptgam
#' out = rptgam(formula = y ~ s(fac,bs="re")+s(x0)+s(x1)+
#' s(x2)+s(x3), data=dat, gam_pars = list(method="REML"),
#'              parallel=TRUE, nboot = 100, nperm = 100,
#'              aic=TRUE, select=TRUE, verbose=TRUE, seed=1,
#'              boot_type = 'all', ci_type = 'all',
#'              case_resample = c(TRUE,FALSE))
#'
#' # bam + discrete method
#' out = rptgam(formula = y ~ s(fac,bs="re")+s(x0)+s(x1)+
#' s(x2)+s(x3), data=dat, bam_pars = list(discrete=TRUE),
#'              parallel=TRUE, nboot = 100, nperm = 100,
#'              aic=TRUE, select=TRUE, verbose=TRUE, seed=1,
#'              boot_type = 'all', ci_type = 'all',
#'              case_resample = c(TRUE,FALSE))
#'
#'
#' @references Nakagawa, S. & Schielzeth, H. (2010) \emph{Repeatability for
#'      Gaussian and non-Gaussian data: a practical guide for biologists}.
#'      Biological Reviews 85: 935-956.
#' @references https://github.com/aloy/lmeresampler
#'
#' @author Eliezer Pickholtz (eyp3@@cornell.edu)
#' @export

rptgam = function(formula = NULL, data = NULL, gamObj = NULL, rterms = NULL,
                   gam_pars = NULL, bam_pars = NULL,
                   nboot = 0, boot_type = 'param', ci = 0.95, ci_type = 'all',
                   case_resample = NULL, case_tries = 100, case_minrows = NULL,
                   nperm = 0, aic = TRUE, select = TRUE,
                   saveboot = TRUE, savepermute = TRUE,
                   seed = NULL, verbose = TRUE,
                   parallel = TRUE, ncores = -1, cl_type = 'PSOCK', logical = FALSE,
                   blas_threads = NULL
                   ) {

  blas_threads1 = FALSE
  if (!is.null(blas_threads) && is.finite(blas_threads)) {
    if (!'RhpcBLASctl' %in% installed.packages()) install.packages('RhpcBLASctl', quiet = TRUE)
    blas_orig = RhpcBLASctl::blas_get_num_procs()
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    blas_threads1 = TRUE
  }

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(RhpcBLASctl::blas_set_num_threads(blas_orig), silent = TRUE)
  })

  t0 = Sys.time()

  # Default cores
  n_cores = 1
  # Mark if clusters are created
  parallel_active = FALSE
  # Mark is parallel should be applied
  parallel1 = FALSE
  cl = NULL

  if (parallel) {
    cluster_out = cluster_tool(ncores = ncores, logical = logical, cl_type = cl_type)
    if ((n_cores <- cluster_out[[1]]) > 1) {
      parallel1 = TRUE
      cluster_type = cluster_out[[2]]
    }
  }

  if (!is.null(gamObj) & (!is.null(formula) | !missing(data))) {
    stop('Only one of gamObj and formula/data can be non-null')
  } else if (is.null(gamObj)) {
      if (is.null(formula) | missing(data)) {
        stop('When gamObj is null, both formula and data should be non-null')
      } else {
          if (identical(is.null(gam_pars), is.null(bam_pars))) {
            stop('When formula and data are provided, exactly one of gam_pars and bam_pars should be non-null')
          } else if (!any(grepl('bs = "re"', attr(terms(formula), 'term.labels')))) {
            stop('Formula should contain at least one random term')
          } else if (!is.null(gam_pars)) {
              if (!isTRUE(gam_pars) & any(grep('formula|data|control.nthreads|nthreads|cluster|discrete|G|fit', names(unlist(gam_pars))))) {
                stop('gam_pars should not include the following arguments: formula, data, nthreads, cluster, discrete, G, fit')
              } else {
                if (verbose) cat('\nRunning gam model...')
                if (!isTRUE(gam_pars) & ('control' %in% names(gam_pars))) {
                  gam_pars_ctrl = gam_pars$control
                  gam_pars_ctrl$nthreads = n_cores
                  gam_pars = gam_pars[!grepl('control', names(gam_pars))]
                  } else gam_pars_ctrl = list(nthreads = n_cores)
                t2 = Sys.time()
                gamObj_pre = try(eval(parse(text=paste0('mgcv::gam(',
                                             paste(paste0('data=', deparse(substitute(data)), ''),
                                             paste0('formula=', paste(deparse(formula), collapse= ' '), ''),
                                             if (!isTRUE(gam_pars)) paste0(names(gam_pars), '=', unlist(lapply(gam_pars, function(x) if(is.character(x)) paste0('"', x, '"') else x)), collapse = ','),
                                             paste0('control=', list(gam_pars_ctrl)),
                                             paste0('fit=FALSE'),
                                             sep=','), ')'))))
                if (gamObj_pre$family[[1]] != 'gaussian') stop("gam family is not Gaussian")
                if (gamObj_pre$n.paraPen != 0) stop('paraPen option is not currently supported in rpt_gam')
                rownames(gamObj_pre$mf) = NULL
                gamObj_pre$call = gamObj_pre$cl
                # Need here family NULL (when family 'gaussian' to prevent possible error (bug))
                gamObj_pre$call$family = NULL
                gamObj_pre$isbam = FALSE
                gamObj_pre$call$formula = eval(gamObj_pre$call$formula)
                gamObj_pre$call$select = eval(gamObj_pre$call$select)
                gamObj_pre$call$drop.unused.levels = eval(gamObj_pre$call$drop.unused.levels)
                gamObj = try(update(gamObj_pre, fit = TRUE, G = gamObj_pre))
                t3 = Sys.time()
                if ('try-error' %in% class(gamObj)) stop('gam failed. Quitting')
                gamObj$call$fit = TRUE
                t_out = minsec_time(t2, t3)
                if (verbose) cat(paste0('done (', round(t_out[[1]], 2), t_out[[2]], ')\n'))
          }
        } else if (!is.null(bam_pars)) {
            bam_modify()
              if (!isTRUE(bam_pars) & any(grep('formula|data|control.nthreads|nthreads|cluster|G|fit', names(unlist(bam_pars))))) {
                stop('bam_pars should not include the following arguments: formula, data, cluster, nthreads, G, fit')
              } else {
                if (parallel1 & (isTRUE(bam_pars) || !isTRUE(eval(bam_pars$discrete)))) {
                  if (verbose) cat(paste0('\nSetting up parallel processing (', n_cores, ' cores, ', cluster_type, ')...'))
                  cl = cluster_tool(ncores = n_cores, cl_type = cluster_type, create = TRUE)[[3]]
                  if (verbose) cat('done\n')
                  parallel_active = TRUE
                  parallel::clusterExport(cl, 'bam_modify', envir=environment())
                  parallel::clusterEvalQ(cl, bam_modify())
                  if (blas_threads1) {
                    parallel::clusterExport(cl, 'blas_threads', envir=environment())
                    parallel::clusterEvalQ(cl, {
                        if (!'RhpcBLASctl' %in% installed.packages()) install.packages('RhpcBLASctl', quiet = TRUE)
                        RhpcBLASctl::blas_set_num_threads(blas_threads)
                    })
                  }
                }
                if (!isTRUE(bam_pars) & 'control' %in% names(bam_pars)) {
                  bam_pars_ctrl = bam_pars$control
                  bam_pars_ctrl$nthreads = n_cores
                  bam_pars = bam_pars[!grepl('control', names(bam_pars))]
                } else bam_pars_ctrl = list(nthreads = n_cores)
                if (verbose) cat('\nRunning bam model...')
                t4 = Sys.time()
                gamObj_pre = try(eval(parse(text=paste0('mgcv::bam(',
                                             paste(paste0('data=', deparse(substitute(data)), ''),
                                             paste0('formula=', paste(deparse(formula), collapse= ' '), ''),
                                             if (!isTRUE(bam_pars)) paste0(names(bam_pars), '=', unlist(lapply(bam_pars, function(x) if(is.character(x)) paste0('"', x, '"') else x)), collapse = ','),
                                             if (parallel_active) paste0('cluster=cl'),
                                             paste0('nthreads=', n_cores),
                                             paste0('control=', list(bam_pars_ctrl)),
                                             paste0('fit=FALSE'),
                                             sep=','), ')'))))
                if (gamObj_pre$family[[1]] != 'gaussian') stop("gam family is not Gaussian")
                if (gamObj_pre$n.paraPen != 0) stop('paraPen option is not currently supported in rpt_gam')
                rownames(gamObj_pre$mf) = NULL
                gamObj_pre$call = gamObj_pre$cl
                gamObj_pre$call$family = NULL
                gamObj_pre$call$formula = eval(gamObj_pre$call$formula)
                gamObj_pre$call$discrete = if (grepl('bam()', gamObj_pre$call[1]) & (is.null(gamObj_pre$call$method) || gamObj_pre$call$method == 'fREML')) eval(gamObj_pre$call$discrete) else NULL
                gamObj_pre$call$select = eval(gamObj_pre$call$select)
                gamObj_pre$call$drop.unused.levels = eval(gamObj_pre$call$drop.unused.levels)
                gamObj_pre$isbam = TRUE
                gamObj <- try(update(gamObj_pre, fit = TRUE, G = gamObj_pre))
                t5 = Sys.time()
                if ('try-error' %in% class(gamObj)) stop('bam failed. Quitting')
                gamObj$call$fit = TRUE
                if (!isTRUE(gamObj_pre$call$discrete)) gamObj_pre$Sl = gamObj$Sl
                t_out = minsec_time(t4, t5)
                if (verbose) cat(paste0('done (', round(t_out[[1]], 2), t_out[[2]], ')\n'))
                }}}
                data = NULL
                rm(data)
      } else {
        if (!inherits(gamObj, 'gam')) stop('gamObj is not an object of class gam')
        if (!is.null(bam_pars)) warning('gamObj was passed. bam_pars will be ignored', immediate. = TRUE)
        if (!is.null(gam_pars)) warning('gamObj was passed. gam_pars will be ignored', immediate. = TRUE)
        if (!is.null(gamObj$paraPen)) stop('paraPen option is not currently supported in rpt_gam')
        # Note that other families requore several checks/changes, such as possibly changing the predict(type=) param (with gaussian link is the same as response and deviance ,etc... (as guassian has identity link, which is the response)). maybe also 'simulate'/'fitted' migh require a change? or not..
        if (gamObj$family[[1]] != 'gaussian') stop("gam family is not Gaussian")
        gamObj$call$family = NULL
        # Just in case, to avoid parallel issues:
        # Add 'mgcv' to the call in order to skip package export to clusters
        if (!grepl('mgcv::', gamObj$call[1])) {
          if ((grepl('gam()', gamObj$call[1]))) {
            gamObj$call[1] = as.call(list(quote(mgcv::gam)))
            } else if ((grepl('bam()', gamObj$call[1]))) {
            gamObj$call[1] = as.call(list(quote(mgcv::bam)))
            }
        }
        # Add the eval formula to the call, in case formula was passed from a variable
        gamObj$call$formula = eval(gamObj$call$formula)
        gamObj$call$discrete = if (grepl('bam()', gamObj$call[1]) & gamObj$method == 'fREML') eval(gamObj$call$discrete) else NULL
        gamObj$call$select = eval(gamObj$call$select)
        gamObj$call$drop.unused.levels = eval(gamObj$call$drop.unused.levels)
        gamObj_pre = update(gamObj, fit = FALSE)
        rownames(gamObj_pre$mf) = NULL
        gamObj_pre$call = gamObj_pre$cl
        gamObj_pre$isbam = inherits(gamObj, 'bam')
        # speedup subsequent fits:
        if (!isTRUE(gamObj_pre$call$discrete) & isTRUE(gamObj_pre$isbam)) gamObj_pre$Sl = gamObj$Sl
        if (inherits(gamObj, 'bam')) bam_modify()
      }

  if (!inherits(gamObj, 'bam')) gamObj_pre$call$in.out = list(sp = gamObj$sp, scale = gamObj$sig2)

  if (!verbose) {
    pb_option = getOption('pboptions')$type
    pbapply::pboptions(type = 'none')
  } else pbapply::pboptions(type = 'timer', style=1)

  try({

    if ((parallel1 & !parallel_active) &
      ((is.numeric(nboot) && nboot > 1) ||
       (is.numeric(nperm) && nperm > 1) ||
       isTRUE(aic) || isTRUE(select))) {

      if (verbose) cat(paste0('\nSetting up parallel processing (', n_cores, ' cores, ', cluster_type, ')...'))
      cl = cluster_tool(ncores = n_cores, cl_type = cluster_type, create = TRUE)[[3]]
      if (verbose) cat('done\n')
      parallel_active = TRUE
      if (inherits(gamObj, 'bam')) {
        parallel::clusterExport(cl, 'bam_modify', envir=environment())
        parallel::clusterEvalQ(cl, bam_modify())
      }
      if (blas_threads1) {
        parallel::clusterExport(cl, 'blas_threads', envir=environment())
        parallel::clusterEvalQ(cl, {
            if (!'RhpcBLASctl' %in% installed.packages()) install.packages('RhpcBLASctl', quiet = TRUE)
            RhpcBLASctl::blas_set_num_threads(blas_threads)
        })
      }
    }

    gamObj_ctrl_slow = gamObj_ctrl_fast = gamObj$control
    gamObj_ctrl_slow$nthreads = 1
    gamObj_ctrl_fast$nthreads = n_cores

    gamObj_info = gam_rand_info(gamObj, rterms = rterms, n_cores = n_cores, cl = cl)

    if (length(gamObj_info$randeff_terms_all) == 0) stop('rpt_gam requires at least one random term')
    if (any(duplicated(gamObj_info$randeff_terms_all))) stop('Random terms cannot be specified more than once in the model')
    if (any(unlist(lapply(gamObj$smooth, function(x) if (isTRUE(x[['random']])) x$by !='NA'))) |
        sum(unlist(lapply(gamObj$smooth, '[[', 'term')) %in% gamObj_info$randeff_terms_all) > length(gamObj_info$randeff_labels_all)) {
      # This error includes interaction of random effects that weren't specified in rterms
      stop('rpt_gam does not support gam models containings random effects with an interaction term')
    }
    if (any(!gamObj_info$rterms %in% gamObj_info$randeff_terms_all)) {
               stop(paste0('rterms ', gamObj_info$rterms[which(!gamObj_info$rterms %in% gamObj_info$randeff_terms_all)],
                           ' are not included or specified correctly as random terms in the model.' ))
    }

    rpt_calc = rpt_func(vcomp = gamObj_info$vcomp, rterms = gamObj_info$rterms, rlabels = gamObj_info$rlabels, var_fixed_pred = gamObj_info$var_fixed_pred)

    # Likelihood ratio test p-values for the random terms (better than using LRT or ANOVA between models),
    # See mgcv::summary.gam for more details, and https://stats.stackexchange.com/questions/274151/anova-to-compare-models
    # according to summary.gam:  "See Wood (2013a) for details: "the test is based on a likelihood ratio statistic (with the
    #                            reference distribution appropriate for the null hypothesis on the boundary of the parameter space).""
    # see https://rdrr.io/cran/mgcv/src/R/mgcv.r for code (reTest) for the test and associated p-values (a little complex)
    gamObj_summary = summary(gamObj)

    df_lrt = data.frame(term = gamObj_info$rterms,
                        edf = gamObj_summary[['s.table']][gamObj_info$rlabels, , drop=FALSE][, 1],
                        Ref.df = gamObj_summary[['s.table']][gamObj_info$rlabels, , drop=FALSE][, 2],
                        chi.sq = gamObj_summary[['chi.sq']][gamObj_info$rlabels],
                        F = gamObj_summary[['s.table']][gamObj_info$rlabels, , drop=FALSE][, 3],
                        p_value = round(gamObj_summary[['s.table']][gamObj_info$rlabels, , drop=FALSE][, 4], 8))
    rpt_check = TRUE

  })

  if (!exists('rpt_check') & !is.null(formula)) {
    warning('rpt_gam failed. Returning gam model', immediate. = TRUE)
    return(gamObj)
  }


  if (parallel_active) {
      parallel::clusterExport(cl, c('gamObj_pre', 'gamObj_info', 'gamObj_ctrl_slow', 'gamObj_ctrl_fast', 'seed'), envir=environment())
      parallel::clusterEvalQ(cl, {gamObj_info$rand_mat = do.call(cbind, lapply(gamObj_info$randeff_terms_all,
                                        function(i) model.matrix(reformulate(c(i, -1)), data = gamObj_pre$mf)))
                                  colnames(gamObj_info$rand_mat) = names(gamObj_info$coef_rand)
                                  rownames(gamObj_info$rand_mat) = NULL
                                  })

    # Adding object to prefit gam objects, to speedup refitted gam
        if (gamObj_pre$isbam) {
            parallel::clusterEvalQ(cl, {gamObj_pre$lpmatrix_out = cbind(gamObj_info$fixed_mat, gamObj_info$rand_mat)[, names(gamObj_info$coefficients)]; NULL})
            if (isTRUE(gamObj_pre$call$discrete)) {
                # below line assumes rho=0. otherwise, this optimization will be skipped
                parallel::clusterEvalQ(cl, {gamObj_pre$qrx_R = mgcv:::XWXd(gamObj_pre$Xd, gamObj_pre$w, gamObj_pre$kd, gamObj_pre$ks, gamObj_pre$ts, gamObj_pre$dt,
                                               gamObj_pre$v, gamObj_pre$qc, 1, gamObj_pre$drop, -1, -1, -1); NULL})
            }
        }

      parallel::clusterExport(cl, c('rpt_func', 'gam_rand_info', 'gam_vcomp_mod', 'gam_boot_rpt', 'gam_boot_param_rpt', 'gam_boot_resid_rpt', 'gam_boot_cgr_rpt', 'gam_boot_reb_rpt', 'gam_boot_case_rpt'), envir=environment())
  }

  gamObj_info$rand_mat = do.call(cbind, lapply(gamObj_info$randeff_terms_all, function(i) model.matrix(reformulate(c(i, -1)), data = gamObj_pre$mf)))

  colnames(gamObj_info$rand_mat) = names(gamObj_info$coef_rand)
  rownames(gamObj_info$rand_mat) = NULL

  if (!parallel_active) {
      # Adding object to prefit gam objects, to speedup refitted gam
        if (gamObj_pre$isbam) {
            gamObj_pre$lpmatrix_out = cbind(gamObj_info$fixed_mat, gamObj_info$rand_mat)[, names(gamObj_info$coefficients)]
            if (isTRUE(gamObj_pre$call$discrete)) {
                # below line assumes rho=0. otherwise, this optimization will be skipped
                gamObj_pre$qrx_R = mgcv:::XWXd(gamObj_pre$Xd, gamObj_pre$w, gamObj_pre$kd, gamObj_pre$ks, gamObj_pre$ts, gamObj_pre$dt,
                                               gamObj_pre$v, gamObj_pre$qc, 1, gamObj_pre$drop, -1, -1, -1)
            }
        }
  } else {gamObj_pre = NULL; rm(gamObj_pre)}

  sim_bind = NULL

  if (is.numeric(nboot) && nboot > 1) {

    # case should go first, as the other modify the y variable
    boot_types_all = c('case', 'param', 'resid', 'cgr', 'reb0', 'reb1', 'reb2')

    boot_type_orig = boot_type
    boot_type = if ('all' %in% boot_type) boot_types_all else boot_types_all[boot_types_all %in% boot_type]

    if (any(c('reb0', 'reb1', 'reb2') %in% boot_type) & (length(gamObj_info$randeff_terms_all) > 1)) {
        boot_type = boot_type[!boot_type %in% c('reb0', 'reb1', 'reb2')]
        if (!'all' %in% boot_type_orig) warning('reb bootstrap methods are not supported currently with > 1 random term. Skipping this boot test', immediate. = TRUE)
    }

    if (length(nboot) > 1 & (length(boot_type) != length(nboot))) stop("Mismastch between number of nboot and boot_type items")
    if (length(nboot) == 1 & (length(boot_type) > 1)) nboot = rep(nboot, length(boot_type))
    if (length(boot_type) == 0) stop("No boot types available. Skipping bootstrapping.")

    if (parallel_active) environment(boot_loop) = .GlobalEnv else environment(boot_loop) = environment()

    try({

    out2 = rpt_calc
    sim_bind = list()
    #if ('bca' %in% ci_type && gamObj_info$mod_nrow > nboot) {
    #    warning('BCa requires at least as many replicates as the number of data points. This ci method will be skipped.')
    #    ci_type = ci_type[ci_type != 'bca']
    #}
    ci_type = if ('all' %in% ci_type) c('perc', 'bca') else unique(ci_type[ci_type %in% c('perc', 'bca')])
    #if ('bca' %in% ci_type && gamObj_info$mod_nrow > nboot) ci_type = ci_type[ci_type != 'bca']
    #if (length(ci_type == 0)) ci_type = 'perc'

    if (verbose) cat('\nCalculating bootstraped CIs\n')

    for (b in seq(boot_type)) {

      if (verbose) cat(paste0('\n  bootstrap type: ', boot_type[b], '\n'))

      if (boot_type[b] == 'case') {
        if ((length(case_resample) != (length(gamObj_info$randeff_terms_all) + 1)) || any(is.na(case_resample) || is.null(case_resample))) {
          warning("'case_resample' is not a vector of logicals with the same length as the number of grouping variables (including the row-level group). Skipping this boot test", immediate. = TRUE)
          next
        } else if (gamObj_info$mod_nrow == 1 || all(case_resample == FALSE)) {
          warning("'case_resample' inputs are all FALSE. Skipping this boot test", immediate. = TRUE)
          next
        }
      }

      boot_data = if (boot_type[b] != 'resid') {
         gam_boot_rpt(gamObj_info = gamObj_info, df = gamObj$model, boot_type = boot_type[b],
                      resample = case_resample, tries = case_tries, minrows = case_minrows, verbose = verbose)
      }

      if (parallel_active) parallel::clusterExport(cl, 'boot_data', envir=environment())

      try(sim_bind[[boot_type[b]]] <- boot_loop(cl, nboot[b], boot_type[b]))
      if (length(sim_bind) == 0) next

      if (boot_type[b] == 'reb2') {
        # Scale and adjust estimates after bootstrapping
        sim_bind[['reb2']] = gam_boot_rpt(gamObj_info = gamObj_info, boot_type = 'reb2', nboot = nboot[b], tstar = sim_bind[['reb2']])
      }

      boot_data = mod_boot = gamObj_pre_case = NULL
      rm(boot_data, mod_boot, gamObj_pre_case)
      if (parallel_active) parallel::clusterEvalQ(cl, {boot_data = mod_boot = gamObj_pre_case = NULL
                                                       rm(boot_data, mod_boot, gamObj_pre_case)})

      sim_bind[[boot_type[b]]] = setNames(data.frame(do.call(rbind, sim_bind[[boot_type[b]]])),
                          c(unlist(lapply(gamObj_info$rterms,
                                         function(x) rbind(paste0(x, '_adj'), paste0(x, '_unadj')))),
                             'Fixed_adj',
                             'Fixed_unadj',
                             'Residual_adj',
                             'Residual_unadj'))

      boot_ci = NULL

      out2 = cbind(out2, data.frame(t(rbind(apply(sim_bind[[boot_type[b]]], 2, stats::sd))),
                                    sapply(1:ncol(sim_bind[[boot_type[b]]]), function(x) mean(sim_bind[[boot_type[b]]][, x]) - rpt_calc[x , 2])))
      names(out2)[(ncol(out2)-1):ncol(out2)] = c(paste0('se_', boot_type[b]),
                                                 paste0('bias_', boot_type[b]))

      for (m in ci_type) {
        boot_ci = try(sapply(1:ncol(sim_bind[[boot_type[b]]]), function(x)
                          ci_methods(point = rpt_calc[x , 2], boots = sim_bind[[boot_type[b]]][, x],
                                     ci = ci, ci_type = m)),
                      silent = FALSE)

        if (!'try-error' %in% class(boot_ci)) {
          out2 = cbind(out2, data.frame(apply(boot_ci, 2, diff), t(boot_ci)))
          names(out2)[(ncol(out2)-2):ncol(out2)] = c(paste0('width_', boot_type[b], '_', m),
                                                   paste0('lower_', boot_type[b], '_', m),
                                                   paste0('upper_', boot_type[b], '_', m))
          boot_check = TRUE
        }
      }
    }
    }, silent = FALSE)

    if (!exists('boot_check', inherits = FALSE)) {

      cat('Bootstrapping failed. Skipping this part.\n')

      rpt_calc = setNames(as.data.frame(t(rpt_calc[, -1])), rpt_calc$term)
      rownames(rpt_calc) = 'rpt'

    } else {
      rpt_calc = setNames(as.data.frame(t(out2[, -1])), out2$term)
      rownames(rpt_calc) = names(out2)[-1]
      rpt_calc$ci_type = sapply(sub('.*\\_', '', row.names(rpt_calc)), function(x) if (x %in% c('perc', 'bca')) x else NA)
      rpt_calc$boot_type = sapply(gsub('.*_(.+)', '\\1', gsub('.*_(.+)_.*', '\\1', rownames(rpt_calc))),
                                  function(x) if (x %in% c('param', 'resid', 'case', 'cgr', 'reb0', 'reb1', 'reb2')) x else NA)
      rpt_calc$metric = sub("\\_.*", "", rownames(rpt_calc))
      rpt_calc = rpt_calc[c('rpt',
        sort(rownames(rpt_calc)[startsWith(rownames(rpt_calc), 'se_')]),
        sort(rownames(rpt_calc)[startsWith(rownames(rpt_calc), 'bias_')]),
        sort(rownames(rpt_calc)[startsWith(rownames(rpt_calc), 'width_')]),
        sort(rownames(rpt_calc)[startsWith(rownames(rpt_calc), 'lower_')]),
        sort(rownames(rpt_calc)[startsWith(rownames(rpt_calc), 'upper_')])), ]
      rpt_calc = rpt_calc[c(ncol(rpt_calc), ncol(rpt_calc)-1, ncol(rpt_calc)-2, seq(ncol(rpt_calc)-3))]
    }

  }

  if (!isTRUE(saveboot)) {sim_bind = NULL; rm(sim_bind)}

  r_permute = NULL
  p_permute = NULL

  if (is.numeric(nperm) && nperm > 1) {

    if (parallel_active) environment(permute_loop) = .GlobalEnv else environment(permute_loop) = environment()

    try({

      if (verbose) cat('\nRunning permutations\n')

      p_permute = matrix(nrow = 2, ncol = gamObj_info$n_randeff)
      r_permute = data.frame()[1:nperm, ]

      if (gamObj_info$n_randeff > 1) {
        errors = gamObj$model[[gamObj_info$dep_name]] - gamObj_info$fixed_pred
      } else errors1 = gamObj$model[[gamObj_info$dep_name]] - gamObj_info$fixed_pred

      for (i in 1:gamObj_info$n_randeff) {

        if (verbose & gamObj_info$n_randeff > 1) cat(paste0('  random term ', i, ' out of ', gamObj_info$n_randeff, ' (', gamObj_info$rterms[i], ')\n'))

          # Lee, O. E., & Braun, T. M. (2012). Permutation tests for random effects in linear mixed models. Biometrics, 68(2), 486-493.

            # "two permutation tests, one that is based on the best linear unbiased predictors and one that is based on the restricted likelihood ratio test statistic. Both methods involve weighted residuals, with the weights determined by the among- and within-subject variance components. The null permutation distributions of our statistics are computed by permuting the residuals both within and among subjects and are valid both asymptotically and in small samples."
            # -- rpt_gam performs the former.

           # adapted from their code: http://www-personal.umich.edu/~tombraun/LeeBraun/Permutation%20Test%20for%20Random%20Slope.R
           # and https://rdrr.io/cran/predictmeans/src/R/permlmer.R

          if (gamObj_info$n_randeff > 1) {

            # Estimate the reduced model
            mod_red = if (isTRUE(gamObj$call$discrete)) {
                        update(gamObj, fit = TRUE,
                             formula = reformulate(gamObj_info$righthand_labels[setdiff(seq_along(gamObj_info$righthand_labels), gamObj_info$randeff_idx_all[i])]),
                             data = gamObj$model,
                             nthreads = n_cores,
                             control = if (parallel1) gamObj_ctrl_fast else gamObj_ctrl_slow)
                        } else {
                            update(gamObj, fit = TRUE,
                             formula = reformulate(gamObj_info$righthand_labels[setdiff(seq_along(gamObj_info$righthand_labels), gamObj_info$randeff_idx_all[i])]),
                             data = gamObj$model,
                             in.out = if (!inherits(gamObj, 'bam')) list(sp = gamObj_info$sp[!names(gamObj_info$sp) %in% gamObj_info$randeff_labels_all[i]], scale = gamObj_info$sig2),
                             cluster = cl,
                             nthreads = n_cores,
                             control = if (parallel1) gamObj_ctrl_fast else gamObj_ctrl_slow)
                        }
            Lambda1nc = {function() {
                         Lambda1nc0 = diag(
                                          nrow = sum(gamObj_info$randeff_terms_all_dim),
                                          x = unlist(c(rep(0, gamObj_info$randeff_terms_all_dim[i]),
                                                sapply(seq(gamObj_info$n_randeff - 1), function(x) {
                                                  rep(sqrt(gam_vcomp_mod(mod_red)[gamObj_info$randeff_labels_all[-i][x], 1] /
                                                           gam_vcomp_mod(mod_red)['scale', 1]),
                                                      gamObj_info$randeff_terms_all_dim[gamObj_info$randeff_terms_all[-i][x]]
                                                      )
                                                })
                                               )))
                         Lambda1nc0 %*% Lambda1nc0}}()
            idx1 = (gamObj_info$randeff_terms_all_dim[i] + 1):ncol(gamObj_info$rand_mat)
            idx2 = seq(sum(gamObj_info$randeff_terms_all_dim))[
                                 which(unlist(sapply(1:gamObj_info$n_randeff,
                                       function(z) rep(z, gamObj_info$randeff_terms_all_dim[z]))) != i)]
            # V1N and (especially) Ut0 can be very slow and require a lot of memory to calculate
            V1n = gamObj_info$rand_mat[, idx2] %*%
                  tcrossprod(Lambda1nc[idx1, idx1], gamObj_info$rand_mat[, idx2])
            diag(V1n) = diag(V1n) + 1
            # THIS WILL FAIL WITH VERY LARGE DATASETS
            Ut0 = chol(V1n)
            # weighted errors
            errors1 = backsolve(r = Ut0, x = errors, transpose = TRUE)
            mod_red = Lambda1nc = V1n = NULL
            rm(mod_red, Lambda1nc, V1n)

          } else Ut0 = NULL

        if (parallel_active) parallel::clusterExport(cl, c('Ut0', 'errors1'), envir=environment())

        permutes = permute_loop(cl, nperm)

        errors1 = Ut0 = NULL
        if (parallel_active) parallel::clusterExport(cl, c('Ut0', 'errors1'), envir=environment())

        permutes = do.call(rbind, permutes)

        if (!exists('boot_check', inherits = FALSE)) {
          rpt_calc_i_adj = rpt_calc[grepl(paste0(gamObj_info$rterms[i], '_adj'), rpt_calc$term), 2]
          rpt_calc_i_unadj = rpt_calc[grepl(paste0(gamObj_info$rterms[i], '_unadj'), rpt_calc$term), 2]
        } else {
          rpt_calc_i_adj = rpt_calc[1, grepl(paste0(gamObj_info$rterms[i], '_adj'), colnames(rpt_calc))]
          rpt_calc_i_unadj = rpt_calc[1, grepl(paste0(gamObj_info$rterms[i], '_unadj'), colnames(rpt_calc))]
        }

        r_permute = cbind(r_permute, rbind(c(rpt_calc_i_adj, rpt_calc_i_unadj), permutes))
        p_permute[1, i] = sum(c(rpt_calc_i_adj, permutes[, 1]) >= rpt_calc_i_adj) / nperm
        p_permute[2, i] = sum(c(rpt_calc_i_unadj, permutes[, 2]) >= rpt_calc_i_unadj) / nperm
    }
    p_permute = setNames(cbind.data.frame(gamObj_info$rterms, t(p_permute)),
                    c('term', 'p_value_adj', 'p_value_unadj'))
    colnames(r_permute) = unlist(lapply(gamObj_info$rterms, function(x) paste0(x, c('_adj', '_unadj'))))
    rownames(r_permute) = NULL
    perm_check = TRUE
    })

    if (!exists('perm_check', inherits = FALSE)) cat('Permutation failed. Skipping this part.\n')

  }

  if (!isTRUE(savepermute)) r_permute = NULL

  df_aic = NULL

  if (aic) {

    if (verbose) cat('\nCalculating AICs...')
      try({
        t12 = Sys.time()
        #TODO: should AIC be done in parallel, if gam/discrete? Check which is faster
        if (parallel_active) {
          aic_out = list()
          for (x in 1:gamObj_info$n_randeff) {
              aic_out[[length(aic_out) + 1]] = if (isTRUE(gamObj$call$discrete)) {
                    update(gamObj,
                          formula = reformulate(gamObj_info$righthand_labels[setdiff(seq_along(gamObj_info$righthand_labels), gamObj_info$randeff_idx_all[x])]),
                          data = gamObj$model[setdiff(names(gamObj$model), gamObj_info$rterms[x])],
                          nthreads = n_cores, control = gamObj_ctrl_fast)
                    } else {
                    update(gamObj,
                          formula = reformulate(gamObj_info$righthand_labels[setdiff(seq_along(gamObj_info$righthand_labels), gamObj_info$randeff_idx_all[x])]),
                          data = gamObj$model[setdiff(names(gamObj$model), gamObj_info$rterms[x])],
                          in.out = if (!inherits(gamObj, 'bam')) list(sp = gamObj_info$sp[!names(gamObj_info$sp) %in% gamObj_info$randeff_labels_all[x]], scale = gamObj_info$sig2),
                          cluster = cl, nthreads = n_cores, control = gamObj_ctrl_fast)

                    }
            }
        } else aic_out = pbapply::pblapply(1:gamObj_info$n_randeff,
                                           function(x){
                 if (isTRUE(gamObj$call$discrete)) {
                      update(gamObj,
                             formula = reformulate(gamObj_info$righthand_labels[setdiff(seq_along(gamObj_info$righthand_labels), gamObj_info$randeff_idx_all[x])]),
                             data = gamObj$model[setdiff(names(gamObj$model), gamObj_info$rterms[x])],
                             nthreads = 1, control = gamObj_ctrl_slow)
                 } else {
                    update(gamObj,
                             formula = reformulate(gamObj_info$righthand_labels[setdiff(seq_along(gamObj_info$righthand_labels), gamObj_info$randeff_idx_all[x])]),
                             data = gamObj$model[setdiff(names(gamObj$model), gamObj_info$rterms[x])],
                             in.out = if (!inherits(gamObj, 'bam')) list(sp = gamObj_info$sp[!names(gamObj_info$sp) %in% gamObj_info$randeff_labels_all[x]], scale = gamObj_info$sig2),
                             cluster = NULL, nthreads = 1, control = gamObj_ctrl_slow)
                 }
                  }, cl = NULL)

        aic_bind = do.call(rbind,
                           lapply(1:gamObj_info$n_randeff,
                                  function(x) AIC(aic_out[[x]])))

        aic_bind = rbind(aic_bind, AIC(gamObj))

        df_aic = data.frame(model = c(unlist(lapply(gamObj_info$rterms,
                                         function(x) rbind(paste0(x, '_excluded')))),
                                 'Full_model'),
                            aic = NA,
                            lowest_aic = NA)

        df_aic$aic = aic_bind[,1]
        df_aic$lowest_aic = c(df_aic$aic == min(aic_bind[,1]))
        aic_check = TRUE

   })
    if (!exists('aic_check', inherits = FALSE)) {
      cat('AIC failed. Skipping this part.\n')
    } else {
      t13 = Sys.time()
      t_out = minsec_time(t12, t13)
      if (verbose) cat(paste0('done (', round(t_out[[1]], 2), t_out[[2]], ')\n'))
    }
  }

  df_select = NULL

  if (select) {

    if (verbose) cat('\nRunning "select = TRUE" model...')
    try({

      t14 = Sys.time()
      select_old = lapply(gamObj_info$rlabels, function(x) c(gamObj_summary[['s.table']][x, ], gamObj_summary[['chi.sq']][x]))
      g1_pre_select = if (!isTRUE(gamObj$call$discrete)) {
                    update(gamObj, fit = FALSE,
                           select = if (!is.null(gamObj$call$select) && isTRUE(gamObj$call$select)) FALSE else TRUE,
                           control = if (parallel1) gamObj_ctrl_fast else gamObj_ctrl_slow,
                           nthreads = n_cores, cluster = cl)
                      } else {
                        update(gamObj, fit = FALSE,
                           select = if (!is.null(gamObj$call$select) && isTRUE(gamObj$call$select)) FALSE else TRUE,
                           control = if (parallel1) gamObj_ctrl_fast else gamObj_ctrl_slow,
                           nthreads = n_cores)
                      }
      g1_pre_select$call$family = NULL
      g1_pre_select$call = g1_pre_select$cl
      g1_pre_select$lpmatrix_out = if (!parallel_active) gamObj_pre$lpmatrix_out else cbind(gamObj_info$fixed_mat, gamObj_info$rand_mat)[, names(gamObj_info$coefficients)]
      select_new = update(g1_pre_select, G = g1_pre_select, fit = TRUE)

      gamObjNew_summary = summary(select_new)

      select_comp = if (isTRUE(select_new$call$select)) {
        select_new = lapply(gamObj_info$rlabels, function(x) c(gamObjNew_summary[['s.table']][x, ], gamObjNew_summary[['chi.sq']][x]))
        do.call(rbind, (lapply(seq(gamObj_info$rlabels), function(x) rbind(select_old[[x]], select_new[[x]]))))
      } else {
        select_new = lapply(gamObj_info$rlabels, function(x) c(gamObjNew_summary[['s.table']][x, ], gamObjNew_summary[['chi.sq']][x]))
        do.call(rbind, (lapply(seq(gamObj_info$rlabels), function(x) rbind(select_new[[x]], select_old[[x]]))))
      }

      df_select = data.frame(model = unlist(lapply(gamObj_info$rterms,
                                       function(x) rbind(paste0(x, '_unpenalized'), paste0(x, '_penalized')))),
                             edf = NA,
                             Ref.df = NA,
                             chi.sq = NA,
                             `F` = NA,
                             p_value = NA
                             )

      df_select$edf = select_comp[, 1]
      df_select$Ref.df = select_comp[, 2]
      df_select$chi.sq = select_comp[, 5]
      df_select$`F` = select_comp[, 3]
      df_select$p_value = select_comp[, 4]

      select_check = TRUE
    })
    if (!exists('select_check', inherits = FALSE)) {
      cat('SELECT test failed. Skipping this part.\n')
    } else {
      t15 = Sys.time()
      t_out = minsec_time(t14, t15)
      if (verbose) cat(paste0('done (', round(t_out[[1]], 2), t_out[[2]], ')\n'))
    }

  }

  t1 = Sys.time()
  if (verbose) {
    t_out = minsec_time(t0, t1)
    cat(paste('\nTotal run time:', round(t_out[[1]], 2), t_out[[2]], '\n'))
    } else pbapply::pboptions(type = pb_option)

  res = list(rpt = rpt_calc,
              lrt = df_lrt,
              p_permute = p_permute,
              aic = df_aic,
              select = df_select,
              rpt_boot_data = if (isTRUE(saveboot)) sim_bind,
              rpt_permute_data = if (isTRUE(savepermute)) r_permute,
              call = match.call(),
              run_info = list(run_time = as.numeric(difftime(t1, t0, units = 'min')),
                              parallel = if (parallel_active) c(n_cores, cluster_type) else c(1, 'None'),
                              session = sessionInfo()),
              gamObj = gamObj
              )
  res = suppressWarnings(Filter(function(x) !all(is.na(x)), res))
  res = lapply(res, function(x) if(is.data.frame(x)) {rownames(x) <- NULL; x} else x)

  class(res) = 'rptgam'
  return(res)

}
