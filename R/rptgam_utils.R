bam_modify = function() {

  if (suppressWarnings(
        (body(mgcv:::bam)[[52]][[3]][[2]][[3]] == substitute(as.numeric(predict.bam(object, newdata = object$model, block.size = chunk.size, cluster = cluster))) ||
         body(mgcv:::bam)[[52]][[3]][[2]][[3]] == substitute(if (!is.null(G$lpmatrix_out)) c(G$lpmatrix_out %*% object$coefficients) else as.numeric(predict.bam(object, newdata = object$model, block.size = chunk.size, cluster = cluster)))) &

        (body(mgcv:::bam.fit)[[4]][[3]][[6]][[3]][[3]][[4]][[5]] == substitute(X <- w*predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind))) ||
         body(mgcv:::bam.fit)[[4]][[3]][[6]][[3]][[3]][[4]][[5]] == substitute(if (!is.null(G$lpmatrix_out)) X <- w*G$lpmatrix_out[ind, ] else X <- w*predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind)))) &

        (body(mgcv:::bgam.fitd)[[34]][[4]][[4]][[3]][[9]] == substitute(qrx$R <- XWXd(G$Xd, w, G$kd, G$ks, G$ts, G$dt, G$v, G$qc, npt[1], G$drop, ar.stop, ar.row, ar.weight)) ||
         body(mgcv:::bgam.fitd)[[34]][[4]][[4]][[3]][[9]] == substitute(if (!is.null(G$qrx_R) & all(c(ar.stop, ar.row, ar.weight) == -1)) qrx$R <- G$qrx_R else qrx$R <- XWXd(G$Xd, w, G$kd, G$ks, G$ts, G$dt, G$v, G$qc, npt[1], G$drop, ar.stop, ar.row, ar.weight))) &

        (body(mgcv:::ar.qr.up)[[5]][[4]][[5]] == substitute(X <- w * predict(arg$G, newdata = arg$mf[ind, ], type = "lpmatrix", newdata.guaranteed = TRUE, block.size = length(ind))) ||
         body(mgcv:::ar.qr.up)[[5]][[4]][[5]] == substitute(if (!is.null(arg$G$lpmatrix_out)) X <- w * arg$G$lpmatrix_out[ind, ] else X <- w * predict(arg$G, newdata = arg$mf[ind, ], type = "lpmatrix", newdata.guaranteed = TRUE, block.size = length(ind))))
     )) {

    bam = mgcv:::bam
    body(bam)[[52]][[3]][[2]][[3]] = substitute(if (!is.null(G$lpmatrix_out)) c(G$lpmatrix_out %*% object$coefficients) else as.numeric(predict.bam(object, newdata = object$model, block.size = chunk.size, cluster = cluster)))
    bam = compiler::cmpfun(bam)
    environment(bam) = asNamespace('mgcv')
    assignInNamespace('bam', bam, ns='mgcv')

    bam.fit = mgcv:::bam.fit
    body(bam.fit)[[4]][[3]][[6]][[3]][[3]][[4]][[5]] = substitute(if (!is.null(G$lpmatrix_out)) X <- w*G$lpmatrix_out[ind, ] else X <- w*predict(G,newdata=mf[ind,],type="lpmatrix",newdata.guaranteed=TRUE,block.size=length(ind)))
    body(bam.fit)[[6]][[4]][[3]][[2]] = substitute(if (!is.null(G$Sl)) Sl = G$Sl else Sl = Sl.setup(G))
    bam.fit = compiler::cmpfun(bam.fit)
    environment(bam.fit) = asNamespace('mgcv')
    assignInNamespace('bam.fit', bam.fit, ns='mgcv')

    bgam.fitd = mgcv:::bgam.fitd
    body(bgam.fitd)[[34]][[4]][[4]][[3]][[9]] = substitute(if (!is.null(G$qrx_R) & all(c(ar.stop, ar.row, ar.weight) == -1)) qrx$R <- G$qrx_R else qrx$R <- XWXd(G$Xd, w, G$kd, G$ks, G$ts, G$dt, G$v, G$qc, npt[1], G$drop, ar.stop, ar.row, ar.weight))
    bgam.fitd = compiler::cmpfun(bgam.fitd)
    environment(bgam.fitd) = asNamespace('mgcv')
    assignInNamespace('bgam.fitd', bgam.fitd, ns='mgcv')

    ar.qr.up = mgcv:::ar.qr.up
    body(ar.qr.up)[[5]][[4]][[5]] = substitute(if (!is.null(arg$G$lpmatrix_out)) X <- w * arg$G$lpmatrix_out[ind, ] else X <- w * predict(arg$G, newdata = arg$mf[ind, ], type = "lpmatrix", newdata.guaranteed = TRUE, block.size = length(ind)))
    ar.qr.up = compiler::cmpfun(ar.qr.up)
    environment(ar.qr.up) = asNamespace('mgcv')
    assignInNamespace('ar.qr.up', ar.qr.up, ns='mgcv')

    NULL
    } else stop('mgcv version not compatible')
}

boot_loop = function(cl, nboot, boot_type) {
      pbapply::pblapply(1:nboot,
                                         function(x, boot_type) {
                                       if (boot_type != 'case') {
                                        gamObj_pre$y = gamObj_pre$mf[[gamObj_info$dep_name]] = gam_boot_rpt(gamObj_info = gamObj_info,
                                                                                                            seed = if(!is.null(seed)) seed + x,
                                                                                                            pt = boot_data,
                                                                                                            boot_type = boot_type)
                                       } else {
                                            samp_rows = gam_boot_rpt(gamObj_info = gamObj_info, df = gamObj_pre$mf,
                                                                   seed = if(!is.null(seed)) seed + x,
                                                                   pt = boot_data,
                                                                   boot_type = 'case')
                                            gamObj_pre_case = update(gamObj_pre, data = gamObj_pre$mf[samp_rows, ], fit = FALSE)
                                            gamObj_pre_case$lpmatrix_out = if (isTRUE(gamObj_pre$isbam)) {
                                                gamObj_pre$lpmatrix_out[samp_rows, gamObj_pre_case$term.names]
                                                } else cbind(gamObj_info$fixed_mat, gamObj_info$rand_mat)[samp_rows, gamObj_pre_case$term.names]
                                       }
                                       mod_boot = update(gamObj_pre, G = if (boot_type != 'case') gamObj_pre else gamObj_pre_case, fit = TRUE,
                                                         cluster = NULL, nthreads = 1, control= gamObj_ctrl_slow)
                                       if (boot_type != 'reb2') {
                                          rpt_func(vcomp = gam_vcomp_mod(mod_boot), rterms = gamObj_info$rterms, rlabels = gamObj_info$rlabels,
                                                   var_fixed_pred = var({
                                                    if (boot_type != 'case') gamObj_info$fixed_mat else gamObj_pre_case$lpmatrix_out[, -unlist(gamObj_info$rand_all_idx_coef_split)]
                                                    } %*% mod_boot$coefficients[-unlist(gamObj_info$rand_all_idx_coef_split)])[1, 1], simple = TRUE)
                                       } else if (boot_type == 'reb2') {
                                          list(fixef = mod_boot$coefficients[-unlist(gamObj_info$rand_all_idx_coef_split)],
                                               varcomp = gam_vcomp_mod(mod_boot)[c(gamObj_info$rlabels, 'scale'), 1])
                                       }},
                                       boot_type = boot_type
                                       , cl = cl
                                       )
}

permute_loop = function(cl, nperm) {
  pbapply::pbsapply(1:(nperm - 1),
                                     function(x) {
                                      set.seed(if (!is.null(seed)) seed + x)
                                      gamObj_pre$y = gamObj_pre$mf[[gamObj_info$dep_name]] = gamObj_info$fixed_pred +
                                                                                              if (gamObj_info$n_randeff > 1) {
                                                                                                # Permute (weighted) resid, then unweight them
                                                                                                crossprod(Ut0, sample(errors1, replace=FALSE))
                                                                                                # Permute resid
                                                                                              } else sample(errors1, replace=FALSE)
                                      mod_new = update(gamObj_pre, G = gamObj_pre, fit = TRUE,
                                                       cluster = NULL, nthreads = 1, control=gamObj_ctrl_slow)
                                      rpt_func(vcomp = gam_vcomp_mod(mod_new), rterms = gamObj_info$rterms, rlabels = gamObj_info$rlabels,
                                               var_fixed_pred = var(gamObj_info$fixed_mat %*% mod_new$coefficients[-unlist(gamObj_info$rand_all_idx_coef_split)])[1, 1], simple = TRUE)[1:2]
                                     },
                               cl = cl,
                               simplify = FALSE)
}

gam_vcomp_mod = function(x) {
  # Return variance
  # Scaled by default

  # The following code is lifted, nearly verbatim, without permission, from mgcv::gam.vcomp.
  # Author: Simon N. Wood simon.wood@r-project.org

  if (!is.null(x$reml.scale) && is.finite(x$reml.scale))
          scale <- x$reml.scale
      else scale <- x$sig2
  m <- length(x$smooth)
  if (is.null(x$paraPen)) {
      k <- 1
      if (is.null(x$full.sp))
          kf <- -1
      else kf <- 1
  }
  else {
      k <- sum(x$paraPen$sp < 0) + 1
      if (is.null(x$full.sp))
          kf <- -1
      else kf <- length(x$paraPen$full.sp.names) + 1
  }
  idx <- rep("", 0)
  idxi <- rep(0, 0)

  if (m > 0)
      for (i in 1:m) {
          if (!is.null(x$smooth[[i]]$id)) {
            if (x$smooth[[i]]$id %in% idx) {
              ok <- FALSE
            }
            else {
              idx <- c(idx, x$smooth[[i]]$id)
              idxi <- c(idxi, i)
              ok <- TRUE
            }
          }
          else {
            ok <- TRUE
          }
          if (ok) {
            if (length(x$smooth[[i]]$S.scale) != length(x$smooth[[i]]$S))
              warning("S.scale vector doesn't match S list - please report to maintainer")
            for (j in 1:length(x$smooth[[i]]$S.scale)) {
              if (x$smooth[[i]]$sp[j] < 0) {
                x$sp[k] <- x$sp[k]/x$smooth[[i]]$S.scale[j]
                k <- k + 1
                if (kf > 0) {
                  x$full.sp[kf] <- x$full.sp[kf]/x$smooth[[i]]$S.scale[j]
                  kf <- kf + 1
                }
              }
              else {
                x$full.sp[kf] <- x$full.sp[kf]/x$smooth[[i]]$S.scale[j]
                kf <- kf + 1
              }
            }
          }
          else {
            ii <- idxi[idx %in% x$smooth[[i]]$id]
            for (j in 1:length(x$smooth[[ii]]$S.scale)) {
              x$full.sp[kf] <- x$full.sp[kf]/x$smooth[[ii]]$S.scale[j]
              kf <- kf + 1
            }
          }
      }

  vc <- c(scale/x$sp)
  names(vc) <- names(x$sp)
  if (is.null(x$full.sp))
      vc.full <- NULL
  else {
      vc.full <- c(scale/x$full.sp)
      names(vc.full) <- names(x$full.sp)
  }


  if (x$method %in% c("ML", "P-ML", "REML", "P-REML", "fREML") &&
      !is.null(x$outer.info$hess)) {
      if (is.null(x$family$n.theta) || x$family$n.theta <= 0) {
          H <- x$outer.info$hess
      } else {
          ind <- 1:x$family$n.theta
          H <- x$outer.info$hess[-ind, -ind, drop = FALSE]
      }
      scale.est = if (ncol(H) > length(x$sp)) TRUE else FALSE

      if (scale.est) {
        if (is.null(vc.full)) {
        vc <- c(vc, scale)
        names(vc) <- c(names(x$sp), "scale")
        } else {
          vc.full <- c(vc.full, scale)
          names(vc.full) <- c(names(x$full.sp), "scale")
        }
      }

      res <- setNames(data.frame(if (is.null(vc.full)) exp(log(sqrt(vc)))^2 else vc.full), 'variance')
      #rownames(res) <- names(vc)
      return (res)

  } else {
          # Trying both sources of scale, just in case
          scale = c(x$sig2, summary(x)$dispersion)
          return (setNames(data.frame(c(if (is.null(vc.full)) vc else vc.full,
                                        setNames(scale[is.finite(scale)][1], 'scale'))
                                      ), 'variance'))
      }
}

gam_rand_info = function(gamObj, rterms = NULL, vcomp = TRUE, n_cores = 1, cl = NULL) {

    out = list()

    if (vcomp) out$vcomp = gam_vcomp_mod(gamObj)

    out$dep_name = attr(attr(gamObj$terms, 'dataClasses'), 'names')[attr(gamObj$terms, 'response')]
    out$indep_names = names(gamObj$model)[!names(gamObj$model) == out$dep_name]
    out$smooth_labels_all = unlist(lapply(gamObj$smooth, '[[', 'label'))
    out$smooth_terms_all = unlist(lapply(gamObj$smooth, function(i) paste(i$term, collapse = ':')))
    out$randeff_labels_all = unlist(lapply(gamObj$smooth, function(x) if (isTRUE(x[['random']])) x$label))
    out$randeff_terms_all = unlist(lapply(gamObj$smooth, function(x) if (isTRUE(x[['random']])) x$term))
    out$randeff_terms_idx_all = which(out$smooth_terms_all %in% out$randeff_terms_all)
    out$rterms = if (is.null(rterms)) out$randeff_terms_all else rterms
    out$n_randeff = length(out$rterms)
    out$randeff_terms_idx = which(out$smooth_terms_all %in% out$rterms)
    out$rlabels = out$smooth_labels_all[out$randeff_terms_idx]
    out$randeff_idx_all = which(attr(terms(gamObj), 'term.labels') %in% out$rterms)
    out$rand_all_idx_coef_split = lapply(out$randeff_terms_idx_all, function(x) gamObj$smooth[[x]]$first.para:gamObj$smooth[[x]]$last.para)
    out$randeff_terms_all_dim = setNames(unlist(lapply(out$randeff_terms_idx_all, function(x) gamObj$smooth[[x]]$bs.dim)),
                                         out$randeff_terms_all)
    out$coef_fixed = gamObj$coefficients[-unlist(out$rand_all_idx_coef_split)]
    out$coef_rand = gamObj$coefficients[unlist(out$rand_all_idx_coef_split)]
    out$righthand_labels = strsplit(as.character(formula(gamObj))[[3]], '\\+')[[1]]
    out$sigma = sigma(gamObj)
    out$mod_nrow = nrow(gamObj$model)
    out$coefficients = gamObj$coefficients
    out$residuals = gamObj$residuals
    out$sp = gamObj$sp
    out$sig2 = gamObj$sig2
    out$drop.unused.levels = gamObj$call$drop.unused.levels

    out$fixed_mat = if (inherits(gamObj, 'bam')) {
                      # discrete will be skipped if model didn't use this method
                      mgcv::predict.bam(gamObj, type = 'lpmatrix', discrete = TRUE, n.threads = n_cores, cluster = cl)[, names(out$coef_fixed)]
                    } else mgcv::predict.gam(gamObj, type = 'lpmatrix')[, names(out$coef_fixed)]
    rownames(out$fixed_mat) = NULL
    out$fixed_pred = as.vector(out$fixed_mat %*% out$coef_fixed)
    out$var_fixed_pred = var(out$fixed_pred)

    out
}

rpt_func = function(vcomp, rterms, rlabels, var_fixed_pred, simple = FALSE) {

    # Group variance
    var_a = vcomp[row.names(vcomp) %in% rlabels, ]

    # Residual variance
    var_e = vcomp[row.names(vcomp) %in% 'scale', ]

    # Fixed effect variance : "variance in the linear predictor"
    # Note: here we exclude all random terms, including those not specified in rterms
    var_f = var_fixed_pred

    # Denominator variance (adjusted and unadjusted)
    var_p_adj = sum(var_a) + var_e
    var_p_unadj = var_p_adj + var_f

    rpt_calc = list(
                rpt_adj = var_a / var_p_adj,
                rpt_unadj = var_a / var_p_unadj,

                rpt_fixed_adj = var_f / var_p_adj,
                rpt_fixed_unadj = var_f / var_p_unadj,

                rpt_resid_adj = var_e / var_p_adj,
                rpt_resid_unadj = var_e / var_p_unadj
            )
    rpt = c(unlist(lapply(seq(rterms),
                          function(x) rbind(rpt_calc$rpt_adj[x], rpt_calc$rpt_unadj[x]))),
                                  rpt_calc$rpt_fixed_adj,
                                  rpt_calc$rpt_fixed_unadj,
                                  rpt_calc$rpt_resid_adj,
                                  rpt_calc$rpt_resid_unadj)

    return(if (!simple) data.frame(term = c(unlist(lapply(rterms,
                                       function(x) rbind(paste0(x, '_adj'), paste0(x, '_unadj')))),
                           'Fixed_adj',
                           'Fixed_unadj',
                           'Residuals_adj',
                           'Residuals_unadj'),
                                  rpt = rpt
            ) else rpt)
}

gam_boot_rpt = function(gamObj_info, df, boot_type, pt = 1, seed = NULL,
                        resample = NULL, tries = 100, minrows = NULL,
                        verbose = TRUE, nboot = 1, reb_type, tstar = NULL) {

  # Van der Leeden, R., Meijer, E. and Busing F. M. (2008) Resampling multilevel models. In J. de Leeuw and E. Meijer, editors, Handbook of Multilevel Analysis, pages 401–433. New York: Springer.
      # The parametric bootstrap requires the strongest assumptions: the explanatory variables are considered fixed, and both the model (specification) and the distribution(s) are assumed to be correct.
      # The residual bootstrap requires weaker assumptions: apart from considering the explanatory variables as fixed, only the model (specification) is assumed to be correct. This implies, for example, that the residuals are assumed to be homoskedastic.
      # The cases bootstrap, finally, requires minimal assumptions: only the hierarchical dependency in the data is assumed to be specified correctly.

    # (Note that I used the parametric boot method from rptR/lme4, which might be slightly (?) different from the param method referecned above.)

  if (boot_type == 'param') return(gam_boot_param_rpt(gamObj_info = gamObj_info, nboot = nboot, pt = pt, seed = seed))
  if (boot_type == 'cgr') return(gam_boot_cgr_rpt(gamObj_info = gamObj_info, pt = pt, seed = seed))
  if (boot_type == 'resid') return(gam_boot_resid_rpt(gamObj_info = gamObj_info, seed = seed))
  if (boot_type == 'case') return(gam_boot_case_rpt(gamObj_info = gamObj_info, df = df, pt = pt, seed = seed,
                                                    resample = resample, tries = tries, minrows = minrows, verbose = verbose))
  if (boot_type %in% c('reb0', 'reb1', 'reb2')) return(gam_boot_reb_rpt(gamObj_info = gamObj_info, df = df, pt = pt, seed = seed, nboot = nboot,
                                                                        tstar = tstar, reb_type = boot_type))
}

gam_boot_param_rpt = function(gamObj_info, nboot = 1, pt = 1, seed = NULL) {
  # Adapted from lme4::.simulateFun (lme4:::simulate.merMod), using re.form = NA (default option in lme4)
  # "simulates responses from a "merMod" fitted model object, i.e., from the model represented by it", and using re.form=NA (the default) -- "condition on none of the random effects, simulating new values for all of the random effects".

  # This is the method in rptR, which basically uses lme4::simulate, which "simulates responses from a "merMod" fitted model object, i.e., from the model represented by it", and using re.form=NA (the default) -- "condition on none of the random effects, simulating new values for all of the random effects".

  # lme4::bootMer -- "each simulation generates new values of both the “spherical” random effects u and the i.i.d. errors ε, using rnorm() with parameters corresponding to the fitted model x."

  # This method is similar to lmeresampler::parametric_bootstrap.lmerMod

  # pt == 0 runs all code, pt == 1 does first part, pt = object does the iteration

  .resamp_gam.param_rpt_pt1 = function (gamObj_info) {
    #U:
    t(diag(nrow = sum(gamObj_info$randeff_terms_all_dim),
           x = unlist(lapply(1:gamObj_info$n_randeff,
            function (i) rep(sqrt(do.call(`/`,  as.list(gamObj_info$vcomp[row.names(gamObj_info$vcomp) %in% c(gamObj_info$rlabels[i], 'scale'), 1]))),
                gamObj_info$randeff_terms_all_dim[i])))) %*%
      t(gamObj_info$rand_mat))
  }

  .resamp_gam.param_rpt_pt2 = function (gamObj_info, U, nboot) {
    gamObj_info$fixed_pred +
     gamObj_info$sigma *
      ((U %*% matrix(rnorm(sum(gamObj_info$randeff_terms_all_dim) * nboot), ncol = nboot)) +
       matrix(rnorm(gamObj_info$mod_nrow * nboot), ncol = nboot))
  }

  if (class(pt) == 'numeric') {
    if (pt == 1) return(.resamp_gam.param_rpt_pt1(gamObj_info))
    pt = .resamp_gam.param_rpt_pt1(gamObj_info)
  }
  set.seed(seed)
  return(.resamp_gam.param_rpt_pt2(gamObj_info = gamObj_info, U = pt, nboot = nboot))

}

gam_boot_cgr_rpt = function(gamObj_info, pt = 1, seed = NULL) {
  # From #lmeresampler:::cgr_bootstrap.lmerMod:
        # The semi-parametric bootstrap algorithm implemented was outlined by Carpenter, Goldstein and Rasbash (2003).

  .resamp_gam.cgr_rpt_pt1 = function (gamObj_info) {

    Uhat.list = lapply(seq_along(gamObj_info$randeff_terms_all), function(i) {
        u = unname(scale(gamObj_info$coefficients[gamObj_info$rand_all_idx_coef_split[[i]]], scale = FALSE))
        #Uhat:
        u %*% sqrt(gamObj_info$vcomp[gamObj_info$randeff_labels_all[i], 1]) %*% solve(chol((t(u) %*% u) / length(u), pivot = TRUE))

    })
    names(Uhat.list) = gamObj_info$randeff_terms_all
    e = as.numeric(scale(gamObj_info$residuals, scale = FALSE))

    list(Uhat.list = Uhat.list,
         ehat = c(gamObj_info$sigma * e * as.numeric((t(e) %*% e)/length(e))^(-1/2)))

  }

  .resamp_gam.cgr_rpt_pt2 = function (gamObj_info, pt) {

    ustar = lapply(pt$Uhat.list, function(dat) {
        index = sample(x = 1:nrow(dat), size = nrow(dat), replace = TRUE)
        dat[index, ]
    })

    ysim = gamObj_info$fixed_pred +
            #Zbstar.sum:
            Reduce('+',
                  # Zbstar:
                  lapply(1:length(ustar), function(i) {
                      gamObj_info$rand_mat[, names(gamObj_info$coefficients)[gamObj_info$rand_all_idx_coef_split[[i]]]] %*% as.matrix(ustar[[i]])
                      })
                  ) +
            #estar:
            sample(x = pt$ehat, size = length(pt$ehat), replace = TRUE)

    ysim[, 1]

  }

  if (class(pt) == 'numeric') {
    if (pt == 1) return(.resamp_gam.cgr_rpt_pt1(gamObj_info))
    pt = .resamp_gam.cgr_rpt_pt1(gamObj_info)

  }

  set.seed(seed)
  return(.resamp_gam.cgr_rpt_pt2(gamObj_info,
                                 pt = pt))

}

gam_boot_resid_rpt = function(gamObj_info, seed = NULL) {
  # From #lmeresampler:::resid_bootstrap.lmerMod:
        # The residual bootstrap resamples the residual quantities from the fitted linear mixed-effects model in order to generate bootstrap resamples. That is, a random sample, drawn with replacement, is taken from the estimated error terms and the EBLUPS (at each level) and the random samples are combined into bootstrap samples via the fitted model equation.

  set.seed(seed)
  bstar = lapply(gamObj_info$rand_all_idx_coef_split, function(x) {
        J = length(gamObj_info$coefficients[x])
        #bstar:
        gamObj_info$coefficients[x][
                 #bstar.index:
                 sample(x = seq_len(J), size = J, replace = TRUE)
                 ]
  })

  ysim = gamObj_info$fixed_pred +
          #Zbstar.sum:
          Reduce('+',
                 #Zbstar:
                 lapply(1:length(bstar), function(i) {
                    gamObj_info$rand_mat[, names(gamObj_info$coefficients)[gamObj_info$rand_all_idx_coef_split[[i]]]] %*% as.matrix(bstar[[i]])
                  })) +
          #estar:
          sample(x = gamObj_info$residuals, size = gamObj_info$mod_nrow, replace = TRUE)

  ysim[, 1]
}

gam_boot_case_rpt = function(gamObj_info, df, resample, pt = 1, minrows = NULL, tries = 100, seed = NULL, verbose = TRUE) {

  # From #lmeresampler:::case_bootstrap.lmerMod:
        # The cases bootstrap is a fully nonparametric bootstrap that resamples the data with respect to the clusters in order to generate bootstrap samples. Depending on the nature of the data, the resampling can be done only for the higher-level cluster(s), only at the observation-level within a cluster, or at all levels. See Van der Leeden et al. (2008) for a nice discussion of this decision.
        # To resample a given level of the model, the corresponding entry in the logical vector specified in the resample parameter must be set to true.
        # Van der Leeden, R., Meijer, E. and Busing F. M. (2008) Resampling multilevel models. In J. de Leeuw and E. Meijer, editors, Handbook of Multilevel Analysis, pages 401–433. New York: Springer.

      # If only one random term, denoting the individuals ("level 2"), and level 1 is repeated measures, then should probably only resample level 2, but not within indiv (according to p. 415)

  # resample: A logical vector specifying whether each level of the model should be resampled in the cases bootsrap. The levels should be specified from the highest level (largest cluster) of the hierarchy to the lowest (observation-level); for example for students within a school, specify the school level first, then the student level.
  # Should be 1 more than the number of random terms (the extra one being the row unit)

  .resamp_gam.case_rpt_pt1 = function(gamObj_info, df, resample, verbose) {

    # Cluster sorts the random terms alphabetically, and then from the smallest to largest number of factors
    # '.row' is just to indicate the additional row unit (i.e., the lowest group is each row)
    clusters = sort(gamObj_info$randeff_terms_all)[order(gamObj_info$randeff_terms_all_dim[order(gamObj_info$randeff_terms_all)])]

    # We add 1 to clusters, to account for the row-level unit (i.e., the lowest group is each row)
    if ((length(clusters) + 1) != length(resample))
          stop("'resample' is not the same length as the number of grouping variables. Please specify whether to resample the data at each level of grouping (including row-level group).")
    if (verbose) cat('    case resamples (from highest to lowest level): \n    ', paste(c(clusters, '.row'), resample, sep = ':', collapse = ', '), '\n')

    list(resample = resample,
         clusters = clusters,
         # Including dep_name, even if not numeric, as it always needs to have > 1 unique values for gam model to run
         non_numeric = c(gamObj_info$dep_name,
                         names(which(sapply(df[gamObj_info$indep_names], class, simplify = TRUE) != 'numeric'))),
         # gam model requires df with nrows >= number of coefficients. Count how many coefs needed regardless of case boot results
                      # Count all non-random coefs
         min_norand = length(gamObj_info$coefficients) - length(unlist(gamObj_info$rand_all_idx_coef_split)) +
                       # If not resampling the first random term, add on its coefs, along with those of other random terms down the hierarchy up to (and not including) the first resampled term (as case boot will contain all unique values from those terms)
                       # The remaining coefs (if any) will be added in pt2
                       if (isFALSE(resample[1])) sum(gamObj_info$randeff_terms_all_dim[clusters[1:(min(which(resample)) - 1)]]) else 0,
          u_clusters = sapply(1:(gamObj_info$n_randeff + 1), function(x)
              if (x == 1) {
                  if (resample[x]) {
                      unique_c = if (!isFALSE(gamObj_info$drop.unused.levels) &
                                     is.factor(df[[clusters[1]]])) levels(df[[clusters[1]]])
                                 else unique(df[[clusters[1]]])
                      if (length(unique_c) < 3 && !'numeric' %in% class(df[[clusters[1]]])) stop('term "', clusters[1], '" contains < 3 unique values and cannot be resampled')
                      if (length(unique_c) == 1 && 'numeric' %in% class(df[[clusters[1]]])) stop('term "', clusters[1], '" contains 1 unique value and cannot be resampled')
                      # Get row index per level
                      else sapply(unique_c, function(y) which(df[[clusters[1]]] %in% y), simplify = FALSE)
                  } else NULL
              } else if (x == (gamObj_info$n_randeff + 1)) {
                  if (resample[x]) {
                    split(1:gamObj_info$mod_nrow, do.call(paste, df[clusters]), drop=TRUE)
                  } else NULL
              } else {
                if (resample[x]) {
                    unique_c = if (!isFALSE(gamObj_info$drop.unused.levels) &
                                   is.factor(df[[clusters[x]]])) levels(df[[clusters[x]]])
                               else unique(df[[clusters[x]]])
                    if (length(unique_c) < 3 && !'numeric' %in% class(df[[clusters[x]]])) stop('term "', clusters[x], '" contains < 3 unique values and cannot be resampled')
                    if (length(unique_c) == 1 && 'numeric' %in% class(df[[clusters[x]]])) stop('term "', clusters[x], '" contains 1 unique value and cannot be resampled')
                    sapply(split(cbind.data.frame(df[[clusters[x]]], 1:gamObj_info$mod_nrow),
                                 do.call(paste, df[clusters[1:(x-1)]]), drop=TRUE),
                           function(z) split(z[[2]], z[[1]], drop = TRUE), simplify = FALSE)
                } else NULL
              }, simplify = FALSE))
  }

  .resamp_gam.case_rpt_pt2 = function(gamObj_info, df, pt, minrows, tries) {

    if (all(pt$resample == FALSE)) stop("'resample' is all FALSE.")
    if (any(!is.finite(tries))) tries = c(100, 100)
    if (length(tries) == 1) tries = rep(tries, 2)
    tries = sapply(tries, function(m) max(m, 1))
    minrows = if (!is.finite(minrows) || is.null(minrows)) 0 else minrows
    for (try_d in 1:tries[1]) {
        samp_rows = NULL
        # Keep track of total unique rand term values that weren't included in pt1
        min_rand = 0

        for (i in 1:(gamObj_info$n_randeff + 1)) {

            if (i == 1 & pt$resample[1]) {
                # Sample with replacement the highest level group, and get all rows for selected units. Data for highest-level group can be repeated. Column and group structure is maintained.
                non_numeric0 = setdiff(pt$non_numeric, pt$clusters[1])
                for (try_i in 1:tries[2]) {
                    samp = sample(names(pt$u_clusters[[1]]), replace = TRUE)
                    min_rand = length(unique(samp))
                    if (min_rand > 1 || ('numeric' %in% class(df[[pt$clusters[1]]]))) {
                        samp_rows = unlist(pt$u_clusters[[1]][samp], use.names = FALSE)
                        if (length(non_numeric0) == 0 ||
                            !any(sapply(df[non_numeric0], function(x) length(unique(x[samp_rows])) == 1))) break
                    }
                    if (try_i == tries[2]) stop('tries limit reached (', pt$clusters[1], ')')
                }

            } else if (i == (gamObj_info$n_randeff + 1) & pt$resample[i]) {
                    # Sample with replacement from the lowest-level group (i.e., rows). Maintains number of the random higher-level units (and the total number of rows in the df), with some of the row-units being repeated
                non_numeric0 = setdiff(pt$non_numeric, pt$clusters)
                for (try_i in 1:tries[2]) {
                    samp_rows = unlist(sapply(
                                              {if (!is.null(samp_rows)) {
                                                  # Subset each sublist to samp_rows, and remove empty sublists (if any)
                                                  Filter(length, lapply(pt$u_clusters[[gamObj_info$n_randeff + 1]], function(w) samp_rows[samp_rows %in% w]))
                                              } else pt$u_clusters[[gamObj_info$n_randeff + 1]]},
                                              function(x) sample(x, replace = TRUE)), use.names = FALSE)
                    if (!any(sapply(df[non_numeric0], function(x) length(unique(x[samp_rows])) == 1))) break
                    if (try_i == tries[2]) stop('tries limit reached (row-level)')
                }

            } else if (pt$resample[i]) {
                # Within each combination of higher-level terms, sample with replacement lower-level term i, and get all rows for term i. Data for term i can be repeated (and varying df size). Column and group structure is maintained.
                u_cluster_i = if (!is.null(samp_rows)) {
                      # Subset each sublist (and sub-sublist) to samp_rows, and remove empty (sub-)sublists (if any)
                      Filter(length, lapply(pt$u_clusters[[i]], function(z)
                                                       Filter(length,
                                                               lapply(z, function(w)
                                                                 samp_rows[samp_rows %in% w]))))
                } else pt$u_clusters[[i]]
                non_numeric0 = setdiff(pt$non_numeric, pt$clusters[1:(i - 1)])
                non_numeric1 = setdiff(non_numeric0, pt$clusters[i])
                      for (try_i in 1:tries[2]) {
                            samp = sapply(u_cluster_i, function(x) sample(names(x), replace = TRUE), simplify = FALSE)
                            min_rand0 = length(unique(unlist(samp, recursive = FALSE, use.names = FALSE)))
                            if (('numeric' %in% class(df[[pt$clusters[i]]])) || length(unique(unlist(samp, use.names = FALSE))) > 1) {
                                samp_rows = unlist(sapply(seq(u_cluster_i), function(x) u_cluster_i[[x]][samp[[x]]]), use.names = FALSE)
                                if ((length(non_numeric1) == 0) ||
                                    !any(sapply(df[non_numeric1], function(x) length(unique(x[samp_rows])) == 1))) {
                                    min_rand = min_rand + min_rand0
                                    break
                                }
                            }
                            if (try_i == tries[2]) stop('tries limit reached (', pt$clusters[i], ')')
                      }
            }
        }

        # Subset df to those non-resampled terms (after the first resampled term, as those before it are included from pt1) and to samp_rows, and, if the previous condition exists, sum unique value for those rendom terms
        remaining_rand = pt$clusters[which(!pt$resample[(min(which(pt$resample))+1):gamObj_info$n_randeff])]
        if (length(remaining_rand) > 1) {
            min_rand = min_rand + sum(sapply(df[remaining_rand],
                                               function(x) length(unique(x[samp_rows])), USE.NAMES = FALSE))
        }

        if (length(samp_rows) >= max(minrows, pt$min_norand + min_rand)) break
        if (try_d == tries[1]) stop('tries limit reached (df size)')
    }

    return(samp_rows)
  }

  if (class(pt) == 'numeric') {
    if (pt == 1) return(.resamp_gam.case_rpt_pt1(gamObj_info = gamObj_info, df = df, resample = resample, verbose = verbose))
    pt = .resamp_gam.case_rpt_pt1(gamObj_info = gamObj_info, df = df, resample = resample, verbose = verbose)
  }

  set.seed(seed)
  return(.resamp_gam.case_rpt_pt2(gamObj_info = gamObj_info, df = df, pt = pt, minrows = minrows, tries = tries))

}

gam_boot_reb_rpt = function(gamObj_info, df, reb_type, tstar = NULL, nboot = 1, pt = 1, seed = NULL) {
  # Based on:
  # Raymond Chambers & Hukum Chandra (2013) A Random Effect Block Bootstrap for Clustered Data, Journal of Computational and Graphical Statistics, 22:2, 452-470, DOI: 10.1080/10618600.2012.681216
    # "Our approach is semiparametric, in the sense that the marginal model is generated parametrically within the bootstrap while the dependence structure of the model residuals is generated nonparametrically. Furthermore, the proposed bootstrap is simple to implement and seems free of both the distribution and the dependence assumptions of the parametric bootstrap, with its main assumption being that the marginal model is correctly specified"
    # "This approach is semiparametric in the sense that although the marginal bootstrap model is based on the parametric fit to the sample data, the dependence structure in the model residuals is generated nonparametrically."

  # As implemented in lmeresampler::reb_bootstrap -- "Generate random effect block (REB) bootstrap replicates of a statistic for a two-level nested linear mixed-effects model."
  # See the help page for more info

  #lmeresampler:::reb_bootstrap.lmerMod
  # which uses: lmeresampler:::.resample.reb

  # Method currently only available for one random effect (as per lmeresampler and the original article)

  # Adapted to return rpt boots (for 'reb2') and/or bootstraped response vectors

  .resamp_gam.reb_rpt_pt1 = function(gamObj_info, df, reb_type) {

    if (length(gamObj_info$randeff_labels_all) != 1) stop ('REB method currently requires models with exactly one random effect')

    model.mresid = df[[gamObj_info$dep_name]] - gamObj_info$fixed_pred
    u = {function() {t_rand = t(gamObj_info$rand_mat)
                     solve(t_rand %*% gamObj_info$rand_mat) %*% t_rand %*% model.mresid}} ()
    e = model.mresid - gamObj_info$rand_mat %*% u
    if (reb_type == 'reb1') {
      Uhat = scale(u %*% sqrt(gamObj_info$vcomp[gamObj_info$randeff_labels_all, 1]) %*% solve(chol((t(u) %*% u) / nrow(u), pivot = TRUE)),
                   scale = FALSE)
    } else Uhat = u

    list(Uhat = Uhat,
         J = nrow(Uhat),
         model.mresid = model.mresid)

  }

  .resamp_gam.reb_rpt_pt2 = function(gamObj_info, pt) {

    ysim = gamObj_info$fixed_pred +
            (gamObj_info$rand_mat %*%
                #ustar:
                pt$Uhat[
                        #ustar.index:
                        sample(x = seq_len(pt$J), size = pt$J, replace = TRUE),
                       ]
            ) +
            #estar:
            sample(x = pt$model.mresid,
                   size = gamObj_info$mod_nrow,
                   replace = TRUE)

    ysim[, 1]

  }

  .resamp_gam.reb_rpt_reb2post = function(gamObj_info, tstar, nboot) {

    t0 = c(beta = gamObj_info$coef_fixed, sigma = c(gamObj_info$vcomp[gamObj_info$randeff_labels_all, 1], gamObj_info$sigma^2))
    vcs = lapply(tstar, function(x) x$varcomp)
    Sb = log(do.call('rbind', vcs))
    Mb = matrix(rep(apply(Sb, 2, mean), times = nboot),
                nrow = nboot, byrow = TRUE)
    CovSb = cov(Sb)
    SdSb = sqrt(diag(CovSb))
    Db = matrix(rep(SdSb, times = nboot),
                nrow = nboot, byrow = TRUE)
    EW = eigen(solve(CovSb), symmetric = TRUE)
    Whalf = EW$vectors %*% diag(sqrt(EW$values))
    Sbmod = (Sb - Mb) %*% Whalf
    Sbmod = Sbmod * Db
    Lb = exp(Mb + Sbmod)
    tstar = lapply(tstar, unlist)
    tstar = do.call('cbind', tstar)
    idx = 1:length(gamObj_info$coef_fixed)
    fe.star = tstar[idx, ]
    fe.adj = sweep(fe.star, MARGIN = 1, STATS = gamObj_info$coef_fixed - rowMeans(fe.star,
        na.rm = TRUE), FUN = "+")
    vc.star = tstar[-idx, ]
    vc.adj = sweep(vc.star, MARGIN = 1, STATS = t0[-idx]/rowMeans(vc.star,
        na.rm = TRUE), FUN = "*")
    tstar = rbind(fe.adj, vc.adj)

    tstar = lapply(1:nboot, function(x) {
              var_a = tstar['varcomp1', x]
              var_e = tstar['varcomp2', x]
              var_f = var(gamObj_info$fixed_mat %*% tstar[idx, x])[1, 1]
              var_p_adj = var_a + var_e
              var_p_unadj = var_p_adj + var_f

              unname(c(
                rpt_adj = var_a / var_p_adj,
                rpt_unadj = var_a / var_p_unadj,

                rpt_fixed_adj = var_f / var_p_adj,
                rpt_fixed_unadj = var_f / var_p_unadj,

                rpt_resid_adj = var_e / var_p_adj,
                rpt_resid_unadj = var_e / var_p_unadj))
            })
    return(tstar)

  }

  if (is.null(tstar)) {
    if (class(pt) == 'numeric') {
      if (pt == 1) return(.resamp_gam.reb_rpt_pt1(gamObj_info = gamObj_info, df = df, reb_type = reb_type))
      pt = .resamp_gam.reb_rpt_pt1(gamObj_info = gamObj_info, df = df, reb_type = reb_type)
    }

    set.seed(seed)
    return(.resamp_gam.reb_rpt_pt2(gamObj_info = gamObj_info, pt = pt))

  } else {

  return(.resamp_gam.reb_rpt_reb2post(gamObj_info = gamObj_info, tstar = tstar, nboot = nboot))

  }
}

ci_methods = function(point, boots, ci_type = 'perc', ci = 0.95) {

      bca_ci5 = function (conf, boots, point) {
        # As in coxed::bca, but using point estimate for bias (w) and for acceleration

          R <- length(boots)
          zalpha = qnorm((1 + c(-conf, conf))/2)
          # Following is according to coxed::bca (which is based on DiCiccio, T. J. and B. Efron. (1996). Bootstrap Confidence Intervals. Statistical Science. 11(3): 189–212. https://doi.org/10.1214/ss/1032280214)
          w = qnorm(length(boots[boots < point]) / R)
          L = (R - 1) * (point - boots)
          adj.alpha = pnorm(w + (w + zalpha)/(1 - (sum(L^3)/(6 * sum(L^2)^1.5)) * (w + zalpha)))
          # coxed::bca uses default quantile type (7), but leaving as 1 to align with original implemntation
          quantile(boots, prob = adj.alpha, type = 1)
      }

    if (ci_type == 'perc') {
      quantile(boots, probs = c((1 - ci) / 2, 1 - (1 - ci) / 2), na.rm = TRUE)
    } else if (ci_type == 'bca') {
      bca_ci5(point = point, boots = boots, conf = ci)
    }

}

bca_sample = function(A, ci = 0.95, perc_dev = 10, prob = 0.95, B1 = NULL, boots = NULL, point = NULL) {
  # Andrews, D. W., & Buchinsky, M. (2002). On the number of bootstrap repetitions for BC a confidence intervals. Econometric Theory, 18(4), 962-984.
  # paper recommends doing this > 1 (twice) for possibly better accuracy

    # The accuracy obtained by a given choice of B is random, because the bootstrap simulations are random.
    # To determine an appropriate value of B, we specify a bound on the percentage deviation, denoted pdb, and we require that the actual percentage deviation is less than this bound with a specified probability, 1 - r, close to one.
    # The three-step method takes pdb and r as given and specifies a data-dependent method of determining a value of B, denoted B*, such that the desired level of accuracy is achieved.
    # For example, one might take (pdb, r) = (10, .05). In this case, the three-step method determines a value B* such that the percentage deviation of the upper and lower confidence interval lengths is less than 10% each with approximate probability .95.

  # A: acceleration factor
  # B1: initial boot number
  # point: sample statistic (original data)
  # perc_dev: number in (0, 100) of max percentage deviation of each of the upper and lower confidence interval lengths. The lower, the more boostraps are required.
  # prob: number in (0,1), indicating the approx probability that perc_dev holds

  a = (1 - ci)/2
  # a should be between 0.01, 0.99, according to paper
  if (a > 0.99 | a < 0.01) stop('a should be in [0.01, 0.99]')
  r = 1 - prob
  pdb = perc_dev

  if (is.null(B1)) {
    B1 =
        # according to formula:
        ceiling(
        # according to paper examples:
        #floor(
        10000 *
        (a * (1 - a) -
        2 * a * dnorm(qnorm(a)) / dnorm(0) +
        (dnorm(qnorm(a))^2) / (dnorm(0)^2)) *
        (qnorm(1 - r/2)^2) /
        ((qnorm(a) * dnorm(qnorm(a)) * pdb)^2)
        )
    return(B1)
  }
  # initial bias correction:
  Z1 = qnorm(sum(boots < point) / B1)
  a1l = max(pnorm(Z1 + ((Z1 + qnorm(a)) / (1 - A * (Z1 + qnorm(a))))), 0.01)
  a1u = min(pnorm(Z1 + ((Z1 + qnorm(1 - a)) / (1 - A * (Z1 + qnorm(1 - a))))), 0.99)
  v1l = max(floor((B1 + 1) * a1l), 1)
  v1u = min(B1, ceiling((B1 + 1) * a1u))
  C_func = function(a) {((1.5 * (qnorm(1 - a/2)^2) * (dnorm(qnorm(1 - a))^2)) /
                       (2 * (qnorm(1 - a)^2) + 1)) ^
                      (1/3)}
  m1l = ceiling(C_func(a1l) * (B1^(2/3)))
  m1u = ceiling(C_func(1 - a1u) * (B1^(2/3)))
  B2l =
    ceiling(
  10000 * (a * (1 - a) - 2 * a * dnorm(qnorm(a)) / dnorm(0) + (dnorm(qnorm(a))^2) / (dnorm(0)^2)) *
  (qnorm(1 - r / 2)^2) *
  #((B1 / (2 * m1l))^2) * ((boots[v1l + m1l] - boots[v1l - m1l])^2) / (((point - boots[v1l]) * pdb)^2)
      ((B1 / (2 * m1l))^2) * ((boots[min(B1, v1l + m1l)] - boots[max(1, v1l - m1l)])^2) / (((point - boots[v1l]) * pdb)^2)
        )
  B2u =
      ceiling(
    10000 * (a * (1 - a) - 2 * a * dnorm(qnorm(a)) / dnorm(0) + (dnorm(qnorm(a))^2) / (dnorm(0)^2)) *
    (qnorm(1 - r / 2)^2) *
    #((B1 / (2 * m1u))^2) * ((boots[v1u + m1u] - boots[v1u - m1u])^2) / (((boots[v1u] - point) * pdb)^2)
      ((B1 / (2 * m1u))^2) * ((boots[min(B1, v1u + m1u)] - boots[max(1, v1u - m1u)])^2) / (((boots[v1u] - point) * pdb)^2)
          )

  B = max(B1, B2l, B2u)
  B

}

cluster_tool = function(ncores, logical = TRUE, cl_type = 'AUTO', create = FALSE, vanilla = FALSE) {
  n_cores = parallel::detectCores(logical = logical)
  n_cores = if (ncores >= 0) {
              min(n_cores, max(1, ncores))
            } else if (ncores < 0) {
              max(1, n_cores + (ncores + 1))
            }
  cl_type1 = NULL
  cl = NULL
  if (n_cores > 1) {
    cl_type1 = if (cl_type == 'PSOCK' | (cl_type == 'AUTO' & .Platform$OS.type == 'windows')) 'PSOCK'
                   else if (cl_type == 'FORK' | (cl_type == 'AUTO' & .Platform$OS.type != 'windows')) 'FORK'
                   else stop('cl_type is not one of "AUTO", "FORK", or "PSOCK"')
    if (create) {
           cl = parallel::makeCluster(n_cores, type = cl_type1, useXDR = FALSE,
                                 rscript_args = if (vanilla & cl_type1 == 'PSOCK') "--vanilla")
           # Quick test (when 'AUTO') to make sure FORK works (e.g., fails on mac with Apple BLAS)
           if (cl_type == 'AUTO' & cl_type1 == 'FORK') {
             if (class(try(mgcv::bam(travel~s(Rail,bs="re"),
                                     data=do.call("rbind", lapply(1:1000, function(x) Rail)),method="fREML",
                                     cluster=cl))) == 'try-error') {
              warning('FORK cluster failed. Switching to PSOCK.', immediate. = TRUE)
              try(parallel::stopCluster(cl), silent = TRUE)
              cl = parallel::makeCluster(n_cores, type = 'PSOCK', useXDR = FALSE,
                                 rscript_args = if (vanilla) "--vanilla")
              cl_type1 == 'PSOCK'
             }
           }
    }
  }
  return(list(n_cores, cl_type1, cl))
}

minsec_time = function(t0, t1, lim = 120) {
  t = as.numeric(difftime(t1, t0, units = 'sec'))
  if (t > lim) list(t/60, 'min') else list(t, 'sec')
}

