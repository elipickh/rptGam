#' Print an rptgam object
#'
#' @param x An rptgam object.
#' @param digits Integer indicating the number of digits to round. Default is 2.
#' @param ... None used currently.
#' 
#' @return 
#' The repeatabilities (rpt) (adjusted and unadjusted), LRT-test p-values, permuations p-values (adjusted and unadjusted), AIC values, penalized and unpenalized model p-values
#' 
#' @author Eliezer Pickholtz (eyp3@@cornell.edu)
#' @export
print.rptgam <- function(x, digits = 2, ...) {
  cat('Repeatability estimates (adj, unadj) for a random effects GAM model\n\n')
  for (i in 1:nrow(out$lrt)) {
    cat('', as.character(out$lrt$term)[i], '\n')
    cat('   rpt: ', round(out$rpt[i, i * 2 + 2], digits), ', ', 
                   round(out$rpt[i, i * 2 + 3], digits), '\n', sep = '')
    cat('   lrt p-value: ', round(out$lrt[i, 6], digits), '\n', sep = '')
    if (!is.null(out$p_permute)) cat('   perm p-value (n = ', nrow(out$rpt_permute_data), '): ', round(out$p_permute[i, 2], digits), ', ', 
                                                                 round(out$p_permute[i, 3], digits), '\n', sep = '')
    if (!is.null(out$aic)) cat('   AIC: ', round(out$aic[nrow(out$aic), 2], digits), ' w/term, ', 
                                          round(out$aic[i, 2], digits), ' wout/term\n', sep = '')
    if (!is.null(out$select)) cat('   Penalized model p-value: ', round(out$select[i * 2, 6], digits), ' penalized, ', 
                                          round(out$select[i * 2 - 1, 6], digits), ' unpenalized\n', sep = '')
  }
  cat(' Fixed', '\n')
  cat('\n   rpt: ', round(out$rpt[i, 6], digits), ', ', round(out$rpt[i, 7], digits), '\n', sep = '')
  cat(' Residuls', '\n')
  cat('\n   rpt: ', round(out$rpt[i, 8], digits), ', ', round(out$rpt[i, 9], digits), '\n', sep = '')
}
  
  