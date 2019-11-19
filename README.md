
# rptGam

## Overview

**rptGam** is a package in R that estimates the repeatabilities of one
or more random rerms in generalized additive models (GAM), and
specifically gam and bam models from
[mgcv](https://cran.r-project.org/web/packages/mgcv/index.html).
**rptGam** was inspired by, and has similar functionality to,
[rptR](https://cran.r-project.org/web/packages/rptR/index.html), which
calculates repeatabilities from
[lme4](https://cran.r-project.org/web/packages/lme4/index.html) linear
mixed-effects models. See rptR’s
[vignette](https://cran.r-project.org/web/packages/rptR/vignettes/rptR.html)
and
[paper](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12797/full)
for an overview og repeatabilities.

Notes: \* **rptGam** currenly supports only Gaussian models \* rptgam
accepts any valid mgcv model containing at least one random term.
Interactions are allowed as long as none of the random terms interact
with any other term (which would invalidate the random terms’
repatability estimate) \* The code in **rptGam** includes a small ‘hack’
to a few mgcv functions to enable faster refits of mgcv models, which
allows for higher bootstrap and permutation repeats.

## Installation

``` r
# To install the latest version from Github:
# install.packages("devtools")
devtools::install_github("elipickh/rptGam")
```

## Usage

``` r
library(rptGam)

# Simulate data with one continuous predictor, one factor predictor and two grouped (factor) predictors
set.seed(1)
df <- mgcv::gamSim(1, n = 200, scale = 2, verbose = FALSE)[, 1:2]
fac1 <- sample(1:3, 200, replace=TRUE)
fac2 <- sample(1:5, 200, replace=TRUE)
b <- rnorm(20) * .5
df$x1 <- sample(1:4, 200, replace = TRUE)
df$y <- df$y + b[fac1] + b[fac2] + df$x1
df$fac1 <- as.factor(fac1)
df$fac2 <- as.factor(fac2)
df$x1 <- as.factor(df$x1)

str(df)
```

    ## 'data.frame':    200 obs. of  5 variables:
    ##  $ y   : num  13.73 18.38 4.87 3.68 4.6 ...
    ##  $ x0  : num  0.266 0.372 0.573 0.908 0.202 ...
    ##  $ x1  : Factor w/ 4 levels "1","2","3","4": 4 4 2 3 4 4 3 4 4 2 ...
    ##  $ fac1: Factor w/ 3 levels "1","2","3": 1 3 1 2 3 3 2 2 2 3 ...
    ##  $ fac2: Factor w/ 5 levels "1","2","3","4",..: 5 1 4 4 2 5 2 3 1 1 ...

The *rptgam* function can take as input either a model output from
mgcv::gam (or mgcv::bam) or a GAM formula. If using the latter option,
set either *gam\_pars* or *bam\_pars* to TRUE (to run a gam or a bam
model, respectively). Alternatively, *gam\_pars*/*bam\_pars* can be a
list of argumetns which will be passed to mgcv when running the model.

``` r
# Option 1: Create a model in mgcv and pass it to rptgam
mod <- mgcv::gam(y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, bs = "re"), 
                 data = df, method = "REML")
out <- rptgam(gamObj = mod)
```

    ## 
    ## Setting up parallel processing (4 cores, PSOCK)...done
    ## 
    ## Calculating AICs...done (0.06sec)
    ## 
    ## Running "select = TRUE" model...done (0.04sec)
    ## 
    ## Total run time: 2.34 sec

``` r
# Option 2: run the model in rptgam
out <- rptgam(formula = y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, bs = "re"), 
              data = df, gam_pars = TRUE)
```

    ## 
    ## Running gam model...done (0.04sec)
    ## 
    ## Setting up parallel processing (4 cores, PSOCK)...done
    ## 
    ## Calculating AICs...done (0.03sec)
    ## 
    ## Running "select = TRUE" model...done (0.03sec)
    ## 
    ## Total run time: 2.46 sec

Note that setting *gam\_pars* to TRUE will use the *default* settings in
mgcv::gam. The default *method* setting in mgcv::gam is ‘GCV.Cp’, which
might not be the best method for estimating random effect models. So
here we change the method to
REML:

``` r
out <- rptgam(formula = y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, bs = "re"), 
              data = df, gam_pars = list(method = 'REML'), verbose = FALSE)

# Alternativley, we could use mgcv::bam, which is often much faster for large datasets 
# (and uses 'fREML as its default *method* setting). Using DISCRETE = TRUE can also help 
# speed up computation (see ?mgcv::bam for more details about *bam* and *discrete*)
#out <- rptgam(formula = y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, bs = "re"), 
# data = df, bam_pars = list(discrete = TRUE))
```

Re-running rptgam with bootstraps and permutations *nboot* and *nperm*
should probably increased much higher than
100.

``` r
out <- rptgam(formula = y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, bs = "re"), 
              data = df, gam_pars = TRUE, parallel = TRUE, nboot = 100, nperm = 100, 
              aic = TRUE, select = TRUE, verbose = TRUE, seed = 1, boot_type = 'all', 
              ci_type = 'all', case_resample = c(TRUE, TRUE, FALSE))
```

    ## 
    ## Running gam model...done (0.03sec)
    ## 
    ## Setting up parallel processing (4 cores, PSOCK)...done
    ## 
    ## Calculating bootstraped CIs
    ## 
    ##   bootstrap type: case
    ##     case resamples (from highest to lowest level): 
    ##      fac1:TRUE, fac2:TRUE, .row:FALSE 
    ## 
    ##   bootstrap type: param
    ## 
    ##   bootstrap type: resid
    ## 
    ##   bootstrap type: cgr
    ## 
    ## Running permutations
    ##   random term 1 out of 2 (fac1)
    ##   random term 2 out of 2 (fac2)
    ## 
    ## Calculating AICs...done (0.03sec)
    ## 
    ## Running "select = TRUE" model...done (0.03sec)
    ## 
    ## Total run time: 5.5 sec

The complete gam object output from mgcv can be accessed using *$gamObj*

``` r
summary(out$gamObj)
```

    ## 
    ## Family: gaussian 
    ## Link function: identity 
    ## 
    ## Formula:
    ## y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, bs = "re")
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  9.01420    0.68198  13.218   <2e-16 ***
    ## x12          0.06464    0.80394   0.080    0.936    
    ## x13          0.08660    0.79215   0.109    0.913    
    ## x14          1.27938    0.83259   1.537    0.126    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##               edf Ref.df     F p-value  
    ## s(x0)   2.056e+00  2.572 2.310  0.0711 .
    ## s(fac1) 1.896e-10  3.000 0.000  0.8993  
    ## s(fac2) 2.198e+00  4.000 1.233  0.0653 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0545   Deviance explained =  8.9%
    ## GCV = 15.018  Scale est. = 14.398    n = 200

Access various outputs from the rptgam object:

  - LRT:

<!-- end list -->

``` r
out$lrt
```

    ##   term          edf Ref.df       chi.sq            F    p_value
    ## 1 fac1 1.896478e-10      3 2.014631e-11 6.715438e-12 0.89934186
    ## 2 fac2 2.198301e+00      4 4.933994e+00 1.233498e+00 0.06529449

  - Bootstrap: includes se, interval width, and ci for various bootstrap
    methods

<!-- end list -->

``` r
print(out$rpt, digits = 2)
```

    ##    metric boot_type ci_type fac1_adj fac1_unadj fac2_adj fac2_unadj Fixed_adj
    ## 1     rpt      <NA>    <NA>  1.5e-12    1.4e-12  3.1e-02    2.9e-02   4.5e-02
    ## 2      se      case    <NA>  1.4e-01    1.4e-01  8.8e-02    5.6e-02   2.5e+01
    ## 3      se       cgr    <NA>  9.7e-03    9.2e-03  3.3e-02    3.2e-02   3.3e-02
    ## 4      se     param    <NA>  9.7e-03    8.9e-03  3.4e-02    3.1e-02   3.3e-02
    ## 5      se     resid    <NA>  1.0e-02    9.4e-03  2.6e-02    2.4e-02   3.5e-02
    ## 6    bias      case    <NA>  3.1e-02    2.6e-02  4.3e-02   -2.7e-03   9.7e+00
    ## 7    bias       cgr    <NA>  4.6e-03    4.4e-03  2.6e-03    2.1e-03   1.9e-02
    ## 8    bias     param    <NA>  4.5e-03    4.2e-03 -8.9e-04   -1.2e-03   2.0e-02
    ## 9    bias     resid    <NA>  4.7e-03    4.4e-03 -9.9e-03   -9.8e-03   2.4e-02
    ## 10  width      case     bca  1.1e-12    7.2e-04  8.6e-02    5.2e-01   1.2e-01
    ## 11  width      case    perc  2.4e-01    1.6e-01  3.1e-01    7.7e-02   5.1e+01
    ## 12  width       cgr     bca  2.4e-03    4.4e-03  1.5e-01    1.5e-01   1.0e-01
    ## 13  width       cgr    perc  3.6e-02    3.3e-02  1.1e-01    1.0e-01   1.3e-01
    ## 14  width     param     bca  6.7e-12    3.7e-04  1.5e-01    1.3e-01   9.0e-02
    ## 15  width     param    perc  3.3e-02    3.0e-02  1.1e-01    1.0e-01   1.1e-01
    ## 16  width     resid     bca  5.6e-03    6.6e-03  1.2e-01    1.1e-01   8.2e-02
    ## 17  width     resid    perc  3.8e-02    3.6e-02  8.7e-02    8.3e-02   1.3e-01
    ## 18  lower      case     bca  5.6e-13    2.0e-14  1.2e-12    3.6e-04   8.2e-07
    ## 19  lower      case    perc  1.4e-12    5.2e-14  3.0e-12    1.0e-12   2.8e-02
    ## 20  lower       cgr     bca  7.2e-14    6.8e-14  1.5e-12    2.1e-12   4.3e-03
    ## 21  lower       cgr    perc  7.4e-14    7.1e-14  1.2e-12    1.1e-12   1.1e-02
    ## 22  lower     param     bca  6.7e-14    6.6e-14  2.5e-12    2.5e-12   1.3e-02
    ## 23  lower     param    perc  1.3e-13    1.3e-13  1.0e-12    9.5e-13   1.4e-02
    ## 24  lower     resid     bca  2.7e-14    2.7e-14  5.7e-12    1.3e-11   6.1e-03
    ## 25  lower     resid    perc  7.4e-14    6.9e-14  8.3e-13    7.9e-13   1.3e-02
    ## 26  upper      case     bca  1.7e-12    7.2e-04  8.6e-02    5.2e-01   1.2e-01
    ## 27  upper      case    perc  2.4e-01    1.6e-01  3.1e-01    7.7e-02   5.1e+01
    ## 28  upper       cgr     bca  2.4e-03    4.4e-03  1.5e-01    1.5e-01   1.1e-01
    ## 29  upper       cgr    perc  3.6e-02    3.3e-02  1.1e-01    1.0e-01   1.4e-01
    ## 30  upper     param     bca  6.7e-12    3.7e-04  1.5e-01    1.3e-01   1.0e-01
    ## 31  upper     param    perc  3.3e-02    3.0e-02  1.1e-01    1.0e-01   1.3e-01
    ## 32  upper     resid     bca  5.6e-03    6.6e-03  1.2e-01    1.1e-01   8.8e-02
    ## 33  upper     resid    perc  3.8e-02    3.6e-02  8.7e-02    8.3e-02   1.4e-01
    ##    Fixed_unadj Residuals_adj Residuals_unadj
    ## 1      4.3e-02        0.9694          0.9281
    ## 2      3.4e-01        0.1693          0.3196
    ## 3      2.9e-02        0.0351          0.0407
    ## 4      2.9e-02        0.0338          0.0403
    ## 5      3.0e-02        0.0281          0.0374
    ## 6      5.1e-01       -0.0745         -0.5319
    ## 7      1.6e-02       -0.0072         -0.0224
    ## 8      1.7e-02       -0.0036         -0.0202
    ## 9      2.0e-02        0.0051         -0.0151
    ## 10     1.1e-01        0.0705          0.0675
    ## 11     9.5e-01        0.6825          0.9384
    ## 12     9.5e-02        0.1193          0.1379
    ## 13     1.1e-01        0.1165          0.1520
    ## 14     8.0e-02        0.1261          0.1083
    ## 15     1.0e-01        0.1115          0.1445
    ## 16     7.5e-02        0.1264          0.1389
    ## 17     1.1e-01        0.1004          0.1371
    ## 18     8.2e-07        0.9295          0.9037
    ## 19     2.7e-02        0.3175          0.0065
    ## 20     4.3e-03        0.8807          0.8516
    ## 21     1.1e-02        0.8835          0.8223
    ## 22     1.3e-02        0.8739          0.8727
    ## 23     1.4e-02        0.8885          0.8258
    ## 24     6.0e-03        0.8736          0.8502
    ## 25     1.2e-02        0.8996          0.8411
    ## 26     1.1e-01        1.0000          0.9712
    ## 27     9.8e-01        1.0000          0.9449
    ## 28     9.9e-02        1.0000          0.9895
    ## 29     1.2e-01        1.0000          0.9743
    ## 30     9.3e-02        1.0000          0.9811
    ## 31     1.1e-01        1.0000          0.9703
    ## 32     8.1e-02        1.0000          0.9890
    ## 33     1.3e-01        1.0000          0.9782

  - AIC:

<!-- end list -->

``` r
out$aic
```

    ##           model      aic lowest_aic
    ## 1 fac1_excluded 1111.070       TRUE
    ## 2 fac2_excluded 1114.059      FALSE
    ## 3    Full_model 1111.070      FALSE

  - Penalized
    model:

<!-- end list -->

``` r
out$select
```

    ##              model          edf Ref.df       chi.sq            F    p_value
    ## 1 fac1_unpenalized 7.011647e-11      3 7.884324e-12 2.628108e-12 0.89375234
    ## 2   fac1_penalized 1.896478e-10      3 2.014631e-11 6.715438e-12 0.89934186
    ## 3 fac2_unpenalized 2.204277e+00      4 4.952247e+00 1.238062e+00 0.06515100
    ## 4   fac2_penalized 2.198301e+00      4 4.933994e+00 1.233498e+00 0.06529449

  - Bootstrapped data is in $rpt\_boot\_data. For example:

<!-- end list -->

``` r
# Histogram of parametric bootstraps for fac1 (adj):
hist(out$rpt_boot_data$param$fac1_adj)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Histogram of case (nonparamteric) bootstraps for fac2 (adj):
hist(out$rpt_boot_data$case$fac2_adj)
```

![](README_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

  - Permutations data is in $rpt\_permute\_data. Note that the
    permutation data includes the *original* rpt estimate as one of its
    entries

<!-- end list -->

``` r
hist(out$rpt_permute_data$fac1_adj)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# Compare the above histogram with the estimated rpt:
out$rpt$fac1_adj[1]
```

    ## [1] 1.45106e-12

``` r
# Can also be accessed here:
#out$rpt_permute_data$fac1_adj[1]
```

  - The call to rptgam:

<!-- end list -->

``` r
out$call
```

    ## rptgam(formula = y ~ x1 + s(x0) + s(fac1, bs = "re") + s(fac2, 
    ##     bs = "re"), data = df, gam_pars = TRUE, nboot = 100, boot_type = "all", 
    ##     ci_type = "all", case_resample = c(TRUE, TRUE, FALSE), nperm = 100, 
    ##     aic = TRUE, select = TRUE, seed = 1, verbose = TRUE, parallel = TRUE)

  - Run info:

<!-- end list -->

``` r
out$run_info
```

    ## $run_time
    ## [1] 0.09164462
    ## 
    ## $parallel
    ## [1] "4"     "PSOCK"
    ## 
    ## $session
    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS:   /opt/intel/compilers_and_libraries_2019.4.233/mac/mkl/lib/libmkl_rt.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] rptGam_0.1.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3           lattice_0.20-38      RhpcBLASctl_0.18-205
    ##  [4] digest_0.6.22        grid_3.6.1           nlme_3.1-142        
    ##  [7] magrittr_1.5         evaluate_0.14        pbapply_1.4-2       
    ## [10] rlang_0.4.1          stringi_1.4.3        Matrix_1.2-17       
    ## [13] rmarkdown_1.17       splines_3.6.1        tools_3.6.1         
    ## [16] stringr_1.4.0        parallel_3.6.1       xfun_0.11           
    ## [19] yaml_2.2.0           compiler_3.6.1       mgcv_1.8-31         
    ## [22] htmltools_0.4.0      knitr_1.26
