Assembling scan results
================
Frederick Boehm
11/25/2018

## Overview

Our goal here is to assemble the scan results for run 3103, which is
part of the chromosome 19 power study. Recall that we split each 2d scan
(of 1000 by 1000 markers) into 25 condor jobs. We now need to assemble
the 25 results files per scan into a single R object.

## Diagnostics

First, let’s ensure that we have the correct number of files. 25 jobs
failed. This was due to the fact that each of 80 expression traits was
paired with Asah2 gene expression, and Asah2 expression itself was among
the 80… so all 25 of its jobs failed.

So, we should have 79 \* 25 jobs.

``` bash
ls ../results/pvl-run3104/*.txt | wc -l
```

    ##     1975

We see that there are 1975 files, as
    needed.

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0     ✔ purrr   0.2.5
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.8
    ## ✔ tidyr   0.8.2     ✔ stringr 1.3.1
    ## ✔ readr   1.2.1     ✔ forcats 0.3.0

    ## ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

## Read results files

``` r
res <- "../results/pvl-run3104"
dir(res) -> fns
results <- list()
for (i in 1:1975){
  results[[i]] <- read.table(file.path(res, fns[i]), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(filename = fns[i])
}
```

``` r
results2 <- results %>%
  bind_rows()
```

``` r
foo <- str_split(fns, pattern = "_") %>%
  purrr::map_chr(.f = function(x)x[3])
results2$gene1 <- rep(foo, each = 40000)
```

## “Create 79 tibbles” by gene1 values

``` r
library(qtl2pleio)
```

``` r
gene1_79 <- unique(results2$gene1)
lrt <- numeric(length = 79)
for (i in 1:79){
  lrt[i] <- results2 %>%
    filter(gene1 == gene1_79[i]) %>%
    select(1:3) %>%
    calc_lrt_tib()
}
(foo <- tibble(gene1 = gene1_79, gene2 = "ENSMUSG00000024887", lrt = lrt))
```

    ## # A tibble: 79 x 3
    ##    gene1              gene2                lrt
    ##    <chr>              <chr>              <dbl>
    ##  1 ENSMUSG00000009378 ENSMUSG00000024887  17.6
    ##  2 ENSMUSG00000024766 ENSMUSG00000024887  38.4
    ##  3 ENSMUSG00000040565 ENSMUSG00000024887  37.3
    ##  4 ENSMUSG00000041731 ENSMUSG00000024887  44.0
    ##  5 ENSMUSG00000044026 ENSMUSG00000024887  58.1
    ##  6 ENSMUSG00000046138 ENSMUSG00000024887  13.3
    ##  7 ENSMUSG00000046324 ENSMUSG00000024887  26.9
    ##  8 ENSMUSG00000047298 ENSMUSG00000024887  29.6
    ##  9 ENSMUSG00000048120 ENSMUSG00000024887  44.9
    ## 10 ENSMUSG00000048612 ENSMUSG00000024887  64.3
    ## # ... with 69 more rows

``` r
write_csv(x = foo, path = "2018-11-26_lrt-tibble.csv") # these all have gene2 = Asah2
```

## Run 3105

``` bash
ls ../results/pvl-run3105/*.txt | wc -l
```

    ##     1975

We see that there are 1975 files, as needed.

### Read results files

``` r
res <- "../results/pvl-run3105"
dir(res) -> fns
results <- list()
for (i in 1:1975){
  results[[i]] <- read.table(file.path(res, fns[i]), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(filename = fns[i])
}
```

``` r
results2 <- results %>%
  bind_rows()
```

``` r
foo <- str_split(fns, pattern = "_") %>%
  purrr::map_chr(.f = function(x)x[3])
results2$gene1 <- rep(foo, each = 40000)
```

### “Create 79 tibbles” by gene1 values

First, we need to get the Ensembl id for the anchor gene.

``` r
colnames(readRDS("../data-to-condor/expr_anchors.rds"))[1] -> anchor
```

``` r
gene1_79 <- unique(results2$gene1)
lrt <- numeric(length = 79)
for (i in 1:79){
  lrt[i] <- results2 %>%
    filter(gene1 == gene1_79[i]) %>%
    select(1:3) %>%
    calc_lrt_tib()
}
(foo <- tibble(gene1 = gene1_79, gene2 = anchor, lrt = lrt))
```

    ## # A tibble: 79 x 3
    ##    gene1              gene2                lrt
    ##    <chr>              <chr>              <dbl>
    ##  1 ENSMUSG00000009378 ENSMUSG00000024766  8.98
    ##  2 ENSMUSG00000040565 ENSMUSG00000024766 21.4 
    ##  3 ENSMUSG00000041731 ENSMUSG00000024766 40.0 
    ##  4 ENSMUSG00000044026 ENSMUSG00000024766 57.6 
    ##  5 ENSMUSG00000046138 ENSMUSG00000024766 16.5 
    ##  6 ENSMUSG00000046324 ENSMUSG00000024766 31.9 
    ##  7 ENSMUSG00000047298 ENSMUSG00000024766 35.0 
    ##  8 ENSMUSG00000048120 ENSMUSG00000024766 41.6 
    ##  9 ENSMUSG00000048612 ENSMUSG00000024766 52.0 
    ## 10 ENSMUSG00000049670 ENSMUSG00000024766 39.9 
    ## # ... with 69 more rows

``` r
write_csv(x = foo, path = "2018-12-03_lrt-tibble-run3105.csv") 
```

## Run 3106

``` bash
ls ../results/pvl-run3106/*.txt | wc -l
```

    ##     1975

We see that there are 1975 files, as needed.

### Read results files

``` r
res <- "../results/pvl-run3106"
dir(res) -> fns
results <- list()
for (i in 1:1975){
  results[[i]] <- read.table(file.path(res, fns[i]), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(filename = fns[i])
}
```

``` r
results2 <- results %>%
  bind_rows()
```

``` r
foo <- str_split(fns, pattern = "_") %>%
  purrr::map_chr(.f = function(x)x[3])
results2$gene1 <- rep(foo, each = 40000)
```

### “Create 79 tibbles” by gene1 values

First, we need to get the Ensembl id for the anchor gene.

``` r
colnames(readRDS("../data-to-condor/expr_anchors.rds"))[2] -> anchor
```

``` r
gene1_79 <- unique(results2$gene1)
lrt <- numeric(length = 79)
for (i in 1:79){
  lrt[i] <- results2 %>%
    filter(gene1 == gene1_79[i]) %>%
    select(1:3) %>%
    calc_lrt_tib()
}
(foo <- tibble(gene1 = gene1_79, gene2 = anchor, lrt = lrt))
```

    ## # A tibble: 79 x 3
    ##    gene1              gene2                lrt
    ##    <chr>              <chr>              <dbl>
    ##  1 ENSMUSG00000009378 ENSMUSG00000087303  8.50
    ##  2 ENSMUSG00000024766 ENSMUSG00000087303  3.07
    ##  3 ENSMUSG00000040565 ENSMUSG00000087303 20.5 
    ##  4 ENSMUSG00000041731 ENSMUSG00000087303 40.6 
    ##  5 ENSMUSG00000044026 ENSMUSG00000087303 58.3 
    ##  6 ENSMUSG00000046138 ENSMUSG00000087303 17.1 
    ##  7 ENSMUSG00000046324 ENSMUSG00000087303 30.9 
    ##  8 ENSMUSG00000047298 ENSMUSG00000087303 36.0 
    ##  9 ENSMUSG00000048120 ENSMUSG00000087303 39.8 
    ## 10 ENSMUSG00000048612 ENSMUSG00000087303 53.9 
    ## # ... with 69 more rows

``` r
write_csv(x = foo, path = "2018-12-03_lrt-tibble-run3106.csv") 
```

## Run 3107

``` bash
ls ../results/pvl-run3107/*.txt | wc -l
```

    ##     1975

We see that there are 1975 files, as needed.

### Read results files

``` r
res <- "../results/pvl-run3107"
dir(res) -> fns
results <- list()
for (i in 1:1975){
  results[[i]] <- read.table(file.path(res, fns[i]), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(filename = fns[i])
}
```

``` r
results2 <- results %>%
  bind_rows()
```

``` r
foo <- str_split(fns, pattern = "_") %>%
  purrr::map_chr(.f = function(x)x[3])
results2$gene1 <- rep(foo, each = 40000)
```

### “Create 79 tibbles” by gene1 values

First, we need to get the Ensembl id for the anchor gene.

``` r
colnames(readRDS("../data-to-condor/expr_anchors.rds"))[3] -> anchor
```

``` r
gene1_79 <- unique(results2$gene1)
lrt <- numeric(length = 79)
for (i in 1:79){
  lrt[i] <- results2 %>%
    filter(gene1 == gene1_79[i]) %>%
    select(1:3) %>%
    calc_lrt_tib()
}
(foo <- tibble(gene1 = gene1_79, gene2 = anchor, lrt = lrt))
```

    ## # A tibble: 79 x 3
    ##    gene1              gene2                 lrt
    ##    <chr>              <chr>               <dbl>
    ##  1 ENSMUSG00000009378 ENSMUSG00000089952 40.6  
    ##  2 ENSMUSG00000024766 ENSMUSG00000089952 74.2  
    ##  3 ENSMUSG00000040565 ENSMUSG00000089952 51.2  
    ##  4 ENSMUSG00000041731 ENSMUSG00000089952 42.3  
    ##  5 ENSMUSG00000044026 ENSMUSG00000089952 66.8  
    ##  6 ENSMUSG00000046138 ENSMUSG00000089952  0.813
    ##  7 ENSMUSG00000046324 ENSMUSG00000089952 10.1  
    ##  8 ENSMUSG00000047298 ENSMUSG00000089952 15.7  
    ##  9 ENSMUSG00000048120 ENSMUSG00000089952 47.0  
    ## 10 ENSMUSG00000048612 ENSMUSG00000089952 73.1  
    ## # ... with 69 more rows

``` r
write_csv(x = foo, path = "2018-12-03_lrt-tibble-run3107.csv") 
```

``` r
devtools::session_info()
```

    ## ─ Session info ──────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 3.5.1 (2018-07-02)
    ##  os       macOS  10.14.1              
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2018-12-03                  
    ## 
    ## ─ Packages ──────────────────────────────────────────────────────────────
    ##  package     * version    date       lib source                           
    ##  assertthat    0.2.0      2017-04-11 [1] CRAN (R 3.5.0)                   
    ##  backports     1.1.2      2017-12-13 [1] CRAN (R 3.5.0)                   
    ##  base64enc     0.1-3      2015-07-28 [1] CRAN (R 3.5.0)                   
    ##  bindr         0.1.1      2018-03-13 [1] CRAN (R 3.5.0)                   
    ##  bindrcpp    * 0.2.2      2018-03-29 [1] CRAN (R 3.5.0)                   
    ##  broom         0.5.0      2018-07-17 [1] CRAN (R 3.5.0)                   
    ##  callr         3.0.0      2018-08-24 [1] CRAN (R 3.5.0)                   
    ##  cellranger    1.1.0      2016-07-27 [1] CRAN (R 3.5.0)                   
    ##  cli           1.0.1      2018-09-25 [1] CRAN (R 3.5.0)                   
    ##  colorspace    1.3-2      2016-12-14 [1] CRAN (R 3.5.0)                   
    ##  crayon        1.3.4      2017-09-16 [1] CRAN (R 3.5.0)                   
    ##  desc          1.2.0      2018-05-01 [1] CRAN (R 3.5.0)                   
    ##  devtools      2.0.1      2018-10-26 [1] CRAN (R 3.5.1)                   
    ##  digest        0.6.18     2018-10-10 [1] CRAN (R 3.5.0)                   
    ##  dplyr       * 0.7.8      2018-11-10 [1] CRAN (R 3.5.0)                   
    ##  evaluate      0.12       2018-10-09 [1] CRAN (R 3.5.0)                   
    ##  fansi         0.4.0      2018-10-05 [1] CRAN (R 3.5.0)                   
    ##  forcats     * 0.3.0      2018-02-19 [1] CRAN (R 3.5.0)                   
    ##  fs            1.2.6      2018-08-23 [1] CRAN (R 3.5.0)                   
    ##  ggplot2     * 3.1.0      2018-10-25 [1] CRAN (R 3.5.0)                   
    ##  glue          1.3.0      2018-07-17 [1] CRAN (R 3.5.0)                   
    ##  gtable        0.2.0      2016-02-26 [1] CRAN (R 3.5.0)                   
    ##  haven         1.1.2      2018-06-27 [1] CRAN (R 3.5.0)                   
    ##  hms           0.4.2      2018-03-10 [1] CRAN (R 3.5.0)                   
    ##  htmltools     0.3.6      2017-04-28 [1] CRAN (R 3.5.0)                   
    ##  httr          1.3.1      2017-08-20 [1] CRAN (R 3.5.0)                   
    ##  jsonlite      1.5        2017-06-01 [1] CRAN (R 3.5.0)                   
    ##  knitr         1.20       2018-02-20 [1] CRAN (R 3.5.0)                   
    ##  lattice       0.20-35    2017-03-25 [1] CRAN (R 3.5.1)                   
    ##  lazyeval      0.2.1      2017-10-29 [1] CRAN (R 3.5.0)                   
    ##  lubridate     1.7.4      2018-04-11 [1] CRAN (R 3.5.0)                   
    ##  magrittr      1.5        2014-11-22 [1] CRAN (R 3.5.0)                   
    ##  memoise       1.1.0      2017-04-21 [1] CRAN (R 3.5.0)                   
    ##  modelr        0.1.2      2018-05-11 [1] CRAN (R 3.5.0)                   
    ##  munsell       0.5.0      2018-06-12 [1] CRAN (R 3.5.0)                   
    ##  nlme          3.1-137    2018-04-07 [1] CRAN (R 3.5.1)                   
    ##  pillar        1.3.0      2018-07-14 [1] CRAN (R 3.5.0)                   
    ##  pkgbuild      1.0.2      2018-10-16 [1] CRAN (R 3.5.0)                   
    ##  pkgconfig     2.0.2      2018-08-16 [1] CRAN (R 3.5.0)                   
    ##  pkgload       1.0.2      2018-10-29 [1] CRAN (R 3.5.0)                   
    ##  plyr          1.8.4      2016-06-08 [1] CRAN (R 3.5.0)                   
    ##  prettyunits   1.0.2      2015-07-13 [1] CRAN (R 3.5.0)                   
    ##  processx      3.2.0      2018-08-16 [1] CRAN (R 3.5.0)                   
    ##  ps            1.2.1      2018-11-06 [1] CRAN (R 3.5.0)                   
    ##  purrr       * 0.2.5      2018-05-29 [1] CRAN (R 3.5.0)                   
    ##  qtl2pleio   * 0.1.2.9000 2018-11-25 [1] Github (fboehm/qtl2pleio@2cd6f40)
    ##  R6            2.3.0      2018-10-04 [1] CRAN (R 3.5.0)                   
    ##  Rcpp          1.0.0.1    2018-11-18 [1] Github (RcppCore/Rcpp@4f168e6)   
    ##  readr       * 1.2.1      2018-11-22 [1] CRAN (R 3.5.0)                   
    ##  readxl        1.1.0      2018-04-20 [1] CRAN (R 3.5.0)                   
    ##  remotes       2.0.2      2018-10-30 [1] CRAN (R 3.5.0)                   
    ##  rlang         0.3.0.1    2018-10-25 [1] CRAN (R 3.5.0)                   
    ##  rmarkdown     1.10       2018-06-11 [1] CRAN (R 3.5.0)                   
    ##  rprojroot     1.3-2      2018-01-03 [1] CRAN (R 3.5.0)                   
    ##  rstudioapi    0.8        2018-10-02 [1] CRAN (R 3.5.0)                   
    ##  rvest         0.3.2      2016-06-17 [1] CRAN (R 3.5.0)                   
    ##  scales        1.0.0      2018-08-09 [1] CRAN (R 3.5.0)                   
    ##  sessioninfo   1.1.1      2018-11-05 [1] CRAN (R 3.5.0)                   
    ##  stringi       1.2.4      2018-07-20 [1] CRAN (R 3.5.0)                   
    ##  stringr     * 1.3.1      2018-05-10 [1] CRAN (R 3.5.0)                   
    ##  testthat      2.0.1      2018-10-13 [1] CRAN (R 3.5.0)                   
    ##  tibble      * 1.4.2      2018-01-22 [1] CRAN (R 3.5.0)                   
    ##  tidyr       * 0.8.2      2018-10-28 [1] CRAN (R 3.5.0)                   
    ##  tidyselect    0.2.5      2018-10-11 [1] CRAN (R 3.5.0)                   
    ##  tidyverse   * 1.2.1      2017-11-14 [1] CRAN (R 3.5.0)                   
    ##  usethis       1.4.0      2018-08-14 [1] CRAN (R 3.5.0)                   
    ##  utf8          1.1.4      2018-05-24 [1] CRAN (R 3.5.0)                   
    ##  withr         2.1.2      2018-03-15 [1] CRAN (R 3.5.0)                   
    ##  xml2          1.2.0      2018-01-24 [1] CRAN (R 3.5.0)                   
    ##  yaml          2.2.0      2018-07-25 [1] CRAN (R 3.5.0)                   
    ## 
    ## [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
