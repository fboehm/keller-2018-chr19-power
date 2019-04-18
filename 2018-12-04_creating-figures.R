## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
library(tidyverse)


## ------------------------------------------------------------------------
read_csv("2018-11-26_lrt-tibble.csv") -> r3104
read_csv("2018-12-03_lrt-tibble-run3105.csv") -> r3105
read_csv("2018-12-03_lrt-tibble-run3106.csv") -> r3106
read_csv("2018-12-03_lrt-tibble-run3107.csv") -> r3107
(bind_rows(r3104, r3105, r3106, r3107) -> rr)


## ------------------------------------------------------------------------
load("../data/Attie_DO378_eQTL_viewer_v1.Rdata")


## ------------------------------------------------------------------------
foo <- dataset.islet.rnaseq$annots %>%
  select(gene_id, symbol, start, end, middle) 
rr2 <- foo %>% 
  right_join(rr, by = c("gene_id" = "gene1")) %>%
  rename(gene1_symbol = symbol, gene1 = gene_id, gene1_start = start, gene1_end = end, gene1_middle = middle) %>%
  left_join(foo, by = c("gene2" = "gene_id")) %>%
  rename(gene2_symbol = symbol, gene2_start = start, gene2_end = end, gene2_middle = middle)


## ------------------------------------------------------------------------
bar <- dataset.islet.rnaseq$lod.peaks %>%
  filter(chrom == 19)
rr3 <- bar %>%
  right_join(rr2, by = c("annot.id" = "gene1")) %>%
  rename(gene1 = annot.id, gene1_marker.id = marker.id, gene1_pos = pos, gene1_lod = lod) %>%
  left_join(bar, by = c("gene2" = "annot.id")) %>%
  rename(gene2_marker.id = marker.id, gene2_pos = pos, gene2_lod = lod) %>%
  select(- chrom.x, - chrom.y)


## ------------------------------------------------------------------------
library(plotly)


## ------------------------------------------------------------------------
(pp <- rr3 %>%
  ggplot() + geom_point(aes(x = gene1_middle, y = lrt, colour = gene1_lod)) + geom_vline(aes(xintercept = gene2_middle)) + xlab("Chromosome 19 position") + ylab("Pleiotropy v separate QTL test statistic") + facet_grid(rows = gene2_symbol ~ .)
)



## ------------------------------------------------------------------------
ggplotly(pp)


## ------------------------------------------------------------------------
(nodup <- rr3 %>%
  filter(!duplicated(gene2)) %>%
  arrange(desc(gene2_lod)))


## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  ggplot() + geom_point(aes(x = gene1_middle, y = lrt, colour = gene1_lod)) + geom_vline(aes(xintercept = gene2_middle)) + xlab("Chromosome 19 position") + ylab("Pleiotropy v separate QTL test statistic") + facet_grid(rows = gene2_symbol_factor ~ .)
ggsave(filename = "lrt-v-middle-of-gene.jpg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-middle-of-gene.svg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-middle-of-gene.eps", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-middle-of-gene.pdf", width = 7, height = 7, units = "in" )



## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  ggplot() + geom_point(aes(x = gene1_middle, y = lrt, colour = gene1_lod)) + geom_vline(aes(xintercept = gene2_middle)) + xlab("Chromosome 19 position") + ylab("Pleiotropy v separate QTL test statistic") + facet_grid(rows = gene2_symbol_factor ~ .) + scale_color_gradient(low = "#500707", high = "#ee0909", name = "Trait LOD") 

ggsave(filename = "lrt-v-middle-of-gene-slides.svg", width = 9, height = 7, units = "in" )


## ------------------------------------------------------------------------
(pp <- rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>% 
  filter(gene2_symbol == "Asah2") %>%
  ggplot() + geom_point(aes(x = gene1_middle, y = lrt, colour = gene1_lod, text = paste("gene1_symbol: ", gene1_symbol))) +
  ggtitle("Asah2") + 
  geom_vline(aes(xintercept = gene2_middle)) + 
  xlab("Chromosome 19 position") + 
  ylab("Pleiotropy v separate QTL test statistic") + 
  scale_color_gradient(low = "#500707", high = "#ee0909", name = "Trait LOD")
)
ggplotly(pp)
ggsave(filename = "lrt-v-middle-of-gene-slides-1.svg", width = 9, height = 7, units = "in" )


## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  mutate(interlocus_distance = abs(gene1_middle - gene2_middle)) %>%
  ggplot() + geom_point(aes(x = gene1_lod, y = lrt, color = interlocus_distance)) + 
  xlab("Univariate LOD") + 
  ylab("Pleiotropy v separate QTL test statistic") + 
  scale_color_gradient(name = "Interlocus distance") + 
  facet_grid(rows = gene2_symbol_factor ~ .)
ggsave(filename = "lrt-v-univariate-lod.jpg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-univariate-lod.svg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-univariate-lod-slides.svg", width = 9, height = 7, units = "in" )

ggsave(filename = "lrt-v-univariate-lod.eps", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-univariate-lod.pdf", width = 7, height = 7, units = "in" )



## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  mutate(interlocus_distance = abs(gene1_middle - gene2_middle)) %>%
  filter(gene2_symbol == "Asah2") %>%
  ggplot() + geom_point(aes(x = gene1_lod, y = lrt, color = interlocus_distance)) + 
  xlab("Univariate LOD") + 
  ylab("Pleiotropy v separate QTL test statistic") + 
  scale_color_gradient(name = "Interlocus distance") + 
  ggtitle("Asah2")

ggsave(filename = "lrt-v-univariate-lod-slides-1.svg", width = 9, height = 7, units = "in" )



## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  ggplot() + 
  geom_point(aes(x = gene1_lod, y = lrt)) +
  xlab("Univariate LOD") + 
  ylab("Pleiotropy v separate QTL test statistic") + 
  geom_smooth(aes(x = gene1_lod, y = lrt), span = 1) + 
  facet_grid(rows = gene2_symbol_factor ~ .)
ggsave(filename = "lrt-v-univariate-lod-loess.jpg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-univariate-lod-loess.svg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-univariate-lod-loess.eps", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-univariate-lod-loess.pdf", width = 7, height = 7, units = "in" )



## ------------------------------------------------------------------------
library(qtl2)


## ------------------------------------------------------------------------
dimnames(genoprobs$`19`)[[3]]


## ------------------------------------------------------------------------
readRDS("../data/chr19_expr_10mb.rds") -> expr
colnames(expr)


## ------------------------------------------------------------------------
ann80 <- rr3 %>%
  select(gene1, gene1_marker.id, gene1_symbol) %>%
  filter(!duplicated(gene1_symbol))
fit1_out <- list()
for (i in 1:80){
  pr <- genoprobs$`19`[ , , which(dimnames(genoprobs$`19`)[[3]] == ann80$gene1_marker.id[i])]
  fit1_out[[i]] <- fit1(genoprobs = pr, 
       pheno = expr[, which(colnames(expr) == ann80$gene1[i]), drop = FALSE], 
       kinship = K$`19`, 
       addcovar = dataset.islet.rnaseq$covar, 
       reml = TRUE
       )
}


## ------------------------------------------------------------------------
sapply(FUN = function(x)x$fitted, X = fit1_out) -> fits
colnames(fits)
colnames(fits) <- ann80$gene1
rownames(fits)


## ------------------------------------------------------------------------
cor(fits) -> fitcor


## ------------------------------------------------------------------------
rr3$fit1_corr <- NA
for (i in 1:nrow(rr3)){
  rr3$fit1_corr[i] <- fitcor[which(rownames(fitcor) == rr3$gene1[i]), which(colnames(fitcor) == rr3$gene2[i])]
}


## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  ggplot() + 
  guides(color = guide_colourbar(title = "Univariate LOD")) + 
  #guides(colour = guide_legend(title = "Univariate LOD")) + 
  geom_point(aes(x = abs(fit1_corr), y = lrt, colour = gene1_lod)) + 
#  geom_smooth(aes(x = abs(fit1_corr), y = lrt), span = 1) + 
  xlab("Fitted values (absolute) correlation") + 
  ylab("Pleiotropy v separate QTL test statistic") + 
  facet_grid(rows = gene2_symbol_factor ~ .)
ggsave(filename = "lrt-v-corr.jpg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-corr.svg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-corr.eps", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-corr.pdf", width = 7, height = 7, units = "in" )



## ------------------------------------------------------------------------
rr3 %>%
  mutate(gene2_symbol_factor = factor(gene2_symbol, levels = nodup$gene2_symbol, labels = paste0(nodup$gene2_symbol, " (", round(nodup$gene2_lod, 1), ")"))) %>%
  ggplot() + 
  guides(color = guide_colourbar(title = "Univariate LOD")) + 
  #guides(colour = guide_legend(title = "Univariate LOD")) + 
  geom_point(aes(x = abs(fit1_corr), y = lrt, colour = gene1_lod)) + 
  geom_smooth(aes(x = abs(fit1_corr), y = lrt), span = 1) + 
  xlab("Fitted values (absolute) correlation") + 
  ylab("Pleiotropy v separate QTL test statistic") + 
  facet_grid(rows = gene2_symbol_factor ~ .)
ggsave(filename = "lrt-v-corr-loess.jpg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-corr-loess.svg", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-corr-loess.eps", width = 7, height = 7, units = "in" )
ggsave(filename = "lrt-v-corr-loess.pdf", width = 7, height = 7, units = "in" )


## ------------------------------------------------------------------------
unique(rr2$gene2) -> id4
unique(rr2$gene2_symbol) -> id4_symbols
expr4 <- expr[ , colnames(expr) %in% id4]


## ------------------------------------------------------------------------
out4 <- list()
for (i in 1:4){
  out4[[i]] <- scan1coef(genoprobs = list(genoprobs$`19`), 
            pheno = expr4[ , i, drop = FALSE], 
            kinship = K$`19`, 
            addcovar = dataset.islet.rnaseq$covar, 
            reml = TRUE,
            se = FALSE)
}
out4se <- list()
for (i in 1:4){
  foobar <- scan1coef(genoprobs = list(genoprobs$`19`), 
            pheno = expr4[ , i, drop = FALSE], 
            kinship = K$`19`, 
            addcovar = dataset.islet.rnaseq$covar, 
            reml = TRUE,
            se = TRUE)
  out4se[[i]] <- foobar
}
names(out4se) <- colnames(expr4)
names(out4) <- colnames(expr4)
  


## ------------------------------------------------------------------------
for (i in 1:4){
  if (i == 1){
    leg <- "bottomright"
  } else {
    leg <- NULL
  }
  uni_peak_pos <- nodup %>%
    filter(gene2_symbol == id4_symbols[i]) %>%
    select(gene2_pos) %>%
    unlist()
  epsfn <- paste0("allele_effects_", id4_symbols[i], ".eps")
  pdffn <- paste0("allele_effects_", id4_symbols[i], ".pdf")
  setEPS()
  postscript(epsfn, 
             width = 6, 
             height = 4
             )
  plot_coefCC(out4[[i]][400:1399, ], 
              map = map[19], 
              legend = leg, 
              ylim = c(- 2.5, 2.5), 
              legend_ncol = 4
              )
  text(x = 50, 
       y = 2, 
       labels = id4_symbols[i]
       )
  #abline(v = uni_peak_pos)
  dev.off()
  pdf(pdffn, 
      width = 6, 
      height = 4
      )
  cex_all <- 0.75
  par(cex = cex_all)
  plot_coefCC(out4[[i]][400:1399, ], 
              map = map[19], 
              legend = leg, 
              ylim = c(- 2.5, 3), 
              legend_ncol = 4, 
              )
  text(x = 50, 
       y = 2, 
       labels = id4_symbols[i],
       cex = 1.5 / cex_all
       )
  #abline(v = uni_peak_pos)
  dev.off()
}


## ------------------------------------------------------------------------
#se4 <- list()
#for (i in 1:4){
#  se4[[i]] <- pluck(out4se, i, attr_getter("SE")) # note that we need capital letters for SE
#} # of course, we could use the full power of purrr (and its map functions) to do this loop in one line, I think...
se4 <- out4se %>%
  map(pluck(attr_getter("SE")))


## ------------------------------------------------------------------------
at1 <- attributes(out4se[[1]])
names(at1)
at1[[4]] # direct from attributes of out4se
at1[[4]] == se4[[1]]


## ------------------------------------------------------------------------
library(qtl2pleio)
library(xtable)


## ------------------------------------------------------------------------
# convert map for Chr 19 to a tibble for subsequent merging
map19_tib <- tibble(marker_position = map$`19`, marker_index = 1:length(map$`19`))
#
ann4 <- ann80 %>%
  filter(gene1 %in% colnames(expr4)) %>%
  left_join(dataset.islet.rnaseq$lod.peaks, by = c("gene1" = "annot.id")) %>%
  filter(chrom == 19) %>%
  left_join(map19_tib, by = c("pos" = "marker_position")) %>%
  rename(peak_position = pos, peak_marker_index = marker_index) %>% arrange(gene1)


## ------------------------------------------------------------------------
i <- 1
# be sure that the two objects - out4se (and se4) have the same ordering of genes as ann4!!
# 
get_effects_table <- function(peak_index, ae_obj, se_obj, map)
{
  ae <- get_effects(marker_index = peak_index, 
                  allele_effects_matrix = ae_obj, 
                  map = map
                  )
  se <- get_effects(marker_index = peak_index, 
                  allele_effects_matrix = se_obj, 
                  map = map
                  )
  mytib <- tibble(founder_allele = LETTERS[1:8], ae, se)
  return(mytib)
}


## ------------------------------------------------------------------------
get_effects_table(peak_index = ann4$peak_marker_index[1], ae_obj = out4[[1]], se_obj = se4[[1]], map = map$`19`)%>%
  xtable(caption = ann4$gene1_symbol[1]) %>%
  print.xtable()


## ------------------------------------------------------------------------
get_effects_table(peak_index = ann4$peak_marker_index[2], ae_obj = out4[[2]], se_obj = se4[[2]], map = map$`19`) %>%
  xtable(caption = ann4$gene1_symbol[2]) %>%
  print.xtable()


## ------------------------------------------------------------------------
get_effects_table(peak_index = ann4$peak_marker_index[3], ae_obj = out4[[3]], se_obj = se4[[3]], map = map$`19`) %>%
  xtable(caption = ann4$gene1_symbol[3]) %>% 
  print.xtable()


## ------------------------------------------------------------------------
get_effects_table(peak_index = ann4$peak_marker_index[4], ae_obj = out4[[4]], se_obj = se4[[4]], map = map$`19`) %>% 
  xtable(caption = ann4$gene1_symbol[4]) %>% 
  print.xtable()


## ---- results = "asis"---------------------------------------------------
nodup %>%
  select(gene2_symbol, gene2_start, gene2_end, gene2_pos, gene2_lod) %>%
  rename(symbol = gene2_symbol, start = gene2_start, end = gene2_end, peak_position = gene2_pos, lod = gene2_lod) %>% 
  xtable() %>%
  print.xtable(include.rownames = FALSE)


## ------------------------------------------------------------------------
rr3 %>%
  select(gene1_symbol, gene1_start, gene1_end, gene1_pos, gene1_lod) %>%
  rename(symbol = gene1_symbol, start = gene1_start, end = gene1_end, peak_position = gene1_pos, lod = gene1_lod) %>%
  filter(!duplicated(symbol)) %>%
  filter(!(symbol %in% c("Asah2", "Lipo1", "Lipo2", "4933413C19Rik"))) %>%
  arrange(desc(lod))%>% 
  xtable() %>%
  print.xtable(include.rownames = FALSE, size = "\\tiny")

