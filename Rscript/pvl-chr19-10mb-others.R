## GOAL: analyze - with pleiotropy v separate QTL tests - the 80 traits all paired with one of the traits that
# is NOT Asah2. Note that the code below will run analyses with each of the lead traits paired with itself, but that's ok


##First read in the arguments listed at the command line

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
print(args)
##args is now a list of character vectors
print(args$argname)
proc_num <- as.numeric(args$argname)
print(proc_num)
#trait_indic <- proc_num + 1 # define trait_indic 
trait_indic <- (proc_num %% 80)  + 1 # 80 is the number of columns in expr matrix
lead_indic <- (proc_num %/% 80) + 1
run_num <- as.numeric(args$run_num)
print(run_num)
#(nsnp <- as.numeric(args$nsnp))
#(s1 <- as.numeric(args$s1))

###############

# define traits to examine per Rmd/2018-10-30_allele-effects-chr19.Rmd
# ENSMUSG00000087303 and ENSMUSG00000024766 are the least correlated with the others, while ENSMUSG00000089952

lead_traits <- c("ENSMUSG00000087303", "ENSMUSG00000024766", "ENSMUSG00000089952")
lead_id <- lead_traits[lead_indic]

library(dplyr)
# load expression traits
readRDS("data-to-condor/chr19_expr_10mb.rds") -> locals
trait_id <- colnames(locals)[trait_indic]
# define phenames 
phenames <- c(trait_id, lead_id)

# load chr2 allele probabilities
readRDS("genoprobs_chr19-keller.rds") -> geno # genoprobs_chr19.rds is on SQUID

# load kinship matrix (LOCO, ie, for chromosome 2, ie, doesn't use chr2 data-to-condor)
readRDS("data-to-condor/kinship-chr19-keller.rds") -> kinship

# load covariates
readRDS("data-to-condor/addcovar-keller.rds") -> covar

# load annot
annot <- readRDS("data-to-condor/annot.rds")

# create matrix of two expression traits
pheno <- cbind(locals[ , trait_indic, drop = FALSE], locals[ , colnames(locals) %in% lead_id, drop = FALSE])
rownames(pheno) <- rownames(locals)
# get only shared subjects with no missing pheno or covariates

# remove subjects with missing data

id2keep <- rownames(locals)
gg <- geno
gg2 <- gg[rownames(gg) %in% id2keep, , , drop = FALSE]
kk <- kinship
kk2 <- kk[rownames(kk) %in% id2keep, colnames(kk) %in% id2keep, drop = FALSE]
cc2 <- covar[rownames(covar) %in% id2keep, , drop = FALSE]


# define s1 and nsnp
lead_index <- annot %>%
  filter(gene_id == lead_id) %>%
  select(marker_index) %>%
  unlist()
trait_index <- annot %>%
  filter(gene_id == trait_id) %>%
  select(marker_index) %>%
  unlist()
s1 <- min(lead_index, trait_index) - 15
# nsnp <- abs(lead_index - trait_index) + 1 + 15 + 15
nsnp <- 5 # for trial runs only!

# two-dimensional scan
library(qtl2pleio)
s_out <- scan_pvl(probs = gg2,
                  pheno = pheno,
                  kinship = kk2,
#                  addcovar = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
                  covariates = cc2[ , -5], # need to remove column 5 because we have no mice from wave 5
                  start_snp = s1,
                  n_snp = nsnp, max_iter = 10000
)


fn_out <- paste0("pvl-run", run_num, "_", proc_num, "_", paste(phenames, collapse = "_"), ".txt")
write.table(s_out, fn_out, quote = FALSE)
devtools::session_info()
q("no")
