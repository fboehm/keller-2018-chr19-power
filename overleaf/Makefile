# Latex Makefile using latexmk
# Modified by Dogukan Cagatay <dcagatay@gmail.com>
# Originally from : http://tex.stackexchange.com/a/40759
#
# Change only the variable below to the name of the main tex file.
PROJNAME=main
LATEXMK=latexmk
RMD_DIR=../Rmd
RMD_FILE=2018-12-04_creating-figures.Rmd
HTML_FILE=2018-12-04_creating-figures.html
# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.
.PHONY: all

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: $(PROJNAME).pdf $(RMD_DIR)/$(HTML_FILE)

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# "raw2tex" and "dat2tex" are just placeholders for whatever custom steps
# you might have.

#%.tex: %.raw
#	./raw2tex $< > $@

#%.tex: %.dat
#	./dat2tex $< > $@

# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interactive=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

$(RMD_DIR)/lrt-v-middle-of-gene.pdf: $(RMD_DIR)/$(RMD_FILE) $(RMD_DIR)/2018-11-26_lrt-tibble.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3105.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3106.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3107.csv ../data/Attie_DO378_eQTL_viewer_v1.Rdata ../data/chr19_expr_10mb.rds
	cd $(RMD_DIR); Rscript -e 'rmarkdown::render("$(RMD_FILE)")'	

$(RMD_DIR)/lrt-v-univariate-lod.pdf: $(RMD_DIR)/$(RMD_FILE) $(RMD_DIR)/2018-11-26_lrt-tibble.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3105.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3106.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3107.csv ../data/Attie_DO378_eQTL_viewer_v1.Rdata ../data/chr19_expr_10mb.rds
	cd $(RMD_DIR); Rscript -e 'rmarkdown::render("$(RMD_FILE)")'	

$(RMD_DIR)/lrt-v-corr.pdf: $(RMD_DIR)/$(RMD_FILE) $(RMD_DIR)/2018-11-26_lrt-tibble.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3105.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3106.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3107.csv ../data/Attie_DO378_eQTL_viewer_v1.Rdata ../data/chr19_expr_10mb.rds
	cd $(RMD_DIR); Rscript -e 'rmarkdown::render("$(RMD_FILE)")'



$(RMD_DIR)/lrt-v-middle-of-gene.eps: $(RMD_DIR)/$(RMD_FILE) $(RMD_DIR)/2018-11-26_lrt-tibble.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3105.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3106.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3107.csv ../data/Attie_DO378_eQTL_viewer_v1.Rdata ../data/chr19_expr_10mb.rds
	Rscript -e 'rmarkdown::render("$<")'	

$(RMD_DIR)/lrt-v-univariate-lod.eps: $(RMD_DIR)/$(RMD_FILE) $(RMD_DIR)/2018-11-26_lrt-tibble.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3105.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3106.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3107.csv ../data/Attie_DO378_eQTL_viewer_v1.Rdata ../data/chr19_expr_10mb.rds
	Rscript -e 'rmarkdown::render("$<")'	

$(RMD_DIR)/lrt-v-corr.eps: $(RMD_DIR)/$(RMD_FILE) $(RMD_DIR)/2018-11-26_lrt-tibble.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3105.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3106.csv $(RMD_DIR)/2018-12-03_lrt-tibble-run3107.csv ../data/Attie_DO378_eQTL_viewer_v1.Rdata ../data/chr19_expr_10mb.rds
	Rscript -e 'rmarkdown::render("$<")'

$(PROJNAME).pdf: $(PROJNAME).tex 
	latexmk -pdf -use-make -pvc $<


$(RMD_DIR)/$(HTML_FILE): $(RMD_DIR)/$(RMD_FILE)
	Rscript -e 'rmarkdown::render("$<")'


cleanall:
	latexmk -C

clean:
	latexmk -c
