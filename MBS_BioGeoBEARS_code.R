
################################################################################

# Code originally written and copyrighted by Nicholas J. Matzke (GPL-3 license)
# http://cran.r-project.org/web/licenses/GPL-3 / http://phylo.wikidot.com/biogeobears#script
# with modifications originally by Jen Mandel and Polla Rogrigues 
# (polla.rodrigues16@gmail.com) and finally compiled by Marisa Blake Szubryt (2026-08-27)

# Set working directory and clear environment
setwd('~/PhD/Heterotheca/Code/RStudio_phylogenies_BioGeoBears')
rm(list=ls()); graphics.off(); gc()

# Install packages and load libraries
install.packages('ape')
install.packages('cladoRcpp')
install.packages('devtools')
install.packages('FD')
install.packages('fdrtool')
install.packages('gdata')
install.packages('GenSA')
install.packages('minqa')
install.packages('MultinomialCI')
install.packages('~/GitHub/Matzke_R_binaries/MultinomialCI_1.2.tar.gz', repos = NULL, type = 'source')
install.packages('optimx')
install.packages('pak')
install.packages('parallel')
install.packages('phangorn')
install.packages('phylobase')
install.packages('phytools')
install.packages('plotrix')
install.packages('Rcpp')
install.packages('rexpokit')
install.packages('snow')
install.packages('spam')
install.packages('SparseM')
install.packages('statmod')
devtools::install_github(repo='nmatzke/BioGeoBEARS')
pak::pak(pkg='nmatzke/BioGeoBEARS', upgrade=FALSE, dependencies=FALSE)
library(ape)
library(BioGeoBEARS)
library(cladoRcpp)
library(devtools)
library(FD)
library(fdrtool)
library(gdata)
library(GenSA)
library(minqa)
library(MultinomialCI)
library(optimx)
library(pak)
library(parallel)
library(phangorn)
library(phylobase)
library(phytools)
library(plotrix)
library(Rcpp)
library(rexpokit)
library(snow)
library(spam)
library(SparseM)
library(statmod)

################################################################################
# DEC analysis (Dispersal–Extinction–Cladogenesis)

# Locate extdata directory
extdata_dir = np(system.file('extdata', package='BioGeoBEARS')); print(extdata_dir)
# [1] 'C:\\Users\\maris\\AppData\\Local\\R\\win-library\\4.5\\BioGeoBEARS\\extdata'
list.files(extdata_dir)

# Load in tree file
trfn <- 'lithophyte_tree_calibrated_newick.tre'
tr <- read.tree(trfn); tr <- ladderize(tr, right=FALSE); plot(tr, cex=0.5)
title('Example phylogeny'); axisPhylo(); mtext('Millions of years ago (Ma)', side=1, line=2)

# Read in data file, set tip ranges, and set the max range size
geogfn <- 'Heterotheca_BioGeoBears_CFP_5states.txt'
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn); print(tipranges)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 2

# Initialize general BioGeoBEARS model parameters
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$on_NaN_error <- -1e50
BioGeoBEARS_run_object$speedup <- TRUE
BioGeoBEARS_run_object$use_optimx <- TRUE
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE 

# Run DEC analysis 
runslow <- TRUE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDEC <- bears_optim_run(BioGeoBEARS_run_object)
save(resDEC, file="resDEC.RData")

# Prepare to plot results (add title, DEC output, scripts, states)
analysis_titletxt <- 'BioGeoBEARS DEC'
results_object <- resDEC
scriptdir <- np(system.file('/extdata/a_scripts', package='BioGeoBEARS'))

# Plot results without ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
     addl_params=list('j'), plotwhat='text', label.offset=0.05, tipcex=0.4, 
     statecex=0.5, splitcex=0.00, titlecex=0.8, plotsplits=FALSE, 
     cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Plot results with ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list('j'), 
     plotwhat='pie', label.offset=0.05, tipcex=0.4, statecex=0.5, 
     splitcex=0.0, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, 
     include_null_range=TRUE, tr=tr, tipranges=tipranges)

################################################################################
# DEC+J analysis (Dispersal–Extinction–Cladogenesis with dispersal-mediated speciation)

# Establish d, e, and j (free) parameters
jstart = 0.0001
dstart = resDEC$outputs@params_table['d','est']
estart = resDEC$outputs@params_table['e','est']
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['d','init'] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['d','est'] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['e','init'] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['e','est'] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','type'] = 'free'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','init'] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','est'] = jstart

# Run DEC+J model
runslow <- TRUE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDECj <- bears_optim_run(BioGeoBEARS_run_object)
save(resDECj, file="resDECj.RData")

# Prepare to plot results (add title, DEC output, scripts, states)
analysis_titletxt ='BioGeoBEARS DEC+J'
results_object = resDECj
scriptdir = np(system.file('extdata/a_scripts', package='BioGeoBEARS'))

# Plot results without ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
     addl_params=list('j'), plotwhat='text', label.offset=0.05, tipcex=0.4, 
     statecex=0.5, splitcex=0.00, titlecex=0.8, plotsplits=FALSE, 
     cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Plot results with ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list('j'), 
    plotwhat='pie', label.offset=0.05, tipcex=0.4, statecex=0.5, splitcex=0.0, 
    titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, 
    include_null_range=TRUE, tr=tr, tipranges=tipranges)

################################################################################
# DIVALIKE analysis (like DEC but allows for vicariance, L form of Ronquist (1997))

# Nullify earlier jstart parameter to avoid carry-over from DEC+J to DIVALIKE
jstart = NULL
dstart = NULL
estart = NULL

# Load in tree file
trfn <- 'lithophyte_tree_calibrated_newick.tre'
tr <- read.tree(trfn); tr <- ladderize(tr, right=FALSE); plot(tr, cex=0.5)
title('Example phylogeny'); axisPhylo(); mtext('Millions of years ago (Ma)', side=1, line=2)

# Read in data file, set tip ranges, and set the max range size
geogfn <- 'Heterotheca_BioGeoBears_CFP_5states.txt'
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn); print(tipranges)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 2

# Initialize general BioGeoBEARS model parameters
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$on_NaN_error <- -1e50
BioGeoBEARS_run_object$speedup <- TRUE
BioGeoBEARS_run_object$use_optimx <- TRUE
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE 

# Set up DIVALIKE model parameters ('mx01v' allows for equiprobable vicariance events)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','init'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','est'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['ysv','type'] = '2-j'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['ys','type'] = 'ysv*1/2'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['y','type'] = 'ysv*1/2'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','type'] = 'ysv*1/2'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01v','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01v','init'] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01v','est'] = 0.5

# Set parameters for running DIVALIKE analysis
runslow <- TRUE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDIVALIKE <- bears_optim_run(BioGeoBEARS_run_object)
save(resDIVALIKE, file="resDIVALIKE.RData")

# Plot DIVALIKE analysis (set title, states, etc.)
analysis_titletxt ='BioGeoBEARS DIVALIKE'
results_object = resDIVALIKE
scriptdir = np(system.file('extdata/a_scripts', package='BioGeoBEARS'))

# Plot results without ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
     addl_params=list('j'), plotwhat='text', label.offset=0.05, tipcex=0.4, 
     statecex=0.5, splitcex=0.00, titlecex=0.8, plotsplits=FALSE, 
     cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Plot results with ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list('j'), 
     plotwhat='pie', label.offset=0.05, tipcex=0.4, statecex=0.5, 
     splitcex=0.0, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, 
     include_null_range=TRUE, tr=tr, tipranges=tipranges)

################################################################################
# DIVALIKE+J analysis (DIVALIKE with dispersal-mediated speciation)

# Set parameters for analysis (including jump dispersal/founder-event speciation)
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','type'] = 'free'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','init'] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','est'] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','min'] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','max'] = 1.99999

# Run DIVALIKE+J analysis
runslow <- TRUE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
resDIVALIKEj <- bears_optim_run(BioGeoBEARS_run_object)
save(resDIVALIKEj, file="resDIVALIKEj.RData")

# Plot DIVALIKE+J analysis (set title, states, etc.)
analysis_titletxt ='BioGeoBEARS DIVALIKE+J'
results_object = resDIVALIKEj
scriptdir = np(system.file('extdata/a_scripts', package='BioGeoBEARS'))

# Plot results without ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
     addl_params=list('j'), plotwhat='text', label.offset=0.05, tipcex=0.4, 
     statecex=0.5, splitcex=0.00, titlecex=0.8, plotsplits=FALSE, 
     cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Plot results with ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list('j'), 
     plotwhat='pie', label.offset=0.05, tipcex=0.4, statecex=0.5, splitcex=0.0, 
     titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, 
     include_null_range=TRUE, tr=tr, tipranges=tipranges)

################################################################################
# BAYAREALIKE analysis (simplified version derived from BayArea from Landis et al. (2013))
# Includes cladogenesis assumption unlike DEC and DIVALIKE

# Nullify earlier jstart parameter to avoid carry-over from DIVALIKE+J to BAYAREALIKE
jstart = NULL

# Load in tree file
trfn <- 'lithophyte_tree_calibrated_newick.tre'
tr <- read.tree(trfn); tr <- ladderize(tr, right=FALSE); plot(tr, cex=0.5)
title('Example phylogeny'); axisPhylo(); mtext('Millions of years ago (Ma)', side=1, line=2)

# Read in data file, set tip ranges, and set the max range size
geogfn <- 'Heterotheca_BioGeoBears_CFP_5states.txt'
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn); print(tipranges)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 2

# Initialize general BioGeoBEARS model parameters
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$on_NaN_error <- -1e50
BioGeoBEARS_run_object$speedup <- TRUE
BioGeoBEARS_run_object$use_optimx <- TRUE
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE 

# Set BAYAREALIKE model (no subset sympatry, vicariance, or dispersal-derived speciation)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','init'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','est'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','init'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','est'] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['ysv','type'] = '1-j'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['ys','type'] = 'ysv*1/1'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['y','type'] = '1-j'

# Allow for only sympatric/range-copying (y) events
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01y','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01y','init'] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01y','est'] = 0.9999

# Run BAYAREALIKE analysis
runslow <- TRUE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
resBAYAREALIKE <- bears_optim_run(BioGeoBEARS_run_object)
save(resBAYAREALIKE, file="resBAYAREALIKE.RData")

# Set parameters for plotting BAYAREALIKE model
analysis_titletxt ='BioGeoBEARS BAYAREALIKE'
results_object = resBAYAREALIKE
scriptdir = np(system.file('extdata/a_scripts', package='BioGeoBEARS'))

# Plot results without ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
     addl_params=list('j'), plotwhat='text', label.offset=0.05, tipcex=0.4, 
     statecex=0.5, splitcex=0.00, titlecex=0.8, plotsplits=FALSE, 
     cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Plot results with ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list('j'), 
     plotwhat='pie', label.offset=0.05, tipcex=0.4, statecex=0.5, splitcex=0.0, 
     titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, 
     include_null_range=TRUE, tr=tr, tipranges=tipranges)

################################################################################
# BAYAREALIKE+J analysis (BAYAREALIKE with dispersal-mediated speciation)
# Includes cladogenesis assumption unlike DEC and DIVALIKE

# Set BAYAREALIKE+J model parameters
jstart = 0.0001
dstart = resBAYAREALIKE$outputs@params_table['d','est']
estart = resBAYAREALIKE$outputs@params_table['e','est']
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['d','init'] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['d','est'] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['e','init'] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['e','est'] = estart

# Do not permit subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','init'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['s','est'] = 0.0

# Do not permit vicariance events
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','init'] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['v','est'] = 0.0

# Do allow for jump dispersal speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','type'] = 'free'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','init'] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','est'] = jstart

# Set the max of 'j' as 1 rather than 3 for BAYAREALIKE+J
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','max'] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['ysv','type'] = '1-j'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['ys','type'] = 'ysv*1/1'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['y','type'] = '1-j'

# Allow only sympatric events  where both descendants always have the same size as the ancestor
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01y','type'] = 'fixed'
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01y','init'] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['mx01y','est'] = 0.9999

# Set parameters to avoid program crashing on Windows machines
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['d','min'] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['d','max'] = 4.9999999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['e','min'] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['e','max'] = 4.9999999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','min'] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table['j','max'] = 0.99999

# Run BAYAREALIKE+J analysis
runslow <- TRUE
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
resBAYAREALIKEj <- bears_optim_run(BioGeoBEARS_run_object)
save(resBAYAREALIKEj, file="resBAYAREALIKEj.RData")

# Plot BAYAREALIKE+J states
analysis_titletxt ='BioGeoBEARS BAYAREALIKE+J'
results_object = resBAYAREALIKEj
scriptdir = np(system.file('extdata/a_scripts', package='BioGeoBEARS'))

# Plot results without ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, 
     addl_params=list('j'), plotwhat='text', label.offset=0.05, tipcex=0.4, 
     statecex=0.5, splitcex=0.00, titlecex=0.8, plotsplits=FALSE, 
     cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Plot results with ancestral pie charts
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list('j'), 
     plotwhat='pie', label.offset=0.05, tipcex=0.4, statecex=0.5, splitcex=0.0, 
     titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, 
     include_null_range=TRUE, tr=tr, tipranges=tipranges)

################################################################################
# Compare summary statistics across models (DEC, DIVALIKE, BAYAREALIKE, w/ or w/o +j)
# helpful reading: http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears

# Generate empty tables to hold results
restable = NULL
teststable = NULL

# Extract DEC and DECj statistics
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
numparams1 = 3; numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2); stats

# DEC = null model for Likelihood Ratio Test (LRT) / DEC+J = alternative model for LRT
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat='table', addl_params=c('j'), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat='table', addl_params=c('j'), paramsstr_digits=4)

# H0 for LRT: two models confer the same likelihood on the data (http://www.brianomeara.info/tutorials/aic)
# Note that LRT p-values indicate whether or not the null should be rejected
# (null = DEC and DEC+J confer equal likelihoods); a low p-value means they are different
rbind(res2, res1); tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1); print(restable)
teststable = rbind(teststable, tmp_tests); print(teststable)

################################################################################
# Extract DIVALIKE and DIVALIKEj statistics
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
numparams1 = 3; numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2); stats

# DIVALIKE = null model for LRT  / DIVALIKE+J = alternative model for LRT
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat='table', addl_params=c('j'), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat='table', addl_params=c('j'), paramsstr_digits=4)

# H0 for LRT: two models confer the same likelihood on the data (http://www.brianomeara.info/tutorials/aic)
# Note that LRT p-values indicate whether or not the null should be rejected
# (null = DIVALIKE and DIVALIKE+J confer equal likelihoods); a low p-value means they are different
rbind(res2, res1)
conditional_format_table(stats)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1); print(restable)
teststable = rbind(teststable, tmp_tests); print(teststable)

################################################################################
# Extract BAYAREALIKE and BAYAREALIKEj statistics
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
numparams1 = 3; numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2); stats

# BAYAREALIKE = null model for LRT  / BAYAREALIKE+J = alternative model for LRT
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat='table', addl_params=c('j'), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat='table', addl_params=c('j'), paramsstr_digits=4)

# H0 for LRT: two models confer the same likelihood on the data (http://www.brianomeara.info/tutorials/aic)
# Note that LRT p-values indicate whether or not the null should be rejected
# (null = BAYAREALIKE and BAYAREALIKE+J confer equal likelihoods); a low p-value means they are different
rbind(res2, res1)
conditional_format_table(stats)
tmp_tests = conditional_format_table(stats)
restable = rbind(restable, res2, res1); print(restable)
teststable = rbind(teststable, tmp_tests); print(teststable)

################################################################################
# Assemble results table
teststable$alt = c('DEC+J', 'DIVALIKE+J', 'BAYAREALIKE+J')
teststable$null = c('DEC', 'DIVALIKE', 'BAYAREALIKE')
row.names(restable) = c('DEC', 'DEC+J', 'DIVALIKE', 'DIVALIKE+J', 'BAYAREALIKE', 'BAYAREALIKE+J')
restable = put_jcol_after_ecol(restable); print(restable)
print(teststable)

# Save restable, teststable, and text files
save(restable, file='restable_v1.Rdata'); load(file='restable_v1.Rdata')
save(teststable, file='teststable_v1.Rdata'); load(file='teststable_v1.Rdata')
write.table(restable, file='restable.txt', quote=FALSE, sep='\t')
write.table(unlist_df(teststable), file='teststable.txt', quote=FALSE, sep='\t')

# Evaluate model weights for all models (for AIC and AICc comparisons)
restable2 = restable

# Evaluate based on AIC
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use='AIC')
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike); restable_AIC_rellike

# Evaluate based on AICc (factors in sample sizes)
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use='AICc')
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike); restable_AICc_rellike

# Combine AIC and AICc tables
restable_AIC_rellike$AICc <- restable_AICc_rellike$AICc
restable_AIC_rellike$AICc_wt <- restable_AICc_rellike$AICc_wt; print(restable_AIC_rellike)

# Save AIC and AICc comparison tables
write.table(restable_AIC_rellike, file='restable_AIC-AICc_rellike.txt', quote=FALSE, sep='\t')
write.table(conditional_format_table(restable_AIC_rellike), file='restable_AIC-AICc_rellike_formatted.txt', quote=FALSE, sep='\t')

# The optimal model should have the lowest AIC and/or AICc values/scores 
# and higher (or less negative) natural log likelihood values (if the number of parameters are the same)
