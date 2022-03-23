library("sbcdata")

## parallelisation on HPC
library("future")
library("future.batchtools")

## helper functions
source("00_rus_ranger_functions.R")

## parallelisation plan
login <- future::tweak(
    future::cluster,
    workers = "brain", ## HPC UMG
    rscript = "/opt/r-4.1.1/bin/Rscript",
    homogeneous = FALSE
)

resources <- list(
    workers = 32L,
    partition = c("batch", "snowball"),
    ncpu = 36,
    mcpu = 2L * 1024L,
    walltime = 48L * 60L * 60L # seconds
)

slurm <- future::tweak(
    future.batchtools::batchtools_slurm,
    template = "slurm_batchtools.tmpl",
    workers = resources$workers,
    resources = resources
)

login2 <- future::tweak(multicore, workers = 5)
node <- sequential

hostname <- Sys.info()[["nodename"]]

if (hostname == "ape") {
    future::plan(list(login, slurm, node))
} else if (grepl("login", hostname)) {
    future::plan(list(slurm, node))
} else {
    ## ranger automatically uses all local cores that's why we don't want
    ## additional parallelisation
    future::plan(sequential)
}
## END parallelisation plan

ukl <- exclude_entries(subset(sbcdata, Center == "Leipzig"))
ukl <- ukl[!ukl$Excluded,]

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

searchspace <- expand.grid(
    mtry = c(2, 3, 4),
    num.trees = c(500, 750, 1000, 1250)
)

#searchspace <- expand.grid(
#    mtry = c(2, 3),
#    num.trees = c(500, 750)
#)

## nested cross validation for hyperparam search
nouterfolds = 5
ninnerfolds = 5 
## number of repeated cross validation in inner loop
nrepcv = 10

set.seed(20220317)
res <- nrcv_ranger(
    x = ukl[, xvar, with = FALSE], y = ifelse(ukl$Diagnosis == "Sepsis", 1, 2),
    searchspace = searchspace,
    nouterfolds = nouterfolds, ninnerfolds = ninnerfolds, nrepcv = nrepcv
)
saveRDS(
    res,
    file = paste0("res-", nouterfolds, "-", ninnerfolds, "-", nrepcv, ".RDS")
)
