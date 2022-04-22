library("sbcdata")
library("rusranger")

## don't exclude by time
ukl <- exclude_entries(
    subset(sbcdata, Center == "Leipzig"), time = c(-Inf, Inf)
)
ukl <- ukl[!ukl$Excluded,]
umg <- exclude_entries(
    subset(sbcdata, Center == "Greifswald"), time = c(-Inf, Inf)
)
umg <- umg[!umg$Excluded,]
mimic <- exclude_entries(
    import_mimic("../../inst/intdata/mimic-iv-1.0/"), time = c(-Inf, Inf)
)
mimic <- mimic[!mimic$Excluded,]

ukl$Sex <- as.integer(ukl$Sex == "M")
umg$Sex <- as.integer(umg$Sex == "M")
mimic$Sex <- as.integer(mimic$Sex == "M")

ukl$Diagnosis <- as.integer(ukl$Diagnosis == "Sepsis")
umg$Diagnosis <- as.integer(umg$Diagnosis == "Sepsis")
mimic$Diagnosis <- as.integer(mimic$Diagnosis == "Sepsis")

set.seed(20220419)

train <- subset(ukl, Set == "Training")
validation <- list(
    UKL = subset(ukl, Set == "Validation"), UMG = umg, MIMIC = mimic
)

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

mtry <- 3
num.trees <- 1000
replace <- TRUE # doing bootstrap, shorter runtime, b/c much lower samplesize/samplefraction
scw <- 10L

caseweights <- function(y, y0, scw = 1) {
    isControl <- !y
    isSepsis <- as.logical(y)
    isSepsisControl <- as.logical(y0) & !as.logical(y)
    ## control = 1
    w <- rep_len(1, length(y))
    ## sepsis cases that are now control
    w[isSepsisControl] <- scw
    ## sepsis cases
    w[isSepsis] <- (sum(isControl) + sum(isSepsisControl) * scw) / sum(isSepsis)
    w
}

cwranger <- function(d, scw = 1) {
    ranger(
        x = as.data.frame(d[, xvar, with = FALSE]), y = d$Diagnosis,
        probability = TRUE, min.node.size = 10,
        mtry = mtry, num.trees = num.trees,
        case.weights = caseweights(d$Diagnosis, d$Diagnosis0, scw),
        sample.fraction = rusranger:::.samplefraction(d$Diagnosis),
        replace = replace
    )
}

train$Diagnosis0 <- train$Diagnosis
train$Diagnosis <- as.integer(
    train$Diagnosis & sbcdata:::.is_time_range(train, range = c(12, 0) * 3600)
)

rcw <- cwranger(train, scw = scw)
validation <- do.call(rbind, validation)
validation$RangerProbability <-
    predict(rcw, as.data.frame(validation[, xvar, with = FALSE]))$prediction[, 2L]

saveRDS(validation, file = paste0("ValidationScw", scw, ".RDS"))
