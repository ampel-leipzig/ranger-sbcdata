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

ukl$Sex <- as.integer(ukl$Sex == "M")
umg$Sex <- as.integer(umg$Sex == "M")

ukl$Diagnosis <- as.integer(ukl$Diagnosis == "Sepsis")
umg$Diagnosis <- as.integer(umg$Diagnosis == "Sepsis")

set.seed(20220419)
#nsub <- 1e4
#ukl <- ukl[sample(nrow(ukl), nsub),]
#umg <- umg[sample(nrow(umg), nsub),]

train <- subset(ukl, Set == "Training")
validation <- list(
    UKL = subset(ukl, Set == "Validation"), UMG = umg
)

## classic approach (6-72) h
train$Diagnosis0 <- as.integer(
    train$Diagnosis & sbcdata:::.is_time_range(train, range = c(72, 6) * 3600)
)
train$Excluded0 <- !train$Diagnosis0 & train$Diagnosis

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

mtry <- 3
num.trees <- 1000
replace <- TRUE # doing bootstrap, shorter runtime, b/c much lower samplesize/samplefraction
scw <- as.integer(c(1, 2, 4, 10, 20, 50, 100))
nrep <- 3

auc <- function(fit, x, y) {
    p <- prediction(predict(fit, as.data.frame(x))$predictions[, 2L], y)
    performance(p, "auc")@y.values[[1L]]
}

rowMedians <- function(x) {mode(x) <- "numeric"; apply(x, 1, median)}

#' Upsampling
#'
#' Upsampling former Sepsis now control cases using Daniel's approach.
#'
#' @param y current y
#' @param y0 original y (without any timerange)
#' @param scw sepsis control weighting factor
#' @return `numeric`, case weights
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

# Baseline
auc.baseline <- rowMedians(do.call(cbind, lapply(seq_len(nrep), function(r) {
    d <- train
    d <- d[!d$Excluded0,]
    rus <- rusranger(
        x = as.data.frame(d[, xvar, with = FALSE]), y = d$Diagnosis0,
        mtry = mtry, num.trees = num.trees, replace = replace
    )
    sapply(validation, function(v)
        auc(rus, v[, xvar, with = FALSE], v$Diagnosis)
    )
})))

# SCW
train$Diagnosis <- as.integer(
    train$Diagnosis & sbcdata:::.is_time_range(train, range = c(12, 0) * 3600)
)

auc.scw <- do.call(cbind, lapply(scw, function(i) {
    rowMedians(do.call(cbind, lapply(seq_len(nrep), function(r) {
        rcw <- cwranger(train, scw = i)
        lapply(validation, function(v)
            auc(rcw, v[, xvar, with = FALSE], v$Diagnosis)
        )
    })))
}))
colnames(auc.scw) <- scw

plot_auc_trend <- function(baseline, scw, main) {
    col <- palette.colors(2)

    xscw <- as.numeric(names(scw))

    plot(
        NA,
        xlim = c(1L, max(xscw)), ylim = c(0L, 1L),
        axes = FALSE, ann = FALSE
    )
    title(main = main, adj = 0L)
    title(ylab = "AUC", adj = 1L)
    title(xlab = "SCW", adj = 1L)
    axis(1, lwd.ticks = 0L, col = "#808080")
    axis(2, lwd.ticks = 0L, col = "#808080")
    abline(h = 0.5, lty = "dashed", col = "#808080")

    abline(h = baseline, col = col[1:2])

    lines(xscw, scw, col = col[2])

    legend(
        "bottomright",
        legend = c("RUS", "SCW"),
        col = col,
        lwd = 1,
        bty = "n"
    )
}

for (nm in names(auc.baseline)) {
    png(paste0("auc-scw-trend_", nm, ".png"), width = 1024, height = 1024)
    plot_auc_trend(
        auc.baseline[nm], auc.scw[nm,], main = paste(nm, "AUC trend")
    )
    dev.off()
}
