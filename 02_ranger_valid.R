library("sbcdata")
library("ranger")
library("ROCR")

source("00_rus_ranger_functions.R")

r <- readRDS("res-3-3-3.RDS")

hp <- lapply(r, "[[", "selectedparams")
hp <- do.call(rbind, hp)
hp
hp <- apply(hp, 2L, median)
hp

ukl <- exclude_entries(subset(sbcdata, Center == "Leipzig"))
ukl <- ukl[!ukl$Excluded,]
umg <- exclude_entries(subset(sbcdata, Center == "Greifswald"))
umg <- umg[!umg$Excluded,]
mimic <- exclude_entries(import_mimic("../../inst/intdata/mimic-iv-1.0/"))
mimic <- mimic[!mimic$Excluded,]

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

set.seed(20220322)
rngr <- do.call(
    rus_ranger,
    c(
        list(
            x = ukl[, xvar, with = FALSE],
            y = ifelse(ukl$Diagnosis == "Sepsis", 1, 2)
        ),
        hp
    )
)

pred <- c(
    setNames(
        lapply(
            r, function(rr) {
                prediction(rr$prediction, rr$truth)
            }
        ),
        paste0("UKL", seq_along(r))
    ),
    lapply(
        list(UMG = umg, MIMIC = mimic),
        function(dd) {
            pred <- predict(rngr, dd[, xvar, with = FALSE])$predictions[, 2L]
            prediction(pred, ifelse(dd$Diagnosis == "Sepsis", 1, 2))
        }
    )
)

roc <- lapply(pred, performance, "tpr", "fpr")
auc <- sapply(pred, function(pp)performance(pp, "auc")@y.values[[1L]])

## ROC
pdf("roc.pdf", width = 7, height = 7)
col <- palette.colors(length(roc))
lwd <- c(rep(1, length(r)), 2, 2)
lty <- c(rep(2, length(r)), 1, 1)
plot(NA, xlim = c(0L, 1L), ylim = c(0L, 1L), axes = FALSE, ann = FALSE)
title(main = "ROC", adj = 0L)
title(ylab = "Sensitivity", adj = 1L)
title(xlab = "1 - Specificity", adj = 1L)
axis(1, lwd.ticks = 0L, col = "#808080")
axis(2, lwd.ticks = 0L, col = "#808080")
abline(0L, 1L, col = "#808080", lty = 2L, lwd = 1L)
for (i in seq(along = roc)) {
    plot(roc[[i]], col = col[i], lwd = lwd[i], lty = lty[i], add = TRUE)
}

lgd <- do.call(
    c,
    lapply(
        seq_along(auc),
        function(i)parse(text =
            paste0("AUC[", names(auc)[i], "] == ", format(auc[i], digits = 3))
        )
    )
)

legend("bottomright", legend = lgd, lwd = lwd, lty = lty, col = col, bty = "n")
dev.off()

## PRC
prc <- lapply(pred, performance, "prec", "rec")
pdf("prc.pdf", width = 7, height = 7)
col <- palette.colors(length(roc))
lwd <- c(rep(1, length(r)), 2, 2)
lty <- c(rep(2, length(r)), 1, 1)
plot(NA, xlim = c(0L, 1L), ylim = c(0.5, 1L), axes = FALSE, ann = FALSE)
title(main = "PRC", adj = 0L)
title(ylab = "Precision", adj = 1L)
title(xlab = "Recall", adj = 1L)
axis(1, lwd.ticks = 0L, col = "#808080")
axis(2, lwd.ticks = 0L, col = "#808080")
for (i in seq(along = prc)) {
    plot(prc[[i]], col = col[i], lwd = lwd[i], lty = lty[i], add = TRUE)
}

legend("bottomright", legend = names(prc), lwd = lwd, lty = lty, col = col, bty = "n")
dev.off()

## CAL
val <- lapply(
    list(UMG = umg, MIMIC = mimic),
    function(dd)
        data.frame(
            predicted = predict(rngr, dd[, xvar, with = FALSE])$predictions[, 2L],
            observed = ifelse(dd$Diagnosis == "Sepsis", 1, 2)
        )
)
val[["UMG"]]$center <- "UMG"
val[["MIMIC"]]$center <- "MIMIC"

val <- do.call(rbind, val)

p <- lapply(
    r,
    function(rr)
        data.frame(
            predicted = rr$prediction,
            observed = rr$truth, center = "UKL"
        )
)

p <- rbind(val, do.call(rbind, p))
rownames(p) <- NULL

m <- 1000
ctpts <- lapply(split(p$predicted, p$center), cutpoints, n = m)
cts <- mapply(
    cut, split(p$predicted, p$center), ctpts, MoreArgs = list(include.lowest = TRUE)
)
ps <- mapply(groupmean, x = split(p$predicted, p$center), f = cts)
os <- mapply(groupmean, x = split(p$observed == 2, p$center), f = cts)

pdf("cal.pdf", width = 7, height = 7)
col <- palette.colors(length(os))
plot(NA, xlim = c(0L, 1L), ylim = c(0L, 1L), axes = FALSE, ann = FALSE)
title(main = "Calibration", adj = 0L)
title(ylab = "Observed", adj = 1L)
title(xlab = "Predicted", adj = 1L)
abline(0L, 1L, col = "#808080", lty = 2L, lwd = 1L)
axis(1, lwd.ticks = 0L, col = "#808080")
axis(2, lwd.ticks = 0L, col = "#808080")

for (i in seq_along(os))
    lines(ps[[i]], os[[i]], col = col[i], type = "b", pch = 19)

legend(
    "bottomright",
    legend = names(os),
    col = col, pch = 19, bty = "n"
)
dev.off()
