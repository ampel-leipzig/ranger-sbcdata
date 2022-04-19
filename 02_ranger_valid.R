library("sbcdata")
library("rusranger")
library("ameld")

r <- readRDS("res-5-5-10.RDS")

hp <- lapply(r, "[[", "selectedparams")
hp <- do.call(rbind, hp)
hp
hp <- apply(hp, 2L, median)
hp

ukl <- exclude_entries(subset(sbcdata, Center == "Leipzig"))
ukl <- ukl[!ukl$Excluded,]
uklt <- subset(ukl, Set == "Training")
uklv <- subset(ukl, Set == "Validation")
umg <- exclude_entries(subset(sbcdata, Center == "Greifswald"))
umg <- umg[!umg$Excluded,]
mimic <- exclude_entries(import_mimic("../../inst/intdata/mimic-iv-1.0/"))
mimic <- mimic[!mimic$Excluded,]

## for development subsample for faster runtime
#set.seed(20220415)
#nsub <- 1e5
#ukl <- ukl[sample(nrow(ukl), nsub),]
#uklt <- uklt[sample(nrow(uklt), nsub),]
#uklv <- uklv[sample(nrow(uklv), min(c(nsub, nrow(uklv)))),]
#umg <- umg[sample(nrow(umg), nsub),]
#mimic <- mimic[sample(nrow(mimic), nsub),]

set.seed(20220325)
uklrandom <- uklt
uklrandom$Diagnosis <- sample(uklrandom$Diagnosis)

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

set.seed(20220324)
rngr <- do.call(
    rusranger,
    c(
        list(
            x = as.data.frame(ukl[, xvar, with = FALSE]),
            y = as.numeric(ukl$Diagnosis == "Sepsis")
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
        list(
            UKLvalidation = uklv, UKLrandom = uklrandom,
            UMG = umg, MIMIC = mimic
        ),
        function(dd) {
            pred <- predict(rngr, as.data.frame(dd[, xvar, with = FALSE]))$predictions[, 2L]
            prediction(pred, as.numeric(dd$Diagnosis == "Sepsis"))
        }
    )
)

roc <- lapply(pred, performance, "tpr", "fpr")
auc <- sapply(pred, function(pp)performance(pp, "auc")@y.values[[1L]])

## ROC
png("roc.png", width = 1024, height = 1024)
col <- palette.colors(length(roc))
lwd <- c(rep(1, length(r)), rep(2, length(pred) - length(r)))
lty <- c(rep(2, length(r)), rep(1, length(pred) - length(r)))
plot(NA, xlim = c(0L, 1L), ylim = c(0L, 1L), axes = FALSE, ann = FALSE)
title(main = "ROC", adj = 0L)
title(ylab = "Sensitivity", adj = 1L)
title(xlab = "1 - Specificity", adj = 1L)
axis(1, lwd.ticks = 0L, col = "#808080")
axis(2, lwd.ticks = 0L, col = "#808080")
abline(0L, 1L, col = "#808080", lty = 2L, lwd = 1L)
for (i in seq(along = roc))
    plot(roc[[i]], col = col[i], lwd = lwd[i], lty = lty[i], add = TRUE)

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
auc <- sapply(pred, function(pp)performance(pp, "aucpr")@y.values[[1L]])
png("prc.png", width = 1024, height = 1024)
plot(NA, xlim = c(0L, 1L), ylim = c(0L, 1L), axes = FALSE, ann = FALSE)
title(main = "PRC", adj = 0L)
title(ylab = "Precision", adj = 1L)
title(xlab = "Recall", adj = 1L)
axis(1, lwd.ticks = 0L, col = "#808080")
axis(2, lwd.ticks = 0L, col = "#808080")
for (i in seq(along = prc))
    plot(prc[[i]], col = col[i], lwd = lwd[i], lty = lty[i], add = TRUE)

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

## CAL
val <- lapply(
    list(UKLvalidation = uklv, UKLrandom = uklrandom, UMG = umg, MIMIC = mimic),
    function(dd)
        data.frame(
            predicted = predict(rngr, as.data.frame(dd[, xvar, with = FALSE]))$predictions[, 2L],
            observed = as.numeric(dd$Diagnosis == "Sepsis")
        )
)
val[["UMG"]]$center <- "UMG"
val[["MIMIC"]]$center <- "MIMIC"
val[["UKLrandom"]]$center <- "UKLrandom"
val[["UKLvalidation"]]$center <- "UKLvalidation"

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
    cut, split(p$predicted, p$center), ctpts,
    MoreArgs = list(include.lowest = TRUE), SIMPLIFY = FALSE
)
ps <- mapply(
    groupmean, x = split(p$predicted, p$center), f = cts, SIMPLIFY = FALSE
)
os <- mapply(
    groupmean, x = split(p$observed, p$center), f = cts, SIMPLIFY = FALSE
)

png("cal.png", width = 1024, height = 1024)
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
    "topleft",
    legend = names(os),
    col = col, pch = 19, bty = "n"
)
dev.off()

## hist
for (cnt in unique(p$center)) {
    png(paste0("hist-", cnt, ".png"), width = 1024, height = 1024)
    par(mfrow = c(2, 1))
    sb <- p[p$center == cnt,]
    hist(sb$predicted[sb$observed == 1], col = 1, xlim = c(0, 1))
    hist(sb$predicted[sb$observed == 0], col = 2, xlim = c(0, 1))
    dev.off()
}
