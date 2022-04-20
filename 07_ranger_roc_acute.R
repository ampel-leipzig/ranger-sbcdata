library("sbcdata")
library("rusranger")
library("ameld")

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

ur <- copy(subset(ukl, Set == "Validation"));
ur$Diagnosis <- sample(ur$Diagnosis)

set.seed(20220420)
train <- copy(subset(ukl, Set == "Training"))
train$Diagnosis <- as.integer(
    train$Diagnosis & sbcdata:::.is_time_range(train, range = c(12, 0) * 3600)
)

validation <- list(
    UKL = copy(subset(ukl, Set == "Validation")),
    UKLrandom = ur, UMG = umg, MIMIC = mimic
)

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

rngr <- rusranger(
    x = as.data.frame(train[, xvar, with = FALSE]),
    y = train$Diagnosis,
    mtry = 3, num.trees = 1000, replace = TRUE
)

pred <- lapply(validation, function(v) {
    pred <- predict(rngr, as.data.frame(v[, xvar, with = FALSE]))$predictions[, 2L]
    prediction(pred, v$Diagnosis)
})

roc <- lapply(pred, performance, "tpr", "fpr")
auc <- sapply(pred, function(pp)performance(pp, "auc")@y.values[[1L]])

## ROC
png("acute-roc.png", width = 1024, height = 1024)
col <- palette.colors(length(roc))
lwd <- rep_len(1, length(roc))
lty <- rep_len(1, length(roc))
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
png("acute-prc.png", width = 1024, height = 1024)
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
val <- lapply(validation, function(v)
    data.frame(
        predicted = predict(rngr, as.data.frame(v[, xvar, with = FALSE]))$predictions[, 2L],
        observed = v$Diagnosis
    )
)
val[["UMG"]]$center <- "UMG"
val[["MIMIC"]]$center <- "MIMIC"
val[["UKLrandom"]]$center <- "UKLrandom"
val[["UKL"]]$center <- "UKL"

val <- do.call(rbind, val)

m <- 1000
ctpts <- lapply(split(val$predicted, val$center), cutpoints, n = m)
cts <- mapply(
    cut, split(val$predicted, val$center), ctpts,
    MoreArgs = list(include.lowest = TRUE), SIMPLIFY = FALSE
)
ps <- mapply(
    groupmean, x = split(val$predicted, val$center), f = cts, SIMPLIFY = FALSE
)
os <- mapply(
    groupmean, x = split(val$observed, val$center), f = cts, SIMPLIFY = FALSE
)

png("acute-cal.png", width = 1024, height = 1024)
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
for (cnt in unique(val$center)) {
    png(paste0("acute-hist-", cnt, ".png"), width = 1024, height = 1024)
    par(mfrow = c(2, 1))
    sb <- val[val$center == cnt,]
    hist(sb$predicted[sb$observed == 1], col = 1, xlim = c(0, 1))
    hist(sb$predicted[sb$observed == 0], col = 2, xlim = c(0, 1))
    dev.off()
}
