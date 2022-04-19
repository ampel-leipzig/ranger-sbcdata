library("sbcdata")
library("rusranger")

ukl <- exclude_entries(subset(sbcdata, Center == "Leipzig"))
ukl <- ukl[!ukl$Excluded,]
umg <- exclude_entries(subset(sbcdata, Center == "Greifswald"))
umg <- umg[!umg$Excluded,]
mimic <- exclude_entries(import_mimic("../../inst/intdata/mimic-iv-1.0/"))
mimic <- mimic[!mimic$Excluded,]

ukl$Sex <- as.numeric(ukl$Sex == "M")
umg$Sex <- as.numeric(umg$Sex == "M")
mimic$Sex <- as.numeric(mimic$Sex == "M")

ukl$Diagnosis <- as.numeric(ukl$Diagnosis == "Sepsis")
umg$Diagnosis <- as.numeric(umg$Diagnosis == "Sepsis")
mimic$Diagnosis <- as.numeric(mimic$Diagnosis == "Sepsis")

set.seed(20220417)
#nsub <- 1e5
#ukl <- ukl[sample(nrow(ukl), nsub),]
#umg <- umg[sample(nrow(umg), nsub),]
#mimic <- mimic[sample(nrow(mimic), nsub),]

train <- subset(ukl, Set == "Training")
validation <- list(
    UKL = subset(ukl, Set == "Validation"), UMG = umg, MIMIC = mimic
)

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

mtry <- 3
num.trees <- 1000
ndups <- as.integer(c(1, 2, 4, 10, 20, 50, 100))
nrep <- 5

auc <- function(fit, x, y) {
    p <- prediction(predict(fit, as.data.frame(x))$predictions[, 2L], y)
    performance(p, "auc")@y.values[[1L]]
}

# baseline
auc.baseline <- do.call(cbind, lapply(seq_len(nrep), function(r) {
    rngr <- ranger(
        x = as.data.frame(train[, xvar, with = FALSE]), y = train$Diagnosis,
        mtry = mtry, num.trees = num.trees, probability = TRUE, min.node.size = 10
    )
    sapply(validation, function(v)
        auc(rngr, v[, xvar, with = FALSE], v$Diagnosis)
    )
}))
mode(auc.baseline) <- "numeric"
auc.baseline <- apply(auc.baseline, 1, median)

# RUS
auc.rus <- do.call(cbind, lapply(seq_len(nrep), function(r) {
    rus <- rusranger(
        x = as.data.frame(train[, xvar, with = FALSE]), y = train$Diagnosis,
        mtry = mtry, num.trees = num.trees
    )
    sapply(validation, function(v)
        auc(rus, v[, xvar, with = FALSE], v$Diagnosis)
    )
}))
mode(auc.rus) <- "numeric"
auc.rus <- apply(auc.rus, 1, median)

# ROS
auc.ros <- do.call(cbind, lapply(ndups, function(i) {
    a <- do.call(cbind, lapply(seq_len(nrep), function(r) {
        ros <- rosranger(
            x = as.data.frame(train[, xvar, with = FALSE]), y = train$Diagnosis,
            mtry = mtry, num.trees = num.trees, ndups = i
        )
        lapply(validation, function(v)
            auc(ros, v[, xvar, with = FALSE], v$Diagnosis)
        )
    }))
    mode(a) <- "numeric"
    apply(a, 1, median)
}))
colnames(auc.ros) <- ndups

# SMOTE
auc.smote <- do.call(cbind, lapply(ndups, function(i) {
    a <- do.call(cbind, lapply(seq_len(nrep), function(r) {
        sm <- smote(
            as.matrix(train[, xvar, with = FALSE]), train$Diagnosis, ndups = i
        )
        rsm <- ranger(
            x = sm$x, y = sm$y,
            probability = TRUE, min.node.size = 10,
            mtry = mtry, num.trees = num.trees
        )
        lapply(validation, function(v)
            auc(rsm, v[, xvar, with = FALSE], v$Diagnosis)
        )
    }))
    mode(a) <- "numeric"
    apply(a, 1, median)
}))
colnames(auc.smote) <- ndups

plot_auc_trend <- function(baseline, rus, ros, smote, main) {
    col <- palette.colors(4)

    xros <- as.numeric(names(ros))
    xsmote <- as.numeric(names(smote))

    plot(
        NA,
        xlim = c(1L, max(xros)), ylim = c(0L, 1L),
        axes = FALSE, ann = FALSE
    )
    title(main = main, adj = 0L)
    title(ylab = "AUC", adj = 1L)
    title(xlab = "x times over-sampled minority class", adj = 1L)
    axis(1, lwd.ticks = 0L, col = "#808080")
    axis(2, lwd.ticks = 0L, col = "#808080")
    abline(h = 0.5, lty = "dashed", col = "#808080")

    abline(h = c(baseline, rus), col = col[1:2])

    lines(xros, ros, col = col[3])
    lines(xsmote, smote, col = col[4])

    legend(
        "bottomright",
        legend = c("default ranger", "RUS", "ROS", "SMOTE"),
        col = col,
        lwd = 1,
        bty = "n"
    )
}

for (nm in names(auc.baseline)) {
    png(paste0("auc-trend_", nm, ".png"), width = 1024, height = 1024)
    plot_auc_trend(
        auc.baseline[nm], auc.rus[nm], auc.ros[nm, ], auc.smote[nm,],
        main = paste(nm, "AUC trend")
    )
    dev.off()
}
