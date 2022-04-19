library("sbcdata")
library("rusranger")

ukl <- exclude_entries(subset(sbcdata, Center == "Leipzig"))
ukl <- ukl[!ukl$Excluded,]

ukl$Sex <- as.numeric(ukl$Sex == "M")

ukl$Diagnosis <- as.numeric(ukl$Diagnosis == "Sepsis")

set.seed(20220417)
#nsub <- 1e5
#ukl <- ukl[sample(nrow(ukl), nsub),]
train <- subset(ukl, Set == "Training")
validation <- subset(ukl, Set == "Validation")

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")

mtry <- 3
num.trees <- 1000

rngr <- rusranger(
    train[, xvar, with = FALSE], train$Diagnosis,
    mtry = mtry, num.trees = num.trees
)

pred <- predict(
    rngr,
    as.data.frame(validation[, xvar, with = FALSE]),
    predict.all = TRUE
)$predictions[, 2L,]

m <- rowMeans(pred)
qu <- apply(pred, 1, quantile, probs = c(0.025, 0.975))
i <- order(m)

png("pred-distribution-UKL.png", width = 1024, height = 512)
col <- palette.colors(3, alpha = 0.25)[2:3]
plot(NA, xlim = c(0L, nrow(pred)), ylim = c(0L, 1L), axes = FALSE, ann = FALSE)
title(main = "Distribution Predictions", adj = 0L)
title(xlab = "ordered predictions", adj = 1L)
axis(1, lwd.ticks = 0L, col = "#808080")
axis(2, lwd.ticks = 0L, col = "#808080")

points(m[i], col = col[1], pch = 19)
points(qu[1, i], col = col[2], pch = 20)
points(qu[2, i], col = col[2], pch = 20)

legend(
    "topleft",
    legend = c("mean prediction", "2.5/97.5 % quantile"),
    bty = "n",
    col = col,
    pch = c(19, 20)
)
dev.off()
