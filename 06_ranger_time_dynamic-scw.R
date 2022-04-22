library("sbcdata")
library("rusranger")
library("beeswarm")

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
#nsub <- 1e5
#ukl <- ukl[sample(nrow(ukl), nsub),]
#umg <- umg[sample(nrow(umg), nsub),]
#mimic <- mimic[sample(nrow(mimic), nsub),]

train <- subset(ukl, Set == "Training")
validation <- list(
    UKL = subset(ukl, Set == "Validation"), UMG = umg, MIMIC = mimic
)

## classic approach (6-72) h
train$Diagnosis0 <- as.integer(
    train$Diagnosis & sbcdata:::.is_time_range(train, range = c(72, 6) * 3600)
)
train$Excluded0 <- !train$Diagnosis0 & train$Diagnosis

xvar <- c("Age", "Sex", "PLT", "RBC", "WBC", "HGB", "MCV")
#timebreaks <-
#    c(0, 6, 12, 24, 48, 72, 7 * 24, 14 * 24, 30 * 24, 365 * 24, Inf) * 3600
#names(timebreaks) <- c(
#    paste("<", c(6, 12, 24), "h"),
#    paste("<", c(2, 3), "Tage"),
#    "< 1 Woche",
#    paste("<", c(2, 4), "Wochen"),
#    "< 1 Jahr", "> 1 Jahr"
#)
timebreaks <-
    c(0, 6, 12, 18, 24, 48, 72, 96, 120, 144, 168, 336, 672, Inf) * 3600
names(timebreaks) <- c(
    paste("<", c(6, 12, 18, 24), "h"),
    paste("<", 2:6, "Tage"),
    paste("<", c(1:2, 4), "Wochen"),
    "> 4 Wochen"
)

mtry <- 3
num.trees <- 1000
replace <- TRUE # doing bootstrap, shorter runtime, b/c much lower samplesize/samplefraction
scw <- as.integer(c(0, 1, 2, 4, 10, 20, 50, 100))
nrep <- 3

auc <- function(p, y) {
    pred <- prediction(p, y)
    performance(pred, "auc")@y.values[[1L]]
}

cut0 <- function(x, breaks) {
    ct <- cut(x, breaks, include.lowest = TRUE)
    levels(ct) <- names(breaks)
    ct
}

.lm <- function(p, time, breaks) {
    g <- split(p, cut0(time, breaks))
    m <- vapply(g, median, NA_real_)
    n <- lengths(g)
    lm(y ~ x, data = data.frame(x = seq_along(m), y = m), weights = n / sum(n))
}

.rlm <- function(p, time, breaks) {
    g <- split(p, cut0(time, breaks))
    m <- vapply(g, median, NA_real_)
    n <- lengths(g)
    MASS::rlm(
        y ~ x, data = data.frame(x = seq_along(m), y = m),
        weights = n, wt.method = "case"
    )
}

.rlm2 <- function(p, time, breaks) {
    g <- split(p, cut0(time, breaks))
    m <- unlist(g)
    n <- rep(seq_along(g), lengths(g))
    MASS::rlm(y ~ x, data = data.frame(x = n, y = m))
}

.rlm3 <- function(p, time, breaks) {
    i <- order(time)
    MASS::rlm(y ~ x, data = data.frame(x = time[i], y = p[i]))
}

bplot <- function(p, time, y, breaks, l, main = "", ...) {
    old.par <- par(mar = c(10, 4, 4, 2) + 0.1, no.readonly = TRUE)
    on.exit(par(old.par))

    col <- c(Control = "#21908C", Sepsis = "#440154")
    bcol <- setNames(paste0(col, "77"), names(col))
    pcol <- setNames(paste0(col, "FF"), names(col))

    d <- data.frame(
        p = p, g = cut0(time, breaks),
        y = as.factor(ifelse(y, "S", "C"))
    )
    bp <- boxplot(
        p ~ y + g, data = d, las = 2, col = bcol, border = pcol, pch = 20,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = FALSE,
        ylim = c(0, 1)
    )
    title(main, adj = 0L)
    title(ylab = "Probability for Sepsis", adj = 1L)
    axis(
        1, lwd.ticks = 0L, col = "#808080",
        at = seq_len(2 * nlevels(d$g)),
        labels = paste0(bp$names, " (", bp$n, ")"),
        las = 2
    )
    axis(2, lwd.ticks = 0L, col = "#808080")
    #beeswarm(
    #    p ~ y + g, data = d, col = pcol, pch = 20, add = TRUE,
    #    cex = 0.75, method = "compactswarm"
    #)
    abline(l$Control, col = pcol["Control"], lwd = 2, lty = 2)
    abline(l$Sepsis, col = pcol["Sepsis"], lwd = 2, lty = 2)
}

alpha <- function(p, time, y, breaks, plot = FALSE, ...) {
    nna <- !is.na(time)
    y <- as.logical(y)
    time <- time[nna]
    p <- p[nna]
    y <- y[nna]
    l <- list(
        Control = .rlm(p[!y], time[!y], breaks = breaks),
        Sepsis = .rlm(p[y], time[y], breaks = breaks)
    )
    a <- c(
        Control = unname(coef(l$Control)[2L]),
        Sepsis = unname(coef(l$Sepsis)[2L])
    )
    ## should be small (negative) for Sepsis and positive for Control
    l$alpha <- unname(
        if (a["Sepsis"] > 0) # penalize positive alpha
            abs(a["Control"]) + abs(a["Sepsis"])
        else
            a["Sepsis"] - a["Control"]
    )
    if (plot)
        bplot(p = p, time = time, y = y, breaks = breaks, l = l, ...)

    print(l)
    l
}

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

dir.create(
        file.path("boxplots", "scw"),
        recursive = TRUE, showWarnings = FALSE
)

# Baseline
a.baseline <- rowMeans(simplify2array(lapply(seq_len(nrep), function(r) {
    d <- copy(train)
    d <- d[!d$Excluded0,]
    rus <- rusranger(
        x = as.data.frame(d[, xvar, with = FALSE]), y = d$Diagnosis0,
        mtry = mtry, num.trees = num.trees, replace = replace
    )
    sapply(validation, function(v) {
        png(
            file.path(
                "boxplots", "scw",
                paste0(v$Center[1], "_baseline-", r, ".png")
            ),
            width = 1024, height = 768
        )
        on.exit(dev.off())
        p <- predict(
            rus, as.data.frame(v[, xvar, with = FALSE])
        )$predictions[, 2L]
        c(
            alpha = alpha(
                p, y = v$Diagnosis, time = v$SecToIcu, breaks = timebreaks,
                plot = TRUE, main = paste(v$Center[1], "Baseline", r
            ))$alpha,
            auc = auc(p, y = v$Diagnosis)
        )
    })
})), dims = 2L)

train$Diagnosis0 <- train$Diagnosis
train$Diagnosis <- as.integer(
    train$Diagnosis & sbcdata:::.is_time_range(train, range = c(12, 0) * 3600)
)

a.scw <- simplify2array(lapply(scw, function(i) {
    rowMeans(simplify2array(lapply(seq_len(nrep), function(r) {
        rcw <- cwranger(train, scw = i)
        sapply(validation, function(v) {
            png(
                file.path(
                    "boxplots", "scw",
                    paste0(v$Center[1], "_scw-", i, "-", r, ".png")
                ),
                width = 1024, height = 768
            )
            on.exit(dev.off())
            p <- predict(
                rcw, as.data.frame(v[, xvar, with = FALSE])
            )$prediction[, 2L]
            c(
                alpha = alpha(
                    p, y = v$Diagnosis, time = v$SecToIcu, breaks = timebreaks,
                    plot = TRUE, main = paste(v$Center[1], i, r
                ))$alpha,
                auc = auc(p, y = v$Diagnosis)
            )
        })
    })), dims = 2L)
}))
dimnames(a.scw)[[3L]] <- scw

knitr::kable(t(a.baseline))

for (i in names(validation))
    print(knitr::kable(t(a.scw[, i, ])))

plot_a_trend <- function(baseline, scw, main) {
    col <- palette.colors(2)
    grd <- c(5, 10, 20, 50, 100)

    old.par <- par(mfrow = c(2, 1))
    on.exit(par(old.par))

    ## alpha
    ylim <- range(c(baseline["alpha"], scw["alpha",]))
    xscw <- as.numeric(colnames(scw))

    par(mar = c(0, 4, 4, 2) + 0.1)
    plot(
        NA,
        xlim = c(0L, max(xscw)), ylim = ylim,
        axes = FALSE, ann = FALSE
    )
    title(main = paste(main, "alpha trend"), adj = 0L)
    title(ylab = "alpha", adj = 1L)
    title(xlab = "SCW", adj = 1L)
    axis(2, lwd.ticks = 0L, col = "#808080")
    abline(v = grd, lwd = 0.75, col = "lightgray", lty = "dotted")

    abline(h = baseline["alpha"], col = col[1])

    lines(xscw, scw["alpha",], col = col[2], type = "b", pch = 20)

    legend(
        "bottomright",
        legend = c("RUS 6-72", "SCW"),
        col = col,
        lwd = 1,
        bty = "n"
    )

    ## auc
    par(mar = c(5, 4, 2, 2) + 0.1)
    plot(
        NA,
        xlim = c(0L, max(xscw)), ylim = c(0.5, 1),
        axes = FALSE, ann = FALSE
    )
    title(main = paste(main, "AUC trend"), adj = 0L)
    title(ylab = "AUC", adj = 1L)
    title(xlab = "SCW", adj = 1L)
    axis(1, lwd.ticks = 0L, col = "#808080")
    axis(2, lwd.ticks = 0L, col = "#808080")
    abline(v = grd, lwd = 0.75, col = "lightgray", lty = "dotted")
    abline(h = 0.5, lty = "dashed", col = "#808080")

    abline(h = baseline["auc"], col = col[1])

    lines(xscw, scw["auc",], col = col[2], type = "b", pch = 20)

    legend(
        "bottomright",
        legend = c("RUS 6-72", "SCW"),
        col = col,
        lwd = 1,
        bty = "n"
    )
}

for (nm in colnames(a.baseline)) {
    png(paste0("alpha-scw-trend_", nm, ".png"), width = 1024, height = 1024)
    plot_a_trend(a.baseline[,nm], a.scw[,nm,], main = nm)
    dev.off()
}
