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
nsub <- 1e5
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
nrep <- 2

auc <- function(fit, x, y) {
    p <- prediction(predict(fit, as.data.frame(x))$predictions[, 2L], y)
    performance(p, "auc")@y.values[[1L]]
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
    col <- c(Control = "#21908C", Sepsis = "#440154")
    bcol <- setNames(paste0(col, "77"), names(col))
    pcol <- setNames(paste0(col, "FF"), names(col))

    d <- data.frame(
        p = p, g = cut0(time, breaks),
        y = as.factor(ifelse(y, "S", "C"))
    )
    bp <- boxplot(
        p ~ y + g, data = d, las = 2, col = bcol, border = pcol, pch = NA,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = FALSE,
        ylim = c(0, 1)
    )
    title(main, adj = 0L)
    title(ylab = "Probability for Sepsis", adj = 1L)
    axis(
        1, lwd.ticks = 0L, col = "#808080",
        at = seq_len(2 * nlevels(d$g)), labels = bp$names, las = 2
    )
    axis(2, lwd.ticks = 0L, col = "#808080")
    beeswarm(
        p ~ y + g, data = d, col = pcol, pch = 20, add = TRUE,
        cex = 0.75, method = "compactswarm"
    )
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

rngr <- function(timerange) {
    d <- copy(train)
    d$Diagnosis <- as.integer(
        d$Diagnosis & sbcdata:::.is_time_range(d, range = timerange * 3600)
    )
    rusranger(
        x = as.data.frame(d[, xvar, with = FALSE]), y = d$Diagnosis,
        mtry = mtry, num.trees = num.trees, replace = replace
    )
}

dir.create(
        file.path("boxplots", "timerange"),
        recursive = TRUE, showWarnings = FALSE
)

# Baseline
a.baseline <- rowMedians(do.call(cbind, lapply(seq_len(nrep), function(r) {
    d <- copy(train)
    d <- d[!d$Excluded0,]
    rus <- rusranger(
        x = as.data.frame(d[, xvar, with = FALSE]), y = d$Diagnosis0,
        mtry = mtry, num.trees = num.trees, replace = replace
    )
    sapply(validation, function(v) {
        png(
            file.path(
                "boxplots", "timerange",
                paste0(v$Center[1], "_baseline-", r, ".png")
            ),
            width = 1024, height = 768
        )
        on.exit(dev.off())
        alpha(
            predict(
                rus,
                as.data.frame(v[, xvar, with = FALSE])
            )$predictions[, 2L],
            y = v$Diagnosis, time = v$SecToIcu, breaks = timebreaks,
            plot = TRUE, main = paste(v$Center[1], "Baseline", r)
        )$alpha
    })
})))

# time ranges
timeranges <- rbind.data.frame(
    c(0, 6),
    c(0, 12),
    c(0, 18),
    c(0, 24),
    c(0, 48),
    c(0, 72),
    c(6, 12),
    c(6, 24),
    c(6, 48),
    c(6, 72),
    c(6, 96),
    c(12, 8760) # 12 h - 365
)
names(timeranges) <- c("from", "to")

a.time <- do.call(cbind, lapply(seq_len(nrow(timeranges)), function(i) {
    rowMedians(do.call(cbind, lapply(seq_len(nrep), function(r) {
        rng <- rngr(unlist(timeranges[i,]))
        lapply(validation, function(v) {
            tr <- paste(unlist(timeranges[i,]), collapse = "-")
            png(
                file.path(
                    "boxplots", "timerange",
                    paste0(v$Center[1], "_", tr, "-", r, ".png")
                ),
                width = 1024, height = 768
            )
            on.exit(dev.off())
            alpha(
                predict(rng,
                    as.data.frame(v[, xvar, with = FALSE])
                )$prediction[, 2L],
                y = v$Diagnosis, time = v$SecToIcu, breaks = timebreaks,
                plot = TRUE, main = paste(v$Center[1], tr, r)
            )$alpha
        })
    })))
}))
colnames(a.time) <- apply(timeranges, 1, paste, collapse = ":")

knitr::kable(t(a.baseline))
#
#
# |       UKL|       UMG|     MIMIC|
# |---------:|---------:|---------:|
# | 0.0295118| 0.0171677| 0.0130184|

knitr::kable(t(a.time))
#
#
# |        |        UKL|       UMG|      MIMIC|
# |:-------|----------:|---------:|----------:|
# |0:6     | -0.0127910| 0.0024242| -0.0050949|
# |0:12    | -0.0122632| 0.0019504| -0.0080606|
# |0:18    | -0.0114314| 0.0019058| -0.0080640|
# |0:24    |  0.0004718| 0.0048308| -0.0083883|
# |0:48    |  0.0101419| 0.0038331| -0.0042198|
# |0:72    |  0.0108929| 0.0094506| -0.0065397|
# |6:12    |  0.0134906| 0.0055807|  0.0001330|
# |6:24    |  0.0249579| 0.0144781|  0.0012085|
# |6:48    |  0.0275609| 0.0219951|  0.0143872|
# |6:72    |  0.0295745| 0.0162717|  0.0109591|
# |6:96    |  0.0316370| 0.0166541|  0.0135355|
# |12:8760 |  0.0473754| 0.0357091|  0.0266228|
