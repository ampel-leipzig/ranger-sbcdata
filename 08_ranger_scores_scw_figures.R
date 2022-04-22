nejm <- c(Control = "#6488c1", Sepsis = "#c8513c")
vdis <- c(Control = "#21908C", Sepsis = "#440154")
cdis <- setNames(rev(viridisLite::cividis(3)[-3]), c("Control", "Sepsis"))
col <- nejm


#' Add Transparency
#'
#' Add Alpha channel to hex encoded rgb colors.
#'
#' @param x `character`, color hex code.
#' @param alpha `numeric`, 0-1.
#' @return `character`, hex encoded color with alpha channel
rgba <- function(x, alpha) {
    if (alpha < 0 || alpha > 1)
        stop("'alpha' has to be between 0 and 1.")
    x <- col2rgb(x, alpha = TRUE) / 255
    x["alpha",] <- alpha
    rgb(x["red",], x["green",], x["blue",], x["alpha",])
}

#' Colorize Title
#'
#' Add colors to title.
#'
#' @param main `character`, main title.
#' @param col.main  `character`, color, has to have the same length as `main`.
#' @param ... further arguments passed to `title`.
#' @examples
#' plot(1:3)
#' cmain(c("FOO", "BAR"), col.main = c(1, 2))
cmain <- function(main, col.main = rep(1, length(main)), ...) {
    n <- length(main)

    if (n != length(col.main))
        stop("'main' and 'col.main' have to have the same length.")

    p <- function(...)paste0(..., collapse = "")

    for (i in seq_along(main)) {
        left <- seq_len(i - 1)
        right <- i + seq_len(n - i)
        title(
            bquote(bold(
                phantom(.(p(main[left])))*.(main[i])*phantom(.(p(main[right])))
            )),
            col.main = col.main[i],
            ...
        )
    }
}
cylab <- function(ylab, col.lab = rep(1, length(ylab)), ...) {
    n <- length(ylab)

    if (n != length(col.lab))
        stop("'ylab' and 'col.lab' have to have the same length.")

    p <- function(...)paste0(..., collapse = "")

    for (i in seq_along(ylab)) {
        left <- seq_len(i - 1)
        right <- i + seq_len(n - i)
        title(
            ylab = bquote(
                phantom(.(p(ylab[left])))*.(ylab[i])*phantom(.(p(ylab[right])))
            ),
            col.lab = col.lab[i],
            ...
        )
    }
}

#' Find ylim
#'
#' Find ylim depending on x$Score.
#'
#' @param x `data.frame`
#' @return `numeric(2)`, ylim
ylim0 <- function(x)c(floor(min(x$Score)), ceiling(max(x$Score)))

#' Filter NA
#'
#' Filter NA from SecToIcu column
#'
#' @param x `data.frame`
#' @return `data.frame` with `SecToIcu != NA`
omitNa <- function(x) {
    x[!is.na(x$SecToIcu),]
}

#' Default boxplot
bp_default <- function(x, ...) {
    boxplot(Score ~ Interval, data = x, ...)
}
bp_default2 <- function(x, ...) {
    boxplot(Score ~ Diagnosis + Interval, data = x, ...)
}

bp_color <- function(x, col = c(Control = 3, Sepsis = 2), main = "", ...) {
    old.par <- par(mar = c(12, 4, 4, 2) + 0.1, no.readonly = TRUE)
    on.exit(par(old.par))

    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))

    bp <- boxplot(
        Score ~ Diagnosis + Interval, data = x,
        col = bcol, border = pcol, pch = 20,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = FALSE,
        ylim = c(0, 1)
    )
    title(main, adj = 0L)
    title(ylab = "Probability for Sepsis", adj = 1L)
    axis(
        1, lwd.ticks = 0L, col = "#808080",
        at = seq_len(2 * nlevels(x$Interval)),
        labels = paste0(
            rep(c("Control", "Sepsis"), nlevels(x$Interval)),
            ", ",
            rep(levels(x$Interval), each = 2),
            " (n = ", bp$n, ")"),
        las = 2
    )
    axis(2, lwd.ticks = 0L, col = "#808080")
}

bp_color2 <- function(x, col = c(Control = 3, Sepsis = 2),
    main, col.main, cex.main = 1.5,
    ylim = c(0, 1), ...) {

    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- sort(c(
        seq(1, by = 3, length.out = nl), seq(2, by = 3, length.out = nl)
    ))

    bp <- boxplot(
        Score ~ Diagnosis + Interval, data = x, at = at,
        col = bcol, border = pcol, pch = 20,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = FALSE,
        ylim = ylim
    )
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    cylab(
        ylab = c("Probability for ", "Sepsis"),
        col.lab = c("black", pcol["Sepsis"]),
        adj = 1L
    )
    axis(
        1, lwd.ticks = 0L, col = "#808080",
        at = at[as.logical(seq_along(at) %% 2L)] + 0.5,
        labels = levels(x$Interval)
    )
    axis(2, lwd.ticks = 0L, col = "#808080", las = 1)
    bp
}

bp_color3 <- function(x, col = c(Control = 3, Sepsis = 2),
    main, col.main, cex.main = 1.5,
    ylim = c(0, 1), ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- sort(c(
        seq(1, by = 3, length.out = nl), seq(2, by = 3, length.out = nl)
    ))
    mid <- at[as.logical(seq_along(at) %% 2L)] + 0.5
    xlim <- c(0, max(at))

    layout(matrix(c(1, 2), nrow = 2), height = c(5, 1))

    ## boxplot
    par(mar = c(0, 10, 4, 2) + 0.1)
    bp <- boxplot(
        Score ~ Diagnosis + Interval, data = x, at = at,
        col = bcol, border = pcol, pch = 20,
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame = FALSE,
        xlim = xlim, ylim = ylim
    )
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    cylab(
        ylab = c("Probability for ", "Sepsis"),
        col.lab = c("black", pcol["Sepsis"]),
        adj = 1L
    )
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)

    ## table
    par(mar = c(4, 10, 0, 2) + 0.1)
    n <- matrix(
        bp$n, ncol = 2,
        dimnames = list(levels(v$Interval), c("Control", "Sepsis"))
    )
    voffset <- 0.25
    plot(NA, xlim = xlim, ylim = c(0, 3), axes = FALSE, ann = FALSE)
    axis(
        side = 2L, at = 2L + voffset,
        labels = bquote(bold(.("Number of Patients"))),
        padj = -1L, las = 1L, tick = FALSE
    )
    for (i in seq_len(ncol(n)))
        axis(
            side = 2L, at = 2L - i + voffset, labels = colnames(n)[i],
            padj = -1L, las = 1L, col.axis = pcol[i], tick = FALSE
        )
    text(x = mid, y = rep(1:0, each = nl) + voffset, pos = 3, labels = n)

    title(xlab = "Time before ICU Admission", adj = 1L)
    axis(
        1, lwd.ticks = 0L, col = "#808080",
        at = mid, labels = levels(x$Interval)
    )
}

bp_median <- function(x, col = c(Control = 3, Sepsis = 2),
                       main, col.main, cex.main = 1.5,
                       ylim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- seq_len(nl)
    xlim <- range(at) + c(-0.5, 0.5)
    width <- rep_len(0.4, length(at))

    .bar <- function(x, y, w, lwd = 3 * par("lwd"), ...) {
        segments(x - w, y, x + w, y, lty = 1, lwd = lwd, lend = 1, ...)
    }

    prob <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                do.call(rbind, lapply(
                    split(xx$Score, xx$Interval),
                    quantile, probs = c(0.25, 0.5, 0.75)
                ))
        )
    )

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Time before ICU Admission", adj = 1L)
    cylab(
        ylab = c("Probability for ", "Sepsis"),
        col.lab = c("black", pcol["Sepsis"]),
        adj = 1L
    )
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)
    lbls <- if (reverse) rev(levels(x$Interval)) else levels(x$Interval)
    axis(
        1L, lwd.ticks = 0L, col = "#808080", at = at, labels = lbls
    )
    for (i in seq_along(pcol)) {
        y <- if (reverse) rev(prob[, 2, i]) else prob[, 2, i]
        lines(
            at, y, col = pcol[i], type = "b", pch = 19, cex = 1.5, lwd = 2
        )
    }
}

bp_median2 <- function(x, col = c(Control = 3, Sepsis = 2),
                       main, col.main, cex.main = 1.5,
                       ylim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- seq_len(nl)
    xlim <- range(at) + c(-0.5, 0.5)
    width <- rep_len(0.4, length(at))

    .bar <- function(x, y, w, lwd = 3 * par("lwd"), ...) {
        segments(x - w, y, x + w, y, lty = 1, lwd = lwd, lend = 1, ...)
    }

    prob <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                do.call(rbind, lapply(
                    split(xx$Score, xx$Interval),
                    quantile, probs = c(0.25, 0.5, 0.75)
                ))
        )
    )

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Time before ICU Admission", adj = 1L)
    cylab(
        ylab = c("Probability for ", "Sepsis"),
        col.lab = c("black", pcol["Sepsis"]),
        adj = 1L
    )
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)
    lbls <- if (reverse) rev(levels(x$Interval)) else levels(x$Interval)
    axis(
        1L, lwd.ticks = 0L, col = "#808080", at = at, labels = lbls
    )

    y <- if (reverse) rev(prob[, 2, ]) else prob[, 2, ]
    pcol <- rep(pcol, each = nl)
    pcol <- if(reverse) rev(pcol) else pcol
    .bar(at, y, width, col = pcol)
    points(rep(at, 2), y, col = pcol, pch = 19, cex = 1.5)
}

bp_median3 <- function(x, col = c(Control = 3, Sepsis = 2),
                       main, col.main, cex.main = 1.5,
                       ylim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    bcol <- setNames(rgba(col, 0.25), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- seq_len(nl)
    xlim <- range(at) + c(-0.5, 0.5)

    .bar <- function(x, y0, y1, width = 0.01, lwd = 8 * par("lwd"), ...) {
        segments(x, y0, x, y1, lty = 1, lwd = lwd, lend = 1, ...)
    }

    prob <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                do.call(rbind, lapply(
                    split(xx$Score, xx$Interval),
                    quantile, probs = c(0.25, 0.5, 0.75)
                ))
        )
    )

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Time before ICU Admission", adj = 1L)
    cylab(
        ylab = c("Probability for ", "Sepsis"),
        col.lab = c("black", pcol["Sepsis"]),
        adj = 1L
    )
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)
    lbls <- if (reverse) rev(levels(x$Interval)) else levels(x$Interval)
    axis(
        1L, lwd.ticks = 0L, col = "#808080", at = at, labels = lbls
    )

    y0 <- if (reverse) rev(prob[, 1, ]) else prob[, 1, ]
    y1 <- if (reverse) rev(prob[, 3, ]) else prob[, 3, ]
    bcol <- rep(bcol, each = nl)
    bcol <- if(reverse) rev(bcol) else bcol
    pcol <- rep(pcol, each = nl)
    pcol <- if(reverse) rev(pcol) else pcol
    .bar(rep(at, 2), y0, y1, col = bcol)
    y <- if (reverse) rev(prob[, 2, ]) else prob[, 2, ]
    points(rep(at, 2), y, col = pcol, pch = 19, cex = 1.5)
}

bp_median_all <- function(x, col = c(Control = 3, Sepsis = 2),
                       main, col.main, cex.main = 1.5,
                       ylim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    pch <- 15:17
    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- seq_len(nl)
    xlim <- range(at) + c(-0.5, 0.5)
    width <- rep_len(0.3, length(at))

    .bar <- function(x, y, w, lwd = 3 * par("lwd"), ...) {
        segments(x - w, y, x + w, y, lty = 1, lwd = lwd, lend = 1, ...)
    }

    medall <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                vapply(split(xx$Score, xx$Interval), median, NA_real_)
        )
    )

    medcenter <- simplify2array(lapply(
        split(x, x$Center),
        function(x2)
            do.call(cbind, lapply(
                split(x2, x2$Diagnosis),
                function(x3)
                    vapply(split(x3$Score, x3$Interval), median, NA_real_)
            ))
    ))

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Time before ICU admission", adj = 0L)
    cylab(
        ylab = c("Probability for ", "Sepsis"),
        col.lab = c("black", pcol["Sepsis"]),
        adj = 1L
    )
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)
    lbls <- if (reverse) rev(levels(x$Interval)) else levels(x$Interval)
    axis(
        1L, lwd.ticks = 0L, col = "#808080", at = at, labels = lbls
    )

    #legend(
    #    "topleft", legend = sort(unique(x$Center)),
    #    pch = pch, col = "#808080", cex = 1.2, bty = "n",
    #    horiz = TRUE
    #)
    y <- if (reverse) rev(medall) else medall
    pcol <- rep(pcol, each = nl)
    pcol <- if(reverse) rev(pcol) else pcol
    .bar(at, y, width, col = pcol)
    text(
        at + width, y, labels = sprintf("%0.2f", y),
        col = pcol, adj = c(-0.25, 0.5)
    )
    #points(rep(at, 2), y, col = pcol, pch = 19, cex = 1.5)
    d <- dim(medcenter)
    bcol <- rep(rep(bcol, each = d[1]), d[3])
    bcol <- if (reverse) rev(bcol) else bcol
    pchr <- rep(pch, each = d[1] * d[2])
    pchr <- if (reverse) rev(pchr) else pchr
    medcenter <- if(reverse) rev(medcenter) else medcenter
    points(rep(at, d[3] * d[2]), medcenter, col = bcol, pch = pchr, cex = 1.2)

    # legend
    centers <- sort(unique(x$Center))
    w <- strwidth(centers) + 0.2
    par(xpd = TRUE)
    x <- c(0, cumsum(w)[-length(w)]) + 0.35
    y <- rep_len(1.04, length(centers))
    text(x, y, labels = centers, adj = c(-0.2, 0.5))
    points(x, y, pch = pch, col = "#808080", cex = 1.2)
}

bp_median_all2 <- function(x, col = c(Control = 3, Sepsis = 2),
                           main, col.main, cex.main = 1.5,
                           ylim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    pch <- c(15, 17, 18, 19)
    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- seq_len(nl)
    xlim <- range(at) + c(-0.5, 0.5)

    medall <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                vapply(split(xx$Score, xx$Interval), median, NA_real_)
        )
    )

    medcenter <- simplify2array(lapply(
        split(x, x$Center),
        function(x2)
            do.call(cbind, lapply(
                split(x2, x2$Diagnosis),
                function(x3)
                    vapply(split(x3$Score, x3$Interval), median, NA_real_)
            ))
    ))

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Time before ICU admission", adj = 0L, col.lab = "#808080")
    title(ylab = "Probability", adj = 1L, col.lab = "#808080")
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)
    lbls <- if (reverse) rev(levels(x$Interval)) else levels(x$Interval)
    axis(
        1L, lwd.ticks = 0L, col = "#808080", at = at, labels = lbls
    )

    y <- medall
    y[] <- if (reverse) rev(medall) else medall
    pcol <- if (reverse) rev(pcol) else pcol
    for (i in seq_len(ncol(medall)))
        lines(at, y[,i], col = pcol[i], type = "b", pch = 20, cex = 2, lwd = 2)
    yl <- y[c(1, nrow(y) + 1)]
    yr <- y[c(nrow(y), 2 * nrow(y))]
    ra <- range(at)
    text(
        ra[1], yl, labels = sprintf("%0.2f", yl),
        col = pcol, adj = c(1.5, 0.5)
    )
    text(
        ra[2], yr, labels = sprintf("%0.2f", yr),
        col = pcol, adj = c(-0.5, 0.5)
    )
    d <- dim(medcenter)
    bcol <- rep(rep(bcol, each = d[1]), d[3])
    bcol <- if (reverse) rev(bcol) else bcol
    pchr <- rep(pch[1:3], each = d[1] * d[2])
    pchr <- if (reverse) rev(pchr) else pchr
    medcenter <- if(reverse) rev(medcenter) else medcenter
    points(rep(at, d[3] * d[2]), medcenter, col = bcol, pch = pchr, cex = 1.2)

    # legend
    centers <- c(sort(unique(x$Center)), "Combined")
    w <- strwidth(centers) + 0.2
    par(xpd = TRUE)
    x <- c(0, cumsum(w)[-length(w)]) + 0.35
    y <- rep_len(1.04, length(centers))
    text(x, y, labels = centers, adj = c(-0.2, 0.5))
    points(x, y, pch = pch, col = "#808080", cex = 1.2)
}

bp_median_all3 <- function(x, col = c(Control = 3, Sepsis = 2),
                           main, col.main, cex.main = 1.5,
                           ylim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    pch <- c(15, 17, 18, 19)
    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    nl <- nlevels(x$Interval)
    at <- seq_len(nl)
    xlim <- range(at) + c(-0.5, 0.5)

    medall <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                vapply(split(xx$Score, xx$Interval), median, NA_real_)
        )
    )

    medcenter <- simplify2array(lapply(
        split(x, x$Center),
        function(x2)
            do.call(cbind, lapply(
                split(x2, x2$Diagnosis),
                function(x3)
                    vapply(split(x3$Score, x3$Interval), median, NA_real_)
            ))
    ))

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Time before ICU admission", adj = 0L, col.lab = "#808080")
    title(ylab = "Probability", adj = 1L, col.lab = "#808080")
    axis(2L, lwd.ticks = 0L, col = "#808080", las = 1L)
    lbls <- if (reverse) rev(levels(x$Interval)) else levels(x$Interval)
    axis(
        1L, lwd.ticks = 0L, col = "#808080", at = at, labels = lbls
    )

    y <- medall
    y[] <- if (reverse) rev(medall) else medall
    pcol <- if (reverse) rev(pcol) else pcol
    for (i in seq_len(ncol(medall)))
        lines(at, y[,i], col = pcol[i], type = "b", pch = 20, cex = 2, lwd = 2)
    yl <- y[c(1, nrow(y) + 1)]
    yr <- y[c(nrow(y), 2 * nrow(y))]
    ra <- range(at)
    text(
        ra[1], yl, labels = sprintf("%0.2f", yl),
        col = pcol, adj = c(1.5, 0.5)
    )
    text(
        ra[2], yr, labels = sprintf("%0.2f", yr),
        col = pcol, adj = c(-0.5, 0.5)
    )
    d <- dim(medcenter)
    for (i in seq_len(d[2]))
        for (j in seq_len(d[3]))
            lines(
                if (reverse) rev(at) else at, medcenter[, i, j],
                col = bcol[i], type = "b", pch = pch[j], cex = 1.2, lwd = 1
            )

    # legend
    centers <- c(sort(unique(x$Center)), "Combined")
    w <- strwidth(centers) + 0.2
    par(xpd = TRUE)
    x <- c(0, cumsum(w)[-length(w)]) + 0.35
    y <- rep_len(1.04, length(centers))
    text(x, y, labels = centers, adj = c(-0.2, 0.5))
    points(x, y, pch = pch, col = "#808080", cex = 1.2)
}

dc_all <- function(x, col = c(Control = 3, Sepsis = 2),
                   main, col.main, cex.main = 1.5,
                   xlim = c(0, 1), reverse = FALSE, ...) {

    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))

    bcol <- setNames(rgba(col, 0.5), names(col))
    pcol <- setNames(col, names(col))
    pch <- c(15, 17, 18, 19)
    cex <- 2
    height = 0.05
    nl <- nlevels(x$Interval)
    lvls <- levels(x$Interval)
    at <- seq_len(nl)
    rat <- rev(at)
    ylim <- range(at) + c(0, 0.5)

    medall <- simplify2array(
        lapply(
            split(x, x$Diagnosis), # could use interaction here but would be harder to construct an array here
            function(xx)
                vapply(split(xx$Score, xx$Interval), median, NA_real_)
        )
    )

    medcenter <- simplify2array(lapply(
        split(x, x$Center),
        function(x2)
            do.call(cbind, lapply(
                split(x2, x2$Diagnosis),
                function(x3)
                    vapply(split(x3$Score, x3$Interval), median, NA_real_)
            ))
    ))

    mai <- par("mai")
    w <- max(strwidth(lvls, "inch")) + 1/8
    if (mai[2L] < w)
        mai[2L] <- mai[4L] + w # taken from dotchart
    par(mai = mai)

    plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, ann = FALSE)
    cmain(main, col.main = col.main, cex.main = cex.main, adj = 0L)
    title(xlab = "Probability", adj = 0L, col.lab = "#808080")
    title(
        ylab = "Time before ICU admission",
        adj = 1L, line = 5L,
        col.lab = "#808080"
    )
    axis(1L, at = c(0, 0.5, 1.0), lwd.ticks = 0L, col = "#808080")
    axis(
        2L, at = if (reverse) rat else at, labels = lvls,
        lty = 0L, lwd.ticks = 0L, col = "#808080", las = 1
    )
    abline(v = 0.5, col = "#808080", lty = "dotted")
    y <- if (reverse) rat else at
    rect(
        xleft = medall[, 1L], xright = medall[, 2L],
        ybottom = y - height, ytop = y + height,
        col = "#E7E7E780", border = 0L
    )
    points(
        medall, rep(y, 2),
        col = rep(pcol, each = nl), pch = pch[4], cex = cex
    )

    d <- dim(medcenter)
    bcol <- rep(rep(bcol, each = d[1]), d[3])
    bcol <- if (reverse) rev(bcol) else bcol
    pchr <- rep(pch[1:3], each = d[1] * d[2])
    pchr <- if (reverse) rev(pchr) else pchr
    medcenter <- if(reverse) rev(medcenter) else medcenter
    points(
        medcenter, rep(at, d[2] * d[3]),
        col = bcol, pch = pchr, cex = 1.25
    )

    # legend
    par(xpd = TRUE)
    centers <- c(sort(unique(x$Center)), "Combined")
    w <- strwidth(centers) + 0.03
    par(xpd = TRUE)
    x <- c(0, cumsum(w)[-length(w)]) - 0.03
    y <- rep_len(max(at) + 0.7, length(centers))
    text(x, y, labels = centers, adj = c(-0.15, 0.5))
    points(x, y, pch = pch, col = "#808080", cex = 1.2)
}

#v <- readRDS("ValidationScw10.RDS")
v$Score <- v$RangerProbability
v$Interval <- cut(
    v$SecToIcu,
    breaks = c(-Inf, 6, 12, 24, 48, 7 * 24, Inf) * 3600,
    labels = c(
        "< 6 h", "6 - 12 h", "12 - 24 h", "1 - 2 d", "2 - 7 d", "> 7 d"
    ),
    include.lowest = TRUE
)
v <- omitNa(v)

ukl <- subset(v, Center == "Leipzig")
ylim <- ylim0(ukl)
control <- subset(ukl, Diagnosis == 0)
sepsis <- subset(ukl, Diagnosis == 1)

pdf("bp.pdf", width = 12, height = 12)
## default

old.par <- par(no.readonly = TRUE)
par(mfcol = c(1, 2))

bp_default(control, main = "Control with ICU Admission", ylim = ylim)
bp_default(sepsis, main = "Sepsis with ICU Admission", ylim = ylim)

par(old.par)

bp_default2(ukl, main = "Sepsis and Control with ICU Admission", ylim = ylim)

bp_color(
    ukl, main = "Sepsis and Control with ICU Admission", ylim = ylim,
    col = nejm
)
bp_color(
    ukl, main = "Sepsis and Control with ICU Admission", ylim = ylim,
    col = vdis
)

bp_color2(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black")
)

bp_color3(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black")
)

bp_median(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black")
)

bp_median(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black"),
    reverse = TRUE
)

bp_median2(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black")
)

bp_median2(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black"),
    reverse = TRUE
)

bp_median3(
    ukl,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black"),
    reverse = TRUE
)

bp_median_all(
    v,
    main = c(
        "Comparison between ", "Sepsis", " and ", "Control",
        " with ICU Admission"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black")
)

bp_median_all(
    v,
    main = c(
        "Probability for ICU admission ", "with", " and ", "without",
        " Sepsis"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black"),
    reverse = TRUE
)

bp_median_all2(
    v,
    main = c(
        "Probability for ICU admission ", "with", " and ", "without",
        " Sepsis"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black"),
    reverse = TRUE
)

bp_median_all3(
    v,
    main = c(
        "Probability for ICU admission ", "with", " and ", "without",
        " Sepsis"
    ),
    col = col,
    col.main = c("black", col["Sepsis"], "black", col["Control"], "black"),
    reverse = TRUE
)

dc_all(
    v,
    main = c(
        "Probability for ICU admission ", "without", " and ", "with",
        " Sepsis"
    ),
    col = col,
    col.main = c("black", col["Control"], "black", col["Sepsis"], "black"),
    xlim = ylim
)

dc_all(
    v,
    main = c(
        "Probability for ICU admission ", "without", " and ", "with",
        " Sepsis"
    ),
    col = col,
    col.main = c("black", col["Control"], "black", col["Sepsis"], "black"),
    xlim = ylim,
    reverse = TRUE
)

par(old.par)

dev.off()
