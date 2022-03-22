library("ranger")
library("ROCR")
library("future")
library("future.apply")

## helper functions
## calculate case weights for 50:50 selection
.caseweights <- function(y, replace = FALSE) {
    if (replace)
        ## these case weights work just for replace = TRUE (=> bootstrap)
        c(1 / table(y))[y]
    else {
        ## if we don't using bootstrap we have to ensure that the cases
        ## (minority) class is selected, so we choose an arbitary high weight
        tbl <- table(y)
        r <- if (tbl[1] > tbl[2]) c(1, 1e3) else c(1e3, 1)
        names(r) <- names(tbl)
        r[y]
    }
}

## calculate sample fraction for 50:50 selection (two times the minority class)
.samplefraction <- function(y, replace = FALSE) {
    tbl <- table(y)
    as.vector(2L * if (tbl[1] > tbl[2]) tbl[2] else tbl[1]) / length(y)
}

#' Create balanced CV folds
#'
#' taken from ampel-leipzig/ameld/R/utils-glmnet.R
#'
#' @param y `factor`, classes
#' @param nfolds `integer(1)`, number of folds
#' @return integer(length(y))
#' @noRd
.bfolds <- function(y, nfolds = 3L) {
    grpn <- table(y)

    if (nfolds < 3L)
        stop("'nfolds' has to be >= 3.")
    if (any(nfolds > grpn))
        warning("'nfolds' > than the groups, reducing to minimal group size.")

    nfolds <- min(nfolds, grpn)
    s <- seq_len(nfolds)
    unlist(lapply(grpn, function(n)sample(rep_len(s, n))), use.names = FALSE)
}

#' Groupwise mean
#'
#' This function calculates the groupwise mean.
#'
#' taken from ampel-leipzig/ameld/R/utils-base.R
#'
#' @param x `double`, values
#' @param f `factor`, splitting factor/grouping variable.
#' @param na.rm `logical(1)`, see [mean()] for details.
#' @return Name `double` with mean values per split/group.
#' @export
#' @examples
#' groupmean(1:9, rep(1:3, 3))
groupmean <- function(x, f, na.rm = TRUE) {
    vapply(split(x, f), mean.default, NA_real_, na.rm = na.rm)
}

#' Calculate cutpoints/breaks
#'
#' Calculate cutpoints/breaks for `cut` to have an equal number of
#' patients/observations per group.
#'
#' taken from ampel-leipzig/ameld/R/utils-predict.R
#'
#' @param x `double` predicted survival
#' @param n `integer(1)` number of patients/observation per interval
#' @return `double`, cutpoints
#' @importFrom stats quantile
#' @export
#' @examples
#' x <- seq(0, 1, length.out = 10)
#' cutpoints(x, n = 2)
cutpoints <- function(x, n = 50L) {
    unique(
        quantile(
            c(0L, x, 1L),
            seq(0L, 1L, length.out = (floor(length(x) / n) + 1L))
        )
    )
}

## wrapper around ranger with RUS settings
rus_ranger <- function(x, y, replace = TRUE, ...) {
    ranger(
        x = as.data.frame(x), y = as.factor(y),
        probability = TRUE,
        min.node.size = 10,
        ## RUS
        ### subsample (replace == FALSE) vs bootstrap (replace == TRUE)
        replace = replace,
        ### case.weights for RUS are inverse proportions of the class for 50:50
        case.weights = .caseweights(y, replace = replace),
        ### sample.fraction is two times the minority class for 50:50 sampling
        sample.fraction = .samplefraction(y, replace = replace),
        ## END OF RUS
        ...,
        # set to TRUE for debugging inbag sampling
        keep.inbag = FALSE
    )
    ## debug
    #b <- do.call(cbind, r$inbag.counts)
    #c(
    #   Predictions = sum(is.na(r$prediction)),
    #   Sepsis = mean(colSums(b[y == "Sepsis",])),
    #   Control = mean(colSums(b[y == "Control",]))
    #)
}

# y has to be numeric (otherwise ROCR::prediction won't work)
cv_ranger <- function(x, y, nfolds = 5, replace = FALSE, ...) {
    folds <- .bfolds(y, nfolds = nfolds)
    xl <- split(x, folds)
    yl <- split(y, folds)
    auc <- unlist(future.apply::future_lapply(
        seq_len(nfolds),
        function(i) {
            xtrain <- do.call(rbind, xl[-i])
            xtest <- xl[[i]]
            ytrain <- do.call(c, yl[-i])
            ytest <- yl[[i]]
            rngr <- rus_ranger(x = xtrain, y = ytrain, replace = replace, ...)
            pred <- as.numeric(predict(rngr, xtest)$predictions[, 2L])
            performance(prediction(pred, ytest), measure = "auc")@y.values[[1L]]
        },
        future.seed = TRUE
    ))
    median(auc)
}

# repeated cross-validation
rcv_ranger <- function(x, y, nfolds = 5, nrepcv = 2, replace = FALSE, ...) {
    auc <- unlist(future.apply::future_lapply(
        seq_len(nrepcv),
        function(i)cv_ranger(
            x = x, y = y, nfolds = nfolds, replace = replace, ...
        ),
        future.seed = TRUE
    ))
    setNames(
        quantile(auc, names = FALSE), c("Min", "Q1", "Median", "Q3", "Max")
    )
}

# grid search
gridsearch <- function(x, y, searchspace,
                       nfolds = 5, nrepcv = 2, replace = FALSE, ...) {
    auc <- future.apply::future_lapply(
        seq_len(nrow(searchspace)),
        function(i) {
            do.call(
                rcv_ranger,
                c(
                    list(
                        x = x, y = y,
                        nfolds = nfolds, nrepcv = nrepcv, replace = replace
                    ), list(...), searchspace[i, ]
                )
            )
        },
        future.seed = TRUE
    )
    cbind.data.frame(searchspace, do.call(rbind, auc))
}

# nested cross-validation
nrcv_ranger <- function(x, y, searchspace,
                 nouterfolds = 5, ninnerfolds = 5, nrepcv = 2,
                 replace = FALSE, ...) {

    folds <- .bfolds(y, nfolds = nouterfolds)
    xl <- split(x, folds)
    yl <- split(y, folds)

    nrcv <- future.apply::future_lapply(
        seq_len(nouterfolds),
        function(i) {
            xtrain <- do.call(rbind, xl[-i])
            xtest <- xl[[i]]
            ytrain <- do.call(c, yl[-i])
            ytest <- yl[[i]]

            gs <- gridsearch(
                xtrain, ytrain, searchspace,
                nfolds = ninnerfolds, nrepcv = nrepcv,
                ...
            )

            top <- which.max(gs$Median)
            selparms <-
                gs[top,
                   !colnames(gs) %in% c("Min", "Q1", "Median", "Q3", "Max"),
                   drop = FALSE
                ]

            ## additional call of an already calculated tree, ...
            ## could be avoided if we would store the results of the trees
            rngr <- do.call(
                rus_ranger,
                c(
                  list(
                    x = xtrain,
                    y = ytrain,
                    replace = replace
                  ), list(...), selparms
                )
            )
            pred <- as.numeric(predict(rngr, xtest)$predictions[, 2L])

            list(
                model = rngr,
                indextrain = unlist(split(seq_along(y), folds)[-i]),
                indextest = unlist(split(seq_along(y), folds)[i]),
                prediction = pred,
                truth = ytest,
                performance = performance(
                    prediction(pred, ytest), measure = "auc"
                )@y.values[[1L]],
                selectedparams = selparms,
                gridsearch = gs,
                nouterfolds = nouterfolds,
                ninnerfolds = ninnerfolds,
                nrepcv = nrepcv
            )
        },
        future.seed = TRUE
    )
    nrcv
}
