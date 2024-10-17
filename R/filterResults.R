setOldClass("KFS")
#'
#' @title FilterResults
#'
#' @description Class for estimated Dynamic Gompertz Curve model and contains
#' methods to extract smoothed/filtered estimates of the states, the level of
#' the incidence variable \eqn{y}, and forecasts of \eqn{y}.
#' @references Harvey, A. C. and Kattuman, P. (2021). A Farewell to R:
#' Time Series Models for Tracking and
#' Forecasting Epidemics, Journal of the Royal Society Interface, vol 18(182):
#' 20210179
#'
#' @importFrom xts periodicity last
#' @importFrom magrittr %>%
#' @importFrom methods new
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#' # Print estimation results
#' res$print_estimation_results()
#' # Forecast 7 days ahead from the end of the estimation window
#' res$predict_level(y.cum = gauteng[idx.est], n.ahead = 7,
#'   confidence_level = 0.68)
#' # Forecast 7 days ahead from the model and return filtered states
#' res$predict_all(n.ahead = 7, return.all = TRUE)
#' # Return the filtered growth rate and its components
#' res$get_growth_y(return.components = TRUE)
#' # Return smoothed growth rate of incidence variable and its confidence
#' # interval
#' res$get_gy_ci(smoothed = TRUE, confidence_level = 0.68)
#'
#' @export
#'
FilterResults <- setRefClass(
  "FilterResults",
  fields = list(
    index = "Date",
    output = "KFS"
  ),
  methods = list(
    initialize = function(index, output)
    {
      index <<- index
      output <<- output
    },
    predict_level = function(
      y.cum,
      n.ahead,
      confidence_level,
      sea.on = FALSE,
      return.diff = FALSE)
    {
      "Forecast the cumulated variable or the incidence of it. This function returns
      the forecast of the cumulated variable \\eqn{Y}, or the forecast of the incidence of the cumulated variable, \\eqn{y}. For
      example, in the case of an epidemic, \\eqn{y} might be daily new cases of
      the disease and
       \\eqn{Y} the cumulative number of recorded infections.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{y.cum} The cumulated variable.}
        \\item{\\code{n.ahead} The number of periods ahead you wish to forecast from
        the end of the estimation window.}
        \\item{\\code{confidence_level} The confidence level for the log growth
         rate that should be used to compute
        the forecast intervals of \\eqn{y}.}
        \\item{\\code{return.diff} Logical value indicating whether to return the cumulated variable,
        \\eqn{Y}, or the incidence of it,
        \\eqn{y} (i.e., the first difference of the cumulated variable). Default is
        \\code{FALSE}.}
      }}
      \\subsection{Return Value}{\\code{xts} object containing the point
      forecasts and upper and lower bounds of
      the forecast interval.}"
      model <- output$model
      n <- attr(model, "n")
      p <- attr(model, "p")

      freq <- unclass(periodicity(y.cum))$label
      endtime <- end(model$y) + c(0, n.ahead)
      filtered.out <- .self$predict_all(n.ahead, sea.on = sea.on,
                                        return.all = FALSE)

      # # 1. Extract parameters.
      timespan <- n + 0:n.ahead

      # Calculate g.t as exponent of y.t
      g.t <- exp(filtered.out$y.hat)

      # Forecast dates
      v_dates_end <- seq(last(index(y.cum)), last(
        index(filtered.out$y.hat)), by = freq
      )

      # Construct CI
      Ptt <- filtered.out$P.t.t
      i.level <- grep("level", colnames(model$T))
      ci <- qnorm((1 - confidence_level) / 2) *
        sqrt(Ptt[i.level, i.level, timespan]) %o% c(1, -1)
      ci <- xts(exp(ci), order.by = v_dates_end)

      y.hat <- xts(matrix(NA, nrow = n.ahead + 1, ncol = 3),
                   order.by = v_dates_end)
      y.hat[v_dates_end[1],] <- y.cum[v_dates_end[1]]
      for (i in seq_len(length(v_dates_end[-1]))) {
        date.forecast <- v_dates_end[i + 1]
        date.lag <- date.forecast - 1
        # Update level
        y.hat[date.forecast, 1] <- as.numeric(y.hat[date.lag, 1]) *
          as.numeric(1 + g.t[date.forecast,])

        # Make prediction intervals
        y.hat[date.forecast, 2:3] <- as.numeric(y.hat[date.lag, 1]) *
          as.numeric(1 + g.t[date.forecast,] %*% ci[date.forecast,])
      }

      # Difference output if requested
      d <- if (return.diff) { diff(y.hat[, 1])[-1] } else { (y.hat[, 1])[-1] }

      ci_bounds <- if (return.diff) {
        (y.hat[, 2:3] - as.vector(y.hat[, 1]))[-1] + as.vector(d)
      } else { y.hat[2:dim(y.hat)[1], 2:3] }

      pred <- vector("list", length = p)
      pred[[p]] <- cbind(fit = d, prediction = ci_bounds)
      pred <- lapply(pred, ts, end = endtime, frequency = 1)

      y.hat <- xts(pred[[p]], order.by = v_dates_end[-1])
      names(y.hat)[2:3] <- list('lower', 'upper')

      return(as.xts(y.hat))
    },
    print_estimation_results = function() {
      "Prints a table of estimated parameters in a format ready to paste into
      LaTeX."
      H <- output$model$H[, , 1]
      Q_gamma <- output$model$Q[2, 2, 1]
      Q_seasonal <- output$model$Q[3, 3, 1]

      tbl <- data.frame(
        a = format(H, digits = 3),
        b = format(Q_gamma, digits = 3),
        c = format(Q_seasonal, digits = 3),
        d = format(Q_gamma / H, digits = 4))
      header.names <- c('$\\sigma_\\varepsilon^2$',
                        '$\\sigma_\\gamma^2$',
                        '$\\sigma_{seas}^2$',
                        'q')

      out <- tbl %>%
        kableExtra::kbl(
          caption = "Estimated parameters",
          col.names = header.names,
          format = 'latex',
          booktabs = TRUE,
          escape = FALSE
        ) %>%
        kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria") %>%
        kableExtra::footnote(general = " ")

      return(out)
    },
    predict_all = function(n.ahead, sea.on = FALSE, return.all = FALSE) {
      "Returns forecasts of the incidence variable \\eqn{y}, the state variables
       and the conditional covariance matrix
      for the states.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.ahead} The number of forecasts you wish to create from
        the end of your sample period.}
        \\item{\\code{sea.on} Logical value indicating whether seasonal
        components should be included in the
        state-space model or not. Default is \\code{TRUE}.}
        \\item{\\code{return.all} Logical value indicating whether to return
        all filtered estimates and forecasts
        (\\code{TRUE}) or only the forecasts (\\code{FALSE}). Default is
        \\code{FALSE}.}
      }}
      \\subsection{Return Value}{\\code{xts} object containing the forecast
      (and filtered, where applicable) level
      of \\eqn{y} (\\code{y.hat}), \\eqn{\\delta} (\\code{level.t.t}),
      \\eqn{\\gamma} (\\code{slope.t.t}), vector of states including the
      seasonals where applicable (\\code{a.t.t}) and covariance matrix of all
      states including seasonals where applicable (\\code{P.t.t}).}"
      idx <- index
      model <- output

      new.model <- output$model
      new.model$y <- rbind(
        new.model$y,
        matrix(NA, ncol = ncol(new.model$y), nrow = n.ahead) %>% as.ts()
      )
      attr(new.model, 'n') <- length(output$model$y) + n.ahead %>%
        as.integer()
      model_output <- KFS(new.model)

      if (sea.on == TRUE) {
        y.hat.kfas <- predict(
          output$model, interval = 'prediction',
          n.ahead = n.ahead, level = 0.68, states = 'all')
      } else {
        y.hat.kfas <- predict(
          output$model, interval = 'prediction',
          n.ahead = n.ahead, level = 0.68, states = 'level')
      }

      n <- attr(model$model, "n")
      dates <- seq(idx[1], by = 'day', length.out = (n + n.ahead))

      # Assumes time invariant Z.t
      y.t.t <- output$att %*% drop(output$model$Z)

      y.hat <- xts::xts(
        rbind(y.t.t, y.hat.kfas[, 1] %>% as.matrix()),
        order.by = dates)

      i.level <- grep("level", colnames(model_output$att))
      level.t.t <- xts::xts(model_output$att[, i.level], order.by = dates) %>%
        as.xts()
      i.slope <- grep("slope", colnames(model_output$att))
      slope.t.t <- xts::xts(model_output$att[, i.slope], order.by = dates) %>%
        as.xts()

      if (!return.all) {
        y.hat <- y.hat %>%
          subset(index(.) > tail(idx, 1))
        level.t.t <- level.t.t %>%
          subset(index(.) > tail(idx, 1))
        slope.t.t <- slope.t.t %>%
          subset(index(.) > tail(idx, 1))
      }

      out <- list(
        y.hat = y.hat,
        level.t.t = level.t.t,
        slope.t.t = slope.t.t,
        a.t.t = model_output$att,
        P.t.t = model_output$Ptt
      )
      return(out)
    },
    get_growth_y = function(smoothed = FALSE, return.components = FALSE) {
      "Returns the growth rate of the incidence (\\eqn{y}) of the cumulated
      variable (\\eqn{Y}). Computed as
      \\deqn{g_t = \\exp\\{\\delta_t\\}+\\gamma_t.}
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{smoothed} Logical value indicating whether to use the
        smoothed estimates of \\eqn{\\delta} and \\eqn{\\gamma} to compute the
        growth rate (\\code{TRUE}), or the contemporaneous filtered estimates
        (\\code{FALSE}). Default is \\code{FALSE}.}
        \\item{\\code{return.components} Logical value indicating whether to
        return the estimates of \\eqn{\\delta} and \\eqn{\\gamma} as well as
        the estimates of the growth rate, or just the growth rate. Default is
        \\code{FALSE}.}
      }}
      \\subsection{Return Value}{\\code{xts} object containing
      smoothed/filtered growth rates and components (\\eqn{\\delta} and
      \\eqn{\\gamma}), where applicable.}"
      kfs_out <- output
      idx <- index

      if (smoothed) {
        att <- kfs_out$alphahat
      } else {
        att <- kfs_out$att
      }

      filtered_slope <- xts(att[, "slope"], order.by = idx)
      filtered.level <- xts(att[, "level"], order.by = idx)
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope
      names(gy.t) <- if (smoothed) { "smoothed gy.t" } else { "filtered gy.t" }
      names(g.t) <- if (smoothed) { "smoothed g.t" } else { "filtered g.t" }
      names(filtered_slope) <- if (smoothed) { "smoothed gamma.t" } else {
        "filtered gamma.t" }
      if (return.components) {
        return(list(gy.t, g.t, filtered_slope))
      } else {
        return(gy.t)
      }
    },
    get_gy_ci = function(smoothed = FALSE, confidence_level = 0.68) {
      "Returns the growth rate of the incidence (\\eqn{y}) of the cumulated
      variable (\\eqn{Y}). Computed as
      \\deqn{g_t = \\exp\\{\\delta_t\\}+\\gamma_t.}
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{smoothed} Logical value indicating whether to use the
        smoothed estimates of \\eqn{\\delta} and \\eqn{\\gamma} to compute the
        growth rate (\\code{TRUE}), or the contemporaneous filtered estimates
        (\\code{FALSE}). Default is \\code{FALSE}.}
        \\item{\\code{confidence_level} Confidence level for the confidence
        interval.  Default is \\eqn{0.68}, which is one standard deviation for
        a normally distributed random variable.}
      }}
      \\subsection{Return Value}{\\code{xts} object containing smoothed/filtered
       growth rates and upper and lower bounds for the confidence intervals.}"

      kfs_out <- output
      idx <- index

      if (smoothed) {
        att <- kfs_out$alphahat
      } else {
        att <- kfs_out$att
      }

      filtered_slope <- xts(att[, "slope"], order.by = idx)
      filtered.level <- xts(att[, "level"], order.by = idx)
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope

      idx.slope <- grep("slope", colnames(kfs_out$att))
      ci <- qnorm((1 - confidence_level) / 2) *
        sqrt(kfs_out$Ptt[idx.slope, idx.slope,]) %o% c(1, -1)
      ci_bounds <- as.vector(gy.t) + ci

      pred <- xts(cbind(gy.t, ci_bounds), order.by = idx)
      colnames(pred) <- c("fit","lower","upper")

      return(pred)
    },
    print=function(){
      H <- output$model$H[, , 1]
      Q_gamma <- output$model$Q[2, 2, 1]
      Q_seasonal <- output$model$Q[3, 3, 1]
      out <- .self$output
      cat("Fitted SSModelDynamicGompertz Model Results\n")
      cat("\n")
      cat("  - Model States and Standard Errors\n")
      base::print(out)
      cat("- Variance parameter estimates\n")
      cat("Observation equation noise:",format(H, digits = 4))
      cat("\n")
      cat("State transition equation noise:",format(Q_gamma, digits = 4))
      cat("\n")
      cat("Signal-to-Noise Ratio (q):", format(Q_gamma / H, digits = 4))
      cat("\n")
      cat("Seasonality noise:",format(Q_seasonal, digits = 4))
    }
  )
)
