# Created by: Craig Thamotheram
# Created on: 27/07/2022

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

setOldClass("KFS")

#'
#' @title  Class for dynamic Gompertz curve state space model object.
#'
#' @description Class for dynamic Gompertz curve state space model object.
#'
#' \subsection{Methods}{
#' \code{get_model(y, q = NULL, sea.type = 'trigonometric', sea.period = 7)}
#' Retrieves the model object.
#' \subsection{Parameters}{\itemize{
#'  \item{\code{y} The cumulated variable.}
#'  \item{\code{q} The signal-to-noise ratio (ratio of slope to irregular
#' variance). Defaults to \code{'NULL'}, in which case no signal-to-noise ratio
#' will be imposed. Instead, it will be estimated.}
#'  \item{\code{sea.type} Seasonal type. Options are \code{'trigonometric'} and
#' \code{'none'}. \code{'trigonometric'} will yield a model with a trigonometric
#' seasonal component and \code{'none'} will yield a model with no seasonal
#' component.}
#'  \item{\code{sea.period}  The period of seasonality. For a day-of-the-week
#' effect with daily data, this would be 7. Not required if
#' \code{sea.type = 'none'}.}
#' }}
#' \subsection{Return Value}{\code{KFS} model object.}
#' }
#'
#' @importFrom xts periodicity last
#' @importFrom methods new
#' @importFrom magrittr %>%
#' @importFrom KFAS SSMtrend SSMseasonal
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' @export SSModelDynamicGompertz
#' @exportClass SSModelDynamicGompertz
SSModelDynamicGompertz <- setRefClass(
  "SSModelDynamicGompertz",
  contains = "SSModelBase",
  fields = list(
    reinit.date = "ANY",
    original.results = "ANY",
    use.presample.info = "ANY"
  ),
  methods = list(initialize = function(Y, q = NULL, sea.type = 'trigonometric',
                                       sea.period = 7,reinit.date=NULL, original.results=NULL,
                                       use.presample.info=TRUE)
  {
    "Create an instance of the \\code{SSModelDynGompertzReinit} class.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{Y} The cumulated variable.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance). Defaults to \\code{'NULL'}, in which case no
         signal-to-noise ratio will be imposed. Instead, it will be estimated.}
        \\item{\\code{reinit.date} The reinitialisation date \\eqn{r}. Should
        be specified as an object of class \\code{\"Date\"}. Must be specified
        if \\code{original.results = NULL} and
        \\code{use.pre.sample.info = TRUE}.}
        \\item{\\code{original.results} In place of a reinitialisation date, a
        \\code{KFS} results object can be specified here and the parameters for
         the reinitialisation will be taken from this object. Must be specified
          if \\code{reinit.date = NULL} and \\code{use.pre.sample.info = TRUE}.}
        \\item{\\code{use.presample.info}  Logical value denoting whether or not
         to use information from before the reinitialisation date in the
         reinitialisation procedure. Default is \\code{TRUE}. If \\code{FALSE},
         the model is estimated from scratch from the reinitialisation date and
         no attempt to use information from before the reinitialisation date is
          made.}
      }}
      \\subsection{Usage}{\\code{SSModelDynGompertzReinit$new(y, q = 0.005,
      reinit.date = as.Date(\"2021-05-12\",format = date.format))}}"
    callSuper(Y, q, sea.type, sea.period)
    reinit.date <<- reinit.date
    original.results <<- original.results
    use.presample.info <<- use.presample.info
  },
    get_model = function(
    y,
    q = NULL,
    sea.type = 'trigonometric',
    sea.period = 7
    )
    {
      "Retrieves the model object.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{y} The cumulated variable.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
        variance). Defaults to \\code{'NULL'}, in which case no signal-to-noise
        ratio will be imposed. Instead, it will be estimated.}
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}. \\code{'trigonometric'} will
         yield a model with a trigonometric seasonal component and
         \\code{'none'} will yield a model with no seasonal component.}
        \\item{\\code{sea.period}  The period of seasonality. For a
        day-of-the-week effect with daily data, this would be 7. Not required
        if \\code{sea.type = 'none'}.}
      }}
      \\subsection{Return Value}{\\code{KFS} model object.}"
      if (is.null(.self$reinit.date)){
        model <- .self$get_dynamic_gompertz_model(
          y, q = q, sea.type = sea.type, sea.period = sea.period
        )
        return(model)
      }
      else{
        # 4.1. Index for reinitialisation, t_0
        stopifnot(length(Y[reinit.date]) == 1)
        Y.t.r_0 <- as.numeric(Y[reinit.date - 1])
        
        # 4.2 Reinitialisation:
        #   ln g_t^r = ln g_t + ln (Y_{t-1}/Y_{t-1}^r), where Y_t^r=Y_t-Y_{r_0}.
        idx.dates <- (index(y) > reinit.date)
        lag.Y <- stats::lag(Y)[idx.dates]
        y.reinit <- y[index(y) > reinit.date] + log(lag.Y / (lag.Y - Y.t.r_0))
        
        # 4.3 Run Kalman filter/smoother on new series with non-diffuse prior
        if (use.presample.info) {
          # Either estimate full model here or take results from previous model.
          if (is.null(original.results)) {
            # NB. Restrict sample to t<=r - date of reinitialisation.
            idx.est <- zoo::index(Y) <= reinit.date
            model <- SSModelDynamicGompertz$new(Y = Y[idx.est], q = q)
            res.original <- model$estimate()
            model_output <- output(res.original)
          } else {
            model_output <- output(original.results)
            model_seasonal <- seasonalComp(original.results)
            sea.type <- if (is.null(model_seasonal)) {'none'} else {
              'trigonometric'}
            sea.period <- if (!is.null(model_seasonal)) {
              ncol(att(model_output)) - 1}
          }
          
          # 4.3 Reset slope to 0 and add constant to initial value for level.
          # where reinit.date is t=r
          idx <- which(reinit.date == index(y))
          stopifnot(length(idx) == 1)
          att <- att(model_output)[idx,]
          Ptt <- Ptt(model_output)[, , idx]
          Tt <- drop(matrixKFS(model_output,"T"))
          Rt <- drop(matrixKFS(model_output,"R"))
          Qt <- drop(matrixKFS(model_output,"Q"))
          Ht <- drop(matrixKFS(model_output,"H"))
          
          # a. Take a_{r|r} and P_{r|r} through prediction step to get a_{r+1}
          # and P_{r+1}
          a1 <- Tt %*% att
          P1 <- Tt %*% Ptt %*% t(Tt) + Rt %*% Qt %*% t(Rt)
          
          # b. Set slope to 0 and add correction (\ln(Y_r/y_r) to level.
          a1["slope",] <- 0
          a1["level",] <- a1["level",] + log(Y[idx] / (Y[idx] - Y.t.r_0))
        } else {
          # Don't use presample info
          a1 <- NULL; P1 <- NULL; Qt <- NULL; Ht <- NULL
        }
        out <- .self$get_dynamic_gompertz_model(
          y = y.reinit, q = q, sea.type = sea.type, sea.period = sea.period,
          a1 = a1, P1 = P1, Q = Qt, H = Ht
        )
        out[['index']] <- index(y.reinit)
        return(out)
      }
    },
    summary = function() {
      out <- suppressWarnings(output(.self$estimate()))
      q<-.self$q
      if(is.null(q)){
        qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      }
      reinit<-!is.null(.self$reinit.date)
      dates<-index(.self$Y)
      
      cat("Summary of SSModelDynamicGompertz Model")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("--------------------------------------\n")
      cat("Cumulated Variable:\n")
      base::print(head(.self$Y))
      cat("Signal-to-Noise Ratio (q):", 
          ifelse(is.null(q), paste(signif(qest,3), "(estimated)"), 
                 paste(q, ("(user specified)"))), "\n")
      cat("Model Details:\n")
      cat("  - Model Type: Dynamic Gompertz Curve")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("  - Seasonal Component: ", ifelse(sea.type == 'none', "None", "Trigonometric"), "\n")
      cat("  - Period of Seasonality: ", ifelse(sea.type == 'none', "N/A", sea.period), "\n")
      cat("  - Dataset start date:", format(as.Date(dates[1], origin = "1970-01-01")))
      cat("\n")
      cat("  - Dataset end date:", format(as.Date(tail(dates,1), origin = "1970-01-01")))
      cat("\n")
      if (reinit){
        cat("  - Reinitialization date:",format(as.Date(.self$reinit.date, origin = "1970-01-01")))
        cat("\n")
        cat("  - Use presample info:", .self$use.presample.info)
        cat("\n")
      }
      cat("  - Model States and Standard Errors\n")
      base::print(out)
    },
    print = function() {
      out <- output(.self$estimate()) #KFS object
      if(is.null(.self$q)){
        qest <- matrixKFS(out,"Q")[2, 2, 1]/matrixKFS(out,"H")[, , 1]
      }
      cat("SSModelDynamicGompertz Model")
      if (reinit) {
        cat(" (Reinitialized)")
      }
      cat("\n")
      cat("\n")
      cat("Cumulated Variable:\n")
      base::print(head(.self$Y))
      cat("Number of observations:", length(.self$Y))
      cat("\n")
      cat("Signal-to-Noise Ratio (q):", 
          ifelse(is.null(.self$q), paste(signif(qest,5), "(estimated)"), 
                 paste(.self$q, ("(user specified)"))), "\n")
      cat("Seasonal components?",
          ifelse(is.null(seasonalComp(out)),
                 "No","Yes"),"\n")
      if (!is.null(.self$reinit.date)){
        cat("Reinit date:",format(as.Date(.self$reinit.date, origin = "1970-01-01")))
        cat("\n")
        cat("Use presample info:", .self$use.presample.info)
      }
    },
    plot_diff =function(title=NULL, ylab=NULL){
      plot(diff(.self$Y), main=title, ylab=ylab)
    }
  )
)
