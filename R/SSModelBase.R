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

setOldClass("xts")
#'
#' @title Base class for estimating time-series growth curve models. Classes
#' \code{SSModelDynamicGompertz} and \code{SSModelDynGompertzReinit} refer back
#' to this base class.
#' @importFrom xts periodicity last
#' @importFrom magrittr %>%
#' @importFrom xts xts
#' @importFrom zoo index
#' @importFrom methods new
#' @importFrom KFAS SSModel fitSSM KFS
#' @examples
#' library(tsgc)
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-06")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#' @export SSModelBase
#' @exportClass SSModelBase
SSModelBase <- setRefClass(
  "SSModelBase",
  fields = list(
    Y = "xts",
    q = "ANY"  # No native option for numeric | NULL - see
    # https://stackoverflow.com/questions/24363069/multiple-acceptable-classes-
    # in-reference-class-field-list
  ),
  methods = list(
    initialize = function(Y, q = NULL)
    {
      Y <<- Y
      q <<- q
    },
    get_model = function(y, q = NULL)
    {
      return(stop('Not Implemented Error: Implement in subclass with your model'))
    },
    get_dynamic_gompertz_model = function(
      y,
      q = NULL,
      sea.type = 'trigonometric',
      sea.period = 7,
      a1 = NULL,
      P1 = NULL,
      Q = NULL,
      H = NULL
    )
    {
      "Returns dynamic Gompertz curve model.
    \\subsection{Parameters}{\\itemize{
      \\item{\\code{y} The cumulated variable}
      \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
      variance). Defaults to \\code{'NULL'}, in which case no signal-to-noise
      ratio will be imposed. Instead, it will be estimated.}
      \\item{\\code{sea.type} Seasonal type. Options are \\code{'trigonometric'}
       and \\code{'none'}. \\code{'trigonometric'} will yield a model with a
       trigonometric seasonal component and \\code{'none'} will yield a model
       with no seasonal component.}
      \\item{\\code{sea.period} The period of seasonality. For a day-of-the-week
       effect with daily data, this would be 7. Not required if
       \\code{sea.type = 'none'}.}
      \\item{\\code{a1} Optional parameter specifying the prior mean of the
      states. Defaults to \\code{'NULL'}. Leave as \\code{'NULL'} for a diffuse
      prior (no prior information). If a proper prior is to be specified, both
      \\code{a1} and \\code{P1} must be given.}
      \\item{\\code{P1} Optional parameter specifying the prior mean of the
      states. Defaults to \\code{'NULL'}. Leave as \\code{'NULL'} for a diffuse
       prior (no prior information). If a proper prior is to be specified,
       both \\code{a1} and \\code{P1} must be given.}
      \\item{\\code{Q} Optional parameter specifying the state error variances
      where these are to be imposed rather than estimated. Defaults to
      \\code{'NULL'} which will see the variances estimated.}
      \\item{\\code{H} Optional parameter specifying the irregular variance
      where this is to be imposed rather than estimated. Defaults to
      \\code{'NULL'} which will see the variance estimated.}
    }}
    \\subsection{Description}{
    The dynamic Gompertz with an integrated random walk (IRW) trend is
    \\deqn{\\ln g_{t}=\\delta_{t}+\\varepsilon_{t},  \\;\\;\\;\\;
    \\varepsilon_{t}\\sim NID(0,\\sigma_{\\varepsilon }^{2}), \\;\\;\\;\\;
    t=2,...,T, }
    where \\eqn{Y_t} is the cumulated variable, \\eqn{y_t = \\Delta Y_t},
    \\eqn{\\ln g_{t}=\\ln y_{t}-\\ln Y_{t-1}} and
    \\deqn{\\delta_{t} =\\delta_{t-1}+\\gamma_{t-1},}
    \\deqn{\\gamma_{t} =\\gamma_{t-1}+\\zeta_{t}, \\;\\;\\;\\;
    \\zeta_{t}\\sim NID(0,\\sigma_{\\zeta }^{2}),}
    where the observation disturbances \\eqn{\\varepsilon_{t}}  and slope
    disturbances \\eqn{\\zeta_{t}}, are iid Normal and mutually independent.
    Note that, the larger the signal-to-noise ratio,
    \\eqn{q_{\\zeta }=\\sigma_{\\zeta }^{2}/\\sigma_{\\varepsilon }^{2}},
    the faster the slope changes in response to new observations. Conversely,
    a lower signal-to-noise ratio induces smoothness.

    For the model without seasonal terms (\\code{sea.type = 'none'}) the are
    priors are
    \\deqn{\\begin{pmatrix} \\delta_1 \\ \\gamma_1 \\end{pmatrix}
    \\sim N(a_1,P_1)}.
    The diffuse prior has \\eqn{P_1 = \\kappa I_{2\\times 2}} with
    \\eqn{\\kappa \\to \\infty}. Implementation of the diffuse prior is handled
     by the package \\code{KFAS} (Helske, 2017). Where the model has a seasonal
      component (\\code{sea.type = 'trigonometric'}), the vector of prior means
       \\eqn{a_1} and the prior covariance matrix \\eqn{P_1} are extended
       accordingly.

    See the vignette for details of the variance matrix \\eqn{Q}.
    \\eqn{H = \\sigma^2_{\\varepsilon}}.
    }
    "
      Qt.slope <- if (is.null(Q)) { NA } else { Q[2, 2] }
      Qt.seas <- if (is.null(Q)) { NA } else { Q[3, 3] }
      Ht <- if (is.null(H)) { NA } else { H }

      # 1. Set prior on state as ~ N(a1, P1) if a1 supplied.
      use.prior <- if (!is.null(a1)) { TRUE } else { FALSE }

      if (use.prior) {
        if (sea.type == 'trigonometric') {
          ss_model <- SSModel(
            y ~
              SSMtrend(
                degree = 2,
                Q = list(matrix(0), matrix(Qt.slope)),
                a1 = a1[1:2],
                P1 = P1[1:2, 1:2]
              ) +
                SSMseasonal(
                  period = sea.period,
                  Q = Qt.seas,
                  sea.type = sea.type,
                  a1 = a1[3:dim(a1)[1]],
                  P1 = P1[3:dim(a1)[1], 3:dim(a1)[1]]
                ),
            H = Ht
          )
          n.pars <- 0
        } else if (sea.type == 'none') {
          ss_model <- SSModel(
            y ~
              SSMtrend(
                degree = 2,
                Q = list(matrix(0), matrix(Qt.slope)),
                a1 = a1[1:2],
                P1 = P1[1:2, 1:2]
              ),
            H = Ht
          )
          n.pars <- 0
        } else {
          stop(sprintf("sea.type= '%s' not implemented", sea.type))
        }
      } else {
        if (sea.type == 'trigonometric') {
          ss_model <- SSModel(
            y ~
              SSMtrend(
                degree = 2,
                Q = list(matrix(0), matrix(Qt.slope))
              ) +
                SSMseasonal(
                  period = sea.period,
                  Q = Qt.seas,
                  sea.type = sea.type),
            H = matrix(Ht)
          )
          n.pars <- if (is.null(q)) { 3 } else { 2 }
        } else if (sea.type == 'none') {
          ss_model <- SSModel(
            y ~
              SSMtrend(
                degree = 2,
                Q = list(matrix(0), matrix(Qt.slope))
              ),
            H = matrix(Ht)
          )
          n.pars <- if (is.null(q)) { 2 } else { 1 }
        } else {
          stop(sprintf("sea.type= '%s' not implemented", sea.type))
        }
      }
      out <- list(model = ss_model, inits = rep(0, n.pars))
      return(out)
    },
    update = function(pars, model, q, sea.type) {
      "Update method for Kalman filter to implement the dynamic Gompertz curve
       model.
       A maximum of 3 parameters are used to set the observation noise
       (1 parameter), the transition equation slope and seasonal noise. If q (signal
        to noise ratio) is not null then the slope noise is set using this
        ratio.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{pars} Vector of parameters.}
        \\item{\\code{model} \\code{KFS} model object.}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance).}
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}.}
      }}
      \\subsection{Return Value}{\\code{KFS} model object.}"
      estH <- any(is.na(model$H))
      estQ <- any(is.na(model$Q))
      if ((!estH) & (!estQ)) {
        # If nothign to update then return model
        return(model)
      } else {
        nparQ <- if (sea.type == 'trigonometric') { 1 } else { 0 }
        # 1. Set seasonal noise
        if (estQ) {
          Q <- as.matrix(model$Q[, , 1])
          # Update diagonal elements
          naQd <- which(is.na(diag(Q)))
          Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
          diag(Q)[naQd] <- exp(0.5 * pars[1])
          # Check for off-diagonal elements and raise error if found.
          naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
          if (length(naQnd) > 0) {
            stop("NotImplmentedError: Unexpected off-diaganol element updating")
          }
        }

        # 2. Set observation noise
        H <- as.matrix(model$H[, , 1])
        if (estH) {
          naHd <- which(is.na(diag(H)))
          H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
          diag(H)[naHd] <- exp(0.5 * pars[(nparQ + 1)])
          model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
        }

        # 3. Set slope noise
        # Get index of slope, 1 before the seasonal component.
        model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
        i.slope <- 2
        # Estimate slope if no signal to noise ratio specified.
        if (is.null(q)) {
          Q.slope <- exp(0.5 * pars[(nparQ + 2)])
        } else {
          Q.slope <- crossprod(H[naHd, naHd]) * q
        }
        model$Q[i.slope, i.slope, 1] <- Q.slope
      }
      return(model)
    },
    estimate = function(sea.type = 'trigonometric', sea.period = 7) {
      "Estimates the dynamic Gompertz curve model when applied to an object of
      class \\code{SSModelDynamicGompertz} or \\code{SSModelDynGompertzReinit}.
      \\subsection{Parameters}{\\itemize{
        \\item{\\code{sea.type} Seasonal type. Options are
        \\code{'trigonometric'} and \\code{'none'}. \\code{'trigonometric'} will
         yield a model with a trigonometric seasonal component and
         \\code{'none'} will yield a model with no seasonal component.}
        \\item{\\code{sea.period} The period of seasonality. For a
        day-of-the-week effect with daily data, this would be 7. Not required
        if \\code{sea.type = 'none'}.}
      }}
      \\subsection{Return Value}{An object of class \\code{FilterResults}
      containing the result output for the estimated dynamic Gompertz curve
      model.}
      "
      # 1. Get LDL of cumulative series Y.
      y <- tsgc::df2ldl(Y)

      # 2. Add update / model methods
      updatefn <- purrr::partial(
        .self$update, ... =, q = q, sea.type = sea.type
      )
      model <- .self$get_model(y, q = q, sea.type, sea.period)
      # 2. Estimate via MLE unknown params
      model_fit <- fitSSM(model$model, inits = model$inits, updatefn = updatefn,
                          method = 'BFGS')

      # 3. Run smoother/filter
      model_output <- KFS(model_fit$model)

      # 4. Get truncated index from model if using a reinitialisation in
      # self$get_model
      date.index <- if (!is.null(model$index)) { model$index } else { index(y) }

      results <- FilterResults$new(
        index = date.index,
        output = model_output
      )
      return(results)
    }
  )
)
