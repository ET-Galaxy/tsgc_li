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
  methods = list(
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
      model <- .self$get_dynamic_gompertz_model(
        y, q = q, sea.type = sea.type, sea.period = sea.period
      )
      return(model)
    }
  )
)
