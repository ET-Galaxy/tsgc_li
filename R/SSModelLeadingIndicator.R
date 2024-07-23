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
#' @title  Class for leading indicator state space model object.
#'
#' @description Class for leading indicator state space model object.
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
SSModelLeadingIndicator <- setRefClass(
  "SSModelLeadingIndicator",
  fields = list(
    Y = "xts",
    q = "ANY",
    n.lag = "integer",
    leading_indicator_col= "integer"
  ),
  methods = list(
    initialize = function(Y, n.lag, LeadIndCol=1, q = NULL)
    {
      "Create an instance of the \\code{SSModelLeadingIndicator} class.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{Y} The cumulated leading indiactor variable and target variable.}
        \\item{\\code{start_date} The start date \\eqn{r}. Should
        be specified as an object of class \\code{\"Date\"}.}
        \\item{\\code{end_date} The end date \\eqn{r}. Should
        be specified as an object of class \\code{\"Date\"}. Must be specified
        if \\code{original.results = NULL} and
        \\code{use.pre.sample.info = TRUE}.}
        \\item{\\code{n.lag} Number of days to lag the leading indicator.}
        \\item{\\code{LeadIndCol} The column in \\code{Y} that contains the leading indicator}
        \\item{\\code{q} The signal-to-noise ratio (ratio of slope to irregular
         variance). Defaults to \\code{'NULL'}, in which case no
         signal-to-noise ratio will be imposed. Instead, it will be estimated.}
      }}
      \\subsection{Usage}{\\code{SSModelDynGompertzReinit$new(y, q = 0.005,
      reinit.date = as.Date(\"2021-05-12\",format = date.format))}}"
      reinit.date <<- reinit.date
      original.results <<- original.results
      use.presample.info <<- use.presample.info
      callSuper(Y, q)
    },
    get_model = function(
      y,
      n.lag,
      q = NULL,
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
