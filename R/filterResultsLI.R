setOldClass("KFS")
#'
#' @title FilterResultsLI
#'
#' @description Class for estimated Leading Indicator Gompertz Curve model and
#' contains methods to extract smoothed/filtered estimates of the states, the
#' level of the incidence variable \eqn{y}, and forecasts of \eqn{y}.
#' @references Harvey, A. C. and Kattuman, P. (2021).
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
FilterResultsLI <- setRefClass(
  "FilterResultsLI",
  fields = list(
    data_xts = "ANY",
    fitmod = "KFS",
    n.lag="numeric",
    sea.period="numeric",
    LeadIndCol="numeric"
  ),
  methods = list(
    initialize = function(data_xts, fitmod,n.lag,sea.period,LeadIndCol)
    {
      data_xts <<- data_xts
      fitmod <<- fitmod
      n.lag <<- n.lag
      sea.period <<- sea.period
      LeadIndCol<<-LeadIndCol
    },
    predict_level = function(n.forc=n.lag, confidence.level=0.68){
      "Forecast the cumulated variable or the incidence of it. This function returns
      the forecast of the cumulated variable \\eqn{Y}, or the forecast of the incidence of the cumulated variable, \\eqn{y}. For
      example, in the case of an epidemic, \\eqn{y} might be daily new cases of
      the disease and
       \\eqn{Y} the cumulative number of recorded infections.
       \\subsection{Parameters}{\\itemize{
        \\item{\\code{n.forc} The number of periods ahead you wish to forecast from
        the end of the estimation window.}
        \\item{\\code{confidence_level} The confidence level for the log growth
         rate that should be used to compute
        the forecast intervals of \\eqn{y}.}
        }
      }
      \\subsection{Return Value}{A list object containing n.lag and 2 \\code{xts}
      objects: the point forecasts and upper and lower bounds of the forecast interval
      for trend and forecast with seasonal component.}"
      if (n.forc==1){
        n.forc=2
        unity=TRUE
      } else{
        unity=FALSE
      }
      out<-fitmod
        start_date<-index(data_xts)[1]
        end_date<-tail(index(data_xts),1)
        data_ldl = data_xts[,c("LDLcases","LDLhosp")] %>% na.omit

        data_ldl$LDLcases = lag(as.vector(data_ldl$LDLcases),n.lag)

        data_ldl <- na.omit(data_ldl)

        data_mat = as.matrix(data_ldl)
        # Create forecast data (using fact have some "future" case observations)
        forcdata <- matrix(NA,ncol=2,nrow=max(n.forc,n.lag))
        colnames(forcdata) = colnames(data_mat)
        forcdata[1:n.lag,1] = as.vector(tail(data_xts,n.lag)$LDLcases)

        # Extract estimate of Q from earlier model
        Qf = out$model$Q[,,1]

        # Create forecast model object
        # Same structure as model before but NAs replaced with the parameter estimates from earlier
        if (is.na(sea.period)) {
          forcmodel = SSModel(forcdata ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1),
                              H = out$model$H)
        }
        else{
          forcmodel = SSModel(forcdata ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+SSMseasonal(sea.period,Q = matrix(c(Qf[4,4],0,0,Qf[5,5]),2,2), sea.type='dummy', type='distinct')+SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1),H = out$model$H)
        }

        # Create the forecasts
        # This gives the forecasts of delta
        forcout = predict(out$model,forcmodel,interval=c('prediction'),
                          level=confidence.level,states=c('trend'))

        # Create empty dataframe to put forecasts in
        forecasts <- matrix(NA,ncol=ncol(data_ldl),nrow=max(n.forc,n.lag)) %>%
          as.data.frame()
        colnames(forecasts) = c('Admissions','Cases')

        # Compute forecasts as per (7) in Andrew's Time Series Models for Epidemics paper
        # Confidence intervals computed as per Harvey, Kattuman and Thamotheram 2021 NIESR paper
        forecasts$Cases[1] = tail(as.vector(data_xts$cCases),(n.lag+1))[1]*exp(forcout$LDLcases[1,1])
        forecasts$Cases[2:n.forc] = tail(as.vector(data_xts$cCases),(n.lag+1))[1]*exp(forcout$LDLcases[2:n.forc,1])*cumprod(1+exp(forcout$LDLcases[1:(n.forc-1),1]))

        forecasts$Admissions[1] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout$LDLhosp[1,1])
        forecasts$Admissions[2:n.forc] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout$LDLhosp[2:n.forc,1])*cumprod(1+exp(forcout$LDLhosp[1:(n.forc-1),1]))

        forecasts$Cases.lwr[1] = tail(as.vector(data_xts$cCases),(n.lag+1))[1]*exp(forcout$LDLcases[1,2])
        forecasts$Cases.lwr[2:n.forc] = tail(as.vector(data_xts$cCases),(n.lag+1))[1]*exp(forcout$LDLcases[2:n.forc,2])*cumprod(1+exp(forcout$LDLcases[1:(n.forc-1),2]))
        forecasts$Admissions.lwr[1] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout$LDLhosp[1,2])
        forecasts$Admissions.lwr[2:n.forc] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout$LDLhosp[2:n.forc,2])*cumprod(1+exp(forcout$LDLhosp[1:(n.forc-1),2]))

        forecasts$Cases.upr[1] = tail(as.vector(data_xts$cCases),(n.lag+1))[1]*exp(forcout$LDLcases[1,3])
        forecasts$Cases.upr[2:n.forc] = tail(as.vector(data_xts$cCases),(n.lag+1))[1]*exp(forcout$LDLcases[2:n.forc,3])*cumprod(1+exp(forcout$LDLcases[1:(n.forc-1),3]))
        forecasts$Admissions.upr[1] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout$LDLhosp[1,3])
        forecasts$Admissions.upr[2:n.forc] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout$LDLhosp[2:n.forc,3])*cumprod(1+exp(forcout$LDLhosp[1:(n.forc-1),3]))

        # Round forecasts to nearest whole number
        forecasts = round(forecasts)

        # Put forecasts into a separate dataframe for admissions and cases
        admissions_forecasts = cbind(forecasts$Admissions,forecasts$Admissions.lwr,forecasts$Admissions.upr)
        colnames(admissions_forecasts) = c('forc','lwr','upr')

        cases_forecasts = cbind(forecasts$Cases,forecasts$Cases.lwr,forecasts$Cases.upr)
        colnames(cases_forecasts) = c('forc','lwr','upr')

        #=============  Re-do with seasonal component=========================

        forcout_sea = predict(out$model,forcmodel,interval=c('prediction'),
                              level=confidence.level,states='all')

        # Create empty dataframe to put forecasts in
        forecasts_sea <- matrix(NA,ncol=ncol(data_ldl),nrow=max(n.forc,n.lag)) %>%
          as.data.frame()
        colnames(forecasts_sea) = c('Admissions','Cases')

        # Compute forecasts as per (7) in Andrew's Time Series Models for Epidemics paper
        # Confidence intervals computed as per Harvey, Kattuman and Thamotheram 2021 NIESR paper
        forecasts_sea$Admissions[1] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout_sea$LDLhosp[1,1])
        forecasts_sea$Admissions[2:n.forc] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout_sea$LDLhosp[2:n.forc,1])*cumprod(1+exp(forcout_sea$LDLhosp[1:(n.forc-1),1]))

        forecasts_sea$Admissions.lwr[1] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout_sea$LDLhosp[1,2])
        forecasts_sea$Admissions.lwr[2:n.forc] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout_sea$LDLhosp[2:n.forc,2])*cumprod(1+exp(forcout_sea$LDLhosp[1:(n.forc-1),2]))

        forecasts_sea$Admissions.upr[1] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout_sea$LDLhosp[1,3])
        forecasts_sea$Admissions.upr[2:n.forc] = tail(as.vector(data_xts$cAdmit),1)*exp(forcout_sea$LDLhosp[2:n.forc,3])*cumprod(1+exp(forcout_sea$LDLhosp[1:(n.forc-1),3]))

        # Round forecasts to nearest whole number
        forecasts_sea = cbind(forecasts_sea$Admissions,forecasts_sea$Admissions.lwr,forecasts_sea$Admissions.upr) %>% round()
        colnames(forecasts_sea) = c('forc','lwr','upr')

        startforc = (data_xts %>% index %>% tail(1))+1

        finds = seq(startforc,length.out = n.forc,by='day')
        fadmits = xts(admissions_forecasts[1:n.forc,],finds)
        sea = xts(forecasts_sea[1:n.forc,],finds)

        if (unity){
          return(list(trend=fadmits[1,], seasonal=sea[1,], n.forc=1))
        } else{
          return(list(trend=fadmits, seasonal=sea, n.forc=n.forc))
        }
        },
        print_estimation_results = function() {
          "Prints a table of estimated parameters in a format ready to paste into
      LaTeX."
          output<-fitmod
          H1 <- output$model$H[1, 1, 1]
          H2 <- output$model$H[2, 2, 1]
          Q_gamma <- output$model$Q[2, 2, 1]
          Q_seasonal <- output$model$Q[3, 3, 1]

          tbl <- data.frame(
            a = format(H1, digits = 3),
            b = format(H2, digits = 3),
            c = format(Q_gamma, digits = 3),
            d = format(Q_seasonal, digits = 3))
          header.names <- c('$\\sigma_\\varepsilon1^2$',
                            '$\\sigma_\\varepsilon2^2$',
                            '$\\sigma_{IRW}^2$',
                            '$\\sigma_{trend1}^2$')

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
      kfs_out <- fitmod
      idx <- index(data_xts)

      if (smoothed) {
        att <- kfs_out$alphahat
      } else {
        att <- kfs_out$att
      }

      filtered_slope <- xts(att[, "slope"], order.by = idx[(n.lag+2):length(idx)])
      filtered.level <- xts(att[, "level"], order.by = idx[(n.lag+2):length(idx)])
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

      kfs_out <- fitmod
      idx <- index(data_xts)

      if (smoothed) {
        att <- kfs_out$alphahat
      } else {
        att <- kfs_out$att
      }

      filtered_slope <- xts(att[, "slope"], order.by = idx[(n.lag+2):length(idx)])
      filtered.level <- xts(att[, "level"], order.by = idx[(n.lag+2):length(idx)])
      g.t <- exp(filtered.level)
      gy.t <- g.t + filtered_slope

      idx.slope <- grep("slope", colnames(kfs_out$att))
      ci <- qnorm((1 - confidence_level) / 2) *
        sqrt(kfs_out$Ptt[idx.slope, idx.slope,]) %o% c(1, -1)
      ci_bounds <- as.vector(gy.t) + ci

      pred <- xts(cbind(gy.t, ci_bounds), order.by = idx[(n.lag+2):length(idx)])
      colnames(pred) <- c("fit","lower","upper")

      return(pred)
    }
  )
)
