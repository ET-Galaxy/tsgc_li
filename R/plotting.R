# Created by: Craig Thamotheram
# Created on: 11/02/2022

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


#' @title Plots the forecast of new cases (the difference of the cumulated
#' variable)
#'
#' @description Plots actual values of the difference in the cumulated variable,
#' the forecasts of the cumulated variable (both including and excluding the
#' seasonal component, where a seasonal is specified) and forecast intervals
#' around the forecasts. The forecast intervals are based on the prediction
#' intervals for \eqn{\ln(g_t)}.
#'
#' @param res Results object estimated using the \code{estimate()} method.
#' @param Y (Only required for vanilla growth curve models) Cumulated variable.
#' @param n.ahead Number of forecasts (i.e. number of periods ahead to forecast
#' from end of estimation window).
#' @param confidence.level Width of prediction interval for \eqn{\ln g_t} to
#' use in forecasts of \eqn{y_t = \Delta Y_t}. Default is 0.68, which is
#' approximately one standard deviation for a Normal distribution.
#' @param date_format Date format. Default is \code{'\%Y-\%m-\%d'}.
#' @param title Title for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param plt.start.date First date of actual data (from estimation sample) to
#' plot on graph.\code{NULL} (i.e. plots all data in estimation window) by
#' default.
#' @param series.name Name of the variable you are forecasting for the purposes
#' of a $y$-axis label. E.g. if \code{series.name = "Cases"} the \eqn{y}-axis
#' will show "New Cases".
#'
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel scale_x_date
#' @importFrom ggplot2 geom_ribbon scale_size_manual margin
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom zoo coredata
#' @importFrom utils head tail
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot forecast of new cases 7 days ahead
#' plot_new_cases(res, Y = gauteng[idx.est], n.ahead = 7,
#' confidence.level = 0.68, date_format = "%Y-%m-%d",
#' title = "Forecast new cases", plt.start.date = as.Date("2020-07-13"))
#'
#'
#'
#' @export
plot_new_cases <- function(
  res, n.ahead, Y=NULL, confidence.level = 0.68, date_format = "%Y-%m-%d",
  title=NULL, plt.start.date=NULL, series.name="target variable")
{
  if (class(res)=="FilterResults") {
  if (is.null(Y)) {
      stop("Y is required when res is of class 'FilterResults'")
    } else{
  Date <- Data <- Forecast <- ForecastTrend <- lower <- upper <- NULL
  if (is.null(title)) {title <- ""}
  est.date.index <- res$index %>% as.Date()
  estimation.date.end <- tail(est.date.index, 1)
  y.level.est <- Y[est.date.index]
  if (is.null(plt.start.date)) {plt.start.date <- head(est.date.index, 1)}

  y.hat.diff.final.ci <- res$predict_level(
    y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
    return.diff = TRUE
  )
  y.hat.diff.final <- res$predict_level(
    y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
    sea.on = TRUE, return.diff = TRUE
  )

  tmp.date <- min(estimation.date.end - 4, plt.start.date)
  s <- sprintf("%s/", format(tmp.date, "%Y-%m-%d"))
  d.plot <- cbind(
    diff(y.level.est)[s],
    y.hat.diff.final[, 1],
    y.hat.diff.final.ci[, 1]
  )
  names(d.plot) <- c('Data', 'Forecast', 'ForecastTrend')

  ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                            (as.Date(index(y.hat.diff.final.ci),
                                     format = date_format))))
  colnames(ci) <- c('lower', 'upper', 'date')
  ci[, 'date'] <- as.Date(
    ci[, 'date'], format = date_format, origin = "1970-01-01"
  )

  df_plot <- as.data.frame(d.plot)
  df_plot$Date <- as.Date(rownames(df_plot), format = date_format)

  ggplot2::ggplot(data = df_plot, aes(x = Date)) +
    ggplot2::geom_line(aes(y = Data, color = "Data"), lwd = 0.85) +
    ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
    ggplot2::geom_line(
      aes(y = ForecastTrend, color = "Forecast\nTrend"), lwd = 0.85
    ) +
    ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045")) +
    ggplot2::geom_ribbon(data = ci, aes(x = date, ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
    labs(x = "Date", y = "New Cases", title = title) +
    theme_economist_white(gray_bg = FALSE, base_size = 12) +
    theme(legend.title = element_blank()) +
    theme(
      text = element_text(size = rel(1.1)),
      axis.text = element_text(size = rel(1)),
      axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
      axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
      plot.title = element_text(margin=margin(b=5)),
      plot.caption = element_text(size = rel(1))
    ) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "solid", "solid")) +
    ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
    ggplot2::scale_size_manual(values = c(1, 1, 1))}}

  else if (class(res)=="FilterResultsLI") {
    forecasts<-res$predict_level(n.forc=n.ahead)
    data_xts<-res$data_xts
    n.lag<-res$n.lag
    out<-res$fitmod

    if (is.null(plt.start.date)){plt.start.date <- head(index(data_xts), 1)}
    # add forecasts to plotting dataframe
    fadmits<-forecasts$trend
    sea<-forecasts$seasonal
    sea <- sea[,1]
    fadmits$zero=NA

    # Create smoothed admissions
    lcadmit = lag(as.vector(data_xts$cAdmit)) %>% na.omit()
    smldlh = predict(out$model,states='trend')$LDLhosp %>% exp %>% as.vector
    smadmit = smldlh*lcadmit[(n.lag+1):length(lcadmit)]
    smAdmit = smadmit %>% xts(index(data_xts[(n.lag+1):(length(lcadmit)),])+1)

    #Plot forecast graph
    df_plot<-rbind(data_xts$newAdmit,fadmits$zero)
    df_plot$Smooth<-smAdmit
    df_plot$Forecast<-sea
    df_plot$ForecastTrend<-fadmits$forc
    df_plot<-df_plot%>% subset(index(.) > plt.start.date)
    df_plot=fortify.zoo(df_plot)

    ci<-fortify.zoo(fadmits)

    p2<-ggplot2::ggplot(data = df_plot, aes(x = Index)) +
      ggplot2::geom_line(aes(y = newAdmit, color = "Data"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = Smooth, color = "Smoothed\ndata"),lwd=0.85)+
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
      ggplot2::geom_line(
        aes(y = ForecastTrend, color = "Forecast\nTrend"), lwd = 0.85
      ) +
      ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045","red")) +
      ggplot2::geom_ribbon(data = ci, aes(x = Index, ymin = lwr, ymax = upr),
                           linetype = 0, linewidth = 0, fill = "#AA2045", alpha = 0.1) +
      labs(x = "Date", y = paste("New",series.name), title = title) +
      theme_economist_white(gray_bg = FALSE, base_size = 12) +
      theme(legend.title = element_blank()) +
      theme(
        text = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
        axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
        plot.title = element_text(margin=margin(b=5)),
        plot.caption = element_text(size = rel(1))
      ) +
      ggplot2::scale_linetype_manual(
        values = c("solid", "solid", "solid", "solid")) +
      ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
      ggplot2::scale_size_manual(values = c(1, 1, 1,1))
    return(p2)
  }
}


#' @title Plots forecast and realised values of the log cumulative growth rate
#'
#' @description Plots actual and filtered values of the log cumulative growth
#' rate (\eqn{\ln(g_t)}) in the estimation sample and the forecast and realised
#' log cumulative growth rate out of the estimation sample.
#'
#' @param res Results object estimated using the \code{estimate()} method.
#' @param Y (Only required for leading indicator model.) An xts object
#' containing the cumulative dataset with two columns: the leading
#' indicator and the target variable, including out-of-sample data.
#' @param y.eval (Only required for vanilla growth curve model.) The out-of-sample realisation of the log growth rate of the
#' cumulated variable (i.e. the actual values to which the forecasts should
#' be compared).
#' @param n.ahead The number of time periods ahead from the end of the sample
#' to be forecast. The default is 14.
#' @param plt.start.date Plot start date. Default is \code{NULL} which is the
#' start of the estimation sample.
#' @param title Plot title. Enter as text string.
#' @param caption Plot caption. Enter as text string.
#'
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel margin
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom utils tail
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' idx.eval <- (zoo::index(gauteng) >= as.Date("2020-07-20")) &
#'      zoo::index(gauteng) <= as.Date("2020-07-27")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot forecast and realised log growth rate of cumulative cases
#' plot_forecast(res, y.eval = df2ldl(gauteng[idx.eval]), n.ahead = 7,
#'   title = "Forecast ln(g)", plt.start.date = as.Date("2020-07-13"))
#'
#' @export
plot_forecast <- function(res, Y=NULL, y.eval=NULL, n.ahead = 14,
                          plt.start.date=NULL, title="", caption = "") {
  if (class(res)=="FilterResults"){
    if (is.null(y.eval)){
      stop("y.eval is required when res is of class 'FilterResults'")
    }
    Date <- NULL
    model <- res$output$model
    est.date.index <- res$index

    y <- xts::xts(res$output$model$y %>% as.numeric(), order.by = est.date.index)
    p <- attr(res$output$model, 'p')

    y.hat.all <- res$predict_all(n.ahead, return.all = TRUE)
    y.pred <-  subset(y.hat.all$y.hat,index(y.hat.all$y.hat) > tail(res$index,1))
    filtered.level <- y.hat.all$level

    if (p == 1) {
      EstimationSample <- FilteredLevel <- Forecast <- RealisedData <- NULL
      d <- cbind(y, filtered.level, y.pred,
                 y.eval[index(y.eval)>tail(index(est.date.index),1),])
      if (!is.null(plt.start.date)) { d <- d[index(d) > plt.start.date] }
      d <- d[index(d) <= tail(index(y.pred),1)]
      names(d) <- c(
        'EstimationSample', 'FilteredLevel', 'Forecast', 'RealisedData'
      )

      df_plot <- as.data.frame(d)
      df_plot$Date <- as.Date(rownames(df_plot))

      p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
        ggplot2::geom_line(aes(
          y = EstimationSample, color = "Estimation\nSample"), lwd = 0.85) +
        ggplot2::geom_line(aes(y = FilteredLevel, color = "Filtered\nLevel"),
                           lwd = 0.85) +
        ggplot2::geom_line(aes(y = Forecast, color = "Forecast"), lwd = 0.85) +
        ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData"),
                           lwd = 0.85) +
        ggplot2::scale_color_manual(values = c(1, 2, 3, 'grey')) +
        scale_linetype_manual(
          values = c("solid", "solid", "solid", "dashed")) +
        scale_x_date(labels = scales::date_format("%d %b %y")) +
        labs(x = "Date", y = "Log Growth Rate", caption = caption,
             title = title
        ) +
        theme_economist_white(gray_bg = FALSE) +
        scale_fill_economist() +
        theme(legend.title = element_blank()) +
        theme(
          text = element_text(size = rel(1)),
          axis.text = element_text(size = rel(1)),
          axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
          axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
          plot.title = element_text(margin=margin(b=5)),
          plot.caption = element_text(size = rel(1)),
        )
    } else if (p == 2) {
      g_1 <- g_2 <- delta <- Forecast <- RealisedData <- NULL
      d <- cbind(y, filtered.level, y.pred[,2],
                 y.eval[index(y.eval)>tail(index(est.date.index),1),2])
      d <- d[index(d) <= tail(index(y.pred),1)]
      names(d) <- c('g_1', 'g_2', 'delta', 'Forecast', 'RealisedData')

      df_plot <- as.data.frame(d)
      df_plot$Date <- as.Date(rownames(df_plot))

      p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
        ggplot2::geom_line(aes(y = g_1, color = "g_1")) +
        ggplot2::geom_line(aes(y = g_2, color = "g_2")) +
        ggplot2::geom_line(aes(y = g_2, color = "delta")) +
        ggplot2::geom_line(aes(y = Forecast, color = "Forecast")) +
        ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData")) +
        ggplot2::scale_color_manual(
          values = c(1, 2, 3, 4, 'grey')) +
        ggplot2::scale_linetype_manual(
          values = c("solid", "solid", "solid", "solid", "dashed")
        ) +
        ggplot2::scale_x_date(labels = scales::date_format("%d %b %y")) +
        labs(x = "Date", y = "Log Growth Rate", caption = caption,
             title = title
        ) +
        theme_economist_white(gray_bg = FALSE) +
        scale_fill_economist() +
        theme(legend.title = element_blank()) +
        theme(
          text = element_text(size = rel(1.)),
          axis.text = element_text(size = rel(1)),
          axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
          axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
          plot.title = element_text(margin=margin(b=5)),
          plot.caption = element_text(size = rel(1))
        )
    } else { stop('NotImplemented Error') }

    return(p1)
  }
  else{
    if (is.null(Y)){
      stop("Y is required when res is of class 'FilterResultsLI'")
    }
    out<-res$fitmod
    data_xts<-res$data_xts
    n.lag<-res$n.lag

    old=res$data_xts[,"LDLhosp"]
    old=old[index(old)>head(index(old),1)+n.lag]

    eng_full<-add_daily_ldl(Y)
    eng_full<-eng_full[index(eng_full)>tail(index(old),1),"LDLhosp"]
    actual=eng_full[1:n.lag]

    lcadmit = lag(as.vector(data_xts$cAdmit)) %>% na.omit()
    smldlh = predict(out$model,states='trend')$LDLhosp
    filtered=xts(smldlh,(head(index(data_xts),1)+n.lag)+(1:length(smldlh)))

    start_date<-index(data_xts)[1]
    end_date<-tail(index(data_xts),1)
    data_ldl = data_xts[,c("LDLcases","LDLhosp")] %>% na.omit

    data_ldl$LDLcases = lag(as.vector(data_ldl$LDLcases),n.lag)

    data_ldl <- na.omit(data_ldl)

    data_mat = as.matrix(data_ldl)
    # Create forecast data (using fact have some "future" case observations)
    forcdata <- matrix(NA,ncol=2,nrow=max(n.ahead,n.lag))
    colnames(forcdata) = colnames(data_mat)
    forcdata[1:n.lag,1] = as.vector(tail(data_xts,n.lag)$LDLcases)

    # Extract estimate of Q from earlier model
    Qf = out$model$Q[,,1]
    if (is.na(res$sea.period)) {
      forcmodel = SSModel(forcdata ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1),
                          H = out$model$H)
    }
    else{
      forcmodel = SSModel(forcdata ~ SSMtrend(degree = 2, Q = matrix(c(0,0,0,Qf[2,2]),2,2),type = 'common')+SSMseasonal(res$sea.period,Q = matrix(c(Qf[4,4],0,0,Qf[5,5]),2,2), sea.type='dummy', type='distinct')+SSMtrend(degree = 1, Q = matrix(Qf[3,3]),index=1),H = out$model$H)
    }
    forcout = predict(out$model,forcmodel,interval=c('prediction'),
                      level=0.68,states=c('trend'))

    trend=xts(forcout$LDLhosp[,"fit"],tail(index(old),1)+(1:n.ahead))

    d.plot<-cbind(old,filtered,trend,actual)
    colnames(d.plot)<-c('EstimationSample', 'FilteredLevel', 'Forecast', 'RealisedData')
    if (!is.null(plt.start.date)){
      d.plot<-d.plot[index(d.plot)>=plt.start.date]
      df_plot <- as.data.frame(d.plot)
    } else{
      df_plot <- as.data.frame(d.plot)
    }
    df_plot$Date <- as.Date(rownames(df_plot))

    ggplot2::ggplot(data = df_plot, aes(x = Date)) +
      ggplot2::geom_line(aes(
        y = EstimationSample, color = "Estimation\nSample"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = FilteredLevel, color = "Filtered\nLevel"),
                         lwd = 0.85) +
      ggplot2::geom_line(aes(y = Forecast, color = "Forecast\nTrend"), lwd = 0.85) +
      ggplot2::geom_line(aes(y = RealisedData, color = "Realised\nData"),
                         lwd = 0.85) +
      ggplot2::scale_color_manual(values = c(1, 2, 3, 'grey')) +
      scale_linetype_manual(
        values = c("solid", "solid", "solid", "dashed")) +
      scale_x_date(labels = scales::date_format("%d %b %y")) +
      labs(x = "Date", y = "Log Growth Rate", caption = caption,
           title = title
      ) +
      theme_economist_white(gray_bg = FALSE) +
      scale_fill_economist() +
      theme(legend.title = element_blank()) +
      theme(
        text = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1),margin = margin(r=10)),
        axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
        plot.title = element_text(margin=margin(b=5)),
        plot.caption = element_text(size = rel(1)),
      )
  }
}

#' @title Plots the growth rates and slope of the log cumulative growth rate
#'
#' @description Plots the smoothed/filtered growth rate of the difference in
#' the cumulated variable (\eqn{g_y}), the smoothed/filtered growth rate of the
#' the cumulated variable (\eqn{g}), and the smoothed/filtered slope of
#' \eqn{\ln(g)}, \eqn{\gamma}.
#' Following Harvey and Kattuman (2021), we compute \eqn{g_{y,t}} as
#' \deqn{g_{y,t} = \exp(\delta_t) + \gamma_t.}
#'
#' @param res Results object estimated using the \code{estimate()} method.
#' @param plt.start.date Plot start date. Default is \code{NULL} which is the
#' start of the estimation sample.
#' @param smoothed Logical value indicating whether to used the smoothed
#' estimates of \eqn{\delta} and \eqn{\gamma}. Default is \code{FALSE}, in
#' which case the filtered estimates are returned.
#' @param title Title for plot. Enter as text string. \code{NULL} (i.e. no
#' title) by default.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot filtered gy, g and gamma
#' plot_gy_components(res, plt.start.date = as.Date("2020-07-06"))
#'
#'@importFrom ggplot2 ggplot geom_line labs scale_x_date scale_y_continuous
#'@importFrom ggplot2 waiver
#'@importFrom ggplot2 theme margin scale_color_manual
#'@importFrom dplyr filter
#'@importFrom tidyr pivot_longer
#'@importFrom ggthemes theme_economist_white
#'@importFrom magrittr %>%
#'
#' @export
plot_gy_components <- function(res, plt.start.date = NULL,
                               smoothed = FALSE, title = NULL){
  Date <- Value <- Variable <- NULL
  # Determine plot start date
  if(is.null(plt.start.date)) {
    if (class(res)=="FilterResults"){
      plt.start.date <- res$index[1]
    } else {
      plt.start.date <- index(res$data_xts)[1]
    }
    }

  # Get gy.t, g.t and gamma
  gy.components <- res$get_growth_y(return.components = TRUE, smoothed =
                                      smoothed)
  gy.t <- gy.components[[1]]
  g.t <- gy.components[[2]]
  gamma.t <- gy.components[[3]]

  d <- cbind(gy.t,g.t,gamma.t)
  names(d) <- c('gy.t','g.t','gamma.t')

  df_plot <- as.data.frame(d)
  df_plot$Date <- as.Date(rownames(df_plot))

  df_long <- df_plot %>%
    dplyr::filter(Date >= plt.start.date) %>%
    pivot_longer(cols = c(gy.t, g.t, gamma.t), names_to = "Variable",
                 values_to = "Value")

  p1 <- ggplot(df_long, aes(x = Date, y = Value, color = Variable)) +
    geom_line(lwd=0.85) +
    ggplot2::facet_wrap(~ factor(
      Variable, c("gy.t", "g.t", "gamma.t")), ncol = 1, scales = "free_y") +
    labs(title = title, y=ggplot2::element_blank()) +
    scale_color_manual(values = c("#AA2045","darkgrey","black")) +
    scale_x_date(labels = scales::date_format("%d %b %y")) +
    scale_y_continuous(breaks = waiver(), n.breaks = 4) +
    theme_economist_white(gray_bg = FALSE, base_size = 14) +
    theme(text = element_text(size= rel(1), margin=ggplot2::margin(b=5)),
          axis.title.x = element_text(size = rel(1),margin = margin(t=10)),
          legend.position = "none")

  return(p1)
}

#' @title Plots the growth rates and slope of the log cumulative growth rate
#'
#' @description Plots the smoothed/filtered growth rate of the difference in the
#' cumulated variable (\eqn{g_y}) and the associated confidence intervals.
#'
#' @param res Results object estimated using the \code{estimate()} method.
#' @param plt.start.date Plot start date. Default is \code{NULL} which is the
#' start of the estimation sample.
#' @param smoothed Logical value indicating whether to used the smoothed
#' estimates of \eqn{\delta} and \eqn{\gamma}. Default is \code{FALSE}, in
#' which case the filtered estimates are returned.
#' @param title Title for plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param series.name The name of the series the growth rate is being computed
#' for. E.g. \code{'New cases'}.
#' @param pad.right Numerical value for the amount of time periods of blank
#' space you wish to leave on the right of the graph. Extends the horizontal
#' axis by the given number of time periods.
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot filtered gy, g and gamma
#' plot_gy_ci(res, plt.start.date = as.Date("2020-07-13"))
#'
#' @importFrom ggplot2 ggplot geom_line geom_hline geom_ribbon labs
#' scale_color_manual scale_linetype_manual margin
#' @importFrom ggthemes theme_economist_white
#' @importFrom utils tail
#'
#' @export
plot_gy_ci <- function(res, plt.start.date = NULL, smoothed = FALSE,
                       title = NULL, series.name = NULL, pad.right = NULL){
  Date <- fit <- upper <- lower <- NULL

  # Determine plot start date
  if(is.null(plt.start.date)) {
    if (class(res)=="FilterResults"){
      plt.start.date <- res$index[1]
    } else {
      plt.start.date <- index(res$data_xts)[1]
    }
  }

  # Get confidence intervals to plot
  gy.ci<- res$get_gy_ci(smoothed = smoothed)

  y.lab <- if(is.null(series.name)) { c("Growth rate") } else {
    paste("Growth rate of"," ",series.name,sep="")
  }

  df_plot <- as.data.frame(gy.ci)
  df_plot$Date <- as.Date(rownames(df_plot))

  p1 <- ggplot2::ggplot(df_plot[df_plot$Date>=plt.start.date,], aes(x=Date)) +
    ggplot2::geom_line(aes(y = fit), lwd = 0.85) +
    ggplot2::geom_hline(yintercept=0, linetype="solid",
                        color = "green", linewidth=1)+
    ggplot2::geom_ribbon(aes(ymin = lower, ymax = upper),
                         linetype = 0, linewidth = 0, fill = "#AA2045",
                         alpha = 0.3) +
    ggplot2::scale_color_manual(values = c("black")) +
    geom_hline(
      aes(yintercept = 0.0), linetype = "solid", color = "green", lwd = 1.
    ) +
    labs(title=title, x="Date", y=y.lab) +
    theme_economist_white(gray_bg = FALSE, base_size = 14) +
    theme(
      legend.title = element_blank(),
      text = element_text(size = rel(1.)),
      axis.text = element_text(size = rel(1.)),
      axis.title.y = element_text(
        size = rel(1.),margin = ggplot2::margin(r=10)),
      axis.title.x = element_text(
        size = rel(1.),margin = ggplot2::margin(t=10)),
      plot.caption = element_text(size = rel(1))
    ) +
    theme(panel.grid.major.x = ggplot2::element_line(
      color = "gray50", linewidth = 0.5)) +
    scale_linetype_manual(
      values = c("solid")) +
    scale_x_date(labels = scales::date_format("%d %b %y"))

  if (!is.null(pad.right)) {
    end.date <- tail(index(gy.ci),1)
    p1 <- p1 +
      ggplot2::scale_x_date(
        limits = c(as.Date(plt.start.date), end.date + pad.right))
  }

  return(p1)
}

#' @title Plots the forecast of new cases (the difference of the cumulated
#' variable) over a holdout sample.
#'
#' @description Plots actual values of the difference in the cumulated variable,
#' the forecasts of the cumulated variable (both including and excluding the
#' seasonal component, where a seasonal is specified) and forecast intervals
#' around the forecasts, plus the actual outcomes from the holdout sample. The
#' forecast intervals are based on the prediction intervals for \eqn{\ln(g_t)}.
#' Also reports the mean absolute percentage prediction error over the holdout
#' sample.
#'
#' @param res Results object estimated using the \code{estimate()} method.
#' @param Y (For vanilla growth curve) Values of the cumulated variable to be
#' used in the estimation window.
#' (For leading indicator) An xts object containing the cumulative dataset with
#' two columns: the leading indicator and the target variable, including
#' out-of-sample data.
#' @param Y.eval (Only required for vanilla growth curve estimation.)
#' Values of the cumulated variable to be used in the holdout
#' sample (i.e. to which the forecasts should be compared to).
#' @param confidence.level Width of prediction interval for \eqn{\ln(g_t)} to
#' use in forecasts of \eqn{y_t = \Delta Y_t}. Default is 0.68, which is
#' approximately one standard deviation for a Normal distribution.
#' @param series.name Name of the variable you are forecasting for the purposes
#' of a $y$-axis label. E.g. if \code{series.name = "Cases"} the \eqn{y}-axis
#' will show "New Cases".
#' @param date_format Date format, e.g. \code{'\%Y-\%m-\%d'}, which is the
#' default. Only needed for Dynamic Gompertz model.
#' @param title Title for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no title) by default.
#' @param caption Caption for forecast plot. Enter as text string. \code{NULL}
#' (i.e. no caption) by default.
#'
#' @importFrom xts as.xts
#' @importFrom ggplot2 scale_color_manual scale_linetype_manual aes labs theme
#' @importFrom ggplot2 element_blank element_text rel autoplot scale_x_date
#' @importFrom ggplot2 geom_ribbon scale_size_manual
#' @importFrom ggthemes theme_economist_white scale_fill_economist
#' @importFrom zoo coredata index
#' @importFrom utils tail
#' @importFrom stats na.omit
#'
#' @returns A \code{ggplot2} plot.
#'
#' @examples
#' library(tsgc)
#' data(gauteng,package="tsgc")
#' idx.est <- zoo::index(gauteng) <= as.Date("2020-07-20")
#' idx.eval <- (zoo::index(gauteng) >= as.Date("2020-07-20")) &
#'      zoo::index(gauteng) <= as.Date("2020-07-27")
#'
#' # Specify a model
#' model <- SSModelDynamicGompertz$new(Y = gauteng[idx.est], q = 0.005)
#' # Estimate a specified model
#' res <- model$estimate()
#'
#' # Plot forecasts and outcomes over evaluation period
#' plot_holdout(res = res, Y = gauteng[idx.est], Y.eval = gauteng[idx.eval])
#'
#' @export
plot_holdout <- function(res, Y, Y.eval=NULL, n.ahead=NULL, confidence.level = 0.68,
                         date_format = "%Y-%m-%d", series.name = "target variable",
                         title= NULL, caption = NULL) {
  if (class(res)=="FilterResults"){
    if (is.null(Y.eval)) {
      stop("Y.eval is required when res is of class 'FilterResults'")
    } else{
      Date <- Actual <- Forecast <- ForecastTrend <- lower <- upper <- NULL

      model <- res$output$model
      est.date.index <- res$index

      y.level.est <- Y[est.date.index]

      p <- attr(res$output$model, 'p')
      if(p!=1) { stop('NotImplementedError') }

      n.ahead <- tail(index(Y.eval),1)-tail(index(y.level.est),1)

      y.eval.diff <- diff(Y.eval) %>% na.omit

      est.date.index <- res$index %>% as.Date()
      estimation.date.end <- tail(est.date.index, 1)

      y.hat.diff.final.ci <- res$predict_level(
        y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
        return.diff = TRUE
      )
      y.hat.diff.final <- res$predict_level(
        y.cum = y.level.est, n.ahead = n.ahead, confidence_level = confidence.level,
        sea.on = TRUE,
        return.diff = TRUE
      )

      d <- cbind(
        y.eval.diff[index(y.eval.diff)>estimation.date.end,],
        y.hat.diff.final[, 1],
        y.hat.diff.final.ci[, 1]
      )
      names(d) <- c('Actual', 'Forecast', 'ForecastTrend')

      df_plot <- as.data.frame(d)
      df_plot$Date <- as.Date(rownames(df_plot), format = date_format)

      d.eval <- na.omit(d)
      mape.trend <- 100*(abs(d.eval$Actual - d.eval$`ForecastTrend`)/
                           d.eval$Actual) %>% mean %>% round(2)
      mape.sea <- 100*(abs(d.eval$Actual - d.eval$Forecast)/d.eval$Actual) %>%
        mean %>% round(2)

      ci <- as.data.frame(cbind(zoo::coredata(y.hat.diff.final.ci[, 2:3]),
                                (as.Date(index(y.hat.diff.final.ci),
                                         format = date_format))))
      colnames(ci) <- c('lower', 'upper', 'date')
      ci[, 'date'] <- as.Date(ci[, 'date'], format = date_format,
                              origin = "1970-01-01")

      p1 <- ggplot2::ggplot(data = df_plot, aes(x = Date)) +
        ggplot2::geom_line(aes(y = Actual, color = "Actual"),lwd = 0.85) +
        ggplot2::geom_line(aes(y = Forecast, color = "Forecast"),lwd = 0.85) +
        ggplot2::geom_line(
          aes(y = ForecastTrend, color = "Forecast\nTrend"),lwd = 0.85) +
        ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045")) +
        ggplot2::geom_ribbon(data = ci, aes(x = date, ymin = lower, ymax = upper),
                             linetype = 0, linewidth = 0, fill = "#AA2045",
                             alpha = 0.1) +
        labs(x = "Date", y = paste("New",series.name), title = title,
             subtitle = paste("MAPE: ",mape.sea,"%. Trend MAPE: ",
                              mape.trend,"%.",sep="")) +
        theme_economist_white(gray_bg = FALSE, base_size = 14) +
        theme(legend.title = element_blank()) +
        theme(
          text = element_text(size = rel(1)),
          axis.text = element_text(size = rel(1)),
          axis.title.y = element_text(size = rel(1), margin = margin(r=10)),
          axis.title.x = element_text(size = rel(1), margin = margin(t=10)),
          plot.title = element_text(margin=margin(b=5)),
          plot.subtitle = element_text(
            size = rel(1), hjust=0,  margin = margin(t=3))
        ) +
        scale_linetype_manual(
          values = c("solid", "solid", "solid")) +
        scale_x_date(labels = scales::date_format("%d %b %y")) +
        scale_size_manual(values = c(1, 1.5, 1))
      return(p1)
    }
  }
  else if (class(res)=="FilterResultsLI"){
    if (is.null(n.ahead)) {
      stop("n.ahead is required when res is of class 'FilterResultsLI'")
    } else{
      forecasts<-res$predict_level(n.forc=n.ahead)
      fadmits<-forecasts$trend    #trend
      sea<-forecasts$seasonal     #seasonal

      end_date<-tail(index(res$data_xts),1)
      sea<-sea[,1] #get the forecast column

      future_data<-add_daily_ldl(Y,LeadIndCol = res$LeadIndCol) %>% subset(index(.) > end_date)
      data_validation<-future_data[1:n.ahead, c("cAdmit", "newAdmit")]

      newAdmit_validation<-data_validation[,c("newAdmit")]
      compare<-cbind(newAdmit_validation,fadmits[,1], sea)
      names(compare)<-c("Actual", "ForecastTrend", "Forecast")

      mape.trend <- 100*(abs(compare$Actual - compare$ForecastTrend)/
                           compare$Actual) %>% mean %>% round(4)
      mape.sea <- 100*(abs(compare$Actual - compare$Forecast)/compare$Actual) %>%
        mean %>% round(4)

      ci<-fadmits[,-1]
      colnames(ci) <- c('lower', 'upper')

      df_plot <- as.data.frame(compare)
      df_plot$Date <- as.Date(rownames(df_plot), format="%Y-%m-%d")

      ci_plot <- as.data.frame(ci)
      ci_plot$Date <- as.Date(rownames(ci_plot), format = "%Y-%m-%d")

      p1<-ggplot2::ggplot(data = df_plot, aes(x = Date)) +
        ggplot2::geom_line(aes(y = Actual, color = "Actual"),lwd = 0.85) +
        ggplot2::geom_line(aes(y = Forecast, color = "Forecast"),lwd = 0.85) +
        ggplot2::geom_line(
          aes(y = ForecastTrend, color = "Forecast\nTrend"),lwd = 0.85) +
        ggplot2::scale_color_manual(values = c("black", "grey", "#AA2045")) +
        ggplot2::geom_ribbon(data = ci_plot, aes(x = Date, ymin = lower, ymax = upper),linetype = 0, linewidth = 0, fill = "#AA2045",
                             alpha = 0.1) +
        labs(x = "Date", y = paste("New",series.name), title = title,
             subtitle = paste("MAPE: ",mape.sea,"%. Trend MAPE: ",
                              mape.trend,"%.",sep="")) +
        theme_economist_white(gray_bg = FALSE, base_size = 14) +
        theme(legend.title = element_blank()) +
        theme(
          text = element_text(size = rel(1)),
          axis.text = element_text(size = rel(1)),
          axis.title.y = element_text(size = rel(1), margin = margin(r=10)),
          axis.title.x = element_text(size = rel(1), margin = margin(t=10)),
          plot.title = element_text(margin=margin(b=5)),
          plot.subtitle = element_text(
            size = rel(1), hjust=0,  margin = margin(t=3))
        ) +
        scale_linetype_manual(
          values = c("solid", "solid", "solid")) +
        scale_x_date(labels = scales::date_format("%d %b %y")) +
        scale_size_manual(values = c(1, 1.5, 1))

      return(p1)
    }
  }
}

