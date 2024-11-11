# Created by: Edwin Tang
# Created on: 07/11/2024

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

#' @title Extract output of FilterResults
#
#' @description Accessor method to access the fitted KFS model
#'
#' @param object FilterResults object
#'
#'
#' @export
output<-function(object){
  return(object$output)
}

#' @title Extract number of seasonal components used in KFS
#
#' @description Accessor method to access the number of seasonal component used in KFS object
#'
#' @param object KFS object
#'
#'
#' @export
seasonalComp<-function(object){
  attr(
    object$model$terms, "specials")$SSMseasonal
}

#' @title Extract filtered state estimates used in KFS
#
#' @description Accessor method to access filtered state estimate in KFS object
#'
#' @param object KFS object
#'
#'
#' @export
att<-function(object){
  object$att
}

#' @title Extract error covariance matrix used in KFS
#
#' @description Accessor method to access the non-diffuse part of the error 
#' covariance matrix in KFS object
#'
#' @param object KFS object
#'
#'
#' @export
Ptt<-function(object){
  object$Ptt
}

#' @title Extract error covariance matrix used in KFS
#
#' @description Accessor method to access the non-diffuse part of the error 
#' covariance matrix in KFS object
#'
#' @param object KFS object
#'
#'
#' @export
Ptt<-function(object){
  object$Ptt
}

#' @title Extract matrices used in observation, state and disturbance equation 
#' in KFS object
#
#' @description Accessor method to access the matrices used in observation, 
#' state and disturbance equation in KFS object
#'
#' @param object KFS object
#' @param matrix String, indicating a matrix component of SSModel, e.g. H,T,R,Q
#'
#' @export
matrixKFS<-function(object,matrix){
  object$model[[matrix]]
}