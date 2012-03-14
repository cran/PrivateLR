################################################################################
#
# File:         dplr.r
# RCS:          $Header: $
# Description:  Differentially private logistic regression.
#               Implements
#               Kamalika Chaudhuri, Claire Monteleoni, Anand Sarwate
#               Privacy-preserving Empirical Risk Minimization.
# Author:       Staal Vinterbo
# Created:      Sat Jan 22 22:48:58 2011
# Modified:     Thu Mar 22 10:22:19 2012 (Staal Vinterbo) staal@ball
# Language:     ESS[S]
# Package:      N/A
# Status:       Experimental, Alpha
#
# dplr.r is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# dplr.r is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with dplr.r; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# (c) Copyright 2011-2012, Staal Vinterbo, all rights reserved.
#
################################################################################
#
# Revisions:
#
# Wed Mar 14 14:50:04 2012 (Staal Vinterbo) staal@dink
#  Added s3 interfaces, initialized parameter vector to 0 for optim,
#  added summary
# Fri Dec  9 15:30:26 2011 (Staal Vinterbo) staal@mats
#  Changed optimizer from Nelder Mead to BFGS
# Sun Jan 23 11:25:30 2011 (Staal Vinterbo) staal@ball
#  Changed implementation from 2008 version to 2010 version.
################################################################################

# these two should probably be moved to 'options'
dplr.verbose = 0 # do we want information printed?
dplr.immediate.warning = TRUE # do we want to see warnings immediately?

# Inputs:
#  y a 0,1 vector of labels
#  x a model matrix of covariates with column names
#    each row in x must have norm <= 1 for privacy guarantees to hold
#  lambda is a regularization parameter
#  eps is the differential privacy parameter. If set to Inf, do non-private
#  logistic regression.
# Returns:
#  the output of stats::optim, and additional parameters.
# NOTE: *only* res$par are differentially private.
dplr.fit = function(y, x, lambda=0.01, eps = 1, dplr.verbose=0) {
  
  if (dplr.verbose > 1) cat('dplr: lambda =', lambda, '\n')
  if (length(unique(y)) < 2) {
    warning('dplr: unique outcome!', immediate. = dplr.immediate.warning)
    res = list(converge=1, par=c(), pred = function(x) 0.5)
    return(res)
  }
  
  y[y == 0] = -1 # recode 0,1 labels to -1, 1
  y[y != -1] = 1
  
  n = length(y)
  d = ncol(x) # we count the intercept, recall x is a model matrix
  if (dplr.verbose > 1) cat('dplr: d =', d, '\n')
  
  b = rep(0, d)
  delta = lambda
  eps2 = NA
  if (eps < Inf) {
    if (dplr.verbose > 1) cat('dplr: eps =', eps, '\n')
    # Chaudhuri-Monteleoni-Sarwate function (2010 version)
    c = 1/4
    q = c/(n*lambda)
    eps2 = eps - log(1 + 2*q + q^2)
    if (eps2 > 0) {
      delta = 0
    } else {
      delta = c/(n*(exp(eps/4) - 1)) - lambda
      eps2 = eps/2
    }
    norm = rgamma(1, shape=d, scale=2/eps2)
    b = runif(d, -1, 1)
    b = b * norm/sqrt(b %*% b)    
    if (dplr.verbose > 1) cat('dplr: norm =', norm, '\n')
  }

  
  target = function(w) {
    xy = y * (x %*% w)
    0.5 * delta * w %*% w + (b %*% w)/n + mean(1 + exp(-xy))
  }

  yx = y*x
  targetgr = function(w) {
    v = as.vector(exp(-yx %*% w))
    delta*w + b/n - apply(v * yx, 2, mean)
  }
    
  #res = optim(rep(1,d), target,  control=list(maxit=10000))
  res = optim(rep(0,d), target, targetgr, method='BFGS')

  names(res$par) = colnames(x)
  class(res) = 'dplr'
  
  # all the other things are not really needed but might be
  # nice to have for posterity, and are used in summary.dplr
  
  pred.p = 1/(1 + exp(-(x %*% res$par)))
  
  CI = function(y, p, n1 = sum(y == 1))
    c(CI=(mean(rank(p)[y == 1]) - (n1+1)/2)/(length(y)-n1))

  res$CIndex = CI(y, pred.p)
  if (dplr.verbose > 1) cat('dplr: CI =', res$CIndex, '\n')  
  res$eps = eps
  res$lambda = lambda
  res$bnorm = norm
  res$delta = delta
  res$eps2 = eps2
  res$n = n
  res$d = d
  res$coefficients = res$par
  
  return(res)
}

dplr = function(object, ...) UseMethod("dplr")

# main interface that all others are using
dplr.formula = function(object, data, lambda=0.01, eps=1, verbose=0, ...) {

  mf = model.frame(object, data)
  mm = model.matrix(terms(mf), mf)
  resp = model.response(mf)

  if(is.factor(resp)){ # last factor level against rest
    #resp = (as.numeric(resp) - 1) == 0
    resp = (resp == tail(levels(resp), 1))
  }
  resp = as.numeric(resp) # in case of bool
  res = dplr.fit(resp, mm, lambda=lambda, eps=eps, dplr.verbose=verbose)

  res$terms = terms(mf)
  res$pred = function(newdata) # nicety 
    predict(res, newdata)
  
  return(res)
}

predict.dplr = function(object, data,...) {
  data = data.frame(cbind(dplr.class=0, data))
  mm = model.matrix(object$terms, data)
  as.vector(1/(1 + exp(-(mm %*% object$par))))
}

sumstring = function(model) {
  l = list(
    #'Summary of Differentially Private Logistic Regression Model:',
    paste(as.character(model$terms)[c(2,1,3)], collapse=' '),
    '',
    paste('Privacy level:', model$eps),
    paste('Regularization parameter:', model$lambda),
    paste('Convergence of optimization:', model$convergence),
    paste('C-index:', model$CIndex),
    '',
    paste(format(names(model$coefficients)), collapse='\t'),
    paste(format(model$coefficients),collapse='\t'))
  paste(paste(' ', l), collapse='\n')
}

summary.dplr = function(object,...) {
  l = sumstring(object)
  cat(l, '\n')
  invisible(l)
}

print.dplr = function(x,...) {
  summary(x)
  invisible(x)
}

# other interfaces:

# y is a 0,1 vector of outcomes,
# x is a matrix/data frame of parameters.
# assumes that no column in x is called dprl.class.
dplr.numeric = function(object, x, ...) {
  df = data.frame(cbind(dplr.class=object, x))
  return(dplr(dplr.class ~ ., data=df, ...))
}

# true/false outcomes
dplr.logical = function(object, x, ...) dplr.numeric(object, x, ...)

# factor outcomes
dplr.factor = function(object, x, ...) dplr.numeric(object, x, ...)


# data frame with target
dplr.data.frame = function(object, target=ncol(object), ...){
  fmla = as.formula(paste(names(object[target]), '~ .'))
  dplr(fmla, data=object,...)
}

# matrix with target
dplr.matrix = function(object, target=ncol(object), ...)
  dplr(data.frame(object), target=target, ...)
  
  
  
