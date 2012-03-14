\name{PrivateLR}
\alias{PrivateLR}
\alias{dplr}
\alias{dplr.formula}
\alias{dplr.numeric}
\alias{dplr.factor}
\alias{dplr.logical}
\alias{dplr.data.frame}
\alias{dplr.matrix}
\alias{predict.dplr}
\alias{summary.dplr}
\alias{print.dplr}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Differentially Private Logistic Regression
}
\description{
  Differential privacy is a relatively recent definition of privacy; it
  provides strong privacy with quantitative guarantees. This package
  implements a differentially private version of logistic
  regression.
}
\usage{
dplr(object, \dots)

\method{dplr}{formula}(object, data, lambda=0.01, eps=1, verbose=0, ...)
\method{dplr}{numeric}(object, x, \dots)
\method{dplr}{logical}(object, x, \dots)
\method{dplr}{factor}(object, x, \dots)
\method{dplr}{data.frame}(object, target=ncol(object),\dots)
\method{dplr}{matrix}(object, target=ncol(object),\dots)
\method{predict}{dplr}(object, data, \dots)
\method{summary}{dplr}(object, \dots)
\method{print}{dplr}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    can be given as an object of \code{formula},
    \code{data.frame}, \code{matrix}, or \code{factor}, \code{logical},
    \code{numeric} vector.

    If a \code{data.frame}, \code{matrix} is given, this object contains
    both the dependent variable indexed by \code{target} as well as the
    independent variables, of which all are used. 
    
    In \code{dplr.formula} \code{object} is
    an object of class \code{formula} or an object that can be coerced
    into one.

    If given as a vector, \code{object} contains the values of the dependent
    variable. The vector \code{object} can be of class \code{numeric}, in which
    case it must only contain values 0 and 1, \code{logical} in which
    case it is coerced into numeric by \code{as.numeric(object)}, or be of
    class \code{factor}, in which case it is coerced into numeric by
    encoding the last factor level as 1 and all the other levels as 0.
  }
  \item{data}{
    a data frame or matrix containing the variables in the model described by
    \code{formula}.
  }
  \item{lambda}{
    the regularization parameter.
  }
  \item{eps}{
    the privacy level. The coefficients of the model are computed by a
    method that guarantees \code{eps}-differential privacy.
  }
  \item{verbose}{
    regulates how much information is printed, 0 nothing, 1 a little, 2 more.
  }
  \item{x}{
    In \code{print} \code{x} is an object of class
    \code{"dplr"}, typically returned by \code{dplr}. Otherwise, the
    parameter \code{x} can either be a 
    numeric matrix containing the covariates or dependent variables (one
    per column) corresponding to the dependent variable \code{object},
    or a data frame containing a mix of numeric and factor columns. Any
    factor is internally recoded as contrasts.
  }
  \item{target}{
    the index of the column in \code{data} that contains the target
    values. Default is the last column of \code{data}.
  }
  \item{\dots}{
    \code{verbose}, \code{lambda}, and \code{eps} parameters.
    Not used in \code{summary}, \code{print}, and \code{predict} functions.
  }
}
\details{
  The function \code{dplr} implements logistic regression using the
  differentially private method of Chaudhuri, Monteleoni, and
  Sarwate. This method is based on adding a random offset to
  the L2-regularized objective function for logistic regression.

  The interface is similar but not identical to that of \code{lm}, with
  the addition of the possibility of supplying a data matrix or
  data.frame together with a target column index (defaults to
  \code{ncol(data)}).

  The returned model has a convenience function \code{model$pred} that
  takes a data matrix or data frame to be classified as input.

  The \code{print} function currently prints the summary.
}
\value{
The \code{dplr} function returns a class \code{"dplr"} list object
comprised of the following elements:
\item{par, coefficients}{the coefficients of the logistic model.}
\item{value, counts, convergence, message}{these are as returned by the
  \code{optim} method.} 
\item{CIndex}{the area under the ROC curve (aka., C-Index) of the model
  on its training data.}
\item{eps}{the supplied privacy level.}
\item{lambda}{the regularization parameter}
\item{bnorm, delta, eps2}{See referenced paper for details.}
\item{n}{the number of data points}
\item{d}{the dimensionality of the data points}
\item{terms}{the list output by \code{terms(mf)} where \code{mf} is the
  model frame constructed for the model.}
\item{pred}{a convenience function that is equivalent to
  \code{function(data) predict(model, data)} where \code{model} is the
  list \code{pred} is an element in.
}
}
\references{
  Chaudhuri K., Monteleoni C., and Sarwate, A.
  Differentially Private Empirical Risk Minimization.
  \emph{JMLR}, 2011, 12, 1069-1109
}
\author{
  Staal A. Vinterbo <sav@ucsd.edu>
}
\section{Warning}{
  The privacy level is only guaranteed for the coefficients of the
  model, not for all the other returned values, and also only in the
  case when the input data (potentially after expansion of factors) is
  of norm <= 1.
}
\note{  
  This implementation was in part supported by NIH NLM grant
  7R01LM007273-07 and NIH Roadmap for Medical Research grant U54
  HL108460.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{glm} and 
  \code{predict}
}
\examples{
data(iris)

# the following two are equivalent
# and predict the last factor level
model <- dplr(iris)
model <- dplr(Species ~ ., iris)

# pick a particular factor level
model <- dplr(I(Species == 'virginica') ~ ., iris, eps=2)

# The following is again equivalent to the two first
# examples. Note that we need to remove 'Species' from the
# covariate matrix/data frame, and
# that the class reported by summary will now
# not be 'Species' but 'dplr.class'.
model <- dplr(iris$Species, iris[,-5])

# two equivalent methods to get at the predicted
# probabilities
p <- model$pred(iris)
p <- predict(model, iris)

# print a summary of the model. Note that
# only the coefficients are guaranteed
# to be generated in an eps-differentially
# private manner.
summary(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{privacy}
