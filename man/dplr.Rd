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
\alias{print.summary.dplr}
\alias{print.dplr}
\alias{scaled}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Differentially Private Logistic Regression
}
\description{
  PrivateLR implements two randomized algorithms for estimating
  \eqn{L_2}{L2}-regularized logistic regression coefficients that allow
  specifying the maximal effect a single point change in the training data are
  allowed to have. Specifically, the algorithms take as parameter the
  maximum allowed change in log-likelihood of producing particular
  coefficients resulting from any single training data point substitution.
}
\usage{
dplr(object, \dots)

\method{dplr}{formula}(object, data, lambda=NA, eps=1, verbose=0,
     rp.dim = 0, threshold='fixed', do.scale=FALSE, ...)
\method{dplr}{numeric}(object, x, \dots)
\method{dplr}{logical}(object, x, \dots)
\method{dplr}{factor}(object, x, \dots)
\method{dplr}{data.frame}(object, target=ncol(object),\dots)
\method{dplr}{matrix}(object, target=ncol(object),\dots)
\method{predict}{dplr}(object, data, type = "probabilities", \dots)
\method{summary}{dplr}(object, \dots)
\method{print.summary}{dplr}(x, \dots)
\method{print}{dplr}(x, \dots)

scaled(fml, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    can be given as an object of \code{formula},
    \code{data.frame}, \code{matrix}, or \code{factor}, \code{logical},
    \code{numeric} vector.

    If a \code{data.frame}, \code{matrix} is given, this object contains
    both the dependent variable indexed by \code{target} as well as the
    independent variables, of which all are used. If the dependent
    variable is a factor, the first level is encoded as 0 and all others
    as 1.
    
    In \code{dplr.formula} \code{object} is
    an object of class \code{formula} or an object that can be coerced
    into one.

    If given as a vector, \code{object} contains the values of the dependent
    variable. The vector \code{object} can be of class \code{numeric}, in which
    case it must only contain values 0 and 1, \code{logical} in which
    case it is coerced into numeric by \code{as.numeric(object)}, or be of
    class \code{factor}, in which case it is coerced into numeric by
    encoding the first factor level as 0 and all the other levels as 1.
  }
  \item{data}{
    a data frame or matrix containing the variables in the model described by
    \code{formula}.
  }
  \item{lambda}{
    the regularization parameter. If \code{NA} (default), the smallest regularizer
    lambda such that \code{2 * log(1 + 1/(4* n * lambda)) == eps/10} is used. If
    \code{eps} is 0, then \code{lambda} is set to 0.001.
  }
  \item{eps}{
    the privacy level. The coefficients of the model are computed by a
    method that guarantees \code{eps}-differential privacy. If \code{eps}
    is 0, then non-private regularized logistic regression is performed.
  }
  \item{verbose}{
    regulates how much information is printed, 0 nothing, 1 a little, 2 more.
  }
  \item{rp.dim}{
    if \code{rp.dim} is non-zero, random projection is performed on the data
    before estimating the model parameters. If \code{rp.dim} is positive,
    the projection will be onto \code{rp.dim} dimensions.  If \code{rp.dim} is
    negative, \code{rp.dim} is set to 1/2 * (1/2)^(-2) * log(n). 
    If \code{rp.dim} is larger than the dimensions of the data, a warning is
    given and no projection is performed.
  }
  \item{threshold}{
    \code{dplr} can non-privately estimate the optimal probability  threshold
    for classification by one of two methods: \code{'youden'}, or \code{'topleft'}.
    The method 'youden' computes the threshold that maximizes the Youden J,
    while 'topleft' computes the threshold corresponding to the point on
    the ROC curve that is closest to (0,1). Any other value (default)
    will result in a threshold of 0.5.
  }
  \item{do.scale}{
    The privacy guarantees are for data where the covariate vectors lie
    within the unit ball. If \code{do.scale} is \code{TRUE}, input data will
    be scaled such that the covariate vectors all lie within the unit
    ball.
  }
  \item{type}{
    \code{predict} can yield two types of results. If \code{type} is
    \code{"probabilities"}, then probabilities are returned, otherwise
    predictions of class values are returned using the threshold
    given by the \code{p.tr} element of \code{object}.
  }
  \item{x}{
    In the \code{print} and \code{print.summary}, \code{x} is an object of class
    \code{"dplr"} or \code{summary.dplr}, typically returned by
    \code{dplr} or \code{summary}. Otherwise, the 
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
  \item{fml}{
    A formula that describes the dimensions of the data that should
    be scaled into the unit ball.
  }
  \item{\dots}{
    \code{verbose}, \code{lambda}, and \code{eps} parameters.
    Not used in \code{summary}, \code{print}, and \code{predict}
    functions. In addition, a Boolean argument \code{op} can be given to
    \code{dplr} to select between \emph{objective perturbation}
    (\code{op = TRUE}, the default) and \emph{output perturbation}
    (\code{op = FALSE}).
  }
}
\details{
  The function \code{dplr} implements logistic regression using the
  differentially private methods by Chaudhuri, Monteleoni, and
  Sarwate. 

  The interface is similar but not identical to that of \code{lm}, with
  the addition of the possibility of supplying a data matrix or
  data.frame together with a target column index (defaults to
  \code{ncol(data)}).

  The returned model instance has a convenience function
  \code{model$pred} that takes a data matrix or data frame to be
  classified as input. 

  The \code{print} function currently prints the summary.

  The \code{scaled} function scales data such that covariate vectors
  lie within the unit ball. Note that the response variable is
  put as the last column in the data frame \code{data} returned.
  Also, the response column name might have changed, depending on
  the left side of the formula given.

  \subsection{Methods details}{

  A randomized algorithm \eqn{A}, taking a dataset as input, is said to be
  \eqn{\epsilon}-differentially private if it holds that
  \deqn{|\log(P(A(D) \in S)) - \log(P(A(D') \in S))| \leq \epsilon}{
    | log(P(A(D) in S)) - log(P(A(D') in S)) | \le \epsilon} for any
  pair of datasets \eqn{D,D'} that differ in exactly one element, and any
  set \eqn{S}. We now turn to the algorithms implemented by \code{dplr}. 
    
  Let \eqn{\|v\|}{l2(v)} denote the L2 norm of a vector \eqn{v}, and let
  \deqn{J(w, \lambda) = ALL(w) + \lambda/2  \|w\|^2}{J(w, \lambda) = ALL(w) + \lambda/2 * l2(w)^2} where \eqn{ALL(w)} is
  the average logistic loss over the
  training data of size \eqn{n} and dimension \eqn{d} with labels
  \eqn{y} and covariates \eqn{x}. L2-regularized logistic regression
  computes \deqn{w^* = \arg\min_w J(w, \lambda)}{w^* = argmin_w J(w,
  \lambda)} for a given \eqn{\lambda}.

  The function \code{dplr} implements two approaches to
  \eqn{\epsilon}-differential private L2 regularized logistic regression
  (see the \dots argument \code{op} above). 
  The first is \emph{output perturbation}, where we compute
  \deqn{w' = w^* + 2/(n \lambda \epsilon) b,}{w' = w^* + 2/(n * \lambda * \epsilon) * b,}
  where \eqn{b} is a \eqn{d}-dimensional real vector sampled with
  probability proportional to 
  \eqn{\exp(-\|b\|)}{exp(-l2(b))}.

  The second is \emph{objective perturbation}. Let
  \deqn{F(w, \lambda, \epsilon) = J(w, \lambda) + 2/(\epsilon n)
  b^Tw}{F(w, \lambda, \epsilon) = J(w, \lambda) + 2/(\epsilon * n) *  b \%*\% w}
  where \eqn{n} and \eqn{b} are as above. Let \eqn{c = 0.25} and let
  \eqn{z = 2\log(1 + c/(\lambda n))}{z = 2 * log1p(c/(\lambda * n))}, then if 
  \deqn{\epsilon - z > 0,}{\epsilon - z > 0,} we compute 
  \deqn{w' = \arg\min_w F(w, \lambda, \epsilon - z)}{w' = argmin_w F(w, \lambda, \epsilon - z)} 
  otherwise we compute an \emph{adjusted lambda} version
  \deqn{w' = \arg\min_w F(w, c/(n (exp(\epsilon/4) - 1)),
  \epsilon/2).}{w' = argmin_w F(w, c/(n * (exp(\epsilon/4) - 1)), \epsilon/2).}

  The logistic regression model coefficients \eqn{w'} are
  then \eqn{\epsilon}-differentially private.
  }
}
\value{
The \code{dplr} function returns a class \code{"dplr"} list object
comprised of elements including:
\item{par}{the coefficients of the logistic model.}
\item{coefficients}{same as \code{par}}
\item{value, counts, convergence, message}{these are as returned by the
  \code{optim} method.} 
\item{CIndex}{the area under the ROC curve (aka., C-Index) of the model
  on its training data.}
\item{eps}{the supplied privacy level.}
\item{lambda}{the regularization parameter used}
\item{n}{the number of data points}
\item{d}{the dimensionality of the data points}
\item{pred}{a convenience function such that \code{predict(model, data, ...)}
  is equivalent to \code{model$pred(data,...)}.
}
\item{p.tr}{this is the classification probability threshold.}
\item{did.rp}{TRUE if random projection was performed.}
\item{rp.dim}{if random projection was performed this contains the
 number of dimensions projected onto. Only present if
  random projection was performed.
}
 \item{rp.p}{the projection matrix used for random projection. Only present if
  random projection was performed.
}
\item{scaled}{TRUE if data was scaled by providing \code{do.scale = TRUE}.}
\item{status}{a text string indicating the status of the computations.
  \code{'ok'} means all is well, \code{'adjusted lambda'} means that the
  regularizer was too small and had to be adjusted, and
  \code{'unique.outcomes'} means that the response had only one value,
  resulting in fixed coefficients returned.
}
  
The \code{scaled} function returns a list of the following:
\item{data}{the scaled data frame}
\item{scale}{the scaling factor used.}
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
  case when input data points (potentially after expansion of factors) are
  of L2-norm <= 1. In particular using prediction thresholds
  estimated using data (methods \code{'youden'} and \code{'topleft'}),
  as well as built in scaling of data is not guaranteed.
  Both of these are turned off by default.
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
# and predict Species being any 
# but the first factor level.
model <- dplr(iris)
model <- dplr(Species ~ ., iris)

# pick a particular factor level and privacy level 2
model <- dplr(I(Species != 'setosa') ~ ., iris, eps=2)

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
