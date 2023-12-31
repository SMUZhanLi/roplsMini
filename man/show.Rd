% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opls-methods.R
\name{show,opls-method}
\alias{show,opls-method}
\alias{show.opls}
\title{Show method for 'opls' objects}
\usage{
\S4method{show}{opls}(object)
}
\arguments{
\item{object}{An S4 object of class \code{opls}, created by the \code{opls}
function.}
}
\value{
Invisible.
}
\description{
Displays information about the dataset and the model.
}
\examples{

data(sacurine)
attach(sacurine)
sacurine.plsda <- opls(dataMatrix, sampleMetadata[, "gender"])

show(sacurine.plsda)

detach(sacurine)

}
\author{
Philippe Rinaudo and Etienne Thevenot (CEA)
}
