% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{preserve_cl}
\alias{preserve_cl}
\title{Checks if Ck, Ck' in C(x'(phi))}
\usage{
preserve_cl(cl, cl_phi, k1, k2)
}
\arguments{
\item{cl}{clustering of x}

\item{cl_phi}{clustering of x'(phi)}

\item{k1, k2}{index of clusters involved in the test}
}
\value{
Returns TRUE if Ck, Ck' in C(x'(phi)), and FALSE otherwise
}
\description{
Checks if Ck, Ck' in C(x'(phi))
}
\keyword{internal}
