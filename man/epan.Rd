\name{epan}
\alias{epan}
\title{
A C++ function for numeric calculation of float format.
}
\description{
The C++ function epan() does numeric calculation depending on the absolute value of the input. If the absolute value of the input is less than or equal to 1, a numeric output is returned. Otherwise, 0 is returned.
}
\usage{
epan(input)
}
\arguments{
  \item{input}{
  A numeric input of float type.
}
}
\value{
Returns a float output vector as 0.75*(1-input^2) if the absolute value of the input is less than or equal to 1; 0 otherwise.
}
\keyword{internal}