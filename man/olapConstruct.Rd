\name{olapConstruct}
\alias{olapConstruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Create Data for Discovery-driven Exploration of OLAP Data Cubes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is an implementation in R to realize the methodology covered in the paper "Discovery-driven Exploration of OLAP Data Cubes" by Sunita Sarawagi, Rakesh Agarawal and Nimord Megiddo. This methodology is to identify anomalies in multi-dimensional OLAP data cubes. The author does not own the methodology proposed in this paper. Please contact the author at g.t.tongchuan@gmail.com for any potential violation or removal of the package.
}
\usage{
olapConstruct(data, measure, dimensions, user_agg_function, exception_tau=2.5, one_side_trim, output_dir, output_csv=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  a data frame containing OLAP data cube structure with multiple columns representing dimensions and one column representing measure.
  }
  \item{measure}{
  a character string representing the name of the column in data which should be treated as measure.
  }
  \item{dimensions}{
  a character vector representing the names of the columns in data which should be treated as dimensions.
  }
  \item{user_agg_function}{
  a aggregation function that user can choose from either "mean" or "sum".
  }
  \item{exception_tau}{
  a numeric number representing the threshold for the standard residual to determine exception, with default value as 2.5 corresponding to a probability of 99 pct in the normal distribution
  }
  \item{one_side_trim}{
  a numeric value representing the fraction (0 to 0.5) of observations to be trimmed from each end in taking averages.
  }
  \item{output_dir}{
  a string representing the output directory.
  }
  \item{output_csv}{
  a boolen representing whether to output the result as csv, with default value as FALSE. Note that the result in .rds format will always be output to the output directory.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The result is a dataframe in .rds format with rows representing all cells in each groupby, containing 'SelfExp', 'InExp' and 'PathExp' to aid the user in identifying anomalies. Values calculated in the intermediate steps to derive the final results are also saved in extra columns.
}
\references{
%% ~put references to the literature/web site here ~
Sunita Sarawagi, Rakesh Agrawal, and Nimrod Megiddo. Discovery-driven exploration of OLAP data cubes. Research Report RJ 10102 (91918), IBM Almaden Research Center, San Jose, CA 95120, January 1998. Available from http://www.almaden.ibm.com/cs/quest.
}
\author{
%%  ~~who you are~~
Tongchuan Yu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(olapConstruct)
recoveryData_dimensions = names(recoveryData)[-ncol(recoveryData)]
current_wd = getwd()

olapConstruct (data=recoveryData,
               measure="AvgOfDEF_PRICE",
               dimensions=recoveryData_dimensions,
               user_agg_function="mean",
               exception_tau, one_side_trim=0,
               output_dir=current_wd, output_csv=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
