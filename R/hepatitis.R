#' Survival from Severe Alcoholic Hepatitis
#'
#' @description The data frame \code{hepatitis} is obtained by digitizing the published
#' Kaplan-Meier curves in Nguyen-Khac et al (2011). The method of digitizing is described in
#' Guyot et al. (2012).
#' See \code{\link{intELtest}} and \code{\link{ptwiseELtest}} for the application.
#' @format The \code{hepatitis} is a data frame with 174 observations of 3 variables,
#' and has the following columns:
#' \itemize{
#'    \item \code{time} the survival time
#'    \item \code{censor} the censoring indicator
#'    \item \code{group} the grouping variable
#' }
#' @source Nguyen-Khac et al., "Glucocorticoids plus N-Acetylcysteine
#' in Severe Alcoholic Hepatitis," \emph{The New England Journal of Medicine},
#' Vol. 365, No. 19, pp. 1781-1789 (2011).
#' \url{http://www.nejm.org/doi/full/10.1056/NEJMoa1101214#t=article}
#' @references P. Guyot, A. E. Ades, M. J. N. M. Ouwens, and N. J. Welton,
#' "Enhanced secondary analysis of survival data: reconstructing the data
#' from published Kaplan-Meier survival curves," \emph{BMC Medical Research Methodology},
#' 12(1):9. \url{http://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-12-9}
#' 
#' @seealso \code{\link{intELtest}}, \code{\link{ptwiseELtest}}

"hepatitis"