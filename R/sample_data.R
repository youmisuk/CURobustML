#' @title Sample two-level data for analysis
#'
#' @description Sample two-level data where individuals (e.g., students) are nested with clusters (e.g., schools)
#'
#' @format Sample data contains:
#'  \describe{
#'    \item{id}{cluster (e.g., school) identifier}
#'    \item{Y}{continuous outcome}
#'    \item{Z}{binary treatment (1=treated, 0=control)}
#'    \item{X1}{first individual-level confounder}
#'    \item{X2}{second individual-level confounder}
#'    \item{X3}{third individual-level confounder}
#'    \item{W1}{first cluster-level confounder}
#'    \item{W2}{second cluster-level confounder, considered as omitted/measured}
#'    \item{lps}{true propensity score logit}
#'    \item{ps}{true propensity score}
#'    \item{Y1}{true potential treatment outcome}
#'    \item{Y0}{true potential control outcome}
#'  }
"twolevel_data"

#' @title Sample cross-classified data for analysis
#'
#' @description Sample cross-classified data where individuals (e.g., students) belong to two clusters simultaneously (e.g., schools and neighborhoods)
#'
#' @format Sample data contains:
#'  \describe{
#'    \item{f1id}{factor-1 cluster (e.g., school) identifier}
#'    \item{f1id}{factor-2 cluster (e.g., neighborhood) identifier}
#'    \item{f12id}{combined factor-12 cluster (e.g., interaction of school and neighborhood) identifier}
#'    \item{Y}{continuous outcome}
#'    \item{Z}{binary treatment (1=treated, 0=control)}
#'    \item{X1}{first individual-level confounder}
#'    \item{X2}{second individual-level confounder}
#'    \item{X3}{third individual-level confounder}
#'    \item{W1}{first factor-1 confounder (e.g., school climate)}
#'    \item{W2}{second factor-1 confounder, considered as omitted/measured}
#'    \item{Q1}{first factor-2 confounder (e.g., perceived neighborhood safety)}
#'    \item{Q2}{second factor-2 confounder, considered as omitted/measured}
#'    \item{lps}{true propensity score logit}
#'    \item{ps}{true propensity score}
#'    \item{Y1}{true potential treatment outcome}
#'    \item{Y0}{true potential control outcome}
#'  }
"crossclassified_data"
