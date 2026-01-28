#' Extract BIC for list of MixtClust results
#'
#' @description After obtaining a list of results for various constraints and
#' number of clusters, extract a data frame contaning the BIC for each model.
#'
#' @param res A list of results given as the output to \code{MixtClust()}
#'
#' @return A data frame with the columns "constraint", "nclusters", and "BIC"
#'
#' @examples
#' set.seed(20241029)
#' d <- subset(iris, select = -Species)
#' res <- MixtClust(d, nclusters = 2:3, sigma.constr = "all")
#' BICs <- extract_bic(res)
#'
#' @author Josh Berlinski
#'
#' @export
extract_bic <- function(res) {
  # number of clusters is nested within constraint
  df <- Map(
    function(constr_res) {
      df_inner <- Map(
        function(nc_res) {
          if (!is(nc_res, "MixtClust"))
            stop("Results must be output from `MixtClust()`")

          data.frame(
            constraint = nc_res$sigma.constr,
            nclusters = nc_res$nclusters,
            BIC = nc_res$bic
          )
        },
        constr_res
      )
      Reduce(rbind, df_inner)
    },
    res
  )

  Reduce(rbind, df)
}

