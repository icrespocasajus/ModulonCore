#' network.TILs Data
#'
#' A dataset containing the transcriptional network inferred using SCENIC pipeline on murine TILs transcriptomics.

#'
#' @format A dataframe encoding the network: \describe{ \item{x}{a dataframe of dimension
#'   \code{12490 x 3} with the following colums:\describe{\item{Source}{Transcription factor} \item{Interaction}{Regulatory interaction} \item{Target}{Target gene regulated by the transcription factor}
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' network.TILs
"network.TILs"

#' network.LCMV Data
#'
#' A dataset containing the transcriptional network inferred using SCENIC pipeline on murine LCMV transcriptomics.

#'
#' @format A dataframe encoding the network: \describe{ \item{x}{a dataframe of dimension
#'   \code{16717 x 3} with the following colums:\describe{\item{Source}{Transcription factor} \item{Interaction}{Regulatory interaction} \item{Target}{Target gene regulated by the transcription factor}
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' network.LCMV
"network.LCMV"




#' modulons.TILs Data
#'
#' A dataset containing the modulons inferred using modulon analysis on murine TILs transcriptomics.

#'
#' @format A list of 7 elements/modulons with the modulon constituent elements
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' modulons.TILs
"modulons.TILs"

#' modulons.LCMV Data
#'
#' A dataset containing the modulons inferred using modulon analysis on murine LCMV transcriptomics.

#'
#' @format A list of 12 elements/modulons with the modulon constituent elements
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' modulons.LCMV
"modulons.LCMV"