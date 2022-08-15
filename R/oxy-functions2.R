#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param b PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  redundancy3(c('A','B','C','D','E'), c('A','B','C'))
#'  }
#' }
#' @rdname redundancy3
#' @export 
redundancy3 <- function(a, b) {
  intersection = length(intersect(a, b))
  min = min(length(a),length(b))
  return (intersection/min)
}
