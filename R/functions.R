#' @title Jaccard Distance
#' @description Calculate Jaccard distance between two vectors
#' #' @param a A vector.
#' @param b A vector.
#' @return Jaccard distance calculated as the ratio between the intersection and the union of the two vectors.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  jaccard(c('A','B','C','D','E'), c('A','B'))
#'  }
#' }
#' @rdname jaccard
#' @export 
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
#' @title Regulon Redundancy
#' Calculate redundancy between two regulons
#' @param a A vector.
#' @param b A vector.
#' @return Redundancy calculated as the ratio between the intersection and the minimum length of the two vectors.
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' redundancy(c('A','B','C','D','E'), c('A','B'))
#'  }
#' }
#' @rdname redundancy
#' @export 
redundancy <- function(a, b) {
  intersection = length(intersect(a, b))
  min = min(length(a),length(b))
  return (intersection/min)
}
#' @title Range [0,1]
#' @description Normalize absolute values of a numeric vector
#' @param x A numeric vector.
#' @return Normalized values (x-min(x))/(max(x)-min(x))
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#' range01(c('-50','-5','5','50','100')
#'  }
#' }
#' @rdname range01
#' @export 
range01 <- function(x){(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))}
