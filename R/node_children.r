
#' Node childre from phylo object
#'
#' @param tre A phylo object
#' @param node Numeric ID or tip label
#'
#' @return Numeric vector of IDs
#' @export
node_children <- function(tre, node){
  if(class(tre) != "phylo")
    stop("ERROR: tre must be of class phylo", call. = TRUE)
  if(length(node) != 1)
    stop("ERROR: only one node can be passed")
  
  if(is.character(node)){
    node <- node_label_to_id(tre = tre, node = node)
  }
  
  if(is.numeric(node)){
    edge <- tre$edge[which(tre$edge[,1] == node), 2]
  }else{
    stop("ERROR: node must be numeric or character", call. = TRUE)
  }
  
  return(edge)
}
