#' Node label to numeric id
#' 
#' Internal
#'
#' @param tre 
#' @param node 
#'
#' @return
node_label_to_id <- function(tre, node){
  # if(class(tre) != "phylo")
  #   stop("ERROR: tre must be of class phylo", call. = TRUE)
  # if(length(node) != 1)
  #   stop("ERROR: only one node can be passed")
  
  # Convert node to numeric
  # if(is.character(node)){
  #   node <- which(tre$tip.label == node)
  #   
  #   if(length(node) != 1)
  #     stop("ERROR: node label must correspond to a tip label in tre", call. = TRUE)
  # }else{
  #   stop("ERROR: node must be a character variable", call. = TRUE)
  # }
  node <- which(tre$tip.label == node)
  
  if(length(node) != 1)
    stop("ERROR: node label must correspond to a tip label in tre", call. = TRUE)
  
  
  return(node)
}
