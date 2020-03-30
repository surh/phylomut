#' All descendants from node in phylo object
#'
#' @param tre Phylo object
#' @param node Node ID or tip label.
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
node_descendants <- function(tre, node){
  if(class(tre) != "phylo")
    stop("ERROR: tre must be of class phylo", call. = TRUE)
  if(length(node) != 1)
    stop("ERROR: only one node can be passed")

  if(is.character(node)){
    node <- node_label_to_id(tre = tre, node = node)
  }

  if(is.numeric(node)){
    descendants <- tre$edge[which(tre$edge[,1] == node), 2]
    res <- descendants
    while(length(descendants) > 0){
      descendants <- descendants %>%
        purrr::map(~node_children(tre = tre, node = .)) %>%
        unlist %>%
        unique

      res <- union(res, descendants)
    }

  }else{
    stop("ERROR: node must be numeric or character", call. = TRUE)
  }

  return(res)

}
