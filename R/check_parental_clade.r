#' Check if parental clade of node is fully contained in ids
#' 
#' Internal
#'
#' @param tre 
#' @param node 
#' @param ids 
#'
#' @return
check_parent_clade <- function(tre, node, ids){
  parent <- node_parent(tre = tre, node = node)
  tips <- phangorn::Descendants(tre, parent, "tips")[[1]]
  
  if(all(tips %in% ids)){
    return(parent)
  }else{
    return(0)
  }
}


#' Recursive check until clade is not in ids
#' 
#' Internal
#'
#' @param tre 
#' @param node 
#' @param ids 
#'
#' @return
#' @export
recursive_check_parent_clade <- function(tre, node, ids){
  parent <- check_parent_clade(tre, node, ids)
  # cat("parent==", parent, "\n")
  if(parent){
    res <- recursive_check_parent_clade(tre, parent, ids)
    # cat("res==", res, "\n")
  }else{
    # cat("node==", node, "\n")
    return(node)
  }
  
  if(res){
    return(res)
  }else{
    return(parent)
  }
}
