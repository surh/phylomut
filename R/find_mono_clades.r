#' Find monophyletic clades
#'
#' Takes a phylo object and a vector of tip labels. Identifies
#' all the clades in the phylogeny that have all their tips
#' contained in the given list.
#'
#' @param tre A phylo object
#' @param ids A vector of tip labels or node IDs
#' @param return Either "all" nodes or "outgroup" and tips.
#'
#' @return A list of numeric vectors where each vector is the
#' node IDs that correspond to a clade in the tree where all
#' the tips are contained in ids.
#' @export
find_mono_clades <- function(tre, ids, return = "all"){
  if(class(tre) != "phylo")
    stop("ERROR: tre must be of class phylo", call. = TRUE)

  # Convert ids to numeric IDs
  if(is.character(ids)){
    ids <- ids %>%
      purrr::map_int(~node_label_to_id(tre = tre, node = .))
  }

  Res <- list()
  while(length(ids) > 0){
    node <- ids[1]
    # cat(node, "\n")
    group_ii <- length(Res) + 1

    group_root <- recursive_check_parent_clade(tre, node, ids)

    if(return == "all"){
      Res[[group_ii]] <- c(group_root, node_descendants(tre, group_root))
    }else if(return == "outgroup"){
      Res[[group_ii]] <- list(outgroup = node_parent(tre, group_root),
                              tips = phangorn::Descendants(tre, group_root, "tips")[[1]])
    }else{
      stop("ERROR::", call. = TRUE)
    }


    ids <- setdiff(ids, phangorn::Descendants(tre, group_root, "tips")[[1]])
  }

  return(Res)
}
