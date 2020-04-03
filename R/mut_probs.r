#' Mutation probabilities
#'
#' @param tree
#' @param aln
#' @param anc
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
mut_probs <- function(tree, aln, anc, root_id = NULL, tips = NULL){
  # get root seq
  if(is.null(root_id)){
    root_id <- as.character(castor::find_root(tree))
  }else{
    root_id <- as.character(root_id)
  }
  root_seq <- anc[[as.character(root_id)]]
  colnames(root_seq) <- levels(anc)
  root_seq <- root_seq %>% tibble::as_tibble() %>%
    dplyr::mutate(position = 1:nrow(root_seq))

  # Get tip distribution
  if(is.null(tips)){
    tips <- tree$tip.label
  }else if(is.numeric(tips)){
    tips <- tree$tip.label[tips]
  }else if(!is.character(tips)){
    stop("ERROR: tips must be numeric, character or null.", call. = TRUE)
  }
  Tips <- tips %>%
    purrr::map_dfr(function(label, anc){
      res <- anc[[label]]
      res <- res / rowSums(res)
      res <- res %>% tibble::as_tibble() %>%
        dplyr::mutate(position = 1:nrow(res),
                      name = label)

      return(res)
    }, anc = anc)

  Tip.probs <- Tips %>%
    split(.$position) %>%
    purrr::map_dfr(function(d){
      pos <- d$position %>% unique
      d <- d %>% dplyr::select(-position, -name)

      char_ids <- colnames(d)
      res <- matrix(apply(d, 2, sum) / nrow(d), nrow = 1)
      colnames(res) <- char_ids
      res <- res %>% tibble::as_tibble()
      res$position <- pos

      return(res)
    })

  # Mutation probabilities
  Probs <- 1:nrow(Tip.probs) %>%
    purrr::map_dfr(function(pos, tip_probs, root_probs){
      x <- root_probs %>%
        dplyr::filter(position == pos) %>%
        dplyr::select(-position) %>%
        unlist
      y <- tip_probs %>%
        dplyr::filter(position == pos) %>%
        dplyr::select(-position) %>%
        unlist

      probs <- outer(x,y)
      char_ids <- colnames(probs)
      diag(probs) <- 0
      probs <- colSums(probs)
      total_prob <- sum(probs)
      probs <- matrix(probs, nrow = 1)
      colnames(probs) <- char_ids

      probs <- probs %>% tibble::as_tibble() %>%
        dplyr::mutate(prob = total_prob,
                      position = pos)

      return(probs)
    }, tip_probs = Tip.probs, root_probs = root_seq)

  return(Probs)
}
