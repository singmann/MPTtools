
parse_model_df <- function(model_df) {
  model_trees <- split(model_df, model_df$Tree)
  model_branches <- lapply(model_trees, function(x) split(x, x$Category))
  model_parsed <- vector("list", length(model_trees))
  for (t in seq_along(model_branches)) {
    for (c in seq_along(model_branches[[t]])) {
      model_parsed[[t]][c] <-  parse(text = paste(model_branches[[t]][[c]][["Equation"]], collapse = " + "))
    }
  }
  names(model_parsed) <- names(model_trees)
  attr(model_parsed, "cat_map") <- lapply(model_branches, names)
  model_parsed
}

make_model_equations <- function(model_df, index) {
  model_trees <- split(model_df, model_df$Tree)
  model_branches <- lapply(model_trees, function(x) split(x, x$Category))
  model_out <- vector("character", length(unique(model_df$Category)))
  prob_out <- vector("character", length(model_df$Category))
  cat <- 1
  for (t in seq_along(model_branches)) {
    for (c in seq_along(model_branches[[t]])) {
      # if (index != "") {
      #   cur_branch <- stringr::str_replace_all(model_branches[[t]][[c]][["Equation"]],
      #                          "([:alpha:][:alnum:]?)", paste0("\\1[", index, "]"))
      #   extra <- paste0(index, ", ")
      # } else {
      #   cur_branch <- model_branches[[t]][[c]][["Equation"]]
      #   extra <- ""
      # }
      model_out[cat] <-  paste0(
        paste0("prob", cat, " <- "),
        paste(model_branches[[t]][[c]][["Equation"]], collapse = " + ")
      )
      if (c == 1) {
        prob_out[t] <- paste0("c(prob", cat, "[", index,"]", ",")
      } else if (c == length(model_branches[[t]])) {
        prob_out[t] <- paste0(prob_out[t], "prob", cat, "[", index,"]", ")")
      } else {
        prob_out[t] <- paste0(prob_out[t], "prob", cat, "[", index,"]",",")
      }
      cat <- cat + 1
    }
  }
  cat(paste(model_out, collapse = "\n"), "\n\n", paste(prob_out))
  return(NULL)
}

find.MPT.params <- function(model_list) {
	tmp <- lapply(model_list,all.vars)
	return(unique(sort(unique(unlist(tmp)))))
}


check.MPT.probabilities <- function(model_list){
	tmp.env <- new.env()
	temp.param.names <- find.MPT.params(model_list)
	temp.param.val <- runif(length(temp.param.names))
	temp.branch <- vapply(model_list,length, NA_integer_)
	prob <- rep(NA_real_,length(temp.branch))
	for (i in 1:length(temp.param.val)) {
		assign(temp.param.names[i],temp.param.val[i], envir = tmp.env)
	}
	tmp_res <- vector("list", length(model_list))
	for (i in seq_along(model_list)) {
	  tmp_res[[i]] <- vapply(model_list[[i]], FUN = eval, FUN.VALUE = NA_real_, envir = tmp.env)
	  names(tmp_res[[i]]) <- attr(model_list, "cat_map")[[i]]
	}
	prob <- vapply(tmp_res, sum, NA_real_)
	prob <- round(prob,digits=6)
	return(prob)
}

