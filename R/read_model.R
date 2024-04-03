#' Read MPT model
#'
#'
#'
read_mpt <- function(file, text,
                       trees, categories,
                       type = c("easy", "eqn", "eqn2")) {
  if (!missing(file) && isTRUE(all.equal(type, c("easy", "eqn", "eqn2"))) &&
      (grepl("\\.eqn$", file) || grepl("\\.EQN$", file))) {
    type <- "eqn"
    message("model type auto-detected as 'eqn'")
  } else {
    type <- match.arg(type)
  }
  if (type == "eqn2") {
    ignore_first <- FALSE
    type <- "eqn"
  } else {
    ignore_first <- TRUE
  }
  model_df <- switch(type,
                  easy = .read.MPT.model(file = file, text = text, trees = trees,
                                         categories = categories),
                  eqn = .read.EQN.model(file = file, text = text,
                                        ignore_first = ignore_first)
  )
  # eqn_extra <- switch(type,
  #                 easy = NULL,
  #                 eqn = .get.EQN.model.order(file = file, text = text,
  #                                       ignore_first = ignore_first)
  # )
  #list(model = model, order = eqn_extra)
  return(model_df)
}

.read.MPT.model <- function(file, text, trees, categories) {
	if (!missing(file)) {
	  whole <- readLines(con = file)
	} else if (!missing(text)) {
	  whole <- stringr::str_split_1(text, "\\n")
	} else {
	  stop("no model specified.", call. = FALSE)
	}
	whole <- gsub("#.*", "", whole)
	whole <- stringr::str_replace_all(whole, " ", "")
	max_rows <- (max(stringr::str_count(whole, "\\+"))+1)*length(whole)
	model_df <- data.frame(tree = rep(NA_character_, max_rows),
	                       category = NA_character_, path = NA_character_)
	if (missing(trees)) {
	  trees <- paste0("t", seq_along(whole))
	}
	if (missing(categories)) {
	  categories <- paste0("c", seq_along(whole))
	}

	row <- 0
	c_tree <- 1
	c_cat <- 1
	s.flag <- FALSE
	for (c1 in seq_along(whole)) {
		if (!(grepl("^[[:space:]]*$", whole[c1]))) {
			s.flag <- TRUE
			path_split <- stringr::str_split_1(whole[c1], stringr::fixed("+"))
			path_seq <- seq_along(path_split)
			model_df[["tree"]][path_seq + row] <- trees[c_tree]
			model_df[["category"]][path_seq + row] <- categories[c_cat]
			model_df[path_seq + row, "path"] <- path_split
			#model[[c2]][c3] <- parse(text = whole[c1])[1]
			c_cat <- c_cat + 1
			row <- row + max(path_seq)
			fin <- row
		}
		else {
			if (s.flag == TRUE) c_tree <- c_tree + 1
			s.flag <- FALSE
		}
	}
	return (model_df[1:fin,])
}

.read.EQN.model <- function(file, text, ignore_first) {
	if (!ignore_first) {
	  nrows <- read.table(file = file, nrows = 1, text = text)[1,1]
	} else {
	  nrows <- -1
	}
	tmp.in <- read.table(file = file, skip = 1, stringsAsFactors = FALSE,
	                     nrows = nrows, text = text)
	colnames(tmp.in) <- c("tree", "category", "path")
	tmp.in
}



# .get.EQN.model.order <- function(file, text, ignore_first) {
# 		if (!ignore_first) {
# 	  nrows <- read.table(file = file, nrows = 1, text = text)[1,1]
# 	} else {
# 	  nrows <- -1
# 	}
# 	tmp.in <- read.table(file = file, skip = 1, stringsAsFactors = FALSE,
# 	                     nrows = nrows, text = text)
# 	order_trees <- sort(unique(tmp.in$V1))
# 	tmp.ordered <- tmp.in[order(tmp.in$V1),]
# 	tmp.spl <- split(tmp.ordered, factor(tmp.ordered$V1))
# 	tmp.spl <- unlist(unname(lapply(tmp.spl, function(d.f) sort(unique(d.f[,2])))))
# 	list(order.trees = order_trees, order.categories = tmp.spl)
# }

