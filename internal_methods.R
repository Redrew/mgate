impute_null <- function(md, p, aa) {
  for (exp in md$experiments) {
    if (is.null(md$pos[[p]][[aa]][[exp]]$score)) {
      return(F)
    }
  }
  return(T)
}

miter <- function(md, experiment=F, skip_null=T) {
  iter <- list()
  for (p in names(md$pos)) {
    for (aa in MGATE_AA) {
      if (skip_null & !impute_null(md, p, aa)) next()
      if (experiment) {
        for (exp in md$experiments)
          iter[[length(iter)+1]] <- list(p, aa, exp, length(iter)+1)
      } else iter[[length(iter)+1]] <- list(p, aa, length(iter)+1)
    }
  }
  return(iter)
}

# formats the mutation codes in data table, 
# extracts the position and amino acids of the wildtype and mutations
parse_table <- function(data_table, offset) {
  head <- filter(data_table, !grepl("\\d", code))
  
  data_table <- filter(data_table, grepl("\\d", code) & !grepl(",", code))
  for (amino_acid in MGATE_AA) 
    data_table$code <- gsub(attr(amino_acid, "term"), amino_acid, data_table$code)
  data_table$code <- gsub("p.", "", data_table$code)
  data_table[["pos"]] <- str_extract(data_table$code, "\\d+") %>% as.numeric() + offset
  data_table[["aa"]] <- sapply(data_table$code, function(x) substr(x, nchar(x), nchar(x)))
  data_table[["base_aa"]] <- sapply(data_table$code, function(x) substr(x, 1, 1))
  data_table <- arrange(data_table, pos, aa)
  suppressWarnings(data_table <- bind_rows(head, data_table))
  
  return(data_table)
}

# loads mutation data given a data table in the format of combined scores
load_combined_values <- function(md, experiment, data_table, value_names) {
  md$wildtype[[experiment]][value_names] <- data_table[data_table$code=="_sy", value_names]
  md$wildtype[["is_wildtype"]] <- T
  
  data_table <- filter(data_table, !is.na(pos))
  by(data_table, 1:nrow(data_table), function(row) {
    p <- row[["pos"]] %>% as.character()
    if (is.null(md$pos[[p]][["base_aa"]])) {
      md$pos[[p]][[row[["base_aa"]]]] <- md$wildtype 
      md$pos[[p]][["base_aa"]] <- row["base_aa"]
    }
    md$pos[[p]][[row[["aa"]]]][[experiment]][value_names] <- row[value_names]
    md$pos[[p]][[row[["aa"]]]][["is_wildtype"]] <- F
  })
}

load_shared_values <- function(md, experiment, data_table, value_names) {
  selection <- data_table[1, 1:ncol(data_table)] %>% find_and_empty() %>% find_and_empty("selection")
  value <- data_table[2, 1:ncol(data_table)] %>% find_and_empty()
  replicates <- selection %>% as.character() %>% unique() %>% (function(x) x[x!=""])
  wt <- filter(data_table, code=="_sy")
  data_table <- filter(data_table, !is.na(pos))
  
  for (repl in replicates) {
    md$wildtype[[experiment]][["replicate"]][[repl]][value_names] <- wt[selection==repl] %>% as.numeric()
    md$wildtype[[experiment]][["replicate"]][["names"]] <- replicates
    
    by(data_table, 1:nrow(data_table), function(row) {
      p <- row[["pos"]] %>% as.character()
      md$pos[[p]][[row[["aa"]]]][[experiment]][["replicate"]][[repl]][value_names] <- row[selection==repl] %>% as.numeric()
      md$pos[[p]][[row[["aa"]]]][[experiment]][["replicate"]][["names"]] <- replicates
    })
  }
}

# reads `md` and creates a dataframe with columns for each value
# Args:
#     values: a character list, can include
#         "experiment" 
#         "replicate"
#         "score"
#         "SE"
#         [mfilter name]
#         "is_wildtype"
#         "c_0"
#         "c_1"
#         add / before the value names to access values at a higher level in the 
#         mutation data structure
#         eg. values=c("experiment", "replicate", "scores", "//scores") will both
#             have a column for replicate scores and experiment scores
# By default the data frame will have a "pos", "aa" (amino acid) and "code" (pos+aa) column
form_data_frame <- function(md, values, ignore=NULL, include_wildtype=F, wildtype=F) {
  derivative_apply <- function(f) {
    for (p_aa in miter(md)) {
      p = p_aa[[1]]; aa = p_aa[[2]]
      
      if (sum(sapply(ignore, function(x) md$pos[[p]][[aa]][[x]] == T) %>% unlist()) > 0) next()
      
      c(lvl, rows) %<-% f(md$pos[[p]][[aa]])
      
      print(paste("READING DATA - position:", p, aa, ", at row:", rows[1]))
      
      df[rows, "pos"] <<- rep(p, length(rows))
      df[rows, "aa"] <<- rep(aa, length(rows))
      df[rows, "code"] <<- rep(paste0(p,aa), length(rows))
      
      for (r in which(opts$level==lvl)) 
        df[rows, opts[r,"name"]] <<- rep(md$pos[[p]][[aa]][[opts[r,"value"]]], length(rows))
    }
  }
  wildtype_apply <- function(f) {
    c(lvl, rows) %<-% f(md$wildtype)
    
    df[rows, "pos"] <<- rep(NA, length(rows))
    df[rows, "aa"] <<- rep("wt", length(rows))
    df[rows, "code"] <<- rep("wt", length(rows))
    
    for (r in which(opts$level==lvl)) 
      df[rows, opts[r,"name"]] <<- rep(md$wildtype[[opts[r,"value"]]], length(rows))
  }
  experiment_wrapper <- function(f) {
    return(function(path) {
      total_rows <- c()
      for (exp in md$experiments) {
        c(lvl, rows) %<-% f(path[[exp]])
        df[rows, "experiment"] <<- rep(exp, length(rows))
        total_rows <- c(total_rows, rows)
        for (r in which(opts$level==lvl)) 
          df[total_rows, opts[r,"name"]] <<- rep(path[[exp]][[opts[r,"value"]]], length(total_rows))
      }
      
      return(list(lvl+1, total_rows))
    })
  }
  replicate_wrapper <- function(f) {
    return(function(path) {
      path <- path[["replicate"]]
      total_rows <- c()
      for (repl in path$names) {
        c(lvl, rows) %<-% f(path[[repl]])
        df[rows, "replicate"] <<- rep(repl, length(rows))
        total_rows <- c(total_rows, rows)
      }
      
      for (r in which(opts$level==lvl)) 
        df[total_rows, opts[r,"name"]] <<- rep(path[[opts[r,"value"]]], length(total_rows))
      
      return(list(lvl+1, total_rows))
    })
  }
  core <- function(path) {
    row <- nrow(df) + 1; lvl <- 1
    for (r in which(opts$level==lvl)) 
      df[row, opts[r,"name"]] <<- path[[opts[r,"value"]]]
    return(list(lvl+1, row))
  }
  
  def_values <- c("pos", "aa", "code")
  rec_modes <- intersect(c("experiment", "replicate"), values)
  paths <- exclude(values, rec_modes, def_values)
  opts <- make_df(c("name", "value", "level"))
  for (path in paths) {
    r <- nrow(opts)+1
    opts[[r, "level"]] <- str_count(path, "/") + 1
    opts[[r, "name"]] <- path
    opts[[r, "value"]] <- str_extract(path, "[^/]+")
  }
  cols <- c(def_values, rec_modes, opts$name)
  df <- make_df(cols)
  
  if (wildtype) 
    apply = wildtype_apply
  else
    apply = derivative_apply
  
  if ("replicate" %in% cols) 
    apply(experiment_wrapper(replicate_wrapper(core)))
  else if ("experiment" %in% cols) 
    apply(experiment_wrapper(core))
  else
    apply(core)
  
  if ("experiment" %in% cols) df$experiment <- factor(df$experiment, levels = md$experiments)
  if (!wildtype) {
    df$aa <- factor(df$aa, MGATE_AA)
    df$pos <- factor(df$pos, df$pos %>% unique() %>% sort(decreasing=T))
  }
  
  if (include_wildtype)
    df <- rbind(form_data_frame(md, values, ignore, wildtype=T), df)
  
  return(df)
}

plot_scatter <- function(md, df, mfilter=NULL, highlight=NULL) {
  
  p <- ggplot() + geom_jitter(aes(experiment, y=score, color=df[[highlight]]), df) 
  if (!is.null(mfilter)) {
    filter_limits <- md$mfilter[[mfilter]] %>% data.frame() %>% gather(value=yintercept, key=experiment)
    p <- p + facet_grid(cols=vars(experiment), scales="free_x", switch="x") + geom_hline(aes(yintercept=yintercept), filter_limits)
  }
  return(p)
}

plot_line <- function(md, df, mfilter, highlight=NULL) {
  #filter_limits <- md$mfilter[[mfilter]] %>% data.frame() %>% gather(value=yintercept, key=experiment)
  
  p <- ggplot() + geom_line(aes(x=experiment, y=score, group=code, color=df[[highlight]]), df) #+ 
    #facet_grid(cols=vars(experiment), scale="free") + geom_hline(aes(yintercept=yintercept), filter_limits)
  return(p)
}

plot_violin <- function(md, df, mfilter=NULL, jitter=T) {
  p <- ggplot() + geom_violin(aes(x=experiment, y=score), df) 
  if (!is.null(mfilter)) {
    filter_limits <- md$mfilter[[mfilter]] %>% data.frame() %>% gather(value=yintercept, key=experiment)
    p <- p + facet_grid(cols=vars(experiment), scales="free_x", switch="x") + geom_hline(aes(yintercept=yintercept), filter_limits)
  }
  if (jitter) p <- p + geom_jitter(aes(x=experiment, y=score), df, width=0.2)
  return(p)
}

plot_bar <- function(md, target_df, target_rep_df, mfilter) {
  p <- ggplot() +
  geom_bar(data=target_df, mapping=aes(x=experiment, y=score, fill=code),
           stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(data=target_df, mapping=aes(x=experiment, group=code, ymin=score_min, ymax=score_max),
                width=.2, position=position_dodge(.9)) +
  geom_point(data=target_rep_df, mapping=aes(x=experiment, y=score, fill=code), 
             pch = 21, position=position_dodge(0.9))
  if (!is.null(mfilter)) {
    filter_limits <- md$mfilter[[mfilter]] %>% data.frame() %>% gather(value=yintercept, key=experiment)
    p <- p + facet_grid(cols=vars(experiment), scale="free_x", switch="x") + geom_hline(aes(yintercept=yintercept), filter_limits)
  }
  
  return(p)
}

plot_heatmap <- function(df, mfilter) {
  df["property"] <- sapply(df$aa, function(aa) {
    attr(MGATE_AA[[aa]], "property")
    })
  p <- ggplot(df, mapping=aes(x=aa, y=pos)) + 
    geom_tile(aes(fill=df[[mfilter]]), colour="black") +
    geom_point(data=filter(df, is_wildtype==1), colour="gray10") +
    facet_grid(cols=vars(property), scales="free", space="free_x", switch="x")
  return(p)
}

# Not used
find_cut <- function(df, is_upper, range, sig=.95, step=0.05) {
  cuts_df <- data.frame(matrix(nrow=0, ncol=3)); names(cuts_df) <- c("cut", "violations", "experiment")
  
  sign <- 2 * is_upper - 1
  df$score <- df$score - sign * zscore(sig, sd=df$SE)
  df[["//score"]] <- df[["//score"]] - sign * zscore(sig, sd=df[["//SE"]])
  
  for (cut in seq(from=range[[1]], to=range[[2]], by=step)){
    violations <- 0
    for (r in 1:nrow(df)) {
      violations <- violations + 
        in_range(cut, c(df[[r, "//score"]], df[[r, "score"]]), 0)
    }
    
    cuts_df[nrow(cuts_df) + 1, c("cut", "violations", "experiment")] <- 
      c(cut, violations, df[[1, "experiment"]])
  }
  
  return(cuts_df)
}