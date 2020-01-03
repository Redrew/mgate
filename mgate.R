# External library dependencies
library(tidyverse)
library(rlang)
library(zeallot)
library(gridExtra) 

# External function files
source("utility.R")
source("internal_methods.R")
source("amino_acids.R")

# This script contains the API functions that you would be using to load, gate, 
# and visualise DMS data.

# [[mutation data methods]]

# This function returns a data structure that is capable of holding all your Enrich2 DMS data.
# `md` is passed around as an argument to most mgate functions.
# Args:
#     experiments: a character vector of the names for each experiment
#     offset: an offset added to the positions of all mutations, use to correct
#         discrepencies between Enrich2 positions and biological positions
# 
# After being loaded, information can be retrieved from `md` through navigating to
# the right level
# $experiments: the list of experiment names
# $mfilter: the list of mfilters, used by mgate to filter mutations
# $wildtype: values for the wildtype sequence is also stored at the top level
#   $score
#   $SE
# $pos: values for each mutation are stored by position and amino acid mutation
#   $[number]
#     $[amino acid letter]
#       $[experiment name]: 
#         $score
#         $SE
#         $replicate
#           $[Rep#]
#             $score
#             $SE
mdata_struct <- function(experiments, offset) {
  md <- env(experiments=experiments,
            offset=offset, wildtype=env(), pos=NULL)
  return(md)
}

# [[loader method]]

# This function loads a list of Enrich2 data files into your md data structure.
# Args:
#     files: need to contain the data files corresponding to each experiment in the
#         exact order as the experiments vector used to create your md data structure
#     type: could be 
#       "combined scores" for main_synonomous_scores.tsv data files 
#           with combine replicate values
#       "replicate scores" for main_synonomous_scores_shared.tsv data files 
#           with separate replicate values
#       "replicate counts" for main_synonomous_counts_shared.tsv data files 
#           with initial and final counts for each replicate
load_mfiles<- function(md, files, type, offset=NULL) {
  if (is.null(md$experiments)) stop("Experiments not in mutation data")
  if (is.null(offset) && is.null(md$offset)) stop("Offset is not provided")
  if (length(md$experiments) != length(files)) stop("Experiments and Score files do not match")
  
  if (is.null(offset)) offset <- md$offset
  
  # for each experiment, create a data.frame from the data file and load into md
  for (n in seq_along(md$experiments)) {
    if (grepl("combined", type)) {
      data_table <- read.table(files[n], sep = "\t", skip = 1, header = T, stringsAsFactors = F)
      # change the left most column's name to code so parse_table recognises the mutation column
      names(data_table)[1] <- "code"
      data_table <- parse_table(data_table, offset)
    }
    if (grepl("replicate", type)) {
      data_table <- read.table(files[n], sep = "\t", header = T, stringsAsFactors = F)
      # change the left most column's name to code so parse_table recognises the mutation column
      names(data_table)[1] <- "code"
      data_table <- parse_table(data_table, offset)
    }
    
    if (type == "combined scores") load_combined_values(md, md$experiments[n], data_table, c("SE", "score"))
    if (type == "replicate scores") load_shared_values(md, md$experiments[n], data_table, c("SE", "score"))
    if (type == "replicate counts") load_shared_values(md, md$experiments[n], data_table, c("c_0", "c_1"))
  }
}

# [[filter methods]]

# This function sets the boundaries of a gating mutation filter.
# Creates new mfilter if none exists with the name.
# Gates the mutations in `md` if all boundaries are supplied
# Args:
#     name: a numeric pair vector. name each pair to a specific experiment if you are only updating
#         the boundaries of specific experiments
# Calls:
#     mgate(...)
set_mfilter <- function(md, name, limits=NULL) {
  if (is.null(md$mfilter[[name]])) {
    md$mfilter[[name]] <- as.list(rep(NA, length(md$experiments)))
    names(md$mfilter[[name]]) <- md$experiments
    attr(md$mfilter[[name]], "name") <- name
  }
  if (!is.null(limits)) {
    if (is.null(names(limits))) names(limits) <- names(md$experiments)
    md$mfilter[[name]] <- modifyList(md$mfilter[[name]], limits)
  }
  
  if (!NA %in% md$mfilter[[name]])
    mgate(md, sig=.95, list(md$mfilter[[name]]))
}

# This function gates each mutation by provided list of mfilters.
# Assumes that the mutations are within the boundaries of the mfilter as the
# null hypothesis.
# Args:
#     sig: the significance level for concluding the mutation is not within the mfilter
#     mfilters: the list of mfilters to filter by
# 
# The function gates mutation for each experiment and as a whole.
#     Mutation is filtered if it's combined scores are within the mfilter boundaries
#         in every experiment
#     Given an experiment, a mutation is filtered if it's combined scores are within the  
#         mfilter boundaries 
#     Given an experiment, a mutation replicate is filtered if it's scores are within the  
#         mfilter boundaries
# mfilter results are recorded within the `md` data structure as follows
# ...$[pos]
#       $[amino acid]: [mfilter_name] = # mutation in filter over all experiments
#                      / # experiments
#         $[experiment name]: [mfilter_name] = mutation in filter, 1 if the combined score
#                             for the given experiment is within the filter, else 0
#           $replicates: [mfilter_name] = # replicate mutation in filter over all experiments
#                        / # experiments
#             $[Rep#]: [mfilter_name] = replicate mutation in filter, 1 if the replicate
#                             score for the given experiment is within the filter, else 0
mgate <- function(md, sig=.95, mfilters=md$mfilter) {
  for (p_aa in miter(md, skip_null=T)) { p = p_aa[[1]]; aa = p_aa[[2]];
    
    for (mfilter in mfilters) { exp_sum <- 0;

      for (exp in md$experiments) { repl_sum <- 0; replicates <- md$pos[[p]][[aa]][[exp]]$replicate$names;

        for (repl in replicates) {
          repl_values <- md$pos[[p]][[aa]][[exp]]$replicate[[repl]]
          repl_in <- in_range(x=repl_values[["score"]],
                              range=mfilter[[exp]], 
                              sd=repl_values[["SE"]], 
                              sig=sig)
          md$pos[[p]][[aa]][[exp]]$replicate[[repl]][[attr(mfilter, "name")]] <- repl_in

          repl_sum <- repl_sum + repl_in
        }

        repl_p <- repl_sum / length(replicates)
        md$pos[[p]][[aa]][[exp]]$replicate[[attr(mfilter, "name")]] <- repl_p

        exp_values <- md$pos[[p]][[aa]][[exp]]
        exp_in <- in_range(x=repl_values[["score"]],
                           range=mfilter[[exp]], 
                           sd=repl_values[["SE"]], 
                           sig=sig)
        md$pos[[p]][[aa]][[exp]][[attr(mfilter, "name")]] <- exp_in
        exp_sum <- exp_sum + exp_in
      }

      md$pos[[p]][[aa]][[attr(mfilter, "name")]] <- exp_sum / length(md$experiments)
    }

  }
}

# This function ______
filter_inactive <- function(md, quant, sig=0.95) {
  limits <- list()
  df <- form_data_frame(md, c("experiment", "score", "SE")) %>% filter(aa=="%")
  df$score <- df$score + zscore(sig, sd=df$SE)
  for (exp in md$experiments) {
    exp_df <- filter(df, experiment==exp)
    limits[[exp]] <- c(-Inf, quantile(exp_df[["score"]], quant))
  }
  set_mfilter(md, "filter_inactive", limits)
}

# [[plotter methods]]
# ______
mplot <- function(md, type, mfilter=NULL, ignore=NULL, args=NULL) {
  if (type == "terminated") {
      df <- form_data_frame(md, c("experiment", "score"), ignore)
      p <- ggplot() + geom_jitter(aes(x=experiment, y=score, color=aa=="%"), df) +
        labs(title="Plot of terminated mutations")
  }
  if (type == "filter by experiment") {
      df <- form_data_frame(md, c("experiment", "score", mfilter), ignore)
      p <- plot_scatter(md, df, mfilter=mfilter, highlight=mfilter) +
        labs(title=paste("Mutation scores filtered by", mfilter, "in each experiment"))
  }
  if (type == "violin") {
    df <- form_data_frame(md, c("experiment", "score", mfilter), ignore)
    p <- plot_violin(md, df, mfilter=mfilter) +
      labs(title="Violin plot of scores")
  }
  if (type == "filter") {
    path <- paste0("//", mfilter)
    df <- form_data_frame(md, c("experiment", "score", path), ignore)
    df[[path]] <- df[[path]] == 1
    p <- plot_scatter(md, df, mfilter=mfilter, highlight=path) +
      labs(title=paste("Filtered by", mfilter))
  }
  if (type == "filter line") {
    path <- paste0("//", mfilter)
    df <- form_data_frame(md, c("experiment", "score", path), ignore)
    df[[path]] <- df[[path]] == 1
    p <- plot_line(md, df, mfilter=mfilter, highlight=path) +
      labs(title=paste("Filtered by", mfilter))
  }
  if (type == "filtered") {
    path <- paste0("//", mfilter)
    df <- form_data_frame(md, c("experiment", "score", path), ignore)
    df <- df[df[[path]]==1,] 
    p <- plot_scatter(md, df, mfilter=mfilter, highlight=path) +
      labs(title=paste("Mutations that match", mfilter))
  }
  if (type == "unfiltered") {
    path <- paste0("//", mfilter)
    df <- form_data_frame(md, c("experiment", "score", path, "//is_wildtype"), ignore)
    df <- df[df[[path]]!=1,] 
    p <- plot_scatter(md, df, mfilter=mfilter, highlight=path) +
      labs(title=paste("Mutations not filtered by", mfilter))
  }
  if (type == "near misses") {
    path <- paste0("//", mfilter)
    df <- form_data_frame(md, c("experiment", "score", path), ignore)
    df <- df[df[[path]] < 1 & df[[path]] > 0,]
    p <- plot_scatter(md, df, mfilter=mfilter, highlight=path) +
      labs(title=paste("Mutations similar to", mfilter))
  }
  if (type == "replicate violin") {
    path <- paste0("//", mfilter)
    df <- form_data_frame(md, c("experiment", "replicate", "score", path), ignore)
    p <- plot_violin(md, df, mfilter=mfilter, jitter=F) +
      labs(title="Violin plot of replicate scores")
    p
  }
  if (type == "replicate violation") {
    path <- paste0("/", mfilter)
    df <- form_data_frame(md, c("experiment", "replicate", "score", path), ignore)
    df <- df[df[[path]] < 1 & df[[path]] > 0,]
    p <- plot_scatter(md, df, mfilter=mfilter, highlight=path) +
      labs(title=paste("Replicates misfiltered by", mfilter))
    p
  }
  if (type == "filter heatmap") {
    df <- form_data_frame(md, c(mfilter, "is_wildtype", ignore))
    for (ifilter in ignore)
      df[df[ifilter]==T, mfilter] <- NA
    df[[mfilter]] <- df[[mfilter]] == 1
    
    p <- plot_heatmap(df, mfilter) +
      labs(title=paste("Heatmap of", mfilter, "mutations"))
    p
  }
  if (type == "bar") {
    df <- form_data_frame(md, c("experiment", "score", "SE"), ignore, include_wildtype=T)
    df["score_max"] <- df["score"] + zscore(0.95) * df["SE"]
    df["score_min"] <- df["score"] - zscore(0.95) * df["SE"]
    
    rep_df <- form_data_frame(md, c("experiment", "replicate", "score"), ignore, include_wildtype=T)
    
    target_code <- (unique(df$code) %>% exclude("wt"))[[1]]
    
    target_df <- filter(df, code==target_code | code=="wt")
    target_rep_df <- filter(rep_df, code==target_code | code=="wt")
    
    p <- plot_bar(md,  target_df, target_rep_df, mfilter) +
      labs(title=paste("Scores of", target_code))
  }
  if (type == "all bars") {
    df <- form_data_frame(md, c("experiment", "score", "SE"), ignore, include_wildtype=T)
    df["score_max"] <- df["score"] + zscore(0.95) * df["SE"]
    df["score_min"] <- df["score"] - zscore(0.95) * df["SE"]
    
    rep_df <- form_data_frame(md, c("experiment", "replicate", "score"), ignore, include_wildtype=T)

    plots <- list()
    for (target_code in unique(df$code) %>% exclude("wt")) {
      target_df <- filter(df, code==target_code | code=="wt")
      target_rep_df <- filter(rep_df, code==target_code | code=="wt")
      
      plots[[length(plots) + 1]] <- plot_bar(md,  target_df, target_rep_df, mfilter) +
        labs(title=paste("Scores of", target_code))
    }
    p <- marrangeGrob(plots, ncol=1, nrow=1)
  }
  
  if(is.null(p)) stop("Plot type is invalid")
  
  return(list(p, df))
}

mlist <- function(md, mfilters, ignore) {
  #______
}
