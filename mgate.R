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
#     experiments: a list of the names for each experiment
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
#         exact order as the experiments list used to create your md data structure
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
#     name: the name of the mfilter
#     limits: a numeric pair list. name each pair to a specific experiment if you are only updating
#         the boundaries of specific experiments.
#         Use the following format:
#             limits[[experiment name]][[1]]: lower boundary
#             limits[[experiment name]][[2]]: upper boundary
#         Example:
#             limits <- list()
#             limits[["RPMI"]][[1]] <- -Inf
#             limits[["RPMI"]][[1]] <- 1
#             set_mfilter(md, "filter_inactive", limits)
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
        exp_in <- in_range(x=exp_values[["score"]],
                           range=mfilter[[exp]], 
                           sd=exp_values[["SE"]], 
                           sig=sig)
        md$pos[[p]][[aa]][[exp]][[attr(mfilter, "name")]] <- exp_in
        exp_sum <- exp_sum + exp_in
      }

      md$pos[[p]][[aa]][[attr(mfilter, "name")]] <- exp_sum / length(md$experiments)
    }

  }
}

# This function creates a mfilter for loss of function mutations named "filter_inactive".
# Lower bounds are set to -Inf
# Upper bounds are set to the top quantile of termination mutation's scores + SE
# Filters all mutations that are lower than a top quantile of termination mutations
# Args:
#     quant: the quantile of termination mutations
#     sig: controls whether the function considers the SE when calculating the top quantile
#         if sig is not 0, the function will find the top quantile of score + SE * z-score
filter_inactive <- function(md, quant, sig=.95) {
  limits <- list()
  df <- form_data_frame(md, c("experiment", "score", "SE")) %>% filter(aa=="%")
  df$score <- df$score + df$SE * zscore(sig)
  for (exp in md$experiments) {
    exp_df <- filter(df, experiment==exp)
    limits[[exp]] <- c(-Inf, quantile(exp_df[["score"]], quant))
  }
  set_mfilter(md, "filter_inactive", limits)
}

# This function creates a mfilter for wildtype mutations named "filter_functional".
# Lower bounds are set to wildtype score - SE * range
# Upper bounds are set to wildtype score + SE * range
# Args:
#     range: the standard deviations away from the wildtype score allowed
filter_functional <- function(md, range=1) {
  limits <- list()
  for (exp in md$experiments) {
    limits[[exp]] <- md$wildtype[[exp]]$score + md$wildtype[[exp]]$SE * c(-range, range)
  }
  set_mfilter(md, "filter_functional", limits)
}

# This function creates a mfilter for constitutively activing mutations named 
#     "filter_mutation_activated".
# Lower bounds are set to the upper boundary of "filter_functional"
# Upper bounds are set to Inf
# Args:
#     range: the standard deviations away from the wildtype score allowed
filter_mutation_activated <- function(md) {
  if (is.null(md$mfilter$filter_functional)) 
    stop("filter_functional must be called before this function")
  limits <- list()
  for (exp in md$experiments) {
    limits[[exp]][[1]] <- md$mfilter$filter_functional[[exp]][[2]] 
    limits[[exp]][[2]] <- Inf
  }
  set_mfilter(md, "filter_mutation_activated", limits)
}

# [[plotter methods]]
# This function plots graphs
# Args:
#     type: used to control what type of plot to output, could be
#         "bar" plot a dynamite plot of the mutations provided in args 
#             (args can also be indices, 1, 2, 3…)
#         "bar all" plot dynamite plots of all mutations
#          "violin" plot a violin plot on top of a scatterplot
#         "filter heatmap" plot a heatmap of filtered and unfiltered mutations
#         "filter line" plot a line plot of mutations
#         "termination" plot a scatterplot highlighting termination mutations
#         "filtered" for a scatterplot of only filtered mutations
#         "unfiltered" 
#         "filter by experiment" for a scatterplot of which scores where filtered
#             in individual experiments
#         "replicate violation" for a scatterplot of replicate scores that were filtered
#             differently with the combined score
#         "replicate violin" plot a violin plot for the scores of replicates
#     mfilter: the name of the mfilter used in the plot. used by most plots
#     ignore: a list of mfilters. the functions will only plot mutations not gated
#         by the mfilters
#     args: optional inputs, used by "bar" plot
# Returns: a list
#     [[1]]: index 1 contains the plot
#     [[2]]: index 2 contains the data.frame used to make the plot
mplot <- function(md, type, mfilter=NULL, ignore=NULL, args=NULL) {
  if (type == "termination") {
      df <- form_data_frame(md, c("experiment", "score"), ignore)
      p <- ggplot() + geom_jitter(aes(x=experiment, y=score, color=aa=="%"), df) +
        labs(title="Plot of termination mutations")
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
    
    if (is.numeric(args))
      target_codes <- (unique(df$code) %>% exclude("wt"))[args]
    else
      target_codes = args

    target_df <- filter(df, code %in% c(target_codes, "wt"))
    target_rep_df <- filter(rep_df, code %in% c(target_codes, "wt"))
    
    p <- plot_bar(md,  target_df, target_rep_df, mfilter) +
      labs(title=paste("Scores of", paste(target_codes, collapse=" ")))
  }
  if (type == "bar all") {
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

# This function returns a list of mutations gated by or not gated by the given
# mfilters
# Args:
#     select: list of mfilter names, mutation list will only include mutations
#         gated by these mfilters
#     ignore: list of mfilter names, mutation list will only include mutations
#         not gated by these mfilters
#     values: list of values to also return with the mutations. can be:
#         "score"
#         "SE"
#         [filter name]
#         "is_wildtype"
#         if values is not null then a data.frame of each mutation and the
#         corresponding values will be returned
mlist <- function(md, select=NULL, ignore=NULL, values=NULL) {
  ret <- form_data_frame(md, c(select, values), ignore)
  for (mfilter in select) {
    ret <- ret[ret[[mfilter]]==T,]
  }
  if (is.null(values)) ret <- unique(ret$code)
  return(ret)
}
