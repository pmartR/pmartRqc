#' Function to take raw output of `imd_anova` and create output for `statRes` object
#' 
#' Not really an interesting function.  Just exporting it so Natalie can use in a different function.
#' 
#' @param imd_out Data.frame containing the results of the imd anova call
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @param comparisons the comparisons made
#' @param test_method the test method
#' @param pval_adjust pvalue adjustment method
#' @param pval_thresh p-value threshold value
#' 
#' @return the final object of class statRes
#' @export
#' 
statRes_output <- function(imd_out,omicsData,comparisons,test_method,pval_adjust,pval_thresh){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
   # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attr(omicsData, "group_DF")
  }
  
  #Flags to determine number of significant
  imd_out_flags <- imd_out[[grep("^flag",tolower(names(imd_out)))]]
  flags <- data.matrix(imd_out_flags)
  
  
  #Add biomolecule to P_values and Flags
  meta <- imd_out$Full_results[,1]
  imd_out$Flags <- cbind.data.frame(meta, flags)
  colnames(imd_out$Flags)[1] <- colnames(imd_out$Full_results)[1]
  
  imd_out_pvals <- imd_out[[grep("^p_values",tolower(names(imd_out)))]]
  pvals <- data.matrix(imd_out_pvals)
  imd_out$P_values <- cbind.data.frame(imd_out$Full_results[,1],pvals)
  colnames(imd_out$P_values)[1] <- colnames(imd_out$Full_results)[1]
  
  #Add attributes
  attr(imd_out, "group_DF") <- groupData
  #if(is.null(comparisons)){
  comparisons <- colnames(imd_out_flags)
  #}
  
  # Count positive, negative biomolecules per test
  total_groups <- pmartR::get_group_table(omicsData)
  nsig <- dplyr::bind_rows(lapply(comparisons, function(comp){
    base <- strsplit(comp, "_vs_")[[1]][2]
    nonbase <- strsplit(comp, "_vs_")[[1]][1]
    
    col_base <- grep(paste0("Count_", base), colnames(imd_out$Full_results))
    col_nonbase <- grep(paste0("Count_", nonbase), colnames(imd_out$Full_results))
    col_pval_g <- grep(paste0("P_value_G_", comp), colnames(imd_out$Full_results))
    
    prop_base <- imd_out$Full_results[col_base]/total_groups[base]
    prop_nonbase <- imd_out$Full_results[col_nonbase]/total_groups[nonbase]
    direction <- sign(prop_nonbase - prop_base)
    g_flag <- as.numeric(imd_out$Full_results[col_pval_g] < pval_thresh) * direction
    
    data.frame(
      Comparison = comp,
      Up_anova = sum(imd_out_flags[comp] == 1),
      Down_anova = sum(imd_out_flags[comp] == -1),
      Up_gtest = sum(g_flag == 1, na.rm = TRUE),
      Down_gtest = sum(g_flag == -1, na.rm = TRUE)
      )
  }))
  
  attr(imd_out, "comparisons") <- comparisons
  attr(imd_out, "number_significant") <-  nsig
  attr(imd_out,"statistical_test") <- test_method
  attr(imd_out, "adjustment_method") <- pval_adjust
  attr(imd_out, "pval_thresh") <- pval_thresh
  attr(imd_out, "data_info") <- attr(omicsData, "data_info")
  class(imd_out) <- "statRes"
  return(imd_out)
}

#' statRes class.
#'
#' Class for statistical results returned by fucntions in this package.
#'
#' @name statRes-class
#' @seealso See \code{\link{imd_anova}}
#'
#' @exportClass statRes
setOldClass("statRes")

#' @export
#' @method summary statRes
summary.statRes <- function(x,...){
  cat("Type of test:",attr(x,"statistical_test"),"\n\n")
  cat("Multiple comparison adjustment:",attr(x,"adjustment_method"),"\n\n")
  cat("p-value threshold:",attr(x,"pval_thresh"),"\n\n")
  cat("Number of significant biomolecules by comparison.  Columns specify fold change direction and type of test:\n\n")
  
  table <- attr(x,"number_significant")
  names(table) <- lapply(names(table), function(name){switch(name,
                                                     Up_total = "Total:Positive",
                                                     Down_total = "Total:Negative",
                                                     Up_anova = "Positive:ANOVA",
                                                     Down_anova = "Negative:ANOVA",
                                                     Up_gtest = "Positive:G-test",
                                                     Down_gtest = "Negative:G-test",
                                                     name)
  })
  
  rownames(table) <- NULL
  
  print(table)
  
  return(invisible(list(test_type = attr(x,"statistical_test"), adjustment = attr(x,"adjustment_method"), pval_thresh = attr(x,"pval_thresh"), sig_table = table,
                        comparisons = attr(x, "comparisons"), group_DF = attr(x,"group_DF"))))
}

#' @export
#' @method print statRes
print.statRes <- function(x,...){
  
  #print.default(x,max=6,...)
  xlen <- length(x)
  xnames <- names(x)
  for(i in 1:xlen){
    cat(xnames[i],":\n")
    print(head(x[[1]]))
    x[[1]] <- NULL
    cat("\n")
  }
  
  str(x,no.list=TRUE,give.head=FALSE)
}

#' Plotting function for `statRes` objects
#' 
#' Produces plots that summarize the results contained in a `statRes` object.
#' 
#' @param x `statRes` object to be plotted, usually the result of `imd_anova`
#' @param plot_method defines which test to use for bar plots, options are "anova" and "gtest"; required for "combined" statistics method, defaults to method used for testing for "anova" or "gtest" statistics
#' @param fc_colors vector of length two with character color values interpretable by ggplot. i.e. c("orange", "blue") with the values being used to color negative and positive fold changes respectively
#' @param invert TRUE/FALSE for whether to plot positive and negative fold change sections in opposite directions for the barplot, defaults to TRUE
#' @param interactive TRUE/FALSE for whether to create an interactive plot using plotly
#' @param bw_theme TRUE/FALSE for whether to apply a black-white background theme. Applies only when interactive = FALSE.
#' 
#' @export
#' @method plot statRes
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' #Transform the data
#' 
#' #Group the data by condition
#' myproData <- group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' #Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = myproData)
#' myproData <- applyFilt(filter_object = imdanova_Filt, omicsData = myproData, min_nonmiss_anova=2)
#' 
#' #Implement the IMD ANOVA method and compuate all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
#' anova_res <- imd_anova(omicsData = myproData, test_method = 'anova')
#' plot(anova_res)
#' 
#' imd_res <- imd_anova(omicsData = myproData, test_method = 'gtest')
#' plot(imd_res)
#' 
#' imd_anova_res <- imd_anova(omicsData = myproData, test_method = 'comb', pval_adjust='bon')
#' plot(imd_anova_res,  plot_method = "gtest", bw_theme = TRUE)
#' plot(imd_anova_res, plot_method = "anova", bw_theme = TRUE)
#' 
#' }
#' 
plot.statRes <- function(x, plot_method = NULL, fc_colors = c("red", "green"), invert = TRUE, interactive = FALSE, bw_theme = FALSE, ...){
 .plot.statRes(x = x, plot_method = plot_method, fc_colors = fc_colors, invert = invert, interactive = interactive, bw_theme = bw_theme, ...) 
}

.plot.statRes <- function(x, plot_method = NULL, fc_colors = c("red", "green"), invert = TRUE, interactive = FALSE, bw_theme = FALSE, 
                          title_theme = element_text(size = 14), axis_text_theme = element_text(size = 12), axis_title_theme = element_text(size = 14),
                          legend_title_theme = element_text(size = 10), legend_text_theme = element_text(size = 8), strip_text_theme  = element_text(size = 12),
                          custom_theme = NULL){
  #For now require ggplot2, consider adding base graphics option too
  if(!require(ggplot2)){
    stop("R package 'ggplot2' is required for non-interactive plots, base graphics aren't available yet")
  }
  
  if(interactive && !require(plotly)){
    stop("R package 'plotly' required for interactive plots")
  }
  
  #Most plots are based on "number_significant" data frame so pull it out
  comp_df <- attr(x,"number_significant")
  
  ##--------##
  #Go through the given plot_method and remove any that aren't currently available
  method <- attr(x, "statistical_test")
  if(method != "combined"){
    if(!is.null(plot_method) && tolower(plot_method) != method) warning(
      paste0(
        "Plot method '", 
        plot_method, 
        "' is invalid for this statRes object. Defaulting to valid plot method '", 
        method, "'"))
    plot_method <- method
  } else {
    if(is.null(plot_method)) stop("Please specify the plot_method argument as 'anova' or 'gtest' when plotting 'combined' statRes objects")
    if(!(tolower(plot_method) %in% c("anova","gtest"))) stop(
      paste0("Plot method '", 
             plot_method, 
             "' is invalid. Please specify the plot_method argument as 'anova' or 'gtest'")
    )
  }
  
  # specified theme parameters
  if(!is.null(custom_theme)){
    if(!bw_theme) warning("Setting both bw_theme to TRUE and specifying a custom theme may cause undesirable results")
    if(!inherits(custom_theme, c("theme", "gg"))) stop("custom_theme must be a valid 'theme' object as used in ggplot")
    mytheme = custom_theme
  }
  else mytheme = theme(plot.title = title_theme, axis.title = axis_title_theme, axis.text = axis_text_theme, legend.title = legend_title_theme, legend.text = legend_text_theme, strip.text = strip_text_theme)
  
  #Bar plot
  comp_df_melt <- reshape2::melt(comp_df,id.vars="Comparison",value.name="Count",variable.name="Direction")
  levels(comp_df_melt$Comparison) <- gsub(pattern = "_",replacement = " ",levels(comp_df_melt$Comparison))
  
  #Bar plots side-by-side, both going up
  #all_plts[[1]] <- ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+geom_bar(stat='identity',position='dodge')

  ##Up direction is positive, down direction is negative
  if(invert) comp_df_melt[grep("Down", comp_df_melt$Direction),]$Count <- (-comp_df_melt[grep("Down", comp_df_melt$Direction),]$Count)
  
  # add whichtest, and posneg columns used for plot grouping and label adjustment
  comp_df_melt <- comp_df_melt %>% 
    dplyr::mutate(whichtest = ifelse(grepl("anova", Direction), "anova", ifelse(grepl("gtest", Direction), "gtest", "total")),
                  posneg = ifelse(grepl("Up", Direction), "Positive", "Negative")) %>%
    dplyr::arrange(desc(posneg))
  
  # get only anova or only g-test rows if user did not specify combined
  if(attr(x, "statistical_test") %in% c("anova", "gtest")){
    comp_df_melt <- comp_df_melt %>%
      dplyr::filter(whichtest == attr(x, "statistical_test"))
  } else {
    comp_df_melt <- comp_df_melt %>%
      dplyr::filter(whichtest == plot_method)
  }
  
  if(interactive){
    
    names(fc_colors) <- c("Negative", "Positive")
    
    p <- plotly::plot_ly(
      x = comp_df_melt$Comparison,
      y = comp_df_melt$Count,
      type = "bar",
      color = comp_df_melt$posneg,
      colors = cols
    ) %>% plotly::layout(yaxis = list(title = "Count of Biomolecules"),
                         title = paste0("Number of DE Biomolecules Between Groups - ", plot_method))
    
    return(p)
    
  } else {
    
    p <- ggplot(data=comp_df_melt,aes(Comparison,Count)) + 
      geom_bar(aes(x = whichtest, fill = posneg, group = whichtest),stat='identity') + 
      geom_text(aes(x = whichtest, label = ifelse(abs(Count) > 0, abs(Count), "")), position = position_stack(vjust = 0.5), size = 3) +
      geom_hline(aes(yintercept=0),colour='gray50') +
      scale_fill_manual(values=c(fc_colors[1], fc_colors[2]), labels = c("Negative", "Positive"), name = "Fold Change Sign") +
      facet_wrap(~Comparison) + 
      xlab(NULL) +
      ylab("Count of Biomolecules") +
      ggtitle("Number of DE Biomolecules Between Groups")
    
    if(bw_theme) p <- p + theme_bw()
    
    return(p + mytheme)
    
  }
  
}



#' Volcano plot function for `statRes` objects
#' 
#' Produces volcano plots that summarize the results contained in a `statRes` object.
#' 
#' @param x `statRes` object to be plotted, usually the result of `imd_anova`
#' @param plot_method defines which test to use for volcano plots, options are "anova" and "gtest"; required for "combined" statistics method, defaults to method used for testing for "anova" or "gtest" statistics
#' @param fc_colors vector of length three with character color values interpretable by ggplot. i.e. c("orange","black", "blue") with the values being used to color negative, non-significant, and positive fold changes respectively
#' @param interactive TRUE/FALSE for whether to create an interactive plot using plotly
#' @param bw_theme TRUE/FALSE for whether to apply a black-white background theme. Applies only when interactive = FALSE.
#' @param fc_threshold optional threshold value for minimum absolute fold change estimates in the volcano plot.
#' @param custom_theme optional custom theme object recognized by ggplot. Not applicable for interactive = TRUE.
#' 
#' @export
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' #Transform the data
#' 
#' #Group the data by condition
#' myproData <- group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' #Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = myproData)
#' myproData <- applyFilt(filter_object = imdanova_Filt, omicsData = myproData, min_nonmiss_anova=2)
#' 
#' #Implement the IMD ANOVA method and compuate all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
#' anova_res <- imd_anova(omicsData = myproData, test_method = 'anova')
#' plot_volcano(anova_res)
#' 
#' imd_res <- imd_anova(omicsData = myproData, test_method = 'gtest')
#' plot_volcano(imd_res)
#' 
#' imd_anova_res <- imd_anova(omicsData = myproData, test_method = 'comb', pval_adjust='bon')
#' plot_volcano(imd_anova_res,  plot_method = "gtest", bw_theme = TRUE)
#' plot_volcano(imd_anova_res, plot_method = "anova", bw_theme = TRUE, interactive = TRUE)
#' 
#' }
#' 

plot_volcano <- function(x, 
                         plot_method = NULL, 
                         fc_colors = c("red","black", "green"),
                         interactive = FALSE, 
                         bw_theme = FALSE,
                         fc_threshold = NULL,
                         custom_theme = NULL){
  
  title_theme = element_text(size = 14)
  axis_text_theme = element_text(size = 12) 
  axis_title_theme = element_text(size = 14)
  legend_title_theme = element_text(size = 10) 
  legend_text_theme = element_text(size = 8) 
  strip_text_theme  = element_text(size = 12)
  
  #For now require ggplot2, consider adding base graphics option too
  if(!require(ggplot2)){
    stop("R package 'ggplot2' is required for non-interactive plots, base graphics aren't available yet")
  }
  
  if(interactive && !require(plotly)){
    stop("R package 'plotly' required for interactive plots")
  }
  
  ##--------##
  #Use correct plot_method
  method <- attr(x, "statistical_test")
  if(method != "combined"){
    if(!is.null(plot_method) && tolower(plot_method) != method) warning(
      paste0(
        "Plot method '", 
        plot_method, 
        "' is invalid for this statRes object. Defaulting to valid plot method '", 
        method, "'"))
    plot_method <- method
  } else {
    if(is.null(plot_method)) stop("Please specify the plot_method argument as 'anova' or 'gtest' when plotting 'combined' statRes objects")
    if(!(tolower(plot_method) %in% c("anova","gtest"))) stop(
      paste0("Plot method '", 
             plot_method, 
             "' is invalid. Please specify the plot_method argument as 'anova' or 'gtest'")
    )
  }
  
  # specified theme parameters
  if(!is.null(custom_theme)){
    if(!bw_theme) warning("Setting both bw_theme to TRUE and specifying a custom theme may cause undesirable results")
    if(!inherits(custom_theme, c("theme", "gg"))) stop("custom_theme must be a valid 'theme' object as used in ggplot")
    mytheme = custom_theme
  }
  else mytheme = theme(plot.title = title_theme, axis.title = axis_title_theme, axis.text = axis_text_theme, legend.title = legend_title_theme, legend.text = legend_text_theme, strip.text = strip_text_theme)
  
  ## Generate Plotting dfs
  comparisons <- colnames(x$Flags)[-1]
  total_groups <- pmartR::get_group_table(x)
  
  list_points <- lapply(comparisons, function(comp){
    base <- strsplit(comp, "_vs_")[[1]][2]
    nonbase <- strsplit(comp, "_vs_")[[1]][1]
    
    col_base <- grep(paste0("Count_", base), colnames(x$Full_results))
    col_nonbase <- grep(paste0("Count_", nonbase), colnames(x$Full_results))
    col_p_val_g <- grep(paste0("P_value_G_", comp), colnames(x$Full_results))
    col_p_val_t <- grep(paste0("P_value_T_", comp), colnames(x$Full_results))
    col_fold_change <- grep(paste0("Fold_change_", comp), colnames(x$Full_results))
    
    prop_base <- x$Full_results[col_base]/total_groups[base]
    prop_nonbase <- x$Full_results[col_nonbase]/total_groups[nonbase]
    per_change <- (prop_nonbase - prop_base) * 100

    pval_g <- unlist(x$Full_results[col_p_val_g])
    pval_t <- unlist(x$Full_results[col_p_val_t])
    fc <- unlist(x$Full_results[col_fold_change])
    
    
    df <- data.frame(
      Comparison = comp,
      Biomolecule = x$Full_results[[1]],
      `G-test p-value` = -log10(if(length(pval_g) == 0) NA else pval_g),
      `ANOVA p-value` = -log10(if(length(pval_t) == 0) NA else pval_t),
      `Difference in % present` = unlist(if(nrow(per_change) == 0) NA else per_change),
      `Fold change` = if(length(fc) == 0) NA else fc,
      check.names = FALSE
    )
    row.names(df) <- NULL
    return(df)
  })
  
  ## Plot functions
  plots <- lapply(list_points, function(comp_df){
    comp <- unique(comp_df$Comparison)
    
    if(plot_method == "anova"){
      
      comp_df$DE <- sign(comp_df$`Fold change`)
      comp_df$DE[comp_df$`ANOVA p-value` < -log10(attr(x, "pval_thresh"))] <- 0
      comp_df$DE <- gsub(1, "Positive", 
                         gsub(-1, "Negative", 
                              gsub(0, "Non-Significant", comp_df$DE)))
      
      names(fc_colors) <- c("Negative", "Non-Significant", "Positive")
      
      if(!is.null(fc_threshold)){
        comp_df <- comp_df[abs(comp_df$`Fold change`) > fc_threshold,]
      }
      
      if(interactive){
        
        a <- list(
          text = paste0(comp, " - ANOVA"),
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        )
        
        p <- plotly::plot_ly(
          x = comp_df$`Fold change`,
          y = comp_df$`ANOVA p-value`,
          color = comp_df$DE,
          colors = fc_colors,
          type = "scatter",
          mode = "markers",
          hovertext = comp_df$Biomolecule
        ) %>% plotly::layout(
          yaxis = list(title = "-log10 P-value"),
          xaxis = list(title = paste0("Fold change (", attr(x, "data_info")$data_scale, ")")),
          annotations = a, 
          showlegend = (comp == comparisons[1])
        )
      } else {
        
        p <- ggplot(data=comp_df, aes(x = `Fold change`, y = `ANOVA p-value`, color = DE)) +
          geom_point(position = position_jitter(width = 0.4, height = 0.4), shape = 1) + 
          scale_color_manual(values = fc_colors) +
          labs(color = NULL, subtitle = paste(comp, "- ANOVA"), 
               y = "-log10 p-value", 
               x = paste0("Fold change (", attr(x, "data_info")$data_scale, ")"))
        
        if(bw_theme) p <- p + theme_bw()
        
        p <- p + mytheme

      }
      
    } else {
      
      comp_df$DE <- sign(comp_df$`Difference in % present`)
      comp_df$DE[(comp_df$`G-test p-value` < -log10(attr(x, "pval_thresh")))] <- 0
      comp_df$DE <- gsub(1, "Positive", 
                         gsub(-1, "Negative", 
                              gsub(0, "Non-Significant", comp_df$DE)))
      
      # if(!is.null(fc_threshold)){
      #   comp_df <- comp_df[abs(comp_df$`Fold change`) > fc_threshold,]
      # }
      
      if(interactive){
        
        a <- list(
          text = paste0(comp, " - G-test"),
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        )
        
        hline <- function(y = 0, color = "blue") {
          list(
            type = "line", 
            xref = "paper",
            x0 = -50,
            x1 = 50,
            y0 = y, 
            y1 = y, 
            line = list(color = color, dash = "dot")
          )
        }
        
        comp_df2 <- reshape2::melt(
          table(comp_df[c("Difference in % present", "G-test p-value")]), 
          value.name = "Count")
        
        comp_df2 <- comp_df2[comp_df2$Count != 0,]
        comp_df2 <- comp_df2[comp_df2$`Difference in % present` != 0,]
        
        p <- plotly::plot_ly(
          x = comp_df2$`Difference in % present`,
          y = comp_df2$`G-test p-value`,
          color = comp_df2$Count,
          # colors = fc_colors,
          type = "scatter",
          mode = "markers"#,
          # hovertext = comp_df$Biomolecule
        ) %>% plotly::layout(
          yaxis = list(title = "-log10 P-value"),
          xaxis = list(title = "Difference in % present"),
          annotations = a, 
          showlegend = (comp == comparisons[1]),
          shapes = list(hline(-log10(attr(x, "pval_thresh")))))
        
        comp_df2 <- comp_df[comp_df$DE != "Non-Significant",]
        
        
      } else {
        
        p <- ggplot(data=comp_df, aes(x = `Fold change`, y = `ANOVA p-value`, color = DE)) +
          geom_point(position = position_jitter(width = 0.4, height = 0.4), shape = 1) + 
          scale_color_manual(values = fc_colors) +
          labs(color = NULL, subtitle = paste(comp, "- ANOVA"), 
               y = "-log10 p-value", 
               x = paste0("Fold change (", attr(x, "data_info")$data_scale, ")"))
        
        if(bw_theme) p <- p + theme_bw()
        
        p <- p + mytheme
        
      }
      
      
      
    }
    
  })
  
  #Both the volcano plot and heatmap need a dataframe of fold changes by comparison/biomolecule
    # fold change values for volcano plot
    fc_data <- x$Full_results[,c(1,grep("^Fold_change",colnames(x$Full_results)))]
    colnames(fc_data) <- gsub(pattern = "^Fold_change_",replacement = "",x = colnames(fc_data))
    fc_data <- reshape2::melt(fc_data,id.vars=1,variable.name="Comparison",value.name="Fold_change")

    # fold change flags for coloring
    fc_flags <- x$Flags
    fc_flags <- reshape2::melt(fc_flags,id.vars=1,variable.name="Comparison",value.name="Fold_change_flag") %>%
      dplyr::mutate(Fold_change_flag = as.character(Fold_change_flag))

    # p values for labeling and y axis in anova volcano plot
    p_data <- x$Full_results[c(1,grep("^P_value",colnames(x$Full_results)))]
    pvals <- reshape2::melt(p_data,id.vars=1,variable.name="Comparison",value.name="P_value")

    # grouping column based on test type
    if(attr(x,"statistical_test")=="combined"){
      pvals$Type <- "G-test"
      pvals$Type[grep(pattern="^P_value_T_",x=pvals$Comparison)] <- "ANOVA"
    }else if(attr(x,"statistical_test")=="gtest"){
      pvals$Type <- "G-test"
    }else if(attr(x,"statistical_test")=="anova"){
      pvals$Type <- "ANOVA"
    }

    levels(pvals$Comparison) <- gsub(pattern="^P_value_G_",replacement = "",levels(pvals$Comparison))
    levels(pvals$Comparison) <- gsub(pattern="^P_value_T_",replacement = "",levels(pvals$Comparison))

    volcano <- merge(merge(fc_data,pvals,all=TRUE), fc_flags, all = TRUE)

    # levels of comparison now of the form 'GROUPNAME_X vs GROUPNAME_Y'
    levels(volcano$Comparison) <- gsub(pattern = "_vs_",replacement = " vs ",levels(volcano$Comparison))

    # create counts for gtest plot (number present in each group)
    if(attr(x, "statistical_test") %in% c("gtest", "combined")){
      counts <- x$Full_results[c(1, grep("^Count_", colnames(x$Full_results)))]

      # trim column names so they are just group names
      colnames(counts) <- gsub("^Count_", replacement = "", colnames(counts))

      counts_df <- data.frame()
      for(comp in as.character(unique(volcano$Comparison))){
        #create a vector of the two group names being compared
        groups = strsplit(comp, " vs ")[[1]]
        gsize_1 = nrow(attr(x, "group_DF") %>% dplyr::filter(Group == groups[1]))
        gsize_2 = nrow(attr(x, "group_DF") %>% dplyr::filter(Group == groups[2]))

        # will contain ID column and count column corresponding to the two groups
        temp_df <- counts[c(1,which(colnames(counts) == groups[1]), which(colnames(counts) == groups[2]))]
        temp_df$Comparison <- comp

        # rename the columns to something static so they can be rbind-ed
        colnames(temp_df)[which(colnames(temp_df)%in% groups)] <- c("Count_First_Group", "Count_Second_Group")

        # store proportion of nonmissing to color g-test values in volcano plot
        temp_df$Prop_First_Group <- temp_df$Count_First_Group/gsize_1
        temp_df$Prop_Second_Group <- temp_df$Count_Second_Group/gsize_2

        counts_df <- rbind(counts_df, temp_df)
      }

      # should automatically left join by ID AND Comparison
      suppressWarnings(
        volcano <- volcano %>% dplyr::left_join(counts_df)
      )
    }
  
  #Volcano plot 
    # global jitter parameter and ID column
    jitter <- position_jitter(width = 0.4, height = 0.4)
    idcol <- colnames(volcano)[1]

    if(attr(x, "statistical_test") %in% c("anova", "combined")){
      # color vector which assigns black to gtest flag values (-2, 2)
      cols_anova <- c("-2" = fc_colors[2], "-1" = fc_colors[1], "0" = fc_colors[2], "1" = fc_colors[3], "2" = fc_colors[2])

      # temp data with rows only for ANOVA
      temp_data_anova <- volcano %>% dplyr::filter(Type == "ANOVA") %>%
        dplyr::mutate(Fold_change_flag = dplyr::case_when(Fold_change > 0 & P_value <= attr(x, "pval_thresh") ~  "1",
                                                          Fold_change < 0 & P_value <= attr(x, "pval_thresh") ~ "-1",
                                                          is.na(Fold_change) ~ "0",
                                                          TRUE ~ Fold_change_flag))

      # interactive plots need manual text applied to prepare for ggplotly conversion
      if(interactive){
        p1 <- ggplot(temp_data_anova, aes(Fold_change,-log(P_value,base=10), text = paste("ID:", !!rlang::sym(idcol), "<br>", "Pval:", round(P_value, 4))))
      }
      else p1 <- ggplot(data = temp_data_anova, aes(Fold_change,-log(P_value,base=10)))

      p1 <- p1 +
          geom_point(aes(color = Fold_change_flag), shape = 1)+
          facet_wrap(~Comparison) +
          ylab("-log[10](p-value)")+xlab(sprintf("Fold-change (%s)", attr(x, 'data_info')$data_scale)) +
          scale_color_manual(values = cols_anova, name = "Fold Change",
                             labels = c("Neg(Anova)", "0", "Pos(Anova)"),
                             breaks = c("-1", "0", "1"))

      if(bw_theme) p1 <- p1 + theme_bw()

      p1 <- p1 + mytheme
    }

    if(attr(x, "statistical_test") %in% c("gtest", "combined")){
      # assign black to anova flag values and filter down to G-test rows, reassign significant g-test flags so that anova significance doesn't cover them up
      cols_gtest <- c("-2" = fc_colors[1], "-1" = fc_colors[2], "0" = fc_colors[2], "1" = fc_colors[2], "2" = fc_colors[3])
      temp_data_gtest <- volcano %>%
        dplyr::filter(Type == "G-test") %>%
        dplyr::mutate(Fold_change_flag = dplyr::case_when(Prop_First_Group > Prop_Second_Group & P_value <= attr(x, "pval_thresh") ~  "2",
                                                          Prop_First_Group < Prop_Second_Group & P_value <= attr(x, "pval_thresh") ~ "-2",
                                                          is.na(Fold_change) ~ "0",
                                                          TRUE ~ Fold_change_flag))

      if(interactive){
        p2 <- ggplot(temp_data_gtest, aes(Count_Second_Group, Count_First_Group, text = paste("ID:", !!rlang::sym(idcol), "<br>", "Pval:", round(P_value, 4))))
      }
      else p2 <- ggplot(data=temp_data_gtest, aes(Count_Second_Group, Count_First_Group))

      p2 <- p2 +
        geom_point(data = temp_data_gtest %>% dplyr::filter(Fold_change_flag %in% c(-1,0,1)), aes(color = Fold_change_flag), position = jitter, shape = 1) +
        geom_point(data = temp_data_gtest %>% dplyr::filter(!(Fold_change_flag %in% c(-1,0,1))), aes(color = Fold_change_flag), position = jitter) +
        facet_wrap(~Comparison) +
        geom_vline(xintercept = c(unique(temp_data_gtest$Count_Second_Group) + 0.5, min(temp_data_gtest$Count_Second_Group) - 0.5)) +
        geom_hline(yintercept = c(unique(temp_data_gtest$Count_First_Group) + 0.5, min(temp_data_gtest$Count_First_Group) - 0.5)) +
        ylab("No. present in first group")+xlab("No. present in second group") +
        scale_color_manual(values = cols_gtest, name = "Fold Change",
                           labels = c("Neg(Gtest)", "0", "Pos(Gtest)"),
                           breaks = c("-2", "0", "2")) +
        scale_x_continuous(breaks = unique(temp_data_gtest$Count_Second_Group), labels = as.character(unique(temp_data_gtest$Count_Second_Group))) +
        scale_y_continuous(breaks = unique(temp_data_gtest$Count_First_Group), labels = as.character(unique(temp_data_gtest$Count_First_Group)))

      if(bw_theme) p2 <- p2 + theme_bw()

      p2 <- p2 + mytheme
    }

    # draw a line if threshold specified
    if(!is.null(fc_threshold)){
      p1 <- p1 + geom_vline(aes(xintercept=abs(fc_threshold)), lty = 2) + geom_vline(aes(xintercept = abs(fc_threshold)*(-1)), lty = 2)
    }

    # apply theme_bw() and interactivity if specified
    if(attr(x, "statistical_test") %in% "anova"){
      if(interactive) return(plotly::ggplotly(p1, tooltip = c("text"))) else return(p1)
    }
    if(attr(x, "statistical_test") %in% "gtest"){
      if(interactive) return(plotly::ggplotly(p2, tooltip = c("text"))) else return(p2)
    }
    if(attr(x, "statistical_test") %in% "combined"){

      if(interactive){
        p1 <- p1 %>% plotly::ggplotly(tooltip = c("text"))
        p2 <- p2 %>% plotly::ggplotly(tooltip = c("text"))
        suppressWarnings(
          p <- plotly::subplot(p1, p2, nrows = 2) %>%
            plotly::layout(showlegend = FALSE,
                           title = sprintf("TOP: -log10-pvalue vs %s fold change | BOTTOM:  #Present in each group | Colored by fold change direction",
                                                               attr(x, 'data_info')$data_scale),
                           font = list(size = 10))
        )
        return(p)
      }

      else{
        gridExtra::grid.arrange(p1,p2,nrow = 2)
        return(invisible(list(p1,p2)))
      }
    }
  
}