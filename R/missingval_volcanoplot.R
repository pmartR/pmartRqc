#' Produces a volcano plot from statRes object
#' 
#' This function creates a fold change vs t-Test p-value volcano plot from a statRes object
#'
#' @param statRes is an object of class 'statRes' 
#' @param comparison can either be a character string name of the comparison to plot, or an integer index refering to the comparisons attribute vector
#' @param x_lab character string to be used for x-axis label. Defaults to NULL
#' @param ... further arguments
#' 
#' \tabular{ll}{
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{legend_title} \tab character string to be used for legend_title label. Defaults to NULL \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{vlines} \tab The x coordinate (integer in absolute value) where to draw vertical lines, defaults to NULL\cr
#' \code{pvalue_threshold} \tab numeric value, draws horizontal line at value, defaults to NULL\cr
#' \code{palette} \tab is a character string indicating the name of the RColorBrewer palette to use; "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"\cr
#' \code{x_lab_angle} \tab integer value indicating the angle of x-axis labels \cr
#' \code{coordinate_flip} \tab logical indicates whether to flip cartesian coordinates so that horizontal becomes vertical and vise versa, defaults to false \cr
#'}
#'
#' @return plots ggplot2 object
#' 
#' @details This function takes in a statRes object and plots log2 fold change vs -log10 t-Test P-values for a comparison from attr(statRes, “comparison”), ‘comparison’ parameter could be either a string or an integer refereing to an index of statRes comparisons
#'  
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' 
#' myproData <- group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' imdanova_Filt <- imdanova_filter(omicsData = myproData)
#' myproData <- applyFilt(filter_object = imdanova_Filt, omicsData = myproData, min_nonmiss_anova=2)
#' 
#' statRes_object <- imd_anova(omicsData = myproData, test_method = 'comb', pval_adjust='bon')
#' 
#' missingval_volcanoplot(statRes_object, comparison = 2, palette = "Set1")
#' missingval_volcanoplot(statRes_object, comparison = 3, pvalue_threshold = .005, palette = "Dark2", x_lab_angle = 45)
#' 
#'}
#' 
#' @name missingval_volcanoplot
#' @rdname missingval_volcanoplot
#' @export
#'

missingval_volcanoplot<- function(statRes, comparison, x_lab = NULL, ...) {
  require(ggplot2)
  .missingval_volcanoplot(statRes, comparison, x_lab, ...)
}

.missingval_volcanoplot<- function(statRes, comparison, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_title = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE, vlines = NULL, pvalue_threshold = NULL, palette = "YlOrRd", x_lab_angle = 0, coordinate_flip = FALSE){

#check that statRes object is of 'statRes' class
if(!inherits(statRes, "statRes")) stop("object must be of class 'statRes'")
  
#check that palette is in the list of RColorBrewer palettes
if(!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", 
                    "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", 
                    "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) stop("palette must be one of RColorBrewer palettes")
  
#check that comparison is of correct class, either integer or character
if(!inherits(comparison, c("character", "numeric"))) stop("comparison must be either of 'character' or 'numeric' class")  

#check that pval_threshold attr of statRes in not NULL
if(is.null(attr(statRes,"pval_thresh"))) stop("pval_thresh attribute of statRes was not found")  
  
#plot label checks
if(!is.null(title_plot)) {
  if(!is.character(title_plot)) stop("title_plot must be a character vector")
}
if(!is.null(x_lab)) {
  if(!is.character(x_lab)) stop("x_lab must be a character vector")
}
if(!is.null(y_lab)) {
  if(!is.character(y_lab)) stop("y_lab must be a character vector")
} 
if(!is.null(legend_title)) {
  if(!is.character(legend_title)) stop("legend_title must be a character vector")
}   
if(!(is.numeric(x_lab_angle))) stop("x_lab_angle must be numeric")
  
#make sure comparison is in attr(statRes, "comparisons")  
if(inherits(comparison, "character")){
  if(!(comparison %in% attr(statRes, "comparisons"))) stop(paste(comparison, "was not found in statRes attributes", sep = " ")) 
}

#check that when comparison is numeric it is a valid index  
if(inherits(comparison, "numeric")){
  if(!(comparison %in% c(1:length(attr(statRes, "comparisons"))))) stop("comparison must be valid index of comparisons attr of statRes object")
  
  #now we take comparison and re-assign it a character string instead of a numeric
  comparison = attr(statRes, "comparisons")[comparison]
}    

#lets check the len of comparison
if(length(comparison) == 1){
  
  #forming pvalue column name based on comparison
  pval_col_name = paste("P_value_T", comparison, sep = "_")
  
  #lets check if names(statRes$Full_results) contains pval_col_name 
  if(!(pval_col_name %in% colnames(statRes$Full_results))) stop(paste(pval_col_name, "was not found in statRes object", sep = " "))
  
  #forming fold change column name based on comparison
  fold_change_col_name = paste("Fold_change", comparison, sep ="_") 
  
  #lets check if names(statRes$Full_results) contains fold_change_col_name 
  if(!(fold_change_col_name %in% colnames(statRes$Full_results))) stop(paste(fold_change_col_name, "was not found in statRes object", sep = " "))
  
  #now lets look for pval_col_name in statRes$Full_results
  pval_ind = which(names(statRes$Full_results) %in% pval_col_name)
  
  #looking for fold_change_col_name in statRes$Full_results
  fold_change_ind = which(names(statRes$Full_results) %in% fold_change_col_name)

  #now lets subset the pvalue col and fold change col from Full_results data frame 
  pvalue_data = statRes$Full_results[,pval_ind]
  fold_change_data = statRes$Full_results[,fold_change_ind]
  
  #cbind data
  plotdata = cbind(pvalue_data, fold_change_data)
  plotdata = as.data.frame(plotdata)

  #now we can apply -log10() to the pvalue_data
  plotdata$pvalue_data = -log10(plotdata$pvalue_data)
  #plotdata$fold_change_data = log2(plotdata$fold_change_data) taking log2 of fold_change_data produces many NaN and inf, going to assume data is already log2 scale
  
  
  # make labels #
  xlabel <- ifelse(is.null(x_lab), "log2 Fold Change", x_lab)
  ylabel <- ifelse(is.null(y_lab), "-log10 t-Test P-value", y_lab)
  plot_title <- ifelse(is.null(title_plot), paste("Volcano Plot", comparison, sep = " "), title_plot)
  legendtitle<- ifelse(is.null(legend_title),"-log10 pvalue", legend_title)
  
  #checks for vlines parameter 
   if(!is.null(vlines)){
    if(!(is.numeric(vlines)) | vlines <= 0) stop("'vlines' must be positive non-zero and numeric")
     vert_lines = c(vlines, -vlines)
   }
  
  #checks for pvalue_threshold paramater
  if(!is.null(pvalue_threshold)){
    if(!(is.numeric(pvalue_threshold)) | pvalue_threshold < 0 | pvalue_threshold > 1) stop("'pvalue threshold' must be numeric, between zero and one")
    pval_thresh = -log10(pvalue_threshold)
  }
    
    p <- ggplot(plotdata, aes(x = fold_change_data, y = pvalue_data)) + geom_point(aes(color = pvalue_data)) +
      xlab(xlabel) +
      ylab(ylabel) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = title_size), axis.title.x = element_text(size = x_lab_size), axis.title.y = element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
      scale_color_distiller(palette = palette, name = legendtitle)
    
    if(coordinate_flip == TRUE){
      p = p + coord_flip()
    }
    
    if(bw_theme == TRUE){
      p = p + theme_bw() + theme(plot.title = element_text(size = title_size), axis.title.x = element_text(size = x_lab_size), axis.title.y = element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1))
    }
    if(!(is.null(vlines))){
      p = p + geom_vline(xintercept = vert_lines, color = "red", linetype = "dashed")
    }
    if(!is.null(pvalue_threshold)){
      p = p + geom_hline(yintercept = pval_thresh, color = "green", linetype = "dashed")
    }
  
  return(p)
}
  
}
