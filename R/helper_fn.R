#' Return comparisons of statRes object
#' 
#' This function returns comparisons from statRes or trellData object
#' 
#' @param compObj is an object with the comparison attribute; specifically objects of class 'statRes' and 'trellData' objects derived from 'statRes' objects in \code{\link{format_data}}
#' @return returns a data frame with comparisons and their indices
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' 
#' my_prodata = group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' imdanova_Filt = imdanova_filter(omicsData = my_prodata)
#' my_prodata = applyFilt(filter_object = imdanova_Filt, omicsData = my_prodata, min_nonmiss_anova=2) 
#' imd_anova_res = imd_anova(omicsData = my_prodata, test_method = 'comb', pval_adjust='bon')
#' 
#' result = get_comparisons(imd_anova_res)
#'}
#' @rdname get_comparisons
#' @export
get_comparisons<- function(compObj){

  #check that compObj object is of 'statRes' or 'trellData' class
  if(!inherits(compObj, c("statRes", "trellData"))) stop("object must be of class 'statRes' or 'trellData'")
  
  #check that compObj object is of 'statRes' or 'trellData' class
  if(inherits(compObj, "trellData") && is.null(attr(compObj, "comparisons"))) stop("trellData object did not inherit 'comparisons' attribute; value is NULL")
  
  #pull comparisons attribute
  comp = attr(compObj, "comparisons")
  
  result = data.frame("comparisons" = as.character(comp), "index" = 1:length(comp), stringsAsFactors = FALSE)
  
  return(result)

}

#' Return data_class of statRes or trellData object
#' 
#' This function returns data_class attribute from statRes or trellData object, inherited from the omicsData used in \code{\link{imd_anova}} or \code{\link{format_data}}
#' 
#' @param dcObj an object of class 'statRes' or 'trellData'
#' @return returns the data_class attribute from a 'statRes' or 'trellData' object
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' 
#' my_prodata = group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' imdanova_Filt = imdanova_filter(omicsData = my_prodata)
#' my_prodata = applyFilt(filter_object = imdanova_Filt, omicsData = my_prodata, min_nonmiss_anova=2) 
#' imd_anova_res = imd_anova(omicsData = my_prodata, test_method = 'comb', pval_adjust='bon')
#' 
#' result = get_data_class(imd_anova_res)
#'}
#' @rdname get_data_class
#' @export
get_data_class<- function(dcObj){
  
  #check that compObj object is of 'statRes' class
  if(!inherits(dcObj, c("statRes", "trellData"))) stop("dcObj object must be of class 'statRes' or 'trellData'")
  
  result = attr(dcObj, "data_class")
  
  return(result)
  
}


#' Return check.names attribute of omicsData object
#' 
#' This function returns check.names attribute from omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @return returns logical, either TRUE or FALSE
#' @examples 
#' dontrun{
#' getchecknames(omicsData)
#'}
#' @rdname getchecknames
#' @export
getchecknames<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  result = attr(omicsData, "check.names")
  
  return(result)
}


#' Set check.names attribute of omicsData object
#' 
#' This function sets the check.names attribute of an omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param set_to logical indicating what to set check.names attribute to. Defaults to TRUE.
#' @return updated omicsData object with check.names attribute
#' @examples 
#' dontrun{
#' setchecknames(omicsData, set_to = TRUE)
#'}
#' @rdname setchecknames
#' @export
setchecknames<- function(omicsData, set_to = TRUE){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  #check that set_to is logical
  if(!is.logical(set_to)) stop ("set_to must be of class 'logical' ")
  
  attr(omicsData, "check.names")<- set_to
  
  return(omicsData)
}


#' Get group info of omicsObject object
#' 
#' This function returns the "group_DF" attribute of an omicsObject object
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a data.frame containing omicsObject object group info
#' @examples 
#' dontrun{
#' get_group_info(omicsObject)
#'}
#' @rdname get_group_info
#' @export
get_group_info<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  if(is.null(attr(omicsObject, "group_DF"))) stop("group_designation has not been run for omicsObject")
  
  res = attr(omicsObject, "group_DF")
  
  return(res)
}


#' Get group table
#' 
#' This function returns a table with number of samples per group
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a table containing number of samples per group
#' @examples 
#' dontrun{
#' get_group_table(omicsObject)
#'}
#' @rdname get_group_table
#' @export
get_group_table<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  if(is.null(attr(omicsObject, "group_DF"))) stop("group_designation has not been run for omicsObject")
  
  #should the result be constructed from omicsObject$f_data or from the group_DF attr? if so...
  group = attr(omicsObject, "group_DF")$Group
  
  return(table(group))
}


#' Get data scale
#' 
#' This function returns data scale attribute
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing data scale
#' @examples 
#' dontrun{
#' get_data_scale(omicsObject)
#'}
#' @rdname get_data_scale
#' @export
get_data_scale<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  data_scale = attr(omicsObject, "data_info")$data_scale
  return(data_scale)
}


#' Get data norm
#' 
#' This function returns data norm attribute
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing data norm
#' @examples 
#' dontrun{
#' get_data_norm(omicsObject)
#'}
#' @rdname get_data_norm
#' @export
get_data_norm<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  data_norm = attr(omicsObject, "data_info")$norm_info$is_normalized
  
  return(data_norm)
}

#' Get isobaric norm
#' 
#' This function returns the isobaric norm attribute
#' 
#' @param omicsData an object of the class 'pepData', 'isobaricpepData' or 'proData', usually created by \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}}.
#' @return a character string describing isobaric norm attribute
#' @examples 
#' dontrun{
#' get_isobaric_norm(omicsData)
#'}
#' @rdname get_isobaric_norm
#' @export
get_isobaric_norm<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "isobaricpepData"))) stop("omicsData must be of class 'pepData', 'proData' or 'isobaricpepData'")
  
  isobaric_norm = attr(omicsData, "isobaric_info")$norm_info$is_normalized
  
  return(isobaric_norm)
}

#' Get e_data cname
#' 
#' This function returns e_data cname
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing e_data cname
#' @examples 
#' dontrun{
#' get_edata_cname(omicsObject)
#'}
#' @rdname get_edata_cname
#' @export
get_edata_cname<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  edata_cname = attr(omicsObject, "cnames")$edata_cname
  return(edata_cname)
}

#' Get f_data cname
#' 
#' This function returns f_data cname
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing f_data cname
#' @examples 
#' dontrun{
#' get_fdata_cname(omicsObject)
#'}
#' @rdname get_fdata_cname
#' @export
get_fdata_cname<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  fdata_cname = attr(omicsObject, "cnames")$fdata_cname
  return(fdata_cname)
}

#' Get e_meta cname
#' 
#' This function returns e_meta cname
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing e_meta cname
#' @examples 
#' dontrun{
#' get_emeta_cname(omicsObject)
#'}
#' @rdname get_emeta_cname
#' @export
get_emeta_cname<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nrmData', 'statRes', or 'trellData'")
  
  #check if emeta is null
  if(is.null(omicsObject$e_meta) && 
     inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop(
       "e_meta is NULL in omicsObject, thus emeta_cname is also NULL")
  
  if(is.null(omicsObject$e_meta) && 
     inherits(omicsObject, "statRes")) stop(
       "emeta_cname of input statsRes object is dependent on omicsData used in imd_anova; emeta_cname is NULL")
  
  emeta_cname = attr(omicsObject, "cnames")$emeta_cname
  return(emeta_cname)
}

#' Get Filters
#' 
#' This function returns filters attribute
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @return a vector of filter names
#' @examples 
#' dontrun{
#' get_filters(omicsData)
#'}
#' @rdname get_filters
#' @export
get_filters<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  filters = names(attr(omicsData, "filters"))
  
  if(is.null(filters)) stop("no filters have been applied")
  
  return(filters)
}

#rollup combine_fn functions
combine_fn_mean<- function(x){
  if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}
}

combine_fn_median<- function(x){median(x, na.rm = T)}

