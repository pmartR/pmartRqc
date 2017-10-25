#' Applies rrollup function 
#' 
#' This function applies the rrollup function to a pepData object for each unique protein and returns a proData object. 
#' 
#' @param pepData an omicsData object of class 'pepData'
#' @param parallel logical indicating whether or not to use "doParallel" loop in applying rrollup function. Defaults to TRUE.
#' 
#' @return an omicsData object of class 'proData'
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = rrollup(pepData = pep_object) 
#'}
#' 
#' @rdname rrollup
#' @export

rrollup<- function(pepData, parallel = TRUE){
  
  # check that pepData is of appropraite class #
  if(class(pepData) != "pepData") stop("pepData is not an object of the appropriate class")
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  
  pep_id = attr(pepData, "cnames")$edata_cname
  pro_id = attr(pepData, "cnames")$emeta_cname
  
  pep = data.table(pepData$e_data)
  pro = data.table(pepData$e_meta[,c(pep_id, pro_id)])
  temp = data.table:::merge.data.table(x = pro, y = pep, by = pep_id, all.x = F, all.y = T)
  temp = as.data.frame(temp, check.names=FALSE)[,-which(names(temp)==pep_id)]
  
  #pull protein column from temp and apply unique function
  unique_proteins<- unique(temp[[pro_id]])
  
  if(parallel == TRUE){
    
  final_list<- vector("list", length(unique_proteins))
  
  require(doParallel)
  cores<- detectCores()
  cl<- makeCluster(cores)
  registerDoParallel(cl)
  
  r<-foreach(i=1:length(unique_proteins))%dopar%{
    
    row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
    current_subset<- temp[row_ind,]
    current_subset<- current_subset[,-which(names(temp) == pro_id)]
   
    #### Perform R_Rollup ####
    ## store number of peptides ##
    num_peps = nrow(current_subset)
    
    res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
    ## if only 1 peptide, set the protein value to the peptide ##
    if(num_peps==1){
      protein_val = unlist(current_subset)
    }else{
      ## Step 1: Select Reference Peptide -- peptide with least amount of missing data ##
      na.cnt = apply(is.na(current_subset),1,sum)
      least.na = which(na.cnt == min(na.cnt))
      
      ## If tied, select one with highest median abundance##
      if(length(least.na)>1){
        mds = apply(current_subset,1,median,na.rm=T)[least.na]
        least.na = least.na[which(mds==max(mds))]		
      }
      prot_val = unlist(current_subset[least.na,])
      
      ## Step 2: Ratio all peptides to the reference.  Since the data is on the log scale, this is the difference ##
      scaling_factor = apply(matrix(prot_val, nrow = num_peps, ncol = ncol(current_subset), byrow=T) - current_subset,1,median,na.rm=T)
      
      ## Step 3: Use the median of the ratio as a scaling factor for each peptide ##
      x_scaled = current_subset + matrix(scaling_factor, nrow = num_peps, ncol = ncol(current_subset))
      
      ## Step 4: Set Abundance as Median Peptide Abundance ##
      protein_val = apply(x_scaled, 2, median, na.rm=T)
    }
    
    res[1,] = protein_val
    res<- data.frame(res)
    names(res)<- names(current_subset)
    
    final_list[[i]]<- res
  }
  stopCluster(cl)
  
  final_result<- do.call(rbind, r)
  final_result<- cbind(unique_proteins, final_result)
  names(final_result)[1]<-pro_id
  
  samp_id = attr(pepData, "cnames")$fdata_cname
  data_scale = attr(pepData, "data_info")$data_scale
  data_norm = attr(pepData, "data_info")$data_norm
  
  #subsetting pepData$e_meta by 'unique_proteins' 
  emeta_indices<- match(unique_proteins, pepData$e_meta[[pro_id]])
  
  if(ncol(pepData$e_meta) == 2){
    e_meta = as.data.frame(pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)])
    names(e_meta)<-pro_id
  }else {e_meta = pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)]} 
  
  prodata = as.proData(e_data = data.frame(final_result, check.names=FALSE), f_data = pepData$f_data, e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, data_norm = data_norm)
  
  #updating prodata attributes
  attr(prodata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
  attr(prodata, "data_info")$data_types = attr(pepData, "data_info")$data_types
  attr(prodata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
  
  attr(prodata, "filters")<- attr(pepData, "filters")
  attr(prodata, "group_DF")<- attr(pepData, "group_DF")
  attr(prodata, "imdanova")<- attr(pepData, "imdanova")
  }
  
  #applying rrollup without doParallel
  
  else{
    final_list<- vector("list", length(unique_proteins))
    
    for(i in 1:length(unique_proteins)){
      
      row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
      current_subset<- temp[row_ind,]
      current_subset<- current_subset[,-which(names(temp) == pro_id)]
      
      #### Perform R_Rollup ####
      ## store number of peptides ##
      num_peps = nrow(current_subset)
      
      res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
      ## if only 1 peptide, set the protein value to the peptide ##
      if(num_peps==1){
        protein_val = unlist(current_subset)
      }else{
        ## Step 1: Select Reference Peptide -- peptide with least amount of missing data ##
        na.cnt = apply(is.na(current_subset),1,sum)
        least.na = which(na.cnt == min(na.cnt))
        
        ## If tied, select one with highest median abundance##
        if(length(least.na)>1){
          mds = apply(current_subset,1,median,na.rm=T)[least.na]
          least.na = least.na[which(mds==max(mds))]		
        }
        prot_val = unlist(current_subset[least.na,])
        
        ## Step 2: Ratio all peptides to the reference.  Since the data is on the log scale, this is the difference ##
        scaling_factor = apply(matrix(prot_val, nrow = num_peps, ncol = ncol(current_subset), byrow=T) - current_subset,1,median,na.rm=T)
        
        ## Step 3: Use the median of the ratio as a scaling factor for each peptide ##
        x_scaled = current_subset + matrix(scaling_factor, nrow = num_peps, ncol = ncol(current_subset))
        
        ## Step 4: Set Abundance as Median Peptide Abundance ##
        protein_val = apply(x_scaled, 2, median, na.rm=T)
      }
      
      res[1,] = protein_val
      res<- data.frame(res)
      names(res)<- names(current_subset)
      
      final_list[[i]]<- res
    }
   
    final_result<- do.call(rbind, final_list)
    final_result<- cbind(unique_proteins, final_result)
    names(final_result)[1]<-pro_id
    
    samp_id = attr(pepData, "cnames")$fdata_cname
    data_scale = attr(pepData, "data_info")$data_scale
    data_norm = attr(pepData, "data_info")$data_norm
    
    #subsetting pepData$e_meta by 'unique_proteins' 
    emeta_indices<- match(unique_proteins, pepData$e_meta[[pro_id]])
    
    if(ncol(pepData$e_meta) == 2){
      e_meta = as.data.frame(pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)])
      names(e_meta)<-pro_id
    }else {e_meta = pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)]} 
    
    prodata = as.proData(e_data = data.frame(final_result, check.names=FALSE), f_data = pepData$f_data, e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, data_norm = data_norm)
    
    #updating prodata attributes
    attr(prodata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
    attr(prodata, "data_info")$data_types = attr(pepData, "data_info")$data_types
    attr(prodata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
    
    attr(prodata, "filters")<- attr(pepData, "filters")
    attr(prodata, "group_DF")<- attr(pepData, "group_DF")
    attr(prodata, "imdanova")<- attr(pepData, "imdanova")
  }
  
  return(prodata)
  } 