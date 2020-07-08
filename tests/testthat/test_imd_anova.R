context("Test error handling and output of imd_anova")
library(pmartR)
library(pmartRdata)
library(testthat)
#Trasnform the data
attr(pep_object,"cnames")$fdata_cname <- "SampleID"
mypepData <- edata_transform(omicsData = pmartRdata::pep_object, data_scale = "log2")

#Group the data by condition
mypepData <- group_designation(omicsData = mypepData, main_effects = c("Condition"))

#Apply the IMD ANOVA filter
imdanova_Filt <- imdanova_filter(omicsData = mypepData)
mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)

#More test objects
omicsData <- pro_object
omicsData$f_data$fakegroup <- c(rep("A",3), rep("B",4), rep("C",4))

# Weird groups
# omicsData$f_data$fakegroup <- 1:12
# omicsData$f_data$fakegroup2 <- c(1:8, 1, 9:11)

# fakegroup = 3 groups, Condition = 2 groups
omicsData_fakegroup <- group_designation(omicsData = omicsData, main_effects = c("fakegroup"))
omicsData_conditiongroup <- group_designation(omicsData = omicsData, main_effects = c("Condition"))

# make filters for the two objects
imdanovafilt_fakegroup <- imdanova_filter(omicsData = omicsData_fakegroup)
imdanovafilt_conditiongroup <- imdanova_filter(omicsData = omicsData_conditiongroup)

##---- Check that warnings and errors are produced in various situations --------##

expect_warning(imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='dunnett'))
expect_warning(imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='tukey'))

# test over 3 combinations of parameters
nonmiss_params <- list(list(2, NULL), list(NULL, 3), list(2, 3))
lapply(nonmiss_params, function(params){
  # two objects with different grouping structures
  omicsData_fakegroup <- applyFilt(filter_object = imdanovafilt_fakegroup, omicsData = omicsData_fakegroup, min_nonmiss_anova = params[[1]], min_nonmiss_gtest = params[[2]])
  omicsData_conditiongroup <- applyFilt(filter_object = imdanovafilt_conditiongroup, omicsData = omicsData_conditiongroup, min_nonmiss_anova = params[[1]], min_nonmiss_gtest = params[[2]])
  
  # test on each test type...
  for(test in list("anova", "gtest", "comb")){
    # ...and on both objects for each test type
    for(omicsData in list(omicsData_fakegroup, omicsData_conditiongroup)){
      statRes_obj <- imd_anova(omicsData = omicsData, test_method = test)
      flag_columns <- names(statRes_obj$Full_results)[grep("^Flag_", names(statRes_obj$Full_results))]
      
      # get different sections of the results
      anova_pval_columns <- names(statRes_obj$Full_results)[grep("^P_value_T", names(statRes_obj$Full_results))]
      g_pval_columns <- names(statRes_obj$Full_results)[grep("^P_value_G", names(statRes_obj$Full_results))]
      count_columns <- names(statRes_obj$Full_results)[grep("^Count_", names(statRes_obj$Full_results))]
      foldchange_columns <- names(statRes_obj$Full_results)[grep("^Fold_change_", names(statRes_obj$Full_results))]
      
      flags <- statRes_obj$Full_results[flag_columns]
      anova_pvals <- statRes_obj$Full_results[anova_pval_columns]
      g_pvals <- statRes_obj$Full_results[g_pval_columns]
      counts <- statRes_obj$Full_results[count_columns]
      fold_changes <- statRes_obj$Full_results[foldchange_columns]
      
      group_sizes <- attr(statRes_obj, "group_DF") %>% dplyr::group_by(Group) %>% dplyr::count() %>% purrr::pluck('n')
      ratios = t(t(counts)/group_sizes)
      
      ratio_diff_list = list()
      for(comp in attr(statRes_obj, "comparisons")){
        names = strsplit(comp, '_vs_')[[1]]
        compare_cols = paste0('Count_', names)
        diff = ratios[,compare_cols[1]] - ratios[,compare_cols[2]]
        ratio_diff_list[[length(ratio_diff_list) + 1]] = diff
      }
      ratio_diffs <- do.call(cbind, ratio_diff_list)
      #
      
      numsig <- attr(statRes_obj, "number_significant")
      
      #test num significant consistency and attributes  
      # for(i in 1:nrow(numsig)){
      #   expect_equal(sum(flags[,i] > 0, na.rm = TRUE), numsig[i,]$Up_total)
      #   expect_equal(sum(flags[,i] < 0, na.rm = TRUE), numsig[i,]$Down_total)
      #   expect_equal(sum(flags[,i] == -2, na.rm = TRUE), numsig[i,]$Down_gtest)
      #   expect_equal(sum(flags[,i] == -1, na.rm = TRUE), numsig[i,]$Down_anova)
      #   expect_equal(sum(flags[,i] == 1, na.rm = TRUE), numsig[i,]$Up_anova)
      #   expect_equal(sum(flags[,i] == 2, na.rm = TRUE), numsig[i,]$Up_gtest)
      # }
      
      # correct number of comparisons for number of groups
      expect_equal(length(attr(statRes_obj, "comparisons")), choose(length(unique(attr(statRes_obj, "group_DF")$Group)),2))
      expect_equal(ncol(flags), length(attr(statRes_obj, "comparisons")))
      
      #
      
      # various truth values to compare
      anova_is_signif <- anova_pvals < 0.05
      g_is_signif <- g_pvals < 0.05
      zero_flag <- flags == 0
      g_flag_pos <- flags == 2
      g_flag_neg <- flags == -2
      anova_flag_pos <- flags == 1
      anova_flag_neg <- flags == -1
      
      # depending on test method, flags should have certain values
      if(test == "gtest"){
        expect_true(all(sapply(flags, function(x){x %in% c(-2,0,2)})))
        expect_false(g_is_signif && zero_flag)
        expect_false((g_flag_pos && (ratio_diffs < 0)) | (g_flag_neg && (ratio_diffs > 0)))
      }
      else if(test == "anova"){
        expect_true(all(sapply(flags, function(x){x %in% c(-1,0,1)})))
        expect_false(anova_is_signif && zero_flag)
        expect_false((anova_flag_pos && (fold_changes < 0)) | (anova_flag_neg && (fold_changes > 0)))
      }
      else if(test == "comb"){
        expect_true(all((unlist(flags) %>% unique()) %in% c(-2,-1,0,1,2,NA)))
        expect_false(anova_is_signif && zero_flag)
        expect_false(g_is_signif && (zero_flag | anova_flag_pos | anova_flag_neg))
        expect_false((g_flag_pos && (ratio_diffs < 0)) | (g_flag_neg && (ratio_diffs > 0)))
        expect_false((anova_flag_pos && (fold_changes < 0)) | (anova_flag_neg && (fold_changes > 0)))
      }
    }
  }
})


#Test with really big dataset
# library(OvarianPepdataBPsubset)
# suppressWarnings(tcga_ovarian_pepdata_bp <- as.pepData(e_data = tcga_ovarian_pepdata_bp_subset$e_data, f_data = tcga_ovarian_pepdata_bp_subset$f_data, e_meta = tcga_ovarian_pepdata_bp_subset$e_meta,
#                                       edata_cname = 'Peptide', fdata_cname = "sampleID", emeta_cname = "Protein", data_scale = "log2"))
# suppressWarnings(tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("race")))
# expect_warning(expect_error(ovarian_res_dunnett <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'dunnett', test_method='combined')))
