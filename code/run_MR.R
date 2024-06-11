library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

mr_results_list <- list()
result_list <- list()
r<-data.frame()

for (ef in exposure_files) {
  # Read exposure data
  exp_data <- read_exposure_data(ef, sep = "\t")

  exp_data <- clump_data(exp_data, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1)
  exp_data<-exp_data %>% filter(SNP %in% clump_data$rsid)

  for (of in outcome_files) {
    # Read outcome data
    out_data <- read_outcome_data(of, sep = "\t")
  
    h <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)
    h <- h %>%
      mutate(
        beta.exposure = if_else(effect_allele.exposure != effect_allele.outcome, -beta.exposure, beta.exposure),
        effect_allele.exposure = if_else(effect_allele.exposure != effect_allele.outcome, other_allele.exposure, effect_allele.exposure),
        other_allele.exposure = if_else(effect_allele.exposure != effect_allele.outcome, effect_allele.exposure, other_allele.exposure)
      )
    
    dat <- mr(h, parameters = default_parameters(), method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

    exposure <- dat$exposure
    outcome <- dat$outcome
    method <- dat$method
    nsnp <- dat$nsnp
    b <- dat$b
    se <- dat$se
    pval <- dat$pval
    result_df <- data.frame(exposure, outcome, method, nsnp, b, se, pval)
    r<-rbind(r,result_df)
    result_list[[basename_prefix]] <- result_df

    write.table(
      x = result_df,
      file = output_file,
      sep = ",", 
      col.names = FALSE, 
      row.names = FALSE, 
      append = TRUE 
    )
  }
}