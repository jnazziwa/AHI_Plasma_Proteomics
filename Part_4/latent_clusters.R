evasVL <- readRDS(file = "evasVL_data.rds")
age_data <- evasVL %>% 
  dplyr::select(pid, Age, Cohort, Sex) %>% 
  dplyr::distinct(pid, .keep_all = TRUE)

CD4_cox <- readxl::read_excel("CD4_groups.xlsx")

#merge
variates <- merge(age_data, CD4_cox,
                  by.x = "pid", by.y = "ID", all.x = TRUE) %>% 
  dplyr::mutate(status = ifelse(Group == "Fast", "1", "0"))

#READ hier clusters

dynamics <- readr::read_csv("cluster6.csv")

variates$pid <- sprintf("%03d", variates$pid)

merged_dyn <-  merge(dynamics, variates,
                     by.x = "subjid", by.y = "pid", all.x = TRUE)
merged_dyn$status <- as.numeric(merged_dyn$status)
merged_dyn$cluster6 <- as.character(merged_dyn$cluster6)


library(coxme)

# Example model: Cox proportional hazards with random effect for patient_id
# time: survival time, status: event/censor, cluster: latent cluster, 
# patient_id as random effect
cox_model <- coxme(Surv(time, status) ~ cluster6 + Age + Sex + Cohort + (1 | subjid), data = merged_dyn)

summary(cox_model)

# Extract summary of the model
model_summary <- summary(cox_model)

coefficients <- summary(cox_model)$coefficients
hazard_ratios <- exp(coefficients[, 1]) 

results <- data.frame(
  Variable = rownames(coefficients),
  Coefficient = coefficients[, 1],
  `Standard Error` = coefficients[, 3],
  `z-value` = coefficients[, 4],
  `p-value` = coefficients[, 5],
  `Hazard Ratio (HR)` = hazard_ratios
)

# View the extracted fixed effects
print(results)

