############################################################
## Distributional-state test results of OTUs across dilution treatments
## Based on Figure3A_IndexOfDispersion.R

library(tidyverse)
library(openxlsx)

## 1. If IoDs already exists from Figure 3A script, use it directly
##    Otherwise, read from IndexOfDispersion.txt

if (!exists("IoDs")) {
  IoDs = read.table(
    "IndexOfDispersion.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
}

treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")

## 2. Add derived variables used in Fig. 3A

IoDs_table = IoDs %>%
  mutate(
    Group = factor(Group, levels = treatment_order), 
    ## Raw index of dispersion used for plotting 
    ## The plotted value in Fig. 3A is log10((variance / mean) * occurrence)
    `Index of dispersion` = (fV / fM) * fN,
    
    `log10 index of dispersion` = ifelse(`Index of dispersion` > 0,
                                         log10(`Index of dispersion`),
                                         NA_real_    ),
    
    ## Occurrence and occupancy
    `Occurrence count` = fN,
    
    ## No. replicate samples per treatment
    `No. of replicates` = case_when(Group == "D0" ~ 14,
                                    Group == "D1" ~ 14,
                                    Group == "D3" ~ 20,
                                    Group == "D5" ~ 20,
                                    Group == "D6" ~ 40,
                                    Group == "D7" ~ 40,
                                    TRUE ~ NA_real_    ),
    
    `Treatment occupancy` = `Occurrence count` / `No. of replicates`,
    
    ## Classification used in Fig. 3A
    `Uniformity result` = ifelse(rChisqP <= 0.05, "Non-uniform", "Uniform"    ),
    
    `Normality result` = ifelse(rKsP <= 0.05, "Non-normal", "Normal"    ) )

## add FDR-adjusted q values for transparency
## Fig. 3A classification remains based on nominal P <= 0.05.
IoDs_table = IoDs_table %>%
  group_by(Group) %>%
  mutate(`Chi-square q` = p.adjust(rChisqP, method = "BH"), `KS q` = p.adjust(rKsP, method = "BH")) %>%
  ungroup()


## 3. Sheet 1: Summary by treatment

summary_by_treatment = IoDs_table %>%
  group_by(Group) %>%
  summarise(`No. of replicates` = first(`No. of replicates`),
            `No. OTUs tested` = n(), 
            
            `No. uniform OTUs` = sum(`Uniformity result` == "Uniform", na.rm = TRUE),
            `No. non-uniform OTUs` = sum(`Uniformity result` == "Non-uniform", na.rm = TRUE),
            `% non-uniform OTUs` = `No. non-uniform OTUs` / `No. OTUs tested` * 100,
            
            `No. normal OTUs` = sum(`Normality result` == "Normal", na.rm = TRUE),
            `No. non-normal OTUs` = sum(`Normality result` == "Non-normal", na.rm = TRUE),
            `% non-normal OTUs` = `No. non-normal OTUs` / `No. OTUs tested` * 100,
            
            `Median occurrence count` = median(`Occurrence count`, na.rm = TRUE),
            `Median treatment occupancy` = median(`Treatment occupancy`, na.rm = TRUE),
            `Median index of dispersion` = median(`Index of dispersion`, na.rm = TRUE),
            `Median log10 index of dispersion` = median(`log10 index of dispersion`, na.rm = TRUE),
            
            .groups = "drop"  ) %>%
  arrange(match(Group, treatment_order))


## 4. Sheet 2: OTU-level tests

otu_level_tests = IoDs_table %>%
  transmute(`Dilution treatment` = as.character(Group),
            `OTU ID` = OTU,
            `No. of replicates`,
            `Occurrence count`,
            `Treatment occupancy`,
            
            `Mean abundance` = fM,
            `Abundance variance` = fV,
            `Index of dispersion`,
            `log10 index of dispersion`,
            
            `Chi-square statistic` = rChisqX,
            `Chi-square P` = rChisqP,
            `Chi-square q` = `Chi-square q`,
            `Uniformity result`,
            
            `KS statistic` = rKsD,
            `KS P` = rKsP,
            `KS q` = `KS q`,
            `Normality result`  ) %>%
  arrange(factor(`Dilution treatment`, levels = treatment_order), `OTU ID`  )


## 5. Optional sheet: occurrence-bin summary

occurrence_bin_summary = IoDs_table %>%
  mutate(`Occurrence bin` = case_when(`Treatment occupancy` == 1 ~ "Persistent",
                                      `Treatment occupancy` >= 0.75 & `Treatment occupancy` < 1 ~ "High occupancy", 
                                      `Treatment occupancy` >= 0.25 & `Treatment occupancy` < 0.75 ~ "Intermediate occupancy",
                                      `Treatment occupancy` > 0 & `Treatment occupancy` < 0.25 ~ "Transient",
                                      TRUE ~ "Not detected" )  ) %>%
  group_by(Group, `Occurrence bin`) %>%
  summarise(`No. OTUs` = n(),
            `No. non-uniform OTUs` = sum(`Uniformity result` == "Non-uniform", na.rm = TRUE),
            `% non-uniform OTUs` = `No. non-uniform OTUs` / `No. OTUs` * 100,
            `No. non-normal OTUs` = sum(`Normality result` == "Non-normal", na.rm = TRUE),
            `% non-normal OTUs` = `No. non-normal OTUs` / `No. OTUs` * 100,
            `Median log10 index of dispersion` = median(`log10 index of dispersion`, na.rm = TRUE),
            .groups = "drop"  ) %>%
  arrange(factor(Group, levels = treatment_order), `Occurrence bin`  )

##   Export Excel workbook

wb = createWorkbook()

addWorksheet(wb, "1.Summary_by_treatment")
writeData(wb, "1.Summary_by_treatment", summary_by_treatment)

addWorksheet(wb, "2.OTU_level_tests")
writeData(wb, "2.OTU_level_tests", otu_level_tests)

addWorksheet(wb, "3.Occurrence_bin_summary")
writeData(wb, "3.Occurrence_bin_summary", occurrence_bin_summary)

saveWorkbook(
  wb,
  "Table_S8_distributional_state_tests_underlying_Fig3A.xlsx",
  overwrite = TRUE
)
