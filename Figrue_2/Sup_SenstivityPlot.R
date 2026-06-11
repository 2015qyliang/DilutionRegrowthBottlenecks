
## Sensitivity analysis of OTU intersections underlying Fig. 2E

library(tidyverse)
library(readxl)
library(ggplot2)
library(patchwork)
library(scales)

## 1. Read Table

infile = "Sup_Sensitivity analysis of OTU intersections.xlsx"

eq = read_excel( infile,   sheet = "Equalized_replicate_summary")
occ = read_excel( infile,   sheet = "Occupancy_threshold_summary" )

## 2. Clean column names

clean_names_simple = function(x) {
  x %>%
    str_replace_all("%", "percent") %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

names(eq) = clean_names_simple(names(eq))
names(occ) = clean_names_simple(names(occ))

## Expected columns after cleaning:
## Equalized_replicate_summary:
## analysis, no_replicates_per_treatment, metric, mean, sd, median,
## 2_5percent_quantile, 97_5percent_quantile, observed_value,
## observed_within_95percent_interval

## Occupancy_threshold_summary:
## metric, value, analysis, threshold_type, observed_value,
## relative_to_original

## 3. Select key metrics shown or discussed in Fig. 2E

selected_metrics = c(
  "Total OTUs in D1",
  "Intersection: D1 only",
  "Intersection: D0&D1",
  "Intersection: D0&D1&D3&D5&D6&D7",
  "Intersection: D1&D3&D5&D6&D7",
  "Intersection: D3&D5&D6&D7",
  "Total OTUs in D6",
  "Total OTUs in D7" )

metric_labels = c(
  "Total OTUs in D1" = "Total OTUs in D1",
  "Intersection: D1 only" = "D1-specific",
  "Intersection: D0&D1" = "D0 ∩ D1",
  "Intersection: D0&D1&D3&D5&D6&D7" = "All-treatment shared",
  "Intersection: D1&D3&D5&D6&D7" = "D1 ∩ D3 ∩ D5 ∩ D6 ∩ D7",
  "Intersection: D3&D5&D6&D7" = "D3 ∩ D5 ∩ D6 ∩ D7",
  "Total OTUs in D6" = "Total OTUs in D6",
  "Total OTUs in D7" = "Total OTUs in D7" )

selected_metrics = intersect(selected_metrics, unique(eq$metric))

eq_plot = eq %>%
  filter(metric %in% selected_metrics) %>%
  mutate(
    Metric_label = recode(metric, !!!metric_labels),
    Metric_label = factor(
      Metric_label,
      levels = rev(metric_labels[selected_metrics])
    )
  )

occ_plot = occ %>%
  filter(metric %in% selected_metrics) %>%
  mutate(
    Metric_label = recode(metric, !!!metric_labels),
    Metric_label = factor(
      Metric_label,
      levels = rev(metric_labels[selected_metrics])
    ),
    threshold_type = factor(
      threshold_type,
      levels = c( ">=1 replicate", ">=2 replicates", ">=5% replicates", ">=10% replicates" ),
      labels = c( "\u22651 rep.", "\u22652 reps.", "\u22655% reps.", "\u226510% reps." )
    )
  )

## 4. Panel A: equalized replicate-number sensitivity

pA = ggplot(eq_plot, aes(y = Metric_label)) +
  
  ## 95% interval from random subsampling
  geom_segment(aes(x = `2_5percent_quantile`, xend = `97_5percent_quantile`,yend = Metric_label ),
               linewidth = 0.55,  color = "grey40" ) +
  
  ## mean after equalized subsampling
  geom_point(aes(x = mean), shape = 21, size = 3.0, 
             fill = "#4C78A8", alpha = 0.5,
             color = "black", stroke = 1  ) +
  
  ## observed value from original Fig. 2E
  geom_point(aes(x = observed_value), shape = 2, size = 3.0, color = "orange", stroke = 1  ) +
  
  labs( x = "OTU count", y = NULL, title = "i. Equalized replicate-number sensitivity"  ) +
  
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(), 
        # panel.grid = element_blank(), 
        plot.title = element_text( face = "bold", size = 10, color = "black" ),
        axis.text.x = element_text( color = "black", size = 8 ),
        axis.text.y = element_text( color = "black", size = 10 ),
        axis.title.x = element_text( color = "black", size = 10 ),
        plot.margin = margin(3, 3, 3, 3, unit = "mm")   )

## 5. Panel B: occupancy-threshold sensitivity

pB = ggplot(occ_plot, aes(x = threshold_type,y = Metric_label,fill = relative_to_original )) +
  geom_tile( color = "white", linewidth = 0.35  ) +
  geom_text( aes(label = value), size = 2.5, color = "black"  ) +
  scale_fill_gradient2( low = "#6BAED6", mid = "#F7F7F7", high = "#CB6A4A", 
                        midpoint = 1, name = "Relative to\noriginal"  ) +

  labs( x = "Treatment-level presence threshold", y = NULL, title = "ii. Occupancy-threshold sensitivity"  ) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(), 
        plot.title = element_text( face = "bold", size = 10, color = "black" ),
        axis.text.x = element_text( color = "black", size = 8, angle = 35, hjust = 1 ),
        axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.title.x = element_text( color = "black", size = 10 ),
        legend.title = element_text( size = 8 ),
        legend.text = element_text(size = 7),
        plot.margin = margin(3, 3, 3, 3, unit = "mm")   )

## 6. Combine panels

figS3 = pA + pB +  plot_layout( widths = c(1.35, 1)  )

ggsave( "FigS3_OTU_intersection_sensitivity.pdf", figS3, width = 8, height = 3 )




