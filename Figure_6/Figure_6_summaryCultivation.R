
# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(ggplot2)
library(ggsci)
library(doBy)
library(tidyverse)
  
####################################################################################################
# Figure_6A

culstrains = data.frame(Group = c("D0", "D1", "D3", "D5", "D6", "D7"), 
                        Strains = c(313, 757, 737, 491, 503, 319))

pA = ggplot(culstrains, aes(x = Group, y = Strains, fill = Group)) + 
  geom_col(width = 0.6, alpha = 0.6) +
  scale_fill_manual(values = pal_nejm()(6), 
                    breaks = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')) + 
  labs(title = "The distribution\nof isolated\nbacteria strains", 
       y = "Number of isolated strains", x = NULL) + 
  theme_bw() + 
  theme(panel.background = element_blank(),
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = 'none',
        plot.title = element_text(size = 10, colour = 'black'),
        axis.title = element_text(size = 8, colour = 'black'),
        axis.text = element_text(size = 8, colour = 'black') ) 

ggsave('Figure_6A.pdf', pA, width = 1.8, height = 2.2)

####################################################################################################
# Figure_6C 

allstrains = read.table('AllIsolatedStrainsTaxon.txt', header = T, sep = '\t' )
allstrains = unique(allstrains, by = "Strains ")

table(allstrains$Group)

newdf = data.frame(Group = allstrains$Group, Count = 1, 
                   Family = gsub(';.*', '', gsub('.*ales;', '', allstrains$Taxon)))
famDF = summaryBy(Count ~ Family, data = newdf, FUN = sum)
famDF = famDF[order(famDF$Count.sum, decreasing = T), ]

newdfSub = newdf[newdf$Family %in% famDF$Family[1:15], ]
dat = summaryBy(Count ~ Family + Group, data = newdfSub, FUN = sum)
# write.csv(dat, 'Figure_6C.csv', row.names = F)

treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")
dat = dat %>%
  mutate(Group = factor(Group, levels = treatment_order),
         Family = as.character(Family),
         Count.sum = as.numeric(Count.sum) )

## ---------------------------------------------------------
## 2. Complete missing combinations

plot_dat = dat %>%
  complete(Family, Group = factor(treatment_order, levels = treatment_order),
           fill = list(Count.sum = 0)  )



## ---------------------------------------------------------
## 3. Candidate families from Figure 5

candidate_families = c("Prolixibacteraceae", "Marinifilaceae", "Marinilabiliaceae",
                       "Vibrionaceae", "Shewanellaceae", "Exiguobacteriaceae",
                       "Fusobacteriaceae", "Dethiosulfatibacteraceae", "Anaerolineaceae" )

plot_dat = plot_dat %>%
  group_by(Family) %>%
  mutate(Total_isolates = sum(Count.sum, na.rm = TRUE)  ) %>%
  ungroup() %>%
  mutate(Candidate_responder = Family %in% candidate_families,
         Family_label = ifelse(Candidate_responder,
                               paste0(Family, " *"),
                               Family)  )

write.csv(plot_dat, 'Figure_6C.csv', row.names = F)

## ---------------------------------------------------------
## 4. Keep top families for main figure

top_n = 15

top_families = plot_dat %>%
  distinct(Family, Total_isolates) %>%
  slice_max(Total_isolates, n = top_n) %>%
  pull(Family)

keep_families = union(top_families, candidate_families)

plot_dat_main = plot_dat %>%
  filter(Family %in% keep_families)

## Reorder by total isolate number
family_order_tbl = plot_dat_main %>%
  distinct(Family, Family_label, Total_isolates) %>%
  arrange(Total_isolates)

family_order = family_order_tbl$Family
family_labels = setNames(family_order_tbl$Family_label, family_order_tbl$Family)

plot_dat_main = plot_dat_main %>%
  mutate(Family = factor(Family, levels = family_order)  )

## ---------------------------------------------------------
## 5. Plot settings

single_fill = "#4C78A8"   # muted blue
label_cutoff = 20         # only label larger bubbles

pC = ggplot(plot_dat_main, aes(x = Group, y = Family)) +
  geom_point( data = plot_dat_main %>% filter(Count.sum > 0), aes(size = Count.sum), shape = 16,  
              color = single_fill,  alpha = 0.6 ) +
  geom_text(data = plot_dat_main %>% filter(Count.sum >= 50 ),
            aes(label = Count.sum), size = 2.8, color = "black"  ) +
  
  scale_y_discrete(labels = family_labels) +
  scale_size_area( max_size = 11, breaks = c(1, 10, 50, 100, 300), name = "No. of isolates" ) +
  labs(x = NULL, y = NULL, title = NULL ) +
  theme_bw() +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(3, 3, 3, 3, unit = "mm") ,
        panel.grid = element_blank(), 
        legend.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8, color  = 'black'),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 10, color = 'black', face = 'italic')) +
  guides(size = guide_legend(override.aes = list(fill = "white", color = "black", alpha = 1  )  )  )

ggsave('Figure_6C.pdf', pC, width = 4.6, height = 4.3)

####################################################################################################
# Figure_6D
# 







####################################################################################################



####################################################################################################  



