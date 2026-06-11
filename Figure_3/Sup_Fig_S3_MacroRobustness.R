
## Subsampling robustness of macroecological abundance patterns
## underlying Fig. 3B

library(tidyverse)
library(pracma)
library(openxlsx)
library(ggplot2)
library(patchwork)

## 1. Basic settings
otu_file = "OTUtable.txt"
treatment_order = c("D0", "D1", "D3", "D5", "D6", "D7")

# n_iter = 100 ## for testing
n_iter = 1000

set.seed(123)

## 2. Read OTU table
otudf = read.table( otu_file, header = T, sep = "\t", row.names = 1, check.names = F )

sample_groups = list(D0 = grep("^D0", colnames(otudf), value = T),
                     D1 = grep("^D1", colnames(otudf), value = T),
                     D3 = grep("^D3", colnames(otudf), value = T),
                     D5 = grep("^D5", colnames(otudf), value = T),
                     D6 = grep("^D6", colnames(otudf), value = T),
                     D7 = grep("^D7", colnames(otudf), value = T) )
print(sapply(sample_groups, length))
target_n = min(sapply(sample_groups, length))

## 3. Function: create MacLaw-like long table
make_maclaw = function(otu_table, selected_samples) {
  
  sub_otu = otu_table[, selected_samples, drop = F]
  nreads = colSums(sub_otu)
  
  lawdf = sub_otu %>%
    as.data.frame() %>%
    rownames_to_column("otu_id") %>%
    pivot_longer( cols = -otu_id, names_to = "sample_id", values_to = "count" ) %>%
    filter(count != 0) %>%
    mutate( run_id = sample_id, project_id = substr(sample_id, 1, 2), classification = project_id )
  
  dfcount = data.frame( sample_id = names(nreads), nreads = as.numeric(nreads) )
  
  lawdf = merge(lawdf, dfcount, by = "sample_id", all.x = T)
  lawdf = lawdf[, c("otu_id", "count", "project_id", "sample_id","run_id", "nreads", "classification"  )]
  
  return(lawdf)
}

## 4. Functions from Figure 3B scripts
fun_erf = function(mu, c, m1, m2) {
  sigma = sqrt(-c * m1 + m2 + c * mu - m1 * mu)
  x = (c - mu) / sigma / sqrt(2)
  f = (mu - m1) * erfc(x) + exp(-x^2) * sqrt(2 / pi) * sigma
  return(f)
}

estimate_mean_f_safe = function(c, m1, m2) {
  out = tryCatch({
    mumin = (c * m1 - m2) / (c - m1)
    uniroot(
      fun_erf,
      c = c,
      m1 = m1,
      m2 = m2,
      interval = c(-50, mumin - 0.001),
      tol = 0.0001
    )$root
  }, error = function(e) {
    NA_real_
  })
  return(out)
}

## 5. Function: calculate macroecological parameters
calc_macro_params = function(lawdf) {
  
  summarydata = lawdf %>%
    group_by(project_id, classification) %>%
    summarise(n_of_runs = n_distinct(run_id),
              mean_nreads = mean(nreads),
              .groups = "drop" ) %>%
    filter(n_of_runs > 0) %>%
    mutate(idall = paste(project_id, classification))
  
  shortnames = data.frame( idall = summarydata$idall, sname = summarydata$project_id, scat = summarydata$project_id )
  
  proj = lawdf %>%
    group_by(project_id, classification, otu_id) %>%
    mutate(tf = mean(count / nreads),
           o = n(),
           tvpf = mean((count^2 - count) / nreads^2) ) %>%
    ungroup() %>%
    group_by(project_id, classification) %>%
    mutate(o = o / n_distinct(run_id)) %>%
    mutate( f = o * tf,  vf = o * tvpf ) %>%
    mutate(vf = vf - f^2) %>%
    ungroup() %>%
    mutate(idall = paste(project_id, classification)) %>%
    select(-tf, -tvpf) %>%
    left_join(shortnames, by = "idall")
  
  gamma_pars = proj %>%
    select(idall, sname, otu_id, o, f, vf) %>%
    mutate(cv = sqrt(vf / f^2)) %>%
    distinct() %>%
    mutate( beta = 1 / cv^2, theta = f / beta )
  
  ## ---------- AFD / Gamma shape ----------
  afd_raw = proj %>%
    filter(o > 0.9) %>%
    mutate(relative_abundance = count / nreads) %>%
    group_by(sname, otu_id) %>%
    mutate(otu_mean_relative_abundance = mean(relative_abundance, na.rm = T),
      normalized_fluctuation = relative_abundance / otu_mean_relative_abundance ) %>%
    ungroup() %>%
    filter(is.finite(normalized_fluctuation), normalized_fluctuation > 0 )
  
  afd_param = afd_raw %>%
    group_by(sname) %>%
    summarise(`AFD Gamma shape` = mean(normalized_fluctuation, na.rm = T)^2 / var(normalized_fluctuation, na.rm = T),
              `No. AFD observations` = n(),
              `No. AFD OTUs` = n_distinct(otu_id),
              .groups = "drop"  )
  
  ## ---------- Taylor's law slope ----------
  nbins = 20
  
  taylor_binned = gamma_pars %>%
    filter(f > 0, vf > 0) %>%
    group_by(idall, sname) %>%
    mutate(lf = log(f), dlf = (max(lf) - min(lf)) / nbins ) %>%
    mutate(b = as.integer((lf - min(lf)) / dlf)) %>%
    ungroup() %>%
    group_by(idall, sname, b) %>%
    summarise(vf = mean(vf, na.rm = T), f = mean(f, na.rm = T), .groups = "drop") %>%
    filter(f > 0, vf > 0) %>%
    mutate( log10_mean = log10(f),log10_variance = log10(vf)  )
  
  taylor_param = taylor_binned %>%
    group_by(sname) %>%
    group_modify(~ {
      fit = lm(log10_variance ~ log10_mean, data = .x)
      fit_sum = summary(fit)
      tibble(`Taylor slope` = coef(fit)[2],
             `Taylor intercept` = coef(fit)[1],
             `Taylor R2` = fit_sum$r.squared,
             `No. Taylor bins` = nrow(.x)    )  }) %>%
    ungroup()
  
  ## ---------- MAD / lognormal sigma ----------
  cutoffs = data.frame(
    sname = treatment_order,
    c = rep(-18, length(treatment_order))
  )
  
  mean_pars = gamma_pars %>%
    mutate(lf = log(f)) %>%
    left_join(cutoffs, by = "sname") %>%
    filter(lf > c) %>%
    group_by(idall, sname) %>%
    summarise(c = mean(c),
              m1 = mean(lf),
              m2 = mean(lf^2),
              ns_obs = n_distinct(otu_id),
              nf = sum(f),
              .groups = "drop"    ) %>%
    rowwise() %>%
    mutate(mu = estimate_mean_f_safe(c, m1, m2)) %>%
    mutate(sigma = sqrt(-c * m1 + m2 + c * mu - m1 * mu)) %>%
    ungroup() %>%
    mutate(stot = 2 * ns_obs / erfc((c - mu) / sigma / sqrt(2)) )
  
  mad_param = mean_pars %>%
    transmute(sname,
              `MAD meanlog` = mu,
              `MAD sdlog` = sigma,
              `MAD estimated species total` = stot,
              `No. MAD OTUs` = ns_obs )
  
  out = afd_param %>%
    full_join(taylor_param, by = "sname") %>%
    full_join(mad_param, by = "sname") %>%
    mutate(sname = factor(sname, levels = treatment_order) ) %>%
    arrange(sname)
  
  return(out)
}

## 6. Observed full-data parameters

all_samples = unlist(sample_groups)

law_observed = make_maclaw( otu_table = otudf, selected_samples = all_samples )

observed_params = calc_macro_params(law_observed) %>%
  mutate(Analysis = "Observed full data")

## 7. Equalized replicate-number subsampling

subsampled_list = vector("list", n_iter)

for (i in seq_len(n_iter)) {
  
  if (i %% 10 == 0 || i == n_iter) {
    message("Running iteration ", i, " / ", n_iter)
  }
  
  sampled_samples = unlist(
    lapply(sample_groups, function(x) {
      sample(x, size = target_n, replace = F)
    })
  )
  
  law_i = make_maclaw( otu_table = otudf, selected_samples = sampled_samples)
  
  pars_i = calc_macro_params(law_i) %>%
    mutate( Analysis = "Equalized replicate-number subsampling",
            Iteration = i, `No. replicates per treatment` = target_n )
  
  subsampled_list[[i]] = pars_i
}

subsampled_params = bind_rows(subsampled_list)

## 8. Summarise subsampling results
metrics_to_plot = c( "AFD Gamma shape", "Taylor slope", "MAD sdlog" )

sub_long = subsampled_params %>%
  select(sname, Iteration, all_of(metrics_to_plot)) %>%
  pivot_longer(cols = all_of(metrics_to_plot),
               names_to = "Metric",
               values_to = "Value"   )

obs_long = observed_params %>%
  select(sname, all_of(metrics_to_plot)) %>%
  pivot_longer(cols = all_of(metrics_to_plot),
               names_to = "Metric",
               values_to = "Observed value"   )

summary_sub = sub_long %>%
  group_by(sname, Metric) %>%
  summarise(Mean = mean(Value, na.rm = T),
            SD = sd(Value, na.rm = T),
            Median = median(Value, na.rm = T),
            `2.5% quantile` = quantile(Value, 0.025, na.rm = T),
            `97.5% quantile` = quantile(Value, 0.975, na.rm = T),
            .groups = "drop" ) %>%
  left_join(obs_long, by = c("sname", "Metric")) %>%
  mutate(`Observed within 95% interval` = `Observed value` >= `2.5% quantile` & `Observed value` <= `97.5% quantile`,
    sname = factor(sname, levels = treatment_order),
    Metric = factor(Metric, levels = metrics_to_plot, 
                    labels = c( "AFD: Gamma shape", "Taylor's law: slope", "MAD: lognormal sdlog" ) 
                    )
  )


## 9. Plot Supplementary Fig. S6

pS6 = ggplot(summary_sub, aes(x = sname)) +
  geom_linerange(aes( ymin = `2.5% quantile`, ymax = `97.5% quantile`), color = "grey45", linewidth = 0.55) +
  geom_point( aes(y = Mean), shape = 21, size = 3.0, 
              fill = "#4C78A8", alpha = 0.5,
              color = "black", stroke = 1) +
  geom_point( aes(y = `Observed value`), shape = 2, size = 3.0, color = "orange", stroke = 1) +
  facet_wrap( ~ Metric, scales = "free_y", nrow = 1  ) +
  labs( x = "Dilution treatment", y = "Macroecological parameter value"  ) +
  # theme_classic(base_size = 10) +
  theme_bw() + 
  theme(plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 10),
        # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.margin = margin(3, 3, 3, 3, unit = "mm")   )

pS6

ggsave( "Sup_Fig_S3_MacroRobustness.pdf", pS6, width = 7.2, height = 2.8 )


