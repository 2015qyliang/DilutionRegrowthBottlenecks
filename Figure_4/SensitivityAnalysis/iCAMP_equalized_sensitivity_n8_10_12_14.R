# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##########################################################################
# iCAMP replicate-equalized sensitivity analysis
# Common-n settings: n = 8, 10, 12, 14 for all dilution treatments
# Subsampling iterations: 100

##########################################################################

library(NST)
library(iCAMP)
library(ape)
library(tidyverse)

rm(list = ls())

##########################################################################
# 1. Global settings

set.seed(20260604)

# iCAMP null randomization/bootstrap iterations.
# Keep this identical to the main iCAMP script unless there is a strong reason to change.
rand.time = 100

# Number of random subsampling iterations for each common-n setting.
subsample.time = 100

# Common-n grid. All dilution treatments must be represented at each n.
n.grid = c(8, 10, 12, 14)
# n.grid = c(10)

nworker = 8
memory.G = 15

# If TRUE, icamp.bins() is also run for each subsampling iteration.
# For treatment-level sensitivity, icamp.boot() is usually sufficient.
# If process extraction from icamp.boot() fails in your iCAMP version, set this to TRUE.
run.bin.summary = TRUE

# For speed, bin-level bootstrap is not necessary for this sensitivity.
# The uncertainty comes from repeated subsampling iterations.
bin.boot = FALSE

# Resume finished iterations if their minimal process csv already exists.
resume = TRUE

groups.order = c('D0', 'D1', 'D3', 'D5', 'D6', 'D7')
diluted.groups = c('D1', 'D3', 'D5', 'D6', 'D7')
strong.groups = c('D3', 'D5', 'D6', 'D7')
process.levels = c('HoS', 'HeS', 'DL', 'HD', 'DR')

out.wd = "iCAMP_equalized"
# if(!dir.exists(out.wd)){dir.create(out.wd, recursive = TRUE)}

##########################################################################
# 2. Data input

comm = t(read.table("../OTUtable.txt", header = TRUE, sep = "\t",
                    row.names = 1, as.is = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE))

treat = read.table("../OTUsamplesGroup.txt", header = TRUE, sep = "\t",
                   row.names = 1, stringsAsFactors = FALSE,
                   check.names = FALSE)

tree = read.tree('../OTUFastTree.txt')

taxon = read.table('../OTUTax.txt', header = TRUE, sep = "\t",
                   row.names = 1, stringsAsFactors = FALSE,
                   check.names = FALSE)

# Align sample order between community table and treatment metadata
common.samples = intersect(rownames(comm), rownames(treat))
comm = comm[common.samples, , drop = FALSE]
treat = treat[common.samples, , drop = FALSE]

# Identify grouping column
if("Group" %in% colnames(treat)){
  group.col = "Group"
} else {
  group.col = colnames(treat)[1]
}
treat[[group.col]] = factor(as.character(treat[[group.col]]), levels = groups.order)

# Basic checks
group.size = table(treat[[group.col]])
print(group.size)

if(any(is.na(treat[[group.col]]))){
  stop("Some samples have NA group labels. Please check 02-samplesGroup.txt.")
}

if(any(group.size[groups.order] < max(n.grid))){
  stop("At least one dilution treatment has fewer replicates than max(n.grid).")
}

if(any(rowSums(comm) == 0)){
  warning("Some samples have zero total abundance and will be removed.")
  keep.sample = rowSums(comm) > 0
  comm = comm[keep.sample, , drop = FALSE]
  treat = treat[rownames(comm), , drop = FALSE]
}

# Remove globally zero OTUs
comm = comm[, colSums(comm) > 0, drop = FALSE]

##########################################################################
# 3. Pairwise phylogenetic distance matrix
# Same logic as 01-iCAMPdiy(1).R

pd.wd = "PD"
if(!dir.exists(pd.wd)){dir.create(pd.wd)}

if(!file.exists(paste0(pd.wd, "/pd.desc"))) {
  pd.big = pdist.big(tree = tree, wd = pd.wd,
                     nworker = nworker, memory.G = memory.G)
} else {
  pd.big = list()
  pd.big$tip.label = read.csv(paste0(pd.wd, "/pd.taxon.name.csv"),
                              row.names = 1, stringsAsFactors = FALSE)[, 1]
  pd.big$pd.wd = pd.wd
  pd.big$pd.file = 'pd.desc'
  pd.big$pd.name.file = 'pd.taxon.name.csv'
}

##########################################################################
# 4. Helper functions

normalize_process_name = function(x) {
  x0 = tolower(gsub("[^a-zA-Z]", "", as.character(x)))

  out = rep(NA_character_, length(x0))
  out[grepl("homogeneousselection|^hos$", x0)] = "HoS"
  out[grepl("heterogeneousselection|^hes$", x0)] = "HeS"
  out[grepl("dispersallimitation|^dl$", x0)] = "DL"
  out[grepl("homogenizingdispersal|^hd$", x0)] = "HD"
  out[grepl("drift|^dr$", x0)] = "DR"

  return(out)
}

pick_group_col = function(df, groups.order) {
  group.candidates = c("Group", "group", "Treatment", "treatment",
                       "Treat", "treat", "dilution", "Dilution")
  cand = intersect(group.candidates, colnames(df))
  if(length(cand) > 0){return(cand[1])}

  overlap = sapply(colnames(df), function(nm){
    sum(as.character(df[[nm]]) %in% groups.order, na.rm = TRUE)
  })

  if(length(overlap) > 0 && max(overlap) > 0){
    return(names(which.max(overlap)))
  } else {
    return(NA_character_)
  }
}

pick_process_col = function(df, process.levels) {
  overlap = sapply(colnames(df), function(nm){
    sum(normalize_process_name(df[[nm]]) %in% process.levels, na.rm = TRUE)
  })

  if(length(overlap) > 0 && max(overlap) > 0){
    return(names(which.max(overlap)))
  } else {
    return(NA_character_)
  }
}

pick_value_col = function(df, exclude.cols = character()) {
  value.candidates = c("Contribution", "contribution",
                       "Percentage", "percentage",
                       "Percent", "percent",
                       "Mean", "mean",
                       "Estimate", "estimate",
                       "Value", "value",
                       "importance", "Importance")

  cand = setdiff(intersect(value.candidates, colnames(df)), exclude.cols)
  if(length(cand) > 0){
    for(cc in cand){
      suppressWarnings(vv <- as.numeric(as.character(df[[cc]])))
      if(sum(is.finite(vv)) > 0){return(cc)}
    }
  }

  numeric.cols = setdiff(colnames(df)[sapply(df, is.numeric)], exclude.cols)
  if(length(numeric.cols) > 0){return(numeric.cols[1])}

  return(NA_character_)
}

extract_process_table = function(x, groups.order, process.levels, source.name) {

  if(is.null(x)){return(NULL)}

  df = as.data.frame(x, stringsAsFactors = FALSE)
  if(nrow(df) == 0){return(NULL)}

  df$RowName_iCAMP = rownames(df)

  # Case 1: wide format with process columns
  norm.cols = normalize_process_name(colnames(df))
  p.idx = which(norm.cols %in% process.levels)

  if(length(p.idx) >= 3){
    gcol = pick_group_col(df, groups.order)
    if(is.na(gcol)){
      df$Group = df$RowName_iCAMP
      gcol = "Group"
    }

    long = map_dfr(p.idx, function(j){
      data.frame(Group = as.character(df[[gcol]]),
                 Process = norm.cols[j],
                 Contribution = as.numeric(as.character(df[[j]])),
                 Source = source.name,
                 stringsAsFactors = FALSE)
    })

    long = long %>%
      filter(Group %in% groups.order, Process %in% process.levels,
             is.finite(Contribution))

    if(nrow(long) > 0){return(long)}
  }

  # Case 2: rows are processes and columns are groups
  row.proc = normalize_process_name(df$RowName_iCAMP)
  gcols = intersect(groups.order, colnames(df))

  if(sum(row.proc %in% process.levels, na.rm = TRUE) >= 3 && length(gcols) >= 2){
    long = map_dfr(seq_len(nrow(df)), function(i){
      pp = row.proc[i]
      if(is.na(pp) || !(pp %in% process.levels)){return(NULL)}
      map_dfr(gcols, function(gg){
        data.frame(Group = gg,
                   Process = pp,
                   Contribution = as.numeric(as.character(df[i, gg])),
                   Source = source.name,
                   stringsAsFactors = FALSE)
      })
    })

    long = long %>% filter(is.finite(Contribution))
    if(nrow(long) > 0){return(long)}
  }

  # Case 3: long format with group, process and value columns
  gcol = pick_group_col(df, groups.order)
  pcol = pick_process_col(df, process.levels)
  vcol = pick_value_col(df, exclude.cols = c(gcol, pcol, "RowName_iCAMP"))

  if(!is.na(gcol) && !is.na(pcol) && !is.na(vcol)){
    long = data.frame(Group = as.character(df[[gcol]]),
                      Process = normalize_process_name(df[[pcol]]),
                      Contribution = as.numeric(as.character(df[[vcol]])),
                      Source = source.name,
                      stringsAsFactors = FALSE) %>%
      filter(Group %in% groups.order, Process %in% process.levels,
             is.finite(Contribution))

    if(nrow(long) > 0){return(long)}
  }

  return(NULL)
}

add_metadata_to_process = function(proc, n.sub, iter, iter.id) {
  proc %>%
    mutate(N = n.sub,
           Iteration = iter,
           IterationID = iter.id,
           Group = factor(Group, levels = groups.order),
           Process = factor(Process, levels = process.levels)) %>%
    select(N, Iteration, IterationID, Group, Process, Contribution, Source)
}

run_one_icamp = function(n.sub, iter) {

  iter.id = paste0("n", n.sub, "_iter", sprintf("%03d", iter))
  iter.wd = file.path(out.wd, paste0("n", n.sub), iter.id)
  if(!dir.exists(iter.wd)){dir.create(iter.wd, recursive = TRUE)}

  process.file = file.path(iter.wd, paste0(iter.id, "_ProcessLong.csv"))
  error.file = file.path(iter.wd, paste0(iter.id, "_ERROR.txt"))

  if(resume && file.exists(process.file)){
    message("Resume: ", iter.id)
    return(read.csv(process.file, stringsAsFactors = FALSE))
  }

  sampled.ids = unlist(lapply(groups.order, function(gg){
    ids = rownames(treat)[treat[[group.col]] == gg]
    sample(ids, n.sub, replace = FALSE)
  }), use.names = FALSE)

  comm.sub = comm[sampled.ids, , drop = FALSE]
  treat.sub = treat[sampled.ids, , drop = FALSE]
  treat.sub[[group.col]] = factor(as.character(treat.sub[[group.col]]),
                                  levels = groups.order)

  # Remove OTUs absent from this subsampling set
  comm.sub = comm.sub[, colSums(comm.sub) > 0, drop = FALSE]

  sample.info = data.frame(SampleID = rownames(comm.sub),
                           Group = as.character(treat.sub[[group.col]]),
                           N = n.sub,
                           Iteration = iter,
                           IterationID = iter.id,
                           stringsAsFactors = FALSE)

  write.csv(sample.info,
            file = file.path(iter.wd, paste0(iter.id, "_SampleSet.csv")),
            row.names = FALSE)

  # Main iCAMP engine: identical key parameters to 01-iCAMPdiy(1).R
  icres = icamp.big(comm = comm.sub,
                    pd.desc = pd.big$pd.file,
                    pd.spname = pd.big$tip.label,
                    pd.wd = pd.big$pd.wd,
                    rand = rand.time,
                    tree = tree,
                    prefix = iter.id,
                    output.wd = iter.wd,
                    nworker = nworker,
                    memory.G = memory.G,
                    unit.sum = rowSums(comm.sub),
                    ds = 0.2, pd.cut = NA,
                    sp.check = TRUE,
                    phylo.rand.scale = "within.bin",
                    taxa.rand.scale = "across.all",
                    phylo.metric = "bMPD",
                    sig.index = "Confidence",
                    bin.size.limit = 12,
                    rtree.save = FALSE,
                    detail.save = TRUE,
                    qp.save = FALSE,
                    detail.null = FALSE,
                    ignore.zero = TRUE,
                    correct.special = TRUE,
                    special.method = "depend",
                    ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975,
                    omit.option = "no", meta.ab = NULL)

  saveRDS(icres, file = file.path(iter.wd, paste0(iter.id, "_icres.rds")))

  icamp.result = icres$CbMPDiCBraya

  icboot = icamp.boot(icamp.result = icamp.result,
                      treat = treat.sub,
                      rand.time = rand.time,
                      compare = TRUE,
                      silent = TRUE,
                      between.group = TRUE,
                      ST.estimation = TRUE)

  saveRDS(icboot, file = file.path(iter.wd, paste0(iter.id, "_icboot.rds")))
  write.csv(icboot$summary,
            file = file.path(iter.wd, paste0(iter.id, "_BootSummary.csv")),
            row.names = FALSE)
  write.csv(icboot$compare,
            file = file.path(iter.wd, paste0(iter.id, "_Compare.csv")),
            row.names = FALSE)

  # Try extracting group-level process summaries from icamp.boot output first
  proc = extract_process_table(icboot$summary,
                               groups.order = groups.order,
                               process.levels = process.levels,
                               source.name = "icamp.boot.summary")

  # Optional fallback using icamp.bins()
  if(is.null(proc) && run.bin.summary){

    icbin = icamp.bins(icamp.detail = icres$detail,
                       treat = treat.sub,
                       clas = taxon,
                       silent = TRUE,
                       boot = bin.boot,
                       rand.time = rand.time,
                       between.group = TRUE)

    saveRDS(icbin, file = file.path(iter.wd, paste0(iter.id, "_icbin.rds")))
    write.csv(icbin$Pt,
              file = file.path(iter.wd, paste0(iter.id, "_ProcessImportance_EachGroup.csv")),
              row.names = FALSE)

    proc = extract_process_table(icbin$Pt,
                                 groups.order = groups.order,
                                 process.levels = process.levels,
                                 source.name = "icamp.bins.Pt")
  }

  if(is.null(proc) || nrow(proc) == 0){
    writeLines("Could not automatically extract HoS/HeS/DL/HD/DR from icamp.boot$summary or icamp.bins$Pt. Please inspect the raw RDS/CSV files.",
               con = error.file)
    return(data.frame())
  }

  proc = add_metadata_to_process(proc, n.sub = n.sub, iter = iter,
                                 iter.id = iter.id)

  write.csv(proc, file = process.file, row.names = FALSE)

  return(proc)
}

##########################################################################
# 5. Run subsampling iCAMP sensitivity

all.process = list()
counter = 1

for(nn in n.grid){
  for(ii in seq_len(subsample.time)){
  # for(ii in  seq(63, subsample.time) ){

    message("Running iCAMP equalized sensitivity: n = ", nn,
            ", iteration = ", ii, " / ", subsample.time)

    res = tryCatch({
      run_one_icamp(n.sub = nn, iter = ii)
    }, error = function(e){
      err.path = file.path(out.wd, paste0("ERROR_n", nn, "_iter",
                                          sprintf("%03d", ii), ".txt"))
      writeLines(c(paste0("n = ", nn),
                   paste0("iteration = ", ii),
                   paste0("error = ", conditionMessage(e))),
                 con = err.path)
      return(data.frame())
    })

    if(nrow(res) > 0){
      all.process[[counter]] = res
      counter = counter + 1
    }
  }

  # Save incremental combined output after each n
  tmp = bind_rows(all.process)
  if(nrow(tmp) > 0){
    write.csv(tmp,
              file = file.path(out.wd, "iCAMP_equalized_ProcessLong_incremental.csv"),
              row.names = FALSE)
  }
}

process.long = bind_rows(all.process)

if(nrow(process.long) == 0){
  stop("No process summary was extracted. Please inspect icamp.boot$summary and icamp.bins$Pt formats.")
}

write.csv(process.long,
          file = file.path(out.wd, "iCAMP_equalized_ProcessLong.csv"),
          row.names = FALSE)

##########################################################################
# 6. Process-level summaries

process.iter = process.long %>%
  mutate(Group = factor(as.character(Group), levels = groups.order),
         Process = factor(as.character(Process), levels = process.levels)) %>%
  group_by(N, Iteration, IterationID, Group, Process) %>%
  summarise(Contribution = mean(Contribution, na.rm = TRUE),
            .groups = "drop")

process.summary = process.iter %>%
  group_by(N, Group, Process) %>%
  summarise(N_success = n_distinct(Iteration),
            Mean = mean(Contribution, na.rm = TRUE),
            Median = median(Contribution, na.rm = TRUE),
            SD = sd(Contribution, na.rm = TRUE),
            CI2.5 = quantile(Contribution, 0.025, na.rm = TRUE),
            CI97.5 = quantile(Contribution, 0.975, na.rm = TRUE),
            .groups = "drop")

write.csv(process.summary,
          file = file.path(out.wd, "iCAMP_equalized_ProcessSummary_MeanMedianCI.csv"),
          row.names = FALSE)

##########################################################################
# 7. Deterministic-stochastic regrouping
# Deterministic = HoS + HeS
# Stochastic/dispersal-related = DL + HD + DR

state.wide = process.iter %>%
  select(N, Iteration, IterationID, Group, Process, Contribution) %>%
  pivot_wider(names_from = Process, values_from = Contribution,
              values_fill = 0) %>%
  mutate(Deterministic = HoS + HeS,
         Stochastic_DispersalRelated = DL + HD + DR,
         DetMinusStoch = Deterministic - Stochastic_DispersalRelated,
         HoS_minus_DLDR = HoS - (DL + DR),
         HoS_minus_DR = HoS - DR,
         EndpointState = ifelse(Deterministic > Stochastic_DispersalRelated,
                                "Deterministic-enriched",
                                "Stochastic/dispersal-related-enriched"))

write.csv(state.wide,
          file = file.path(out.wd, "iCAMP_equalized_StateWide_ByIteration.csv"),
          row.names = FALSE)

state.long = state.wide %>%
  select(N, Iteration, IterationID, Group,
         Deterministic, Stochastic_DispersalRelated,
         DetMinusStoch, HoS_minus_DLDR, HoS_minus_DR) %>%
  pivot_longer(cols = c(Deterministic, Stochastic_DispersalRelated,
                        DetMinusStoch, HoS_minus_DLDR, HoS_minus_DR),
               names_to = "Metric", values_to = "Value")

state.summary = state.long %>%
  group_by(N, Group, Metric) %>%
  summarise(N_success = n_distinct(Iteration),
            Mean = mean(Value, na.rm = TRUE),
            Median = median(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            CI2.5 = quantile(Value, 0.025, na.rm = TRUE),
            CI97.5 = quantile(Value, 0.975, na.rm = TRUE),
            .groups = "drop")

write.csv(state.summary,
          file = file.path(out.wd, "iCAMP_equalized_StateSummary_MeanMedianCI.csv"),
          row.names = FALSE)

##########################################################################
# 8. Framework-level trend summaries
# These metrics directly correspond to the manuscript framework:
#   - D0 as stochastic/dispersal-related baseline
#   - diluted treatments as deterministic-enriched endpoints
#   - HoS increase after dilution/regrowth
#   - stochastic components remaining detectable

framework.iter = state.wide %>%
  group_by(N, Iteration, IterationID) %>%
  summarise(
    D0_Deterministic = Deterministic[Group == "D0"],
    D0_Stochastic = Stochastic_DispersalRelated[Group == "D0"],
    D0_HoS = HoS[Group == "D0"],
    D0_DL = DL[Group == "D0"],
    D0_DR = DR[Group == "D0"],

    D0_BelowParity = D0_Stochastic > D0_Deterministic,

    MeanDet_Diluted = mean(Deterministic[Group %in% diluted.groups], na.rm = TRUE),
    MeanStoch_Diluted = mean(Stochastic_DispersalRelated[Group %in% diluted.groups], na.rm = TRUE),
    MeanHoS_Diluted = mean(HoS[Group %in% diluted.groups], na.rm = TRUE),
    MeanDL_Diluted = mean(DL[Group %in% diluted.groups], na.rm = TRUE),
    MeanDR_Diluted = mean(DR[Group %in% diluted.groups], na.rm = TRUE),

    MeanDet_Strong = mean(Deterministic[Group %in% strong.groups], na.rm = TRUE),
    MeanStoch_Strong = mean(Stochastic_DispersalRelated[Group %in% strong.groups], na.rm = TRUE),
    MeanHoS_Strong = mean(HoS[Group %in% strong.groups], na.rm = TRUE),

    Diluted_DetFraction = mean(Deterministic[Group %in% diluted.groups] >
                                 Stochastic_DispersalRelated[Group %in% diluted.groups],
                               na.rm = TRUE),
    Strong_DetFraction = mean(Deterministic[Group %in% strong.groups] >
                               Stochastic_DispersalRelated[Group %in% strong.groups],
                             na.rm = TRUE),

    AllDiluted_DeterministicEnriched = all(Deterministic[Group %in% diluted.groups] >
                                             Stochastic_DispersalRelated[Group %in% diluted.groups]),
    AllStrong_DeterministicEnriched = all(Deterministic[Group %in% strong.groups] >
                                            Stochastic_DispersalRelated[Group %in% strong.groups]),

    DetDiff_DilutedMinusD0 = MeanDet_Diluted - D0_Deterministic,
    DetDiff_StrongMinusD0 = MeanDet_Strong - D0_Deterministic,
    HoSDiff_DilutedMinusD0 = MeanHoS_Diluted - D0_HoS,
    HoSDiff_StrongMinusD0 = MeanHoS_Strong - D0_HoS,

    # Alternative state-space summaries used to support Fig. S3 logic
    Mean_HoS_minus_DLDR_Diluted = mean(HoS_minus_DLDR[Group %in% diluted.groups], na.rm = TRUE),
    Mean_HoS_minus_DR_Diluted = mean(HoS_minus_DR[Group %in% diluted.groups], na.rm = TRUE),
    .groups = "drop"
  )

write.csv(framework.iter,
          file = file.path(out.wd, "iCAMP_equalized_FrameworkTrend_ByIteration.csv"),
          row.names = FALSE)

framework.summary = framework.iter %>%
  group_by(N) %>%
  summarise(
    N_success = n_distinct(Iteration),

    Prob_D0_BelowParity = mean(D0_BelowParity, na.rm = TRUE),
    Prob_AllDiluted_DeterministicEnriched = mean(AllDiluted_DeterministicEnriched, na.rm = TRUE),
    Prob_AllStrong_DeterministicEnriched = mean(AllStrong_DeterministicEnriched, na.rm = TRUE),

    Median_Diluted_DetFraction = median(Diluted_DetFraction, na.rm = TRUE),
    CI2.5_Diluted_DetFraction = quantile(Diluted_DetFraction, 0.025, na.rm = TRUE),
    CI97.5_Diluted_DetFraction = quantile(Diluted_DetFraction, 0.975, na.rm = TRUE),

    Median_DetDiff_DilutedMinusD0 = median(DetDiff_DilutedMinusD0, na.rm = TRUE),
    CI2.5_DetDiff_DilutedMinusD0 = quantile(DetDiff_DilutedMinusD0, 0.025, na.rm = TRUE),
    CI97.5_DetDiff_DilutedMinusD0 = quantile(DetDiff_DilutedMinusD0, 0.975, na.rm = TRUE),

    Median_DetDiff_StrongMinusD0 = median(DetDiff_StrongMinusD0, na.rm = TRUE),
    CI2.5_DetDiff_StrongMinusD0 = quantile(DetDiff_StrongMinusD0, 0.025, na.rm = TRUE),
    CI97.5_DetDiff_StrongMinusD0 = quantile(DetDiff_StrongMinusD0, 0.975, na.rm = TRUE),

    Median_HoSDiff_DilutedMinusD0 = median(HoSDiff_DilutedMinusD0, na.rm = TRUE),
    CI2.5_HoSDiff_DilutedMinusD0 = quantile(HoSDiff_DilutedMinusD0, 0.025, na.rm = TRUE),
    CI97.5_HoSDiff_DilutedMinusD0 = quantile(HoSDiff_DilutedMinusD0, 0.975, na.rm = TRUE),

    Median_HoSDiff_StrongMinusD0 = median(HoSDiff_StrongMinusD0, na.rm = TRUE),
    CI2.5_HoSDiff_StrongMinusD0 = quantile(HoSDiff_StrongMinusD0, 0.025, na.rm = TRUE),
    CI97.5_HoSDiff_StrongMinusD0 = quantile(HoSDiff_StrongMinusD0, 0.975, na.rm = TRUE),

    Median_HoS_minus_DLDR_Diluted = median(Mean_HoS_minus_DLDR_Diluted, na.rm = TRUE),
    Median_HoS_minus_DR_Diluted = median(Mean_HoS_minus_DR_Diluted, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(framework.summary,
          file = file.path(out.wd, "iCAMP_equalized_FrameworkTrend_Summary.csv"),
          row.names = FALSE)

##########################################################################
# 9. Trend consistency relative to n = 14
# n = 14 is the primary all-treatment common-n sensitivity.
# n = 8, 10 and 12 are used to test whether the trend is dependent on one n.

metric.long = state.wide %>%
  select(N, Iteration, Group, HoS, HeS, DL, HD, DR,
         Deterministic, Stochastic_DispersalRelated,
         DetMinusStoch, HoS_minus_DLDR, HoS_minus_DR) %>%
  pivot_longer(cols = c(HoS, HeS, DL, HD, DR,
                        Deterministic, Stochastic_DispersalRelated,
                        DetMinusStoch, HoS_minus_DLDR, HoS_minus_DR),
               names_to = "Metric", values_to = "Value")

metric.median = metric.long %>%
  group_by(N, Group, Metric) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")

ref.n14 = metric.median %>%
  filter(N == 14) %>%
  select(Group, Metric, Median_n14 = Median)

trend.consistency = metric.median %>%
  filter(N != 14) %>%
  inner_join(ref.n14, by = c("Group", "Metric")) %>%
  group_by(N, Metric) %>%
  summarise(SpearmanRho_vs_n14 = suppressWarnings(cor(Median, Median_n14,
                                                       method = "spearman",
                                                       use = "complete.obs")),
            PearsonR_vs_n14 = suppressWarnings(cor(Median, Median_n14,
                                                    method = "pearson",
                                                    use = "complete.obs")),
            MaxAbsDiff_vs_n14 = max(abs(Median - Median_n14), na.rm = TRUE),
            .groups = "drop")

direction.table = metric.median %>%
  group_by(N, Metric) %>%
  summarise(D0 = Median[Group == "D0"],
            MeanDiluted = mean(Median[Group %in% diluted.groups], na.rm = TRUE),
            MeanStrong = mean(Median[Group %in% strong.groups], na.rm = TRUE),
            Direction_DilutedMinusD0 = sign(MeanDiluted - D0),
            Direction_StrongMinusD0 = sign(MeanStrong - D0),
            .groups = "drop")

direction.ref = direction.table %>%
  filter(N == 14) %>%
  select(Metric,
         Direction_DilutedMinusD0_n14 = Direction_DilutedMinusD0,
         Direction_StrongMinusD0_n14 = Direction_StrongMinusD0)

direction.consistency = direction.table %>%
  filter(N != 14) %>%
  inner_join(direction.ref, by = "Metric") %>%
  mutate(SameDirection_DilutedMinusD0 =
           Direction_DilutedMinusD0 == Direction_DilutedMinusD0_n14,
         SameDirection_StrongMinusD0 =
           Direction_StrongMinusD0 == Direction_StrongMinusD0_n14)

trend.consistency = trend.consistency %>%
  left_join(direction.consistency,
            by = c("N", "Metric"))

write.csv(trend.consistency,
          file = file.path(out.wd, "iCAMP_equalized_TrendConsistency_vs_n14.csv"),
          row.names = FALSE)

##########################################################################
# 10. Optional compact plots for checking
# The numerical tables above should be used as the primary supplementary outputs.

pdf(file.path(out.wd, "iCAMP_equalized_n_grid_DeterministicStochastic.pdf"),
    width = 7, height = 4.5)
print(
  ggplot(state.summary %>%
           filter(Metric %in% c("Deterministic", "Stochastic_DispersalRelated")),
         aes(x = Group, y = Median, group = Metric, shape = Metric)) +
    geom_point(size = 2, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = CI2.5, ymax = CI97.5),
                  width = 0.15, position = position_dodge(width = 0.4)) +
    facet_wrap(~ N, nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Contribution",
         title = "iCAMP common-n sensitivity: deterministic vs stochastic/dispersal-related")
)
dev.off()

pdf(file.path(out.wd, "iCAMP_equalized_n_grid_ProcessContributions.pdf"),
    width = 8, height = 5)
print(
  ggplot(process.summary,
         aes(x = Group, y = Median, group = Process, shape = Process)) +
    geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = CI2.5, ymax = CI97.5),
                  width = 0.15, position = position_dodge(width = 0.5)) +
    facet_grid(Process ~ N, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Contribution",
         title = "iCAMP common-n sensitivity: HoS, HeS, DL, HD and DR")
)
dev.off()

##########################################################################
# 11. Session information

sink(file.path(out.wd, "SessionInfo_iCAMP_equalized_sensitivity.txt"))
print(sessionInfo())
sink()

##########################################################################
# End
##########################################################################
