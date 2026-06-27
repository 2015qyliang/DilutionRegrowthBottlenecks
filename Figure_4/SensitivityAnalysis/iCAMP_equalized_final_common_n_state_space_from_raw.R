# coding: utf-8

##########################################################################
# This script is standalone for the final figure:
#   - It starts from the per-iteration raw process files:
#       iCAMP_equalized/n8/n8_iter001/n8_iter001_ProcessLong.csv
#       ...
#       iCAMP_equalized/n14/n14_iter100/n14_iter100_ProcessLong.csv
#
##########################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

##########################################################################
# 1. Settings

out.wd = "iCAMP_equalized"
n.grid = c(8, 10, 12, 14)
subsample.time = 100

groups.order = c("D0", "D1", "D3", "D5", "D6", "D7")
process.levels = c("HoS", "HeS", "DL", "HD", "DR")

if(!dir.exists(out.wd)){
  stop("Cannot find output directory: ", out.wd)
}

##########################################################################
# 2. Helpers

qtl = function(x, p) {
  x = suppressWarnings(as.numeric(x))
  x = x[is.finite(x)]
  if(length(x) == 0){return(NA_real_)}
  unname(quantile(x, p, na.rm = TRUE, names = FALSE))
}

read_process_file = function(file) {
  df = read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required.cols = c("N", "Iteration", "IterationID",
                    "Group", "Process", "Contribution")
  missing.cols = setdiff(required.cols, colnames(df))

  if(length(missing.cols) > 0){
    stop("Missing required columns in ", file, ": ",
         paste(missing.cols, collapse = ", "))
  }

  if(!("Source" %in% colnames(df))){
    df$Source = NA_character_
  }

  data.frame(
    N = as.integer(df$N),
    Iteration = as.integer(df$Iteration),
    IterationID = as.character(df$IterationID),
    Group = as.character(df$Group),
    Process = as.character(df$Process),
    Contribution = suppressWarnings(as.numeric(df$Contribution)),
    Source = as.character(df$Source),
    SourceFile = normalizePath(file, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

process_summary_one = function(df) {
  data.frame(
    N = df$N[1],
    Iteration = df$Iteration[1],
    IterationID = df$IterationID[1],
    Group = df$Group[1],
    Process = df$Process[1],
    Contribution = mean(df$Contribution, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

get_process_value = function(df, process.name) {
  value = df$Contribution[df$Process == process.name]
  if(length(value) == 0){return(0)}
  mean(value, na.rm = TRUE)
}

state_one = function(df) {
  hos = get_process_value(df, "HoS")
  hes = get_process_value(df, "HeS")
  dl = get_process_value(df, "DL")
  hd = get_process_value(df, "HD")
  dr = get_process_value(df, "DR")

  deterministic = hos + hes
  stochastic = dl + hd + dr

  data.frame(
    N = df$N[1],
    Iteration = df$Iteration[1],
    IterationID = df$IterationID[1],
    Group = df$Group[1],
    HoS = hos,
    HeS = hes,
    DL = dl,
    HD = hd,
    DR = dr,
    Deterministic = deterministic,
    Stochastic_DispersalRelated = stochastic,
    DetMinusStoch = deterministic - stochastic,
    stringsAsFactors = FALSE
  )
}

state_summary_one = function(df) {
  det = df$Deterministic
  stoch = df$Stochastic_DispersalRelated

  data.frame(
    N = df$N[1],
    Group = df$Group[1],
    N_success = length(unique(df$Iteration)),
    Deterministic_median = median(det, na.rm = TRUE),
    Deterministic_CI2.5 = qtl(det, 0.025),
    Deterministic_CI97.5 = qtl(det, 0.975),
    Stochastic_median = median(stoch, na.rm = TRUE),
    Stochastic_CI2.5 = qtl(stoch, 0.025),
    Stochastic_CI97.5 = qtl(stoch, 0.975),
    DetMinusStoch_median = median(df$DetMinusStoch, na.rm = TRUE),
    Prob_deterministic_enriched = mean(det > stoch, na.rm = TRUE),
    Prob_stochastic_dominated = mean(stoch > det, na.rm = TRUE),
    Median_state = ifelse(median(det, na.rm = TRUE) >
                            median(stoch, na.rm = TRUE),
                          "Deterministic-enriched",
                          "Stochastic-dominated"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

draw_error_xy = function(x, y, xmin, xmax, ymin, ymax, col) {
  segments(xmin, y, xmax, y, col = col, lwd = 1.35)
  segments(x, ymin, x, ymax, col = col, lwd = 1.35)
}

##########################################################################
# 3. Read per-iteration process files from scratch

process.files = unlist(lapply(n.grid, function(nn){
  n.dir = file.path(out.wd, paste0("n", nn))
  if(!dir.exists(n.dir)){return(character())}
  list.files(n.dir, pattern = "_ProcessLong[.]csv$",
             recursive = TRUE, full.names = TRUE)
}), use.names = FALSE)

process.files = sort(unique(process.files))

if(length(process.files) == 0){
  stop("No *_ProcessLong.csv files found under ", out.wd, "/n*/")
}

message("Reading ", length(process.files), " per-iteration ProcessLong files...")

process.raw = do.call(rbind, lapply(process.files, read_process_file))
rownames(process.raw) = NULL

process.raw = process.raw[
  process.raw$N %in% n.grid &
    process.raw$Group %in% groups.order &
    process.raw$Process %in% process.levels &
    is.finite(process.raw$Contribution),
]

if(nrow(process.raw) == 0){
  stop("No valid HoS/HeS/DL/HD/DR rows were found.")
}

process.raw = unique(process.raw)

##########################################################################
# 4. Rebuild state-space coordinates from process rows

proc.key = interaction(process.raw$N, process.raw$Iteration,
                       process.raw$IterationID,
                       process.raw$Group, process.raw$Process,
                       drop = TRUE, lex.order = TRUE)

process.iter = do.call(rbind, lapply(split(process.raw, proc.key),
                                     process_summary_one))
rownames(process.iter) = NULL

state.key = interaction(process.iter$N, process.iter$Iteration,
                        process.iter$IterationID, process.iter$Group,
                        drop = TRUE, lex.order = TRUE)

state.iter = do.call(rbind, lapply(split(process.iter, state.key), state_one))
state.iter = state.iter[order(state.iter$N, state.iter$Iteration,
                              match(state.iter$Group, groups.order)), ]
rownames(state.iter) = NULL

state.summary = do.call(
  rbind,
  lapply(split(state.iter, list(state.iter$N, state.iter$Group), drop = TRUE),
         state_summary_one)
)

state.summary = state.summary[
  order(state.summary$N, match(state.summary$Group, groups.order)), ]
rownames(state.summary) = NULL

write.csv(state.summary,
          file.path(out.wd,
                    "iCAMP_equalized_FinalCommonN_StateSpace_from_raw_Summary.csv"),
          row.names = FALSE)

iteration.check = aggregate(Iteration ~ N,
                            data = unique(state.iter[c("N", "Iteration")]),
                            FUN = length)
colnames(iteration.check) = c("N", "N_success")
iteration.check = merge(data.frame(N = n.grid), iteration.check,
                        by = "N", all.x = TRUE)
iteration.check$N_success[is.na(iteration.check$N_success)] = 0
iteration.check = iteration.check[order(iteration.check$N), ]
rownames(iteration.check) = NULL

if(any(iteration.check$N_success != subsample.time)){
  warning("Some common-n settings do not have ",
          subsample.time, " successful iterations.")
}

##########################################################################
# 5. ggplot state-space figure

s1.blue = "#0072B2"
s1.orange = "#D55E00"

group.cols = c(D0 = "#A65E43",
               D1 = "#56B4E9",
               D3 = "#E69F00",
               D5 = "#009E73",
               D6 = "#7B3294",
               D7 = "#6A9FB5")
group.shapes = c(D0 = 24, D1 = 21, D3 = 21,
                 D5 = 21, D6 = 21, D7 = 21)

plot.df = state.summary
plot.df$Group = factor(plot.df$Group, levels = groups.order)
plot.df$N_label = factor(paste0("n = ", plot.df$N),
                         levels = paste0("n = ", n.grid))
plot.df$Stochastic_percent = plot.df$Stochastic_median * 100
plot.df$Stochastic_CI2.5_percent = plot.df$Stochastic_CI2.5 * 100
plot.df$Stochastic_CI97.5_percent = plot.df$Stochastic_CI97.5 * 100
plot.df$Deterministic_percent = plot.df$Deterministic_median * 100
plot.df$Deterministic_CI2.5_percent = plot.df$Deterministic_CI2.5 * 100
plot.df$Deterministic_CI97.5_percent = plot.df$Deterministic_CI97.5 * 100

label.dx = c(D0 = 4.5, D1 = 3.2, D3 = -7.0,
             D5 = -7.0, D6 = 3.2, D7 = 3.2)
label.dy = c(D0 = 1.2, D1 = 1.2, D3 = 3.2,
             D5 = 3.2, D6 = 3.2, D7 = 1.2)
plot.df$Label_x = plot.df$Stochastic_percent + label.dx[as.character(plot.df$Group)]
plot.df$Label_y = plot.df$Deterministic_percent + label.dy[as.character(plot.df$Group)]

region.df = do.call(rbind, lapply(levels(plot.df$N_label), function(nn){
  data.frame(
    N_label = nn,
    Region = rep(c("Deterministic-enriched", "Stochastic-dominated"), each = 3),
    x = c(0, 0, 85, 0, 85, 85),
    y = c(0, 85, 85, 0, 0, 85),
    stringsAsFactors = FALSE
  )
}))
region.df$N_label = factor(region.df$N_label, levels = levels(plot.df$N_label))

region.label.df = do.call(rbind, lapply(levels(plot.df$N_label), function(nn){
  data.frame(
    N_label = nn,
    Label = c("Deterministic\nenriched", "Stochastic\ndominated"),
    x = c(17, 58),
    y = c(45, 15),
    LabelColor = c("deterministic_label", "stochastic_label"),
    stringsAsFactors = FALSE
  )
}))
region.label.df$N_label = factor(region.label.df$N_label,
                                 levels = levels(plot.df$N_label))

state.plot = ggplot() +
  geom_polygon(data = region.df,
               aes(x = x, y = y, group = Region, fill = Region),
               alpha = 0.08, color = NA, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", linewidth = 0.55, color = "gray30") +
  geom_segment(data = plot.df,
               aes(x = Stochastic_CI2.5_percent,
                   xend = Stochastic_CI97.5_percent,
                   y = Deterministic_percent,
                   yend = Deterministic_percent,
                   color = Group),
               linewidth = 0.55, show.legend = FALSE) +
  geom_segment(data = plot.df,
               aes(x = Stochastic_percent,
                   xend = Stochastic_percent,
                   y = Deterministic_CI2.5_percent,
                   yend = Deterministic_CI97.5_percent,
                   color = Group),
               linewidth = 0.55, show.legend = FALSE) +
  geom_point(data = plot.df,
             aes(x = Stochastic_percent,
                 y = Deterministic_percent,
                 color = Group, fill = Group, shape = Group),
             size = 3.0, stroke = 0.75) +
  geom_text(data = plot.df,
            aes(x = Label_x, y = Label_y, label = Group),
            size = 3.3, color = "gray15") +
  geom_text(data = region.label.df,
            aes(x = x, y = y, label = Label, color = LabelColor),
            inherit.aes = FALSE, fontface = "italic", size = 3.0,
            show.legend = FALSE) +
  facet_wrap(~ N_label, ncol = 2) +
  coord_fixed(xlim = c(0, 85), ylim = c(0, 85), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 80, 20)) +
  scale_y_continuous(breaks = seq(0, 80, 20)) +
  scale_fill_manual(values = c("Deterministic-enriched" = s1.blue,
                               "Stochastic-dominated" = s1.orange,
                               group.cols),
                    breaks = groups.order) +
  scale_color_manual(values = c(group.cols,
                                deterministic_label = s1.blue,
                                stochastic_label = s1.orange),
                     breaks = groups.order) +
  scale_shape_manual(values = group.shapes, breaks = groups.order) +
  labs(
    title = "Common-n sensitivity judgement of deterministic-enriched endpoint trend",
    x = "Stochastic / dispersal-related contribution\n(DL + HD + DR, %)",
    y = "Deterministic contribution\n(HoS + HeS, %)",
    color = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.55),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "top",
    legend.justification = "center",
    legend.key = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    axis.title = element_text(size = 11),
    axis.text = element_text(color = "black")
  )

ggsave(file.path(out.wd, "iCAMP_equalized_FinalCommonN_StateSpace_from_raw_ggplot2x2.pdf"),
       plot = state.plot, width = 7.6, height = 7.4, units = "in")

ggsave(file.path(out.wd, "iCAMP_equalized_FinalCommonN_StateSpace_from_raw_ggplot2x2.png"),
       plot = state.plot, width = 7.6, height = 7.4, units = "in",
       dpi = 300)

message("Saved: ",
        file.path(out.wd, "iCAMP_equalized_FinalCommonN_StateSpace_from_raw_ggplot2x2.pdf"))
message("Saved: ",
        file.path(out.wd, "iCAMP_equalized_FinalCommonN_StateSpace_from_raw_ggplot2x2.png"))
message("Saved: ",
        file.path(out.wd,
                  "iCAMP_equalized_FinalCommonN_StateSpace_from_raw_Summary.csv"))
print(iteration.check)

##########################################################################
# End
##########################################################################
