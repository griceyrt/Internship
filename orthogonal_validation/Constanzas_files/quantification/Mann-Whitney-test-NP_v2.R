library(ggplot2)
library(dplyr)

# =========================
# DATA
# =========================
# <-- EDIT HERE for new numbers: replace Condition names, Experiment labels, and Values below.
#     Values are listed as: 6 values for the first group, then 6 values for the second group,
#     in the same order as they appear in the "Sumup" sheet of the source Excel file (e.g. Sumup-2).
data <- data.frame(
  Condition = rep(c("WT", "PERKO"), each = 6),                              # <-- EDIT: group names (2 groups, must match `levels=` below)
  Experiment = rep(c("Exp1", "Exp2", "Exp3", "Exp4", "Exp5", "Exp6"), 2),   # <-- EDIT: replicate/experiment labels
  Value = c(
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,   # <-- EDIT: group 1 (WT) values, copy-pasted from the Excel "Sumup" table
    0.9, 0.9, 1.2, 0.8, 0.6, 0.7    # <-- EDIT: group 2 (PERKO) values, copy-pasted from the Excel "Sumup" table
  )
)

# FORCE ORDER (IMPORTANT)
data$Condition <- factor(
  data$Condition,
  levels = c("WT", "PERKO")   # <-- EDIT: must match the Condition names above, in the left-to-right plotting order you want
)

# numeric x positions, closer together than the default 1/2 discrete spacing (bars closer + less side margin)
x_positions <- c(WT = 1, PERKO = 1.6)   # <-- EDIT: shrink/grow the gap between bars by changing 1.6
data$x <- x_positions[as.character(data$Condition)]

# =========================
# SUMMARY (Mean + SEM)
# =========================
summary_data <- data %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Value),
    SEM = sd(Value) / sqrt(n()),
    x = first(x)
  )

# =========================
# MANN-WHITNEY TEST
# =========================
mw_test <- wilcox.test(
  Value ~ Condition,
  data = data,
  exact = FALSE
)

p_value <- mw_test$p.value

print(mw_test)
print(paste("p-value =", p_value))

# significance label (stars)
sig_label <- ifelse(
  p_value < 0.001, "***",
  ifelse(
    p_value < 0.01, "**",
    ifelse(
      p_value < 0.05, "*", "ns"
    )
  )
)

# numeric p-value, formatted for display on the plot (e.g. "0.0022" or "<0.0001")
p_label <- paste0("p = ", format.pval(p_value, digits = 2, eps = 0.0001))  # <-- EDIT: change `digits` for more/fewer decimal places

# combined text shown on the plot: significance stars on top, numeric p-value below (same format as the Chrom script)
plot_label <- paste0(sig_label, "\n", p_label)

# =========================
# Y POSITION FOR BRACKET
# =========================
y_max <- max(data$Value)
y_bracket <- y_max + 0.8   # <-- EDIT: raise/lower this offset if the bracket sits too close to (or far from) the data points
bracket_tick <- 0.1        # <-- EDIT: length of the short vertical lines at each end of the bracket
x1 <- x_positions[["WT"]]
x2 <- x_positions[["PERKO"]]

# =========================
# PLOT
# =========================
p <- ggplot() +

  # bars (mean)
  geom_bar(
    data = summary_data,
    aes(x = x, y = Mean),
    stat = "identity",
    fill = "grey80",
    width = 0.4   # <-- EDIT: bar width (narrowed from 0.6 -> 0.4; smaller number = thinner bars)
  ) +

  # SEM error bars
  geom_errorbar(
    data = summary_data,
    aes(
      x = x,
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
    width = 0.15,   # <-- EDIT: error bar cap width (narrowed to match the thinner bars above)
    linewidth = 0.8
  ) +

  # individual points
  geom_point(
    data = data,
    aes(
      x = x,
      y = Value,
      color = Experiment
    ),
    size = 3,
    position = position_jitter(width = 0.08)
  ) +

  # significance bracket line
  annotate(
    "segment",
    x = x1, xend = x2,
    y = y_bracket,
    yend = y_bracket,
    linewidth = 0.8
  ) +

  # vertical left line
  annotate(
    "segment",
    x = x1, xend = x1,
    y = y_bracket - bracket_tick,   # <-- EDIT: shorten/lengthen via `bracket_tick` above
    yend = y_bracket,
    linewidth = 0.8
  ) +

  # vertical right line
  annotate(
    "segment",
    x = x2, xend = x2,
    y = y_bracket - bracket_tick,   # <-- EDIT: shorten/lengthen via `bracket_tick` above
    yend = y_bracket,
    linewidth = 0.8
  ) +

  # significance stars (top line) + numeric p-value (bottom line), shown together above the bracket
  annotate(
    "text",
    x = mean(c(x1, x2)),
    y = y_bracket + 0.15,   # <-- EDIT: vertical offset of the label above the bracket (smaller than Chrom's since NP's y-axis range is much smaller)
    label = plot_label,     # <-- shows stars + "p = ..." on two lines (built above)
    size = 4,
    lineheight = 0.6   # <-- EDIT: lower value = stars sit closer to the p-value line below them
  ) +

  # x-axis: named breaks at the bar positions, trimmed side margins
  scale_x_continuous(
    breaks = x_positions,
    labels = names(x_positions),
    expand = expansion(add = 0.25)   # <-- EDIT: lower value = less empty space on the left/right sides
  ) +

  # y-axis: whole-number ticks only
  scale_y_continuous(breaks = seq(0, ceiling(y_bracket), by = 1)) +   # <-- EDIT: change `by=` for tick spacing

  labs(
    title = "SRSF3 Nucleoplasm levels\nWT(CT20) vs PERKO",   # <-- EDIT: plot title (use "\n" to wrap onto a 2nd line if it's long)
    x = "",
    y = "Fold change of SRSF3 levels (vs WT)"                 # <-- EDIT: y-axis label
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, lineheight = 1),   # <-- EDIT: shrink/center the title so it fits the narrower plot
    axis.text.x = element_text(size = 16)   # <-- EDIT: font size of the "WT" / "PERKO" labels
  )

print(p)

# =========================
# SAVE FIGURE
# =========================
ggsave(
  "SRSF3_Nucleoplasm-WTvsKO-MannWhitney_v2.png",   # <-- EDIT: output filename
  plot = p,
  width = 5,    # <-- EDIT: overall figure width in inches (reduced from 8 for a smaller plot)
  height = 4,   # <-- EDIT: overall figure height in inches (reduced from 6)
  dpi = 300
)
