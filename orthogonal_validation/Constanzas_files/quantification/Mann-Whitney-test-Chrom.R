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
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,                # <-- EDIT: group 1 (WT) values, copy-pasted from the Excel "Sumup" table
    9.36, 4.29, 1.58, 1.79, 2.36, 1.36            # <-- EDIT: group 2 (PERKO) values, copy-pasted from the Excel "Sumup" table
  )
)

# FORCE ORDER (IMPORTANT)
data$Condition <- factor(
  data$Condition,
  levels = c("WT", "PERKO")   # <-- EDIT: must match the Condition names above, in the left-to-right plotting order you want
)

# =========================
# SUMMARY (Mean + SEM)
# =========================
summary_data <- data %>%
  group_by(Condition) %>%
  summarise(
    Mean = mean(Value),
    SEM = sd(Value) / sqrt(n())
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

# combined text shown on the plot: significance stars on top, numeric p-value below
plot_label <- paste0(sig_label, "\n", p_label)

# =========================
# Y POSITION FOR BRACKET
# =========================
y_max <- max(data$Value)
y_bracket <- y_max + 0.8   # <-- EDIT: raise/lower this offset if the bracket sits too close to (or far from) the data points

# =========================
# PLOT
# =========================
p <- ggplot() +

  # bars (mean)
  geom_bar(
    data = summary_data,
    aes(x = Condition, y = Mean),
    stat = "identity",
    fill = "grey80",
    width = 0.4   # <-- EDIT: bar width (narrowed from 0.6 -> 0.4; smaller number = thinner bars)
  ) +

  # SEM error bars
  geom_errorbar(
    data = summary_data,
    aes(
      x = Condition,
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
      x = Condition,
      y = Value,
      color = Experiment
    ),
    size = 3,
    position = position_jitter(width = 0.08)
  ) +

  # significance bracket line
  annotate(
    "segment",
    x = 1, xend = 2,
    y = y_bracket,
    yend = y_bracket,
    linewidth = 0.8
  ) +

  # vertical left line
  annotate(
    "segment",
    x = 1, xend = 1,
    y = y_bracket - 0.2,
    yend = y_bracket,
    linewidth = 0.8
  ) +

  # vertical right line
  annotate(
    "segment",
    x = 2, xend = 2,
    y = y_bracket - 0.2,
    yend = y_bracket,
    linewidth = 0.8
  ) +

  # significance stars + numeric p-value, shown together above the bracket
  annotate(
    "text",
    x = 1.5,
    y = y_bracket + 0.5,   # <-- EDIT: vertical offset of the label above the bracket (raise if the 2-line label overlaps the bracket)
    label = plot_label,    # <-- shows stars + "p = ..." (built above); was just `sig_label` before
    size = 5,
    lineheight = 0.9
  ) +

  labs(
    title = "SRSF3 Chromatin levels WT(CT20) vs PERKO",   # <-- EDIT: plot title
    x = "",
    y = "Fold change of SRSF3 levels (vs WT)"              # <-- EDIT: y-axis label
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.title = element_blank()
  )

print(p)

# =========================
# SAVE FIGURE
# =========================
ggsave(
  "SRSF3_chrom-WTvsKO-MannWhitney.png",   # <-- EDIT: output filename
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
