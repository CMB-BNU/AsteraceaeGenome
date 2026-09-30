#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
# ===============================================================================
# 脚本名称: 01plot_final_figures.R
# 功能描述: 
#   1. Figure A: 纬度与 psbA 相对表达量拟合图（ANCOVA左上，相关性右上左对齐）
#   2. Figure B: 低/高纬度 psbA 相对表达量精美柱状图 (FC + italic p 值标注)
#   3. Figure C: 基于残差的 psbA vs NLR 偏相关图（拉长版）
#   4. Combined Figure: 布局 (A + B) / C，支持 SVG、PDF 与 PNG 导出
# ===============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(svglite)
})

# 1. 命令行参数解析
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  input_csv <- "00data/normalized_summary_matrix.csv"
  out_dir <- "02plots"
} else {
  input_csv <- args[1]
  out_dir <- args[2]
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 全局绘图主题
theme_custom <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
}

# ===============================================================================
# 绘图函数: Figure A
# ===============================================================================
plot_figA <- function(data) {
  df_long <- data %>%
    dplyr::select(Sample, Population, Latitude, Sum_psbA_pos, Sum_psbA_non_pos) %>%
    pivot_longer(
      cols = c(Sum_psbA_pos, Sum_psbA_non_pos),
      names_to = "Status",
      values_to = "Norm_Exp"
    ) %>%
    mutate(
      Status = ifelse(Status == "Sum_psbA_pos", "Positively Selected", "Non-positively Selected"),
      Status = factor(Status, levels = c("Positively Selected", "Non-positively Selected"))
    )
  
  # ANCOVA 分析计算
  ancova_model <- lm(Norm_Exp ~ Latitude * Status, data = df_long)
  coef_matrix <- summary(ancova_model)$coefficients
  
  slope_term <- grep("Latitude:Status", rownames(coef_matrix), value = TRUE)
  intercept_term <- grep("^Status", rownames(coef_matrix), value = TRUE)
  
  p_slope <- coef_matrix[slope_term[1], "Pr(>|t|)"]
  p_intercept <- coef_matrix[intercept_term[1], "Pr(>|t|)"]
  
  p <- ggplot(df_long, aes(x = Latitude, y = Norm_Exp, color = Status, fill = Status)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("Positively Selected" = "#70CDBE", "Non-positively Selected" = "#AC99D2")) +
    scale_fill_manual(values = c("Positively Selected" = "#70CDBE", "Non-positively Selected" = "#AC99D2")) +
    labs(
      x = "Latitude (°N)",
      y = expression(bold("Relative Expression")),
      color = "psbA Status",
      fill = "psbA Status"
    ) +
    # 左上角 ANCOVA 结果（严格转义字符，修复解析错误）
    annotate("text", x = 23.5, y = 2.85, label = "ANCOVA:", hjust = 0, size = 3.8, fontface = "bold", color = "black") +
    annotate("text", x = 23.5, y = 2.68, label = sprintf("'Slope diff ' ~ italic(p) == %.3f", p_slope), hjust = 0, parse = TRUE, size = 3.8, color = "black") +
    annotate("text", x = 23.5, y = 2.51, label = sprintf("'Intercept diff ' ~ italic(p) == %.3f", p_intercept), hjust = 0, parse = TRUE, size = 3.8, color = "black") +
    
    # 右上角相关性结果（x=38.5 处 hjust=0 实现严格“左对齐”）
    annotate("text", x = 34.5, y = 2.85, label = "'Positively selected: R = 0.440, ' ~ italic(p) == 0.052", hjust = 0, parse = TRUE, size = 3.8, color = "#70CDBE") +
    annotate("text", x = 34.5, y = 2.68, label = "'Non-selected: R = 0.398, ' ~ italic(p) == 0.082", hjust = 0, parse = TRUE, size = 3.8, color = "#AC99D2") +
    
    scale_y_continuous(limits = c(0.6, 2.95), breaks = seq(1.0, 2.5, by = 0.5)) +
    theme_custom() +
    theme(
      legend.position = "top",
      legend.margin = margin(b = -5)
    )
  
  return(p)
}

# ===============================================================================
# 绘图函数: Figure B
# ===============================================================================
plot_figB <- function() {
  df_raw <- data.frame(
    Sample = paste0("S", 1:5),
    Population = factor(c("FK", "FK", "FK", "MDJ", "MDJ"), levels = c("FK", "MDJ")),
    Lat_Group = factor(c("Low Latitude", "Low Latitude", "Low Latitude", 
                         "High Latitude", "High Latitude"),
                       levels = c("High Latitude", "Low Latitude")),
    `Positively Selected` = c(0.24165764, 0.22393111, 0.393748142, 0.737933949, 0.741708825),
    `Non-positively Selected` = c(0.290769989, 0.245422331, 0.352835966, 0.608761574, 0.612056674),
    check.names = FALSE
  )
  
  df_long <- df_raw %>%
    pivot_longer(
      cols = c("Positively Selected", "Non-positively Selected"),
      names_to = "Gene_Class",
      values_to = "Relative_Expression"
    ) %>%
    mutate(
      Gene_Class = factor(Gene_Class, levels = c("Positively Selected", "Non-positively Selected")),
      Fill_Group = factor(
        paste(Gene_Class, Lat_Group, sep = "_"),
        levels = c(
          "Positively Selected_High Latitude",
          "Positively Selected_Low Latitude",
          "Non-positively Selected_High Latitude",
          "Non-positively Selected_Low Latitude"
        )
      )
    )
  
  df_summary <- df_long %>%
    group_by(Gene_Class, Lat_Group, Fill_Group) %>%
    summarise(
      Mean_Expression = mean(Relative_Expression),
      .groups = "drop"
    )
  
  df_bracket <- data.frame(
    Gene_Class = factor(c("Positively Selected", "Non-positively Selected"), 
                        levels = c("Positively Selected", "Non-positively Selected")),
    x1 = 1, 
    x2 = 2,
    y_top = 0.81,
    y_tick = 0.79,
    y_text = 0.87,
    label = c("atop('FC = 2.58', italic(p) == 0.0139)", "atop('FC = 2.06', italic(p) == 0.0096)")
  )
  
  palette_colors <- c(
    "Positively Selected_High Latitude"    = "#70CDBE",
    "Positively Selected_Low Latitude"     = "#C2F0E8",
    "Non-positively Selected_High Latitude"= "#AC99D2",
    "Non-positively Selected_Low Latitude" = "#E4DCF2"
  )
  
  p <- ggplot() +
    geom_bar(
      data = df_summary,
      aes(x = Lat_Group, y = Mean_Expression, fill = Fill_Group),
      stat = "identity",
      width = 0.52,
      color = NA,
      show.legend = FALSE
    ) +
    geom_jitter(
      data = df_long,
      aes(x = Lat_Group, y = Relative_Expression),
      color = "#E45756",
      fill = "#E45756",
      size = 1.3,
      width = 0.06,
      height = 0,
      shape = 21,
      alpha = 0.90
    ) +
    geom_segment(
      data = df_bracket,
      aes(x = x1, xend = x2, y = y_top, yend = y_top),
      linewidth = 0.35,
      color = "black"
    ) +
    geom_segment(
      data = df_bracket,
      aes(x = x1, xend = x1, y = y_top, yend = y_tick),
      linewidth = 0.35,
      color = "black"
    ) +
    geom_segment(
      data = df_bracket,
      aes(x = x2, xend = x2, y = y_top, yend = y_tick),
      linewidth = 0.35,
      color = "black"
    ) +
    geom_text(
      data = df_bracket,
      aes(x = 1.5, y = y_text, label = label),
      size = 3.5,
      parse = TRUE,
      lineheight = 0.85
    ) +
    facet_wrap(~Gene_Class, scales = "fixed") +
    scale_fill_manual(values = palette_colors) +
    scale_y_continuous(
      limits = c(0, 0.98),
      expand = c(0, 0),
      breaks = seq(0, 0.8, by = 0.2)
    ) +
    labs(
      x = NULL,
      y = expression(bold("Relative Expression"))
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 9.5, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, color = "black", margin = margin(r = 8)),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      panel.spacing = unit(1.0, "lines")
    )
  
  return(p)
}

# ===============================================================================
# 绘图函数: Figure C
# ===============================================================================
plot_figC <- function(data) {
  res_psbA <- lm(Sum_psbA_pos ~ Latitude, data = data)$residuals
  res_NLR  <- lm(Sum_NLR ~ Latitude, data = data)$residuals
  
  df_res <- data %>%
    mutate(
      psbA_res = res_psbA,
      NLR_res  = res_NLR
    )
  
  cor_test <- cor.test(df_res$psbA_res, df_res$NLR_res, method = "pearson")
  r_val <- cor_test$estimate[[1]]
  p_val <- cor_test$p.value
  
  p <- ggplot(df_res, aes(x = psbA_res, y = NLR_res)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(aes(color = Latitude), size = 3.5, alpha = 0.9) +
    geom_smooth(method = "lm", color = "black", fill = "grey70", se = TRUE, alpha = 0.2) +
    scale_color_gradient(low = "#8FB4DC", high = "#F5AA61", name = "Latitude (°N)") +
    labs(
      x = "psbA genes Expression (Residuals)",
      y = "NLR genes Expression (Residuals)"
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = sprintf("'Partial R = %.3f, ' ~ italic(p) == %.3f", r_val, p_val), 
             hjust = -0.05, vjust = 1.2, size = 4, parse = TRUE) +
    theme_custom()
  
  return(p)
}

# ===============================================================================
# 主程序运行与导出
# ===============================================================================
df_full <- read.csv(input_csv, stringsAsFactors = FALSE)

# 1. 绘制三张独立图件
p_figA <- plot_figA(df_full)
p_figB <- plot_figB()
p_figC <- plot_figC(df_full)

# 2. 导出单图 (SVG, PDF, PNG)
ggsave(file.path(out_dir, "FigureA_Main.svg"), p_figA, width = 7.0, height = 5.2, device = "svg")
ggsave(file.path(out_dir, "FigureA_Main.pdf"), p_figA, width = 7.0, height = 5.2)

ggsave(file.path(out_dir, "FigureB_Main.svg"), p_figB, width = 5.2, height = 4.2, device = "svg")
ggsave(file.path(out_dir, "FigureB_Main.pdf"), p_figB, width = 5.2, height = 4.2)

ggsave(file.path(out_dir, "FigureC_Main.svg"), p_figC, width = 10.0, height = 4.5, device = "svg")
ggsave(file.path(out_dir, "FigureC_Main.pdf"), p_figC, width = 10.0, height = 4.5)

# 3. 组合图布局
combined_figure <- (p_figA + p_figB) / p_figC +
  plot_layout(heights = c(1, 0.85), widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")

# 4. 导出组合图
ggsave(file.path(out_dir, "Combined_Main_Figure.svg"), combined_figure, width = 11, height = 9.0, device = "svg")
ggsave(file.path(out_dir, "Combined_Main_Figure.pdf"), combined_figure, width = 11, height = 9.0)
ggsave(file.path(out_dir, "Combined_Main_Figure.png"), combined_figure, width = 11, height = 9.0, dpi = 300)

cat("[Complete] 运行成功！已在", out_dir, "输出修正后的 SVG/PDF/PNG 图件。\n")
