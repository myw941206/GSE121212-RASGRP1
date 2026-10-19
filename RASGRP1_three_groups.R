# ============================================================
# GSE121212 原始 read counts 分析并绘制 RASGRP1 表达图
# 健康对照 38例 / 银屑病皮损 28例 / 特应性皮炎皮损 27例
# ============================================================

# 1. 路径与加载包 ----------------------------------------------

ROOT <- "D:/visual/GSE121212_RASGRP1_repro"

if (!dir.exists("D:/")) {
  stop("D盘不可用，请检查。")
}

DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "outputs", "three_groups_from_raw")

dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(ROOT, "R_library"), .libPaths()))

library(edgeR)
library(limma)
library(statmod)
library(digest)

set.seed(20261019)


# 2. 下载并读取原始计数矩阵 ------------------------------------

message("第1步：读取GEO原始read counts")

count_file <- file.path(
  DATA,
  "GSE121212_readcount.txt.gz"
)

if (!file.exists(count_file)) {

  options(timeout = 600)

  download.file(
    url = paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/",
      "GSE121nnn/GSE121212/suppl/",
      "GSE121212_readcount.txt.gz"
    ),
    destfile = count_file,
    mode = "wb"
  )
}

# 检查是否为本次核验过的原始文件
expected_sha256 <- paste0(
  "bdcfd8fb24858b4d11a5bed16f8383ba5",
  "b9d1203f6855104fc38ce4f1a3e6ca6"
)

stopifnot(
  digest(file = count_file, algo = "sha256") ==
    expected_sha256
)

raw <- read.delim(
  gzfile(count_file),
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

gene <- trimws(raw[[1]])
counts <- as.matrix(raw[, -1, drop = FALSE])

stopifnot(
  nrow(counts) == 31364L,
  ncol(counts) == 147L,
  is.numeric(counts),
  !anyNA(counts),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts))
)


# 3. 处理重复基因名称 ------------------------------------------
# 1-Mar、2-Mar各出现两次，排除全部4行歧义记录

duplicated_gene <- duplicated(gene) |
  duplicated(gene, fromLast = TRUE)

stopifnot(
  sum(duplicated_gene) == 4L,
  setequal(
    unique(gene[duplicated_gene]),
    c("1-Mar", "2-Mar")
  )
)

counts <- counts[!duplicated_gene, , drop = FALSE]
rownames(counts) <- gene[!duplicated_gene]

stopifnot(
  nrow(counts) == 31360L,
  sum(rownames(counts) == "RASGRP1") == 1L
)


# 4. 根据原始样本名称建立分组 ----------------------------------

sample_id <- colnames(counts)

meta <- data.frame(
  sample_id = sample_id,
  subject_id = sub(
    "^((CTRL|PSO|AD)_[0-9]+)_.*$",
    "\\1",
    sample_id
  ),
  disease = sub("_.*$", "", sample_id),
  tissue = sub(
    "^(CTRL|PSO|AD)_[0-9]+_",
    "",
    sample_id
  )
)

meta$group <- NA_character_

meta$group[
  meta$disease == "CTRL" &
    meta$tissue == "healthy"
] <- "CTRL"

meta$group[
  meta$disease == "PSO" &
    meta$tissue == "lesional"
] <- "PSO"

# AD的普通皮损与chronic_lesion均计入皮损组
meta$group[
  meta$disease == "AD" &
    meta$tissue %in% c("lesional", "chronic_lesion")
] <- "AD"

meta <- meta[!is.na(meta$group), , drop = FALSE]

meta$group <- factor(
  meta$group,
  levels = c("CTRL", "PSO", "AD")
)

counts <- counts[, meta$sample_id, drop = FALSE]

stopifnot(
  identical(
    as.integer(table(meta$group)),
    c(38L, 28L, 27L)
  ),
  !anyDuplicated(meta$subject_id),
  identical(colnames(counts), meta$sample_id)
)

# 在控制台查看实际分组人数
print(table(meta$group))


# 5. 过滤低表达基因并进行TMM标准化 -----------------------------

message("第2步：过滤低表达基因并进行TMM标准化")

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

y <- DGEList(
  counts = counts,
  genes = data.frame(gene = rownames(counts))
)

keep <- filterByExpr(
  y,
  group = meta$group,
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7
)

y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

if (!"RASGRP1" %in% rownames(y)) {
  stop("RASGRP1未通过表达过滤，不能进行主差异检验。")
}

cat("过滤后保留基因数：", nrow(y), "\n")


# 6. limma-voom差异分析 ---------------------------------------

message("第3步：根据原始计数拟合模型并计算P值")

v <- voom(
  y,
  design,
  normalize.method = "none",
  plot = FALSE
)

contrast_matrix <- makeContrasts(
  PSO_vs_CTRL = PSO - CTRL,
  AD_vs_CTRL  = AD - CTRL,
  levels = design
)

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = FALSE, robust = TRUE)

# 提取RASGRP1的两项检验结果
ras_result <- do.call(
  rbind,
  lapply(colnames(contrast_matrix), function(comparison) {

    result <- topTable(
      fit,
      coef = comparison,
      number = Inf,
      sort.by = "none",
      confint = 0.95
    )

    result <- result[
      result$gene == "RASGRP1",
      ,
      drop = FALSE
    ]

    result$comparison <- comparison
    result
  })
)

rownames(ras_result) <- NULL

# 本图预设两项比较：
# 银屑病vs对照、特应性皮炎vs对照
ras_result$P_Holm <- p.adjust(
  ras_result$P.Value,
  method = "holm"
)

print(
  ras_result[, c(
    "comparison", "logFC", "CI.L", "CI.R",
    "P.Value", "P_Holm"
  )]
)


# 7. 从本次标准化对象中计算绘图表达量 --------------------------

log_cpm <- cpm(
  y,
  log = TRUE,
  prior.count = 0.5
)

plot_data <- meta

plot_data$log2CPM <- as.numeric(
  log_cpm["RASGRP1", meta$sample_id]
)

# 查看刚刚计算出来的样本表达量
head(plot_data)


# 8. 使用R绘图 ------------------------------------------------

message("第4步：根据刚计算的表达量绘图")

draw_plot <- function() {

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  colors <- c("#8D99AE", "#D76A6A", "#4C9EAD")

  values <- split(
    plot_data$log2CPM,
    plot_data$group
  )

  ymax <- max(plot_data$log2CPM)
  ymin <- min(plot_data$log2CPM)
  span <- ymax - ymin

  par(
    mar = c(5.5, 5, 4, 1),
    mgp = c(2.8, 0.8, 0),
    las = 1
  )

  plot(
    NA,
    xlim = c(0.45, 3.55),
    ylim = c(ymin - 0.10 * span, ymax + 0.55 * span),
    xaxt = "n",
    bty = "l",
    xlab = "",
    ylab = "RASGRP1 expression (TMM log2CPM)"
  )

  abline(h = axTicks(2), col = "#EBEDF0")

  boxplot(
    values,
    at = 1:3,
    add = TRUE,
    axes = FALSE,
    col = adjustcolor(colors, alpha.f = 0.35),
    border = colors,
    outline = FALSE,
    boxwex = 0.48,
    lwd = 1.8
  )

  set.seed(20261019)

  for (i in 1:3) {

    points(
      jitter(rep(i, length(values[[i]])), amount = 0.13),
      values[[i]],
      pch = 21,
      bg = adjustcolor(colors[i], alpha.f = 0.85),
      col = "white",
      cex = 0.95
    )
  }

  axis(
    1,
    at = 1:3,
    labels = c(
      "Healthy control\n(n = 38)",
      "Psoriasis, lesional\n(n = 28)",
      "Atopic dermatitis, lesional\n(n = 27)"
    ),
    tick = FALSE,
    cex.axis = 0.9
  )

  title(
    expression(italic(RASGRP1) ~ "expression in GSE121212"),
    line = 2
  )

  add_bracket <- function(x2, height, p_value) {

    h <- 0.025 * span

    segments(1, height - h, 1, height)
    segments(1, height, x2, height)
    segments(x2, height, x2, height - h)

    text(
      (1 + x2) / 2,
      height + 0.045 * span,
      paste0(
        "Holm P = ",
        formatC(p_value, format = "e", digits = 2)
      ),
      cex = 0.85
    )
  }

  p_values <- setNames(
    ras_result$P_Holm,
    ras_result$comparison
  )

  add_bracket(
    2,
    ymax + 0.15 * span,
    p_values["PSO_vs_CTRL"]
  )

  add_bracket(
    3,
    ymax + 0.37 * span,
    p_values["AD_vs_CTRL"]
  )

  mtext(
    "limma-voom; Holm correction for 2 comparisons versus controls",
    side = 1,
    line = 4.3,
    cex = 0.7
  )
}

# 在RStudio右下角Plots窗口显示本次计算生成的图
if (interactive()) {
  draw_plot()
}


# 9. 将R刚生成的图保存到D盘 ------------------------------------

png(
  filename = file.path(
    OUT,
    "RASGRP1_three_groups_from_raw.png"
  ),
  width = 2400,
  height = 1950,
  res = 300
)

tryCatch(
  draw_plot(),
  finally = dev.off()
)

# 保存本次重新计算的数据和统计结果
write.csv(
  plot_data,
  file.path(OUT, "RASGRP1_expression_from_raw.csv"),
  row.names = FALSE
)

write.csv(
  ras_result,
  file.path(OUT, "RASGRP1_statistics_from_raw.csv"),
  row.names = FALSE
)

cat("分析和绘图完成，结果保存在：\n", OUT, "\n")