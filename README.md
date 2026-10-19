# GSE121212 RASGRP1 三组表达分析

从 GEO 原始 read counts 重新计算 RASGRP1 的标准化表达量、差异检验和比较图。

- 健康对照：38 位独立受试者。
- 银屑病皮损：28 位独立受试者。
- 特应性皮炎皮损：27 位独立受试者，包含 chronic_lesion。
- 过滤低表达基因后依次进行 TMM、limma-voom、稳健经验贝叶斯检验。
- 两项预设比较为 PSO vs CTRL、AD vs CTRL，对两项 RASGRP1 P 值作 Holm 校正。
- 排除 1-Mar 和 2-Mar 的全部四行重复歧义记录，剩余 31,360 行。
- 全部下载、依赖、安装包缓存、临时文件和结果位于 D:/visual 下。

## 文件

- RASGRP1_three_groups.R：用户提供的完整 R 代码，原样保留。
- install_dependencies.R：首次使用时安装缺失依赖。
- run_analysis.cmd：设置 D 盘临时目录，然后安装依赖并执行分析。
- verified_package_versions.csv：本次真实运行验证所用的项目库版本。

## 运行

将仓库克隆或下载到 D 盘。编辑 run_analysis.cmd 中的 RSCRIPT，使其指向本机 Rscript.exe。
当前验证机器使用 D:\Program Files\R\R-4.5.1\bin\Rscript.exe。

Windows PowerShell 中进入仓库后执行：

    .\run_analysis.cmd

已装好依赖时，也可以在 RStudio 打开 RASGRP1_three_groups.R，按顺序执行全部代码，
图会显示在 Plots 窗口。

脚本中 ROOT 固定为：

    D:/visual/GSE121212_RASGRP1_repro

原始数据不存在时自动从 GEO 下载。原始 gzip 文件通过 SHA-256 校验后直接读取。
无需读取此前导出的结果 CSV 或已经绘制好的图片。

## 输出

输出目录为：

    D:/visual/GSE121212_RASGRP1_repro/outputs/three_groups_from_raw

- RASGRP1_three_groups_from_raw.png
- RASGRP1_expression_from_raw.csv
- RASGRP1_statistics_from_raw.csv

## 验证

使用用户粘贴的原始脚本，在 Windows / R 4.5.1 / Bioconductor 3.22 上实际运行通过。

| 比较 | log2FC | Holm P（本图两项检验） |
|---|---:|---:|
| PSO vs CTRL | 1.1713289 | 1.546143e-17 |
| AD vs CTRL | 0.9016389 | 1.693991e-12 |

本图校正家族仅包含这两项与健康对照的比较，不包含另一个 AD 皮损/非皮损配对检验。
95% 置信区间为单项比较的 moderated 区间。
当前独立三组模型未调整年龄、性别、批次或细胞组成。

原始数据与安装缓存未上传。该分析结果属于此处明确实现的流程，
不代表已经逐项复现原研究文章的所有方法。

## 数据来源

[GEO GSE121212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121212)

