ROOT <- "D:/visual/GSE121212_RASGRP1_repro"
if (!dir.exists("D:/")) stop("D: drive unavailable.")
for (d in c("R_library", "package_cache", "temp"))
  dir.create(file.path(ROOT, d), recursive = TRUE, showWarnings = FALSE)
if (!startsWith(tolower(normalizePath(tempdir(), winslash = "/")),
                tolower(paste0(ROOT, "/temp/"))))
  stop("Start with run_analysis.cmd to keep installation temporary files on D:.")
.libPaths(c(file.path(ROOT, "R_library"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
cache <- file.path(ROOT, "package_cache")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = .libPaths()[1], destdir = cache)
pkgs <- c("edgeR", "limma", "statmod", "digest")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  BiocManager::install(missing, lib = .libPaths()[1], destdir = cache,
                       ask = FALSE, update = FALSE)
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))

