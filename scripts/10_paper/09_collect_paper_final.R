#!/usr/bin/env Rscript
# Stage 10: copies every final figure and table into results/paper_final/<version>/ under its paper name, with a
# MANIFEST.tsv (source path, size, modification time, md5). The list is config/paper_outputs.tsv. Each source pattern
# must match exactly one file; anything missing or ambiguous is reported and not copied.
# Run from the repository root:  Rscript scripts/10_paper/09_collect_paper_final.R
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
tab <- fread(file.path(cfg$root, "config", "paper_outputs.tsv"), sep = "\t", na.strings = "")
od <- file.path(pa$results_root, "paper_final", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)
old <- list.files(od, full.names = TRUE); unlink(old[!dir.exists(old)])   # start clean each run
man <- rbindlist(lapply(seq_len(nrow(tab)), function(i) {
  pat <- file.path(pa$results_root, gsub("{ver}", ver, tab$source_glob[i], fixed = TRUE))
  hits <- Sys.glob(pat)
  if (grepl("**", pat, fixed = TRUE)) {                       # recursive pattern
    base <- sub("/\\*\\*.*$", "", pat); rx <- glob2rx(basename(pat))
    hits <- list.files(base, pattern = rx, recursive = TRUE, full.names = TRUE)
  }
  hits <- hits[grepl(paste0("/", ver, "/"), hits, fixed = TRUE)]
  st <- if (length(hits) == 1) "copied" else if (!length(hits)) "MISSING" else paste0("AMBIGUOUS (", length(hits), " matches)")
  if (st == "copied") file.copy(hits, file.path(od, tab$paper_name[i]), overwrite = TRUE, copy.date = TRUE)
  data.table(item = tab$item[i], paper_name = tab$paper_name[i], status = st,
             source = if (length(hits) == 1) sub(paste0("^", pa$results_root, "/"), "", hits) else tab$source_glob[i],
             bytes = if (st == "copied") file.size(hits) else NA_real_,
             modified = if (st == "copied") format(file.mtime(hits), "%Y-%m-%d %H:%M") else "",
             md5 = if (st == "copied") unname(tools::md5sum(hits)) else "", note = tab$note[i])
}))
fwrite(man, file.path(od, "MANIFEST.tsv"), sep = "\t")
cat("Callset version:", ver, "| output:", od, "\n\n"); print(man[, .(item, paper_name, status, modified, source)])
cat("\nCopied", sum(man$status == "copied"), "of", nrow(man), "| not copied:", sum(man$status != "copied"), "\n")
