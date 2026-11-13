# Functions used by the call-set scripts (scripts/01_callset).

# Sample ID from the per-sample PennCNV file path,
# e.g. ".../01_00000000_CK00000000G_A01.penncnv.cnv.txt" -> "CK00000000"
extract_ckid <- function(x) {
  pos <- regexpr("CK[0-9]+G(?:-[A-Za-z0-9]+)*", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  hit <- pos > 0
  out[hit] <- sub("G", "", regmatches(x, pos))
  out
}

# Parse one raw PennCNV call file. Returns NULL for an empty file.
parse_penncnv <- function(f) {
  if (file.size(f) == 0) return(NULL)
  d <- read.table(f, header = FALSE, sep = "", stringsAsFactors = FALSE, fill = TRUE,
                  quote = "", comment.char = "", colClasses = "character")
  if (nrow(d) == 0) return(NULL)
  if (ncol(d) < 8) for (i in (ncol(d) + 1):8) d[[paste0("V", i)]] <- NA_character_
  d <- as.data.table(d[, 1:8])
  setnames(d, c("Region", "NumSNP", "Length", "State_CN", "File", "StartSNP", "EndSNP", "Conf"))
  sc <- tstrsplit(d$State_CN, ",", fixed = TRUE, fill = NA_character_, keep = 1:2)
  rg <- tstrsplit(d$Region, ":|-", fill = NA_character_, keep = 1:3)
  num <- function(v) suppressWarnings(as.numeric(v))
  data.table(
    Region     = d$Region,
    NumSNP     = as.integer(num(sub("numsnp=", "", d$NumSNP))),
    Length     = as.integer(num(gsub(",", "", sub("length=", "", d$Length)))),
    File       = d$File,
    StartSNP   = sub("startsnp=", "", d$StartSNP),
    EndSNP     = sub("endsnp=", "", d$EndSNP),
    Conf       = num(sub("conf=", "", d$Conf)),
    State      = sc[[1]],
    CN         = as.integer(num(sub("cn=", "", sc[[2]]))),
    Chromosome = rg[[1]],
    StartPos   = as.integer(num(rg[[2]])),
    EndPos     = as.integer(num(rg[[3]])))
}

# Merge groups. Input must be sorted by key then start. A call joins the current merged
# call when (start - current end) / (end - current start) <= gap_fraction, where the
# current end is the furthest end so far.
assign_merge_groups <- function(key, start, end, gap_fraction) {
  n <- length(start); grp <- integer(n); g <- 0L
  cur_start <- 0; cur_end <- 0
  for (i in seq_len(n)) {
    if (i == 1L || key[i] != key[i - 1L]) {
      g <- g + 1L; cur_start <- start[i]; cur_end <- end[i]
    } else {
      span <- end[i] - cur_start
      if (span <= 0) span <- 1
      if ((start[i] - cur_end) / span <= gap_fraction) {
        if (end[i] > cur_end) cur_end <- end[i]
      } else {
        g <- g + 1L; cur_start <- start[i]; cur_end <- end[i]
      }
    }
    grp[i] <- g
  }
  grp
}

# Merge fragmented calls within sample (File), chromosome and CN_type.
# A merged call gets its own coordinates, Region and Length, StartSNP of the first fragment,
# EndSNP of the last, NumSNP and Conf combined by the configured rules, and State/CN of the
# fragment with most probes.
merge_fragments <- function(d, gap_fraction, numsnp_rule = "sum", conf_rule = "max") {
  stopifnot(numsnp_rule %in% c("sum", "first"), conf_rule %in% c("max", "sum", "first"))
  d <- copy(d)
  setorderv(d, c("File", "Chromosome", "CN_type", "StartPos"))     # stable sort
  key <- paste(d$File, d$Chromosome, d$CN_type, sep = "\r")
  d[, merge_id := assign_merge_groups(key, StartPos, EndPos, gap_fraction)]
  d[, n_frag := .N, by = merge_id]

  single <- d[n_frag == 1L]
  single[, `:=`(n_merged = 1L, mixed_cn = FALSE)]

  multi <- d[n_frag > 1L]
  merged <- multi[, {
    s <- min(StartPos); e <- max(EndPos)
    i_end <- which.max(EndPos); i_main <- which.max(NumSNP)
    list(Region = paste0(Chromosome[1], ":", s, "-", e),
         NumSNP = if (numsnp_rule == "sum") sum(NumSNP) else NumSNP[1],
         Length = as.integer(e - s + 1L),
         File = File[1], StartSNP = StartSNP[1], EndSNP = EndSNP[i_end],
         Conf = switch(conf_rule, max = max(Conf), sum = sum(Conf), first = Conf[1]),
         State = State[i_main], CN = CN[i_main], Chromosome = Chromosome[1],
         StartPos = s, EndPos = e, CN_type = CN_type[1], n_merged = .N,
         mixed_cn = uniqueN(CN) > 1L)
  }, by = merge_id]

  keep_cols <- c("Region", "NumSNP", "Length", "File", "StartSNP", "EndSNP", "Conf", "State", "CN",
                 "Chromosome", "StartPos", "EndPos", "CN_type", "n_merged", "mixed_cn")
  calls <- rbind(single[, ..keep_cols], merged[, ..keep_cols])
  details <- merge(multi[, .(merge_id, fragment_region = Region, NumSNP, Conf, CN)],
                   merged[, .(merge_id, merged_region = Region, n_merged)], by = "merge_id")
  list(calls = calls, details = details)
}

# Exclusion regions: one region per line as chr:start-end (comments '#' and blank lines skipped).
read_regions <- function(f) {
  l <- trimws(gsub("\r", "", readLines(f, warn = FALSE)))
  l <- l[l != "" & !startsWith(l, "#")]
  tok <- sub("[[:space:]].*$", "", l)
  m <- regmatches(tok, regexec("^(chr)?([0-9A-Za-z]+):([0-9,]+)-([0-9,]+)$", tok))
  bad <- lengths(m) == 0
  if (any(bad)) stop("Could not parse exclusion line(s): ", paste(l[bad], collapse = " | "))
  data.table(chr = toupper(vapply(m, `[`, "", 3)),
             start = as.integer(gsub(",", "", vapply(m, `[`, "", 4))),
             end   = as.integer(gsub(",", "", vapply(m, `[`, "", 5))),
             region = tok)
}

# Merge overlapping regions into consensus regions (as PennCNV scan_region.pl does)
condense_regions <- function(r) {
  r <- r[order(chr, start, end)]
  r[, grp := cumsum(c(1L, (start[-1] > cummax(end)[-.N]))), by = chr]
  r[, .(start = min(start), end = max(end), regions = paste(region, collapse = ";")), by = .(chr, grp)][, grp := NULL][]
}

# For each call, the largest fraction of its length inside any single condensed region.
# A call is excluded when this is >= min_frac (scan_region.pl --minqueryfrac).
max_region_fraction <- function(chr, start, end, regions_condensed) {
  q <- data.table(row = seq_along(start), chr = toupper(sub("^chr", "", chr)),
                  start = as.integer(start), end = as.integer(end))
  y <- copy(regions_condensed)
  setkey(y, chr, start, end)
  ov <- foverlaps(q[!is.na(start) & !is.na(end)], y, by.x = c("chr", "start", "end"),
                  type = "any", nomatch = NULL)
  out <- data.table(row = q$row, frac = 0, hit = NA_character_)
  if (nrow(ov)) {
    ov[, frac := (pmin(end, i.end) - pmax(start, i.start) + 1) / (i.end - i.start + 1)]
    best <- ov[order(row, -frac)][!duplicated(row)]
    out[best, on = "row", `:=`(frac = i.frac, hit = i.regions)]
  }
  out
}
