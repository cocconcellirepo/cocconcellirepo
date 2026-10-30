# Runs 02_tCRE_embedding.R on all the timepoints, one after the other
f <- "04_tCRE_SCAFE/scripts/02_tCRE_embedding.R"
for (tp_run in 0:4) {
  l <- readLines(f)
  stopifnot(sum(grepl("^tp <- 0 ", l)) == 1)
  l <- sub("^tp <- 0 .*", paste0("tp <- ", tp_run), l)
  message("##### T", tp_run, " start ", Sys.time())
  source(textConnection(l), echo = TRUE, max.deparse.length = 200, local = FALSE)
  message("##### T", tp_run, " end ", Sys.time())
  rm(list = setdiff(ls(), c("f", "tp_run")))
  invisible(gc())
}
message("##### ALL DONE ", Sys.time())
