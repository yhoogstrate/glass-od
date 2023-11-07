#!/usr/bin/env R

format_subtitle <- function(prefix) {
  prefix <- trimws(prefix)
  
  basename <- basename(rstudioapi::getSourceEditorContext()$path)
  git_rev <- trimws(system("git rev-parse --short HEAD", intern=T))
  cur_exec <- Sys.time()
  
  out <- paste0(basename, "  - git:", git_rev," - ", Sys.time())
  
  if(nchar(prefix) > 0) {
    out <- paste0(prefix, " [",out,"]")
  }

  return(out)
}


