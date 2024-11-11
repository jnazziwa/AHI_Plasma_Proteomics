# An usual 'setup' chunk in R markdown reports

## Chunk hooks
knitr::knit_hooks$set(
  optipng = knitr::hook_optipng   # Ready to activate optipng (Optimize PNG)
)

## Set default chunk options
knitr::opts_chunk$set(
  # echo = FALSE,    # comment out for html.code_folding
  
  cache.path = "../data/interim/",   # cache the results from heavy analyses
  purl = FALSE,      # because most of the code are not to be extracted
  
  results = "asis",   # for HTML
  message = FALSE,    # no message from the code
  warning = FALSE,    # no warning from the code
  
  # # export figures to pdf files
  # dev = c('pdf'), 
  # fig.path = '../reports/figures/',
  
  fig.width  = 8,
  fig.height = 6,
  dpi = 300,
  fig.align = "center",  # pdf, line break between figures
  out.width = if(knitr::is_html_output()) NULL else "70%",  # pdf
  # optipng = "-o1", # reduce file size for 'png' plots
  # optipng = "-o7 -zc9 -zm8 -zs0", # more time but smaller size
  eval.after = "fig.cap"    # add fig. caption within the chunk
)

## Control knitr general options
options(
  # knitr.duplicate.label = "allow",   # avoid error from multiple usage of child
  knitr.kable.NA = "",         # NA is shown as ""
  knitr.table.format = if(knitr::is_html_output()) "html" else "markdown"
)