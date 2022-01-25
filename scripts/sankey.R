#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(plotly))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(readxl))
suppressMessages(library(stringr))

parse_sheet <- function(excel_sheet, excel_file) {
  df <- readxl::read_excel(excel_file, sheet = excel_sheet)
  df <- df[, c('...1', 'dist mean', 'significant')]
  df['source'] <- excel_sheet
  colnames(df) <- c('target', 'value', 'significant', 'source')
  
  return (df)
}

parse_excel <- function(excel_file, only_sig = TRUE, skip = NULL) {
    
    sheets <- excel_sheets(excel_file)[-1] # Skip Dashboard
    dfs <- lapply(sheets, parse_sheet, excel_file=excel_file)
    dfs <- do.call("rbind", dfs)

    if (only_sig) {
      dfs <- dfs[dfs$significant, ]
      dfs$color <- "#cccccc"
    } else {
      dfs$color <- mapvalues(dfs$significant, from=c(TRUE, FALSE), to=c("#4c92c3","#cccccc"))
    }
    
    if (!is.null(skip)) {
      clusters <- str_split(skip, ",", simplify = TRUE)
      dfs <- dfs[!(dfs$source %in% clusters), ]
      dfs <- dfs[!(dfs$target %in% clusters), ]
    }
    
    labels <- unique(c(dfs$target, dfs$source))
    dfs$source_num <- match(dfs$source, labels) -1
    dfs$target_num <- match(dfs$target, labels) -1
    dfs$updated_value <- abs(dfs$value - max(dfs$value))
    
    data <- NULL
    data$df <- dfs
    data$labels <- labels
    
    return (data)
}

parser <- OptionParser()
parser <- add_option(parser, "--excel", type="character", default=NULL, help="CAT excel result", metavar="character")
parser <- add_option(parser, "--output", type="character", default="./figures/", help="output file name [default= %default]", metavar="character")
parser <- add_option(parser, "--sig_only", action="store_true", default=TRUE, help="Plots only significant connections")
parser <- add_option(parser, "--skip", type="character", default=NULL, help="Clusters names to skip separated by comma", metavar="character")
args <- parse_args(parser)

if (is.null(args$excel)) {
  stop("CAT excel results not provided!")
}

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive=TRUE)
}

cat_res <- parse_excel(args$excel, args$sig_only, args$skip)

fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  node = list(
    label = cat_res$labels,
    pad = 15,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = cat_res$df$source_num,
    target = cat_res$df$target_num,
    color = cat_res$df$color,
    value =  cat_res$df$value
  )
) %>% layout(font = list(size = 20))

htmlwidgets::saveWidget(
  as_widget(fig), 
  paste0(args$output, "/", str_replace(basename(args$excel), '.xlsx', ''), '.html')
)

# fig.write_image("surface-plot.svg", engine="kaleido")
# htmlwidgets::saveWidget(fig, paste0(args$ouput, "/index.html"))
