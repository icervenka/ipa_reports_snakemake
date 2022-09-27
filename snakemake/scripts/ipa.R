suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(janitor))
suppressMessages(library(plotly))
suppressMessages(library(dplyr))
rmarkdown::find_pandoc(dir =snakemake@params[["pandoc_path"]])

# functions -------------------------------------------------------------------------------
minmax_scale = function(x) {
  (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}

#TODO deal with end of parsing block to avoid warnings and parsing errors
create_ipa_object = function(filename, parse_strings) {

  conn = file(filename, open = "r")
  flines = readLines(conn)
  parse_empty_lines= grep("^$", flines)

  parse_positions = lapply(names(parse_strings), function(x) {
    start = min(grep(parse_strings[x], flines))
    numrows = min(parse_empty_lines[parse_empty_lines > start]) - start
    return(c(start, numrows-1))
  })
  names(parse_positions) = names(parse_strings)
  ipa_object = vector(mode="list", length=length(parse_strings))

  ipa_object = lapply(names(parse_strings), function(x) {
    df = suppressWarnings(read_delim(filename, delim = "\t", skip = parse_positions[[x]][1], n_max = parse_positions[[x]][2], skip_empty_rows = F))
    df = remove_empty(df, c("rows", "cols"))
    #Filter(function(x) !all(is.na(x)), df)
  })

  names(ipa_object) = names(parse_strings)

  return(ipa_object)
}

extract_canonical = function(ipa_obj, zscore_threshold = 0.01) {
  temp_df = ipa_obj[["canonical"]]

  new_colnames = c("category", "-log10(pvalue)", "zscore", "ratio", "molecules")
  if(length(colnames(temp_df)) != length(new_colnames)) {
    print("Not compatible file specification.")
  }
  colnames(temp_df) = new_colnames

  temp_df$`-log10(pvalue)` = suppressWarnings(as.numeric(temp_df$`-log10(pvalue)`))
  temp_df$zscore = suppressWarnings(as.numeric(temp_df$zscore))
  temp_df$zscore = round(temp_df$zscore, digits = 3)
  temp_df$ratio = suppressWarnings(as.numeric(temp_df$ratio))
  temp_df[is.na(temp_df)] = 0
  temp_df$direction = cut(temp_df$zscore, breaks = c(-Inf, -zscore_threshold, zscore_threshold, Inf), labels = c("negative", "no direction", "positive"))
  temp_df$direction = factor(temp_df$direction, levels = c("no direction", "positive", "negative"))
  temp_df = temp_df[order(temp_df$`-log10(pvalue)`, decreasing = T), ]
  temp_df$category = factor(temp_df$category, levels=temp_df$category)
  temp_df$molecules = gsub(",", ", ", temp_df$molecules)

  return(temp_df)
}

extract_upstream = function(ipa_obj) {
  temp_df = ipa_obj[["upstream"]]
  temp_df = temp_df %>% dplyr::select(-matches("Mechanistic Network"))

  new_colnames = c(analysis = quote(`Analysis`),
                   regulator = quote(`Upstream Regulator`),
                   fc = quote(`Expr Log Ratio`),
                   type = quote(`Molecule Type`),
                   state = quote(`Predicted Activation State`),
                   zscore = quote(`Activation z-score`),
                   bias = quote(`Bias Term`),
                   bias_zscore = quote(`Bias-corrected z-score`),
                   pvalue = quote(`p-value of overlap`),
                   target = quote(`Target molecules in dataset`))

  temp_df = temp_df %>% rename(!!!new_colnames)
  temp_df = transform(temp_df, fc = as.numeric(fc),
                      zscore = as.numeric(zscore),
                      bias = as.numeric(bias),
                      bias_zscore = as.numeric(bias_zscore),
                      pvalue = as.numeric(pvalue))
  temp_df[is.na(temp_df)] = 0
  temp_df$type = factor(temp_df$type, levels=unique(temp_df$type))
  temp_df$target = gsub(",", ", ", temp_df$target)
  temp_df = temp_df %>% dplyr::arrange(pvalue)

  return(temp_df)
}

extract_disease = function(ipa_obj) {
  temp_df = ipa_obj[["disease"]]

  new_colnames = c(category = quote(`Categories`),
                   functions = quote(`Functions`),
                   disease = quote(`Diseases or Functions Annotation`),
                   pvalue = quote(`p-Value`),
                   state = quote(`Predicted Activation State`),
                   zscore = quote(`Activation z-score`),
                   bias_zscore = quote(`Bias-corrected z-score`),
                   molecules = quote(`Molecules`),
                   no_molecules = quote(`# Molecules`))

  temp_df = temp_df %>% rename(!!!new_colnames)
  temp_df = transform(temp_df, zscore = as.numeric(zscore), bias_zscore = as.numeric(bias_zscore),
                      pvalue = as.numeric(pvalue), no_molecules = as.numeric(no_molecules))
  temp_df[is.na(temp_df)] = 0
  temp_df = temp_df[order((temp_df$zscore)), ]
  temp_df$molecules = gsub(",", ", ", temp_df$molecules)
  temp_df = temp_df %>% dplyr::arrange(pvalue)
  return(temp_df)
}

graph_canonical = function(canonical_df, nitems = 20, ylegend = TRUE, select_rows = NULL, exclude_rows = NULL) {
  if(!is.null(select_rows)) {
    canonical_df = canonical_df[select_rows,]
    canonical_df = droplevels(canonical_df)
  }
  if(!is.null(exclude_rows)) {
    canonical_df = canonical_df[-exclude_rows,]
    canonical_df = droplevels(canonical_df)
  }
  p <- ggplot(canonical_df[1:nitems, ], aes(x = category, y = `-log10(pvalue)`, color = direction)) +
    ggtitle("") +
    xlab("") +
    ylab("-log10(P-value)") +
    geom_point(size = 60/nitems) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(canonical_df$category)[1:nitems])) +
    theme_bw() +
    scale_color_manual(name = "z-score", values = volPallete)
  if(ylegend) {
    ggplotly(p + theme(axis.text.y = element_text(size = -0.005*(nitems^2) + 0.05*nitems + 11), legend.position = "left"))
  } else {
    gplotly(p + theme(axis.text.y = element_blank(), axis.title = element_text(size=10), legend.position="none"))
  }

}

graph_upstream = function(df, pval_threshold, zscore_threshold) {

  df$threshold = as.factor(df$pvalue < pval_threshold & abs(df$bias_zscore) > zscore_threshold)

  p = df %>%
    filter(pvalue > 0) %>%
    ggplot(aes(x = bias_zscore, y = -log10(pvalue), color = threshold, text = regulator)) +
    geom_point(size = 3) +
    theme_bw() +
    theme(legend.position='none') +
    scale_color_manual(name = "threshold", values = volPallete)

  ggplotly(p, tooltip = "text")
}

graph_disease = function(df, pval_threshold, zscore_threshold, selected_labels = NULL) {

  df$threshold = as.factor(df$pvalue < pval_threshold & abs(df$bias_zscore) > zscore_threshold)
  p = df %>%
    filter(pvalue > 0) %>%
    ggplot(aes(x = bias_zscore, y = -log10(pvalue), color = threshold, text = paste0("category: ", category, '\n', "disease: ", disease))) +
    geom_point(size = 3) +
    theme_bw() +
    theme(legend.position='none') +
    scale_color_manual(name = "threshold", values = volPallete)

  # p <- figure(data = df %>% filter(pvalue > 0), width = 900, height = 500, legend_location = NULL) %>%
  #   ly_points(bias_zscore, -log10(pvalue),
  #             color = threshold,
  #             hover = list(category, disease),
  #             size=15) %>%
  #   theme_plot(min_border_left = 50, min_border_bottom = 50) %>%
  #   theme_axis(axis_label_standoff = 25)
  ggplotly(p, tooltip = "text")
}

# export = function(ipa_object, filename) {
#   base_filename = file_path_sans_ext(filename)
#   write.table(canonical_df(ipa_object), paste0(base_filename, "_candf.csv"), quote = F, row.names = F, sep = "\t")
#   write.table(disease_df(ipa_object), paste0(base_filename, "_diseasedf.csv"), quote = F, row.names = F, sep = "\t")
#   write.table(upstream_df(ipa_object), paste0(base_filename, "_upstreamdf.csv"), quote = F, row.names = F, sep = "\t")
# }


# run analysis -----------------------------------------------------------------------------
volPallete <- c("#B0B0B0", "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")

filename = snakemake@input[[1]]

parse_strings = list(
  canonical=snakemake@params[["canonical_string"]],
  upstream=snakemake@params[["upstream_string"]],
  disease=snakemake@params[["disease_string"]],
  tox=snakemake@params[["tox_string"]],
  regulator=snakemake@params[["regulator_string"]],
  networks=snakemake@params[["networks_string"]],
  molecules=snakemake@params[["molecules_string"]]
)

# KEEP FOR DEBUGGING FOR NOW
# setwd("~/programming/analysis_pipeline/input/")
# filename = "ipa_alpha1_gfp_nervous.txt"
#canonical_nitems = 20

# parse_strings = list(
#   canonical="Canonical Pathways for",
#   upstream="Upstream Regulators",
#   disease="Diseases and Bio Functions",
#   tox="Tox Functions",
#   regulator="Regulator Effects",
#   networks="Networks",
#   molecules=""
# )

ipa_object = create_ipa_object(filename, parse_strings = parse_strings)

canonical_df = extract_canonical(ipa_object)
canonical_graph = graph_canonical(canonical_df, nitems = snakemake@params[["graph_nitems"]])

upstream_df = extract_upstream(ipa_object)
upstream_graph = graph_upstream(upstream_df,
                                pval_threshold = snakemake@params[["upstream_pval_threshold"]],
                                zscore_threshold = snakemake@params[["upstream_zscore_threshold"]])

disease_df = extract_disease(ipa_object)
disease_graph = graph_disease(disease_df,
                              pval_threshold = snakemake@params[["disease_pval_threshold"]],
                              zscore_threshold = snakemake@params[["disease_zscore_threshold"]])

write.table(canonical_df, snakemake@output[[2]], sep = '\t', quote = F, row.names = F)
write.table(upstream_df, snakemake@output[[3]], sep = '\t', quote = F, row.names = F)
write.table(disease_df, snakemake@output[[4]], sep = '\t', quote = F, row.names = F)

knitr_output_options = list(
  mathjax = NULL,
  self_contained = snakemake@params[["self_contained"]]#,
  # lib_dir = paste0(outdir, "/libs")
)

rmarkdown::render("snakemake/scripts/ipa.Rmd",
                  output_file = paste0("../../", snakemake@output[[1]]),
                  output_format = "html_document",
                  output_options = knitr_output_options)
