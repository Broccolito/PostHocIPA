#' Analyze IPA output
#'
#' @param ipa_input_file A string indicating the file location of IPA input gene list
#' @param expression_matrix_file A string indicating the file location of the expression matrix file
#' @param pathway_file A string indicating the file location of IPA output pathway .xls file
#' @param sample_list A vector indicating columns in the expression matrix that the user wants to keep
#' @param sample_names A vector indicating the user-imposed column names
#' @param group_names A vector indicating the grouping of the samples
#' @param pathway_header A boolean value specifying wether the IPA output pathway .xls file has a meta header
#' @examples
#' IPA_PostHoc(
#' ipa_input_file = "mint_flavor.csv",
#' expression_matrix_file = "ecig_count_matrix",
#' pathway_file = "mint_flavor_pathway.xls",
#'
#' sample_list = c(
#'   "Air_JUUL_1Month_noLPS_Female_S1.results",
#'   "Air_JUUL_1Month_noLPS_Female_S2.results",
#'   "Air_JUUL_1Month_noLPS_Female_S3.results",
#'   "Air_JUUL_1Month_noLPS_Female_S4.results",
#'   "MintJUUL_JUUL_1Month_noLPS_Female_S1.results",
#'   "MintJUUL_JUUL_1Month_noLPS_Female_S2.results",
#'   "MintJUUL_JUUL_1Month_noLPS_Female_S3.results",
#'   "MintJUUL_JUUL_1Month_noLPS_Female_S4.results"
#' ),
#'
#' sample_names = c(
#'   "Air1",
#'   "Air2",
#'   "Air3",
#'   "Air4",
#'   "Mint1",
#'   "Mint2",
#'   "Mint3",
#'   "Mint4"
#' ),
#'
#' group_names = c(
#'   "Air",
#'   "Air",
#'   "Air",
#'   "Air",
#'   "Mint",
#'   "Mint",
#'   "Mint",
#'   "Mint"
#' ),
#'
#' pathway_header = TRUE
#' )
#' @export IPA_PostHoc

IPA_PostHoc = function(

  ipa_input_file = "mint_flavor.csv",
  expression_matrix_file = "ecig_count_matrix",
  pathway_file = "mint_flavor_pathway.xls",

  sample_list = c(
    "Air_JUUL_1Month_noLPS_Female_S1.results",
    "Air_JUUL_1Month_noLPS_Female_S2.results",
    "Air_JUUL_1Month_noLPS_Female_S3.results",
    "Air_JUUL_1Month_noLPS_Female_S4.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S1.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S2.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S3.results",
    "MintJUUL_JUUL_1Month_noLPS_Female_S4.results"
  ),

  sample_names = c(
    "Air1",
    "Air2",
    "Air3",
    "Air4",
    "Mint1",
    "Mint2",
    "Mint3",
    "Mint4"
  ),

  group_names = c(
    "Air",
    "Air",
    "Air",
    "Air",
    "Mint",
    "Mint",
    "Mint",
    "Mint"
  ),

  pathway_header = TRUE

){

  sample_list = c("gene_name", sample_list)
  sample_names = c("gene_name", sample_names)

  # Read pathway file
  pathways = read_excel(pathway_file)
  if(pathway_header == TRUE){
    names(pathways) = pathways[1,]
    pathways = pathways[-1,]
  }

  # Read expression matrix
  expression_matrix = read.delim(expression_matrix_file)
  names(expression_matrix)[1] = "gene_name"
  expression_matrix = as.data.frame(expression_matrix)
  expression_matrix = expression_matrix[,(names(expression_matrix) %in% sample_list)]

  # Read IPA input
  ipa_input = read.csv(ipa_input_file)

  # Annotate expression matrix
  annotated_expression_matrix = ipa_input %>%
    select(gene_name, logFC, PValue, symbol) %>%
    inner_join(expression_matrix, by = "gene_name") %>%
    mutate(symbol = gsub("\\s*\\([^\\)]+\\)","",symbol)) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    mutate(FDR = p.adjust(PValue, method = "fdr")) %>%
    select(gene_name, logFC, PValue, FDR, symbol, everything())

  # Generate pathway list
  pathways_list = pathways %>%
    split(.$`Ingenuity Canonical Pathways`) %>%
    map(function(x){
      tibble(symbol = unlist(strsplit(x$Molecules,split = ","))) %>%
        inner_join(annotated_expression_matrix, by = "symbol") %>%
        arrange(PValue)
    })
  empty_index = unlist(map(pathways_list, function(x){dim(x)[1]}))==0
  pathways_list = pathways_list[!empty_index]

  # Generate significant pathway list
  pathways_significant_list = pathways %>%
    split(.$`Ingenuity Canonical Pathways`) %>%
    map(function(x){
      tibble(symbol = unlist(strsplit(x$Molecules,split = ","))) %>%
        inner_join(annotated_expression_matrix, by = "symbol") %>%
        arrange(PValue) %>%
        filter(FDR <= 0.05)
    })
  significant_empty_index = unlist(map(pathways_significant_list, function(x){dim(x)[1]}))==0
  pathways_significant_list = pathways_significant_list[!significant_empty_index]

  # Define heatmap plotting function
  make_heatmap = function(pathway_df, pathway_name){

    cat(paste0("Making heatmap for ", pathway_name, "...\n"))
    x11()

    if(dim(pathway_df)[1]>=2){
      logcounts = pathway_df
      gene_symbol = logcounts$symbol
      logcounts = logcounts %>%
        select(-(symbol:FDR))
      row_mean = apply(logcounts, MARGIN = 1, mean)
      row_sd = apply(logcounts, MARGIN = 1, sd)

      logcounts_zscore = vector()
      for(i in 1:dim(logcounts)[1]){
        row_zscore = (logcounts[i,] - row_mean[i]) / row_sd[i]
        logcounts_zscore = rbind(logcounts_zscore, row_zscore)
      }
      logcounts_zscore = data.frame(logcounts_zscore)
      logcounts_zscore$gene_name = gene_symbol
      logcounts_zscore$row_mean = row_mean
      logcounts_zscore$row_sd = row_sd
      logcounts_zscore = logcounts_zscore %>%
        select(gene_name, row_mean, row_sd, everything())

      heatmap_source_matrix = logcounts_zscore %>%
        select(-(gene_name:row_sd)) %>%
        as.matrix()
      rownames(heatmap_source_matrix) = logcounts_zscore$gene_name
      heatmap_object = heatmap.2(heatmap_source_matrix, Colv = FALSE, Rowv = TRUE,
                                 col = "green2red", trace = "none")
      heatmap_source_matrix = heatmap_object$carpet
      heatmap_source_matrix = t(heatmap_source_matrix)
      heatmap_source_matrix = heatmap_source_matrix %>%
        as.data.frame() %>%
        mutate(gene_names = rownames(heatmap_source_matrix)) %>%
        select(gene_names, everything())
      rownames(heatmap_source_matrix) = NULL

      names(heatmap_source_matrix) = sample_names

      heatmap_source_processed = vector()
      for(j in 1:(dim(heatmap_source_matrix)[2]-1)){
        for(i in 1:dim(heatmap_source_matrix)[1]){
          heatmap_source_fragment = tibble(gene_names = heatmap_source_matrix[i,1],
                                           sample = names(heatmap_source_matrix)[j+1],
                                           intensity = heatmap_source_matrix[i,j+1])
          heatmap_source_processed = rbind.data.frame(heatmap_source_processed,
                                                      heatmap_source_fragment)
        }
      }

      plt = ggplot(data = heatmap_source_processed,
                   aes(x = sample, y = gene_names, fill = intensity)) +
        geom_tile() +
        scale_fill_gradient2(low="green",mid = "black", high="red",midpoint = 0) +
        ylab("") +
        xlab("") +
        labs(fill = "Corrected Z-Score:  ") +
        ggtitle(label = NULL, subtitle = pathway_name) +
        theme_pubclean() +
        theme(axis.text.y=element_text(face="italic"),
              axis.text.x = element_text(angle = 45,hjust = 1),
              text = element_text(size = 15),
              legend.position="top")
      graphics.off()
      return(plt)

    }else{
      pathway_df$pathway = pathway_name
      pathway_df = pathway_df %>%
        select(pathway, symbol:FDR, everything())
      graphics.off()
      return(pathway_df)
    }
  }

  # Plotting heatmaps for all pathways
  pathways_heatmaps = map2(pathways_list,
                           names(pathways_list),
                           make_heatmap)

  # Plotting heatmaps for all pathways with significant genes only
  pathways_significnat_heatmaps = map2(pathways_significant_list,
                                       names(pathways_significant_list),
                                       make_heatmap)

  # Plot boxplots for significant genes
  gene_boxplots = annotated_expression_matrix %>%
    filter(FDR <= 0.05) %>%
    split(.$symbol) %>%
    map(function(x){
      gene_expression = unlist(select(x, -(gene_name:symbol)))
      gene_name = x$symbol
      d = tibble(group_names, gene_expression)
      ggplot(data = d, aes(x = group_names, y = gene_expression)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(aes(fill = group_names), color = "black", shape = 21, size = 3) +
        labs(fill = "") +
        xlab("") +
        ylab("Relative Expression") +
        ggtitle(label = gene_name) +
        stat_compare_means(method = "t.test",
                           comparisons = list(unique(group_names)),
                           symnum.args = list(cutpoints = c(0, 1),
                                              symbols = c("*")),
                           size = 5) +
        theme_pubr() +
        theme(text = element_text(size = 15),
              plot.title = element_text(face = "italic"))
    })

  result = list(
    pathways_heatmaps = pathways_heatmaps,
    pathways_significnat_heatmaps = pathways_significnat_heatmaps,
    gene_boxplots = gene_boxplots
  )

  return(result)

}


