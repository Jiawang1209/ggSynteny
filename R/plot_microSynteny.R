#' Plot micro-synteny for a selected block
#'
#' @param data A `ggsynteny_data` object.
#' @param block_id Character block identifier to plot.
#' @param flank_genes Number of genes to include on each side of the selected block.
#' @param top_y Y position of the upper genome track.
#' @param bottom_y Y position of the lower genome track.
#' @param track_height Height of the gene tracks.
#' @param link_alpha Alpha for anchor links.
#' @param gene_alpha Alpha for gene arrows.
#' @param gene_color_by One of `"anchor_status"` or `"genome"`.
#' @param color_by One of `"block_id"` or `"genome_pair"`.
#' @param anchor_gene_color Color for anchor genes when `gene_color_by = "anchor_status"`.
#' @param flank_gene_color Color for non-anchor flank genes when `gene_color_by = "anchor_status"`.
#' @param highlight_block Logical, whether to highlight the selected block region on each genome.
#' @param block_fill Fill color for highlighted block regions.
#' @param block_fill_alpha Alpha for highlighted block regions.
#' @param show_features Logical, whether to draw exon/CDS features when available.
#' @param feature_types Character vector of feature types to draw when available.
#' @param transcript_mode One of `"longest_transcript"`, `"longest_cds"`, `"canonical"`,
#'   `"all"`, or `"custom"`.
#' @param transcript_ids Optional named character vector mapping `gene_id` to `transcript_id`
#'   when `transcript_mode = "custom"`.
#' @param show_gene_label Logical, whether to label genes participating in anchors.
#'
#' @returns A `ggplot` object.
#' @export
#'
#' @examples NULL
plot_microSynteny <- function(data,
                              block_id,
                              flank_genes = 1,
                              top_y = 1,
                              bottom_y = 0,
                              track_height = 0.16,
                              link_alpha = 0.35,
                              gene_alpha = 0.9,
                              gene_color_by = c("anchor_status", "genome"),
                              color_by = c("block_id", "genome_pair"),
                              anchor_gene_color = "#d95f02",
                              flank_gene_color = "grey35",
                              highlight_block = TRUE,
                              block_fill = "#fee8c8",
                              block_fill_alpha = 0.5,
                              show_features = TRUE,
                              feature_types = c("exon", "CDS"),
                              transcript_mode = "longest_transcript",
                              transcript_ids = NULL,
                              show_gene_label = FALSE) {
  data <- validate_ggsynteny_data(data)
  gene_color_by <- rlang::arg_match(gene_color_by)
  color_by <- rlang::arg_match(color_by)
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  micro_data <- prepare_micro_synteny_data(
    data = data,
    block_id = block_id,
    flank_genes = flank_genes,
    transcript_mode = transcript_mode,
    transcript_ids = transcript_ids
  )

  genes_plot <- micro_data$genes |>
    dplyr::mutate(
      y = dplyr::if_else(.data$genome_id == micro_data$genomes$genome_id[[1]], top_y, bottom_y),
      segment_start = dplyr::if_else(.data$strand == "-", .data$plot_end, .data$plot_start),
      segment_end = dplyr::if_else(.data$strand == "-", .data$plot_start, .data$plot_end),
      xmin = pmin(.data$plot_start, .data$plot_end),
      xmax = pmax(.data$plot_start, .data$plot_end),
      xmid = (.data$xmin + .data$xmax) / 2,
      yend = .data$y,
      ymin = .data$y - track_height / 2,
      ymax = .data$y + track_height / 2,
      gene_color_value = if (gene_color_by == "anchor_status") {
        dplyr::if_else(.data$is_anchor_gene, "anchor", "flank")
      } else {
        .data$genome_id
      }
    )

  block_regions <- micro_data$block_regions |>
    dplyr::mutate(
      y = dplyr::if_else(.data$genome_id == micro_data$genomes$genome_id[[1]], top_y, bottom_y),
      ymin = .data$y - track_height / 1.6,
      ymax = .data$y + track_height / 1.6
    )

  anchor_plot <- micro_data$anchors |>
    dplyr::mutate(
      x = .data$plot_mid_a,
      y = top_y,
      group = .data$anchor_id,
      color_value = if (color_by == "block_id") .data$block_id else .data$genome_pair
    ) |>
    dplyr::select(anchor_id, color_value, x, y, group) |>
    dplyr::bind_rows(
      micro_data$anchors |>
        dplyr::mutate(
          x = (.data$plot_mid_a + .data$plot_mid_b) / 2,
          y = (top_y + bottom_y) / 2 + 0.25,
          group = .data$anchor_id,
          color_value = if (color_by == "block_id") .data$block_id else .data$genome_pair
        ) |>
        dplyr::select(anchor_id, color_value, x, y, group),
      micro_data$anchors |>
        dplyr::mutate(
          x = .data$plot_mid_b,
          y = bottom_y,
          group = .data$anchor_id,
          color_value = if (color_by == "block_id") .data$block_id else .data$genome_pair
        ) |>
        dplyr::select(anchor_id, color_value, x, y, group)
    )

  feature_plot <- micro_data$features |>
    dplyr::filter(.data$feature_type %in% .env$feature_types) |>
    dplyr::mutate(
      y = dplyr::if_else(.data$genome_id == micro_data$genomes$genome_id[[1]], top_y, bottom_y),
      ymin = .data$y - track_height / 3.2,
      ymax = .data$y + track_height / 3.2
    )

  p <- ggplot2::ggplot()

  if (highlight_block) {
    p <- p +
      ggplot2::geom_rect(
        data = block_regions,
        ggplot2::aes(
          xmin = .data$plot_start,
          xmax = .data$plot_end,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        inherit.aes = FALSE,
        fill = block_fill,
        alpha = block_fill_alpha,
        color = NA
      )
  }

  p <- p +
    ggplot2::geom_rect(
      data = genes_plot,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax,
        fill = .data$gene_color_value
      ),
      color = NA,
      alpha = gene_alpha
    ) +
    ggplot2::geom_segment(
      data = genes_plot,
      ggplot2::aes(
        x = .data$segment_start,
        xend = .data$segment_end,
        y = .data$y,
        yend = .data$y
      ),
      linewidth = 0.9,
      alpha = gene_alpha,
      color = "grey10",
      arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed")
    ) +
    ggforce::geom_bezier(
      data = anchor_plot,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$group,
        color = .data$color_value
      ),
      alpha = link_alpha,
      linewidth = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::labs(x = NULL, y = NULL, color = color_by, fill = "Genes") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  if (show_features && nrow(feature_plot) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = feature_plot,
        ggplot2::aes(
          xmin = .data$plot_start,
          xmax = .data$plot_end,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        inherit.aes = FALSE,
        fill = "grey10",
        alpha = 0.9,
        color = NA
      )
  }

  values <- sort(unique(anchor_plot$color_value))
  p <- p +
    ggplot2::scale_color_manual(
      values = stats::setNames(grDevices::hcl.colors(length(values), "Dynamic"), values)
    )

  p <- add_micro_gene_color_scale(
    plot_obj = p,
    genes_plot = genes_plot,
    gene_color_by = gene_color_by,
    anchor_gene_color = anchor_gene_color,
    flank_gene_color = flank_gene_color
  )

  p <- p +
    ggplot2::geom_text(
      data = micro_data$genomes |>
        dplyr::mutate(
          x = min(genes_plot$plot_start),
          y = c(top_y + 0.18, bottom_y - 0.18)
        ),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$genome_label),
      hjust = 0,
      fontface = "bold",
      size = 4
    )

  if (show_gene_label) {
    p <- p +
      ggplot2::geom_text(
        data = genes_plot |>
          dplyr::filter(.data$is_anchor_gene),
        ggplot2::aes(
          x = .data$xmid,
          y = dplyr::if_else(.data$genome_id == micro_data$genomes$genome_id[[1]], top_y + track_height, bottom_y - track_height),
          label = dplyr::coalesce(.data$gene_name, .data$gene_id)
        ),
        size = 2.7
      )
  }

  p
}

#' Prepare block-level data for micro-synteny plotting
#'
#' @param data A `ggsynteny_data` object.
#' @param block_id Character block identifier.
#' @param flank_genes Number of genes to include on each side of the block.
#' @param transcript_mode One of `"longest_transcript"`, `"longest_cds"`, `"canonical"`,
#'   `"all"`, or `"custom"`.
#' @param transcript_ids Optional named character vector mapping `gene_id` to `transcript_id`
#'   when `transcript_mode = "custom"`.
#'
#' @returns A list with `genomes`, `genes`, `anchors`, `features`, and `block`.
#' @export
#'
#' @examples NULL
prepare_micro_synteny_data <- function(data,
                                       block_id,
                                       flank_genes = 1,
                                       transcript_mode = "longest_transcript",
                                       transcript_ids = NULL) {
  data <- validate_ggsynteny_data(data)
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  block_tbl <- data$blocks |>
    dplyr::filter(.data$block_id == .env$block_id)

  if (nrow(block_tbl) != 1) {
    stop("`block_id` must identify exactly one block.", call. = FALSE)
  }

  anchor_tbl <- data$anchors |>
    dplyr::filter(.data$block_id == .env$block_id) |>
    dplyr::mutate(genome_pair = stringr::str_c(.data$genome_a, .data$genome_b, sep = " vs "))

  genes_a <- subset_block_genes(
    genes = data$genes,
    genome_id = block_tbl$genome_a[[1]],
    chr_id = block_tbl$chr_a[[1]],
    start = block_tbl$start_a[[1]],
    end = block_tbl$end_a[[1]],
    flank_genes = flank_genes
  )

  genes_b <- subset_block_genes(
    genes = data$genes,
    genome_id = block_tbl$genome_b[[1]],
    chr_id = block_tbl$chr_b[[1]],
    start = block_tbl$start_b[[1]],
    end = block_tbl$end_b[[1]],
    flank_genes = flank_genes
  )

  genes_plot <- dplyr::bind_rows(genes_a, genes_b)

  genes_plot <- genes_plot |>
    dplyr::group_by(.data$genome_id) |>
    dplyr::mutate(
      region_start = min(.data$start),
      plot_start = .data$start - .data$region_start,
      plot_end = .data$end - .data$region_start,
      plot_mid = (.data$plot_start + .data$plot_end) / 2
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      is_anchor_gene = .data$gene_id %in% c(anchor_tbl$gene_a, anchor_tbl$gene_b)
    )

  genes_a_lookup <- genes_plot |>
    dplyr::filter(.data$genome_id == block_tbl$genome_a[[1]]) |>
    dplyr::select(gene_a = gene_id, plot_mid_a = plot_mid)

  genes_b_lookup <- genes_plot |>
    dplyr::filter(.data$genome_id == block_tbl$genome_b[[1]]) |>
    dplyr::select(gene_b = gene_id, plot_mid_b = plot_mid)

  anchor_plot <- anchor_tbl |>
    dplyr::left_join(genes_a_lookup, by = "gene_a") |>
    dplyr::left_join(genes_b_lookup, by = "gene_b")

  block_regions <- tibble::tibble(
    genome_id = c(block_tbl$genome_a[[1]], block_tbl$genome_b[[1]]),
    plot_start = c(
      block_tbl$start_a[[1]] - min(genes_a$start),
      block_tbl$start_b[[1]] - min(genes_b$start)
    ),
    plot_end = c(
      block_tbl$end_a[[1]] - min(genes_a$start),
      block_tbl$end_b[[1]] - min(genes_b$start)
    )
  )

  feature_tbl <- build_micro_feature_table(
    features = if ("features" %in% names(data)) data$features else tibble::tibble(),
    genes_plot = genes_plot,
    transcripts = if ("transcripts" %in% names(data)) data$transcripts else tibble::tibble(),
    transcript_mode = transcript_mode,
    transcript_ids = transcript_ids
  )

  list(
    genomes = data$genomes |>
      dplyr::filter(.data$genome_id %in% c(block_tbl$genome_a[[1]], block_tbl$genome_b[[1]])) |>
      dplyr::arrange(.data$order),
    genes = genes_plot,
    anchors = anchor_plot,
    features = feature_tbl,
    block = block_tbl,
    block_regions = block_regions
  )
}

subset_block_genes <- function(genes, genome_id, chr_id, start, end, flank_genes = 1) {
  genes_chr <- genes |>
    dplyr::filter(.data$genome_id == .env$genome_id, .data$chr_id == .env$chr_id) |>
    dplyr::arrange(.data$gene_order)

  block_genes <- genes_chr |>
    dplyr::filter(.data$start <= .env$end, .data$end >= .env$start)

  if (nrow(block_genes) == 0) {
    stop("No genes overlap the selected block region.", call. = FALSE)
  }

  lower_order <- max(min(block_genes$gene_order) - flank_genes, 1)
  upper_order <- min(max(block_genes$gene_order) + flank_genes, max(genes_chr$gene_order))

  genes_chr |>
    dplyr::filter(.data$gene_order >= lower_order, .data$gene_order <= upper_order)
}

add_micro_gene_color_scale <- function(plot_obj,
                                       genes_plot,
                                       gene_color_by,
                                       anchor_gene_color,
                                       flank_gene_color) {
  if (gene_color_by == "anchor_status") {
    return(
      plot_obj +
        ggplot2::scale_fill_manual(
          values = c(anchor = anchor_gene_color, flank = flank_gene_color)
        )
    )
  }

  values <- sort(unique(genes_plot$gene_color_value))
  plot_obj +
    ggplot2::scale_fill_manual(
      values = stats::setNames(grDevices::hcl.colors(length(values), "Set 2"), values)
    )
}

build_micro_feature_table <- function(features,
                                      genes_plot,
                                      transcripts = tibble::tibble(),
                                      transcript_mode = "longest_transcript",
                                      transcript_ids = NULL) {
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  if (nrow(features) == 0) {
    return(empty_micro_feature_table())
  }

  selected_transcripts <- resolve_micro_transcripts(
    genes_plot = genes_plot,
    transcripts = transcripts,
    transcript_mode = transcript_mode,
    transcript_ids = transcript_ids
  )

  genes_lookup <- genes_plot |>
    dplyr::select(
      genome_id,
      gene_id,
      region_start,
      gene_plot_start = plot_start,
      gene_plot_end = plot_end
    )

  features |>
    dplyr::semi_join(selected_transcripts, by = c("gene_id", "transcript_id")) |>
    dplyr::inner_join(genes_lookup, by = c("genome_id", "gene_id")) |>
    dplyr::mutate(
      plot_start = .data$start - .data$region_start,
      plot_end = .data$end - .data$region_start
    ) |>
    dplyr::filter(
      .data$plot_end >= .data$gene_plot_start,
      .data$plot_start <= .data$gene_plot_end
    )
}

empty_micro_feature_table <- function() {
  tibble::tibble(
    genome_id = character(),
    gene_id = character(),
    transcript_id = character(),
    chr_id = character(),
    feature_type = character(),
    start = numeric(),
    end = numeric(),
    strand = character(),
    region_start = numeric(),
    gene_plot_start = numeric(),
    gene_plot_end = numeric(),
    plot_start = numeric(),
    plot_end = numeric()
  )
}

resolve_micro_transcripts <- function(genes_plot,
                                      transcripts,
                                      transcript_mode = "longest_transcript",
                                      transcript_ids = NULL) {
  transcript_mode <- normalize_transcript_mode(transcript_mode)

  if (nrow(transcripts) == 0) {
    return(tibble::tibble(gene_id = character(), transcript_id = character()))
  }

  transcripts <- transcripts |>
    dplyr::semi_join(genes_plot |> dplyr::select(gene_id), by = "gene_id")

  select_transcripts(
    transcripts = transcripts,
    transcript_mode = transcript_mode,
    transcript_ids = transcript_ids
  ) |>
    dplyr::select(gene_id, transcript_id)
}
