#' Plot genome-level synteny
#'
#' @param data Output of [prepare_mcscanx_data()].
#' @param min_pairs Minimum number of anchor pairs required to draw a block.
#' @param min_block_width Minimum genomic width required on both genomes for a block to be drawn.
#' @param genome_order Optional character vector giving the plotting order of genome IDs.
#' @param chr_order Optional named list of chromosome orders. Names should match genome IDs.
#' @param color_by One of `"n_anchors"`, `"block_id"`, or `"genome_pair"`.
#' @param chr_alpha Alpha for chromosome bars.
#' @param block_alpha Alpha for synteny links.
#' @param chr_color Fill color for chromosome bars.
#' @param top_y Y position of the upper genome.
#' @param bottom_y Y position of the lower genome.
#' @param bar_height Height of chromosome bars.
#' @param curvature_height Height of bezier control points.
#' @param show_chr_label Logical, whether to draw chromosome labels.
#' @param show_genome_label Logical, whether to draw genome labels.
#'
#' @returns A `ggplot` object.
#' @export
#'
#' @examples NULL
plot_Synteny <- function(data,
                         min_pairs = 5,
                         min_block_width = 0,
                         genome_order = NULL,
                         chr_order = NULL,
                         color_by = c("n_anchors", "block_id", "genome_pair"),
                         chr_alpha = 0.9,
                         block_alpha = 0.35,
                         chr_color = "grey30",
                         top_y = 1,
                         bottom_y = 0,
                         bar_height = 0.08,
                         curvature_height = 0.35,
                         show_chr_label = TRUE,
                         show_genome_label = TRUE) {
  data <- validate_ggsynteny_data(data)
  color_by <- rlang::arg_match(color_by)

  genomes <- resolve_genome_order(data$genomes, genome_order)
  if (nrow(genomes) != 2) {
    stop("plot_Synteny() currently supports exactly two genomes.", call. = FALSE)
  }

  chromosomes <- apply_chr_order(data$chromosomes, chr_order)
  chromosomes <- recompute_plot_offsets(chromosomes, data$metadata$gap_size %||% 1e6)

  genes <- map_genes_to_linear(data$genes, chromosomes)
  anchors <- rebuild_anchor_linear_positions(data$anchors, genes)
  blocks <- build_blocks_from_anchors(anchors) |>
    compute_block_layout(chromosomes) |>
    dplyr::filter(
      .data$n_anchors >= min_pairs,
      (.data$end_a - .data$start_a) >= min_block_width,
      (.data$end_b - .data$start_b) >= min_block_width
    )

  if (nrow(blocks) == 0) {
    stop("No synteny blocks passed the current filters.", call. = FALSE)
  }

  chromosome_tracks <- chromosomes |>
    dplyr::left_join(genomes |> dplyr::select(genome_id, plot_order), by = "genome_id") |>
    dplyr::mutate(
      y = dplyr::if_else(.data$plot_order == 1L, top_y, bottom_y),
      ymin = .data$y - bar_height / 2,
      ymax = .data$y + bar_height / 2
    )

  blocks <- blocks |>
    dplyr::mutate(
      genome_pair = stringr::str_c(.data$genome_a, .data$genome_b, sep = " vs "),
      color_value = if (color_by == "n_anchors") {
        as.character(.data$n_anchors)
      } else if (color_by == "block_id") {
        .data$block_id
      } else {
        .data$genome_pair
      }
    )

  block_paths <- blocks |>
    dplyr::mutate(
      x = .data$global_mid_a,
      y = top_y,
      group = .data$block_id
    ) |>
    dplyr::select(block_id, n_anchors, color_value, x, y, group) |>
    dplyr::bind_rows(
      blocks |>
        dplyr::mutate(
          x = (.data$global_mid_a + .data$global_mid_b) / 2,
          y = (top_y + bottom_y) / 2 + curvature_height,
          group = .data$block_id
        ) |>
        dplyr::select(block_id, n_anchors, color_value, x, y, group),
      blocks |>
        dplyr::mutate(
          x = .data$global_mid_b,
          y = bottom_y,
          group = .data$block_id
        ) |>
        dplyr::select(block_id, n_anchors, color_value, x, y, group)
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = chromosome_tracks,
      ggplot2::aes(
        xmin = .data$global_start,
        xmax = .data$global_end,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      fill = chr_color,
      alpha = chr_alpha
    ) +
    ggforce::geom_bezier(
      data = block_paths,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        group = .data$group,
        color = .data$color_value
      ),
      alpha = block_alpha,
      linewidth = 0.7,
      show.legend = TRUE
    ) +
    ggplot2::labs(x = NULL, y = NULL, color = color_by) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  p <- add_block_color_scale(p, blocks, color_by)

  if (show_chr_label) {
    label_df <- chromosome_tracks |>
      dplyr::mutate(
        label_x = (.data$global_start + .data$global_end) / 2,
        label_y = dplyr::if_else(.data$plot_order == 1L, top_y + 0.1, bottom_y - 0.1)
      )

    p <- p +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$chr_label),
        size = 3
      )
  }

  if (show_genome_label) {
    genome_label_df <- genomes |>
      dplyr::mutate(
        x = min(chromosome_tracks$global_start),
        y = dplyr::if_else(.data$plot_order == 1L, top_y + 0.22, bottom_y - 0.22)
      )

    p <- p +
      ggplot2::geom_text(
        data = genome_label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$genome_label),
        hjust = 0,
        fontface = "bold",
        size = 4
      )
  }

  p
}

resolve_genome_order <- function(genomes, genome_order = NULL) {
  genomes <- genomes |>
    dplyr::arrange(.data$order)

  if (is.null(genome_order)) {
    return(dplyr::mutate(genomes, plot_order = dplyr::row_number()))
  }

  if (!setequal(genome_order, genomes$genome_id)) {
    stop("`genome_order` must contain exactly the genome IDs in `data$genomes`.", call. = FALSE)
  }

  genomes |>
    dplyr::mutate(plot_order = match(.data$genome_id, genome_order)) |>
    dplyr::arrange(.data$plot_order)
}

apply_chr_order <- function(chromosomes, chr_order = NULL) {
  if (is.null(chr_order)) {
    return(chromosomes)
  }

  if (!is.list(chr_order) || is.null(names(chr_order))) {
    stop("`chr_order` must be a named list keyed by genome ID.", call. = FALSE)
  }

  chromosomes |>
    dplyr::group_by(.data$genome_id) |>
    dplyr::group_modify(function(.x, .y) {
      genome_id <- .y$genome_id[[1]]
      requested_order <- chr_order[[genome_id]]

      if (is.null(requested_order)) {
        return(.x)
      }

      if (!setequal(requested_order, .x$chr_id)) {
        stop(
          "Chromosome order for genome `", genome_id, "` must contain exactly its chromosome IDs.",
          call. = FALSE
        )
      }

      .x |>
        dplyr::mutate(chr_order = match(.data$chr_id, requested_order)) |>
        dplyr::arrange(.data$chr_order)
    }) |>
    dplyr::ungroup()
}

recompute_plot_offsets <- function(chromosomes, gap_size) {
  chromosomes |>
    dplyr::group_by(.data$genome_id) |>
    dplyr::arrange(.data$chr_order, .by_group = TRUE) |>
    dplyr::mutate(
      offset = dplyr::lag(cumsum(.data$length + gap_size), default = 0),
      global_start = .data$offset,
      global_end = .data$offset + .data$length
    ) |>
    dplyr::ungroup()
}

rebuild_anchor_linear_positions <- function(anchors, genes) {
  genes_a <- genes |>
    dplyr::select(
      genome_a = genome_id,
      gene_a = gene_id,
      chr_a = chr_id,
      start_a = start,
      end_a = end,
      global_start_a = global_start,
      global_end_a = global_end
    )

  genes_b <- genes |>
    dplyr::select(
      genome_b = genome_id,
      gene_b = gene_id,
      chr_b = chr_id,
      start_b = start,
      end_b = end,
      global_start_b = global_start,
      global_end_b = global_end
    )

  anchors |>
    dplyr::select(-dplyr::any_of(c(
      "chr_a", "start_a", "end_a", "global_start_a", "global_end_a",
      "chr_b", "start_b", "end_b", "global_start_b", "global_end_b"
    ))) |>
    dplyr::left_join(genes_a, by = c("genome_a", "gene_a")) |>
    dplyr::left_join(genes_b, by = c("genome_b", "gene_b"))
}

add_block_color_scale <- function(plot_obj, blocks, color_by) {
  if (color_by == "n_anchors") {
    values <- sort(unique(blocks$n_anchors))
    return(
      plot_obj +
        ggplot2::scale_color_manual(
          values = grDevices::hcl.colors(length(values), "Blues 3"),
          breaks = as.character(values),
          labels = values
        )
    )
  }

  values <- sort(unique(blocks$color_value))
  plot_obj +
    ggplot2::scale_color_manual(
      values = stats::setNames(grDevices::hcl.colors(length(values), "Dynamic"), values)
    )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
