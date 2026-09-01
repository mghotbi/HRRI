#' Cross-domain association or graph-topology summary
#'
#' @description Retains the legacy function name. Correlation magnitude is
#' association, not electron-transfer encounter probability or measured alpha.
#' Network summaries concern unweighted topology, not biochemical connectivity.
#' @param res RRI result; graph method uses meta$graph.
#' @param method cross_domain_magnitude or network.
#' @param per_group Compute association by group.
#' @param group_cols Required grouping columns when per_group=TRUE.
#' @return Score, association coefficients and method/provenance information.
#' @importFrom igraph vcount ecount edge_attr is_directed simplify global_efficiency transitivity centr_degree
#' @examples
#' \dontrun{
#'   sim <- simulate_redox_holobiont(seed = 1)
#'   res <- rri_pipeline(soil = sim$Eh_stability, plant = sim$ROS_flux)
#'   rri_connectivity_score(res)
#' }
#' @export
rri_connectivity_score <- function(res, method = c("cross_domain_magnitude", "network"),
  per_group = FALSE, group_cols = NULL) {
  method <- match.arg(method)
  if (method == "network") {
    if (per_group) stop("Network per_group requires separately supplied group networks.")
    if (!requireNamespace("igraph", quietly = TRUE)) stop("network requires igraph.")
    g <- res$meta$graph
    if (!inherits(g, "igraph")) stop("A single igraph object is required; no correlation fallback is applied.")
    if (igraph::vcount(g) < 2L || igraph::ecount(g) == 0L)
      stop("Network requires at least two vertices and one edge.")
    gw <- igraph::edge_attr(g, "weight")
    if (igraph::is_directed(g) || (!is.null(gw) &&
        (any(!is.finite(gw)) || any(gw < 0))))
      stop("Use an undirected graph without negative weights.")
    g <- igraph::simplify(g)
    # All metrics use the same unweighted topology. No stochastic community detection.
    eff <- igraph::global_efficiency(g, directed = FALSE, weights = NA)
    clus <- igraph::transitivity(g, type = "average", isolates = "zero")
    cent <- igraph::centr_degree(g)$centralization
    metrics <- c(global_efficiency = eff, clustering = clus, degree_centralization = cent)
    return(list(connectivity_score = eff, method_used = "unweighted_global_efficiency",
      pairwise_abs_cor = NULL, network_metrics = metrics,
      interpretation = "Graph topology only; not physical electron-transfer connectivity"))
  }
  df <- as.data.frame(res$row_scores)
  dom <- intersect(c("Physio", "Soil", "Micro"), names(df))
  pairs <- list(Physio_Soil = c("Physio", "Soil"), Physio_Micro = c("Physio", "Micro"),
                Soil_Micro = c("Soil", "Micro"))
  calc <- function(d) vapply(pairs, function(p) {
    if (!all(p %in% dom)) return(NA_real_)
    abs(.rri_cor(d[[p[1]]], d[[p[2]]], min_n = 4L))
  }, numeric(1))
  if (!per_group) {
    r <- calc(df)
    return(list(connectivity_score = .rri_mean(r), method_used = "cross_domain_magnitude",
      pairwise_abs_cor = r, n_pairs_used = sum(is.finite(r)), network_metrics = NULL,
      interpretation = "Descriptive association; shared forcing and repeated measures can explain correlation"))
  }
  if (!length(group_cols)) stop("Provide group_cols for per_group=TRUE.")
  key <- .rri_key(df, group_cols)
  rows <- lapply(unique(key), function(k) {
    d <- df[key == k, , drop = FALSE]; r <- calc(d)
    cbind(d[1, group_cols, drop = FALSE], data.frame(connectivity_score = .rri_mean(r),
      n_pairs_used = sum(is.finite(r))), as.data.frame(as.list(r)))
  })
  list(connectivity_score = do.call(rbind, rows), method_used = "cross_domain_magnitude_per_group",
    pairwise_abs_cor = NULL, network_metrics = NULL)
}
