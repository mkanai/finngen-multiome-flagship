library(dplyr)
source(here::here("R/const.R"))

df.gex.n_peers = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/optimize_n_peer/integrated_gex_batch1_5.fgid.qc.mean.n_peer.txt"
)
df.atac.n_peers = rgsutil::read_gsfile(
  "gs://expansion_areas/multiome/batch1_5/optimize_n_peer/integrated_atac_batch1_5.fgid.sum.n_peer.txt"
)

df.n_peers = dplyr::bind_rows(
  dplyr::transmute(
    df.gex.n_peers,
    QTL = "eQTL",
    type = "PEER",
    cell_type = cell_type,
    n_factors = local_greedy
  ),
  dplyr::transmute(
    df.atac.n_peers,
    QTL = "caQTL",
    type = "PCA",
    cell_type = cell_type,
    n_factors = local_greedy
  )
)

export_table(df.n_peers, "tables/ST34_n_peers.tsv", "ST34")
