"""
Analisis jaringan PPI (STRING) untuk 4 protein target hasil intersection
SwissTargetPrediction x GeneCards: DPP4, MMP9, MMP2, MIF, ditambah satu
lapis ketetanggaan.

Tujuan: memeringkat protein mana yang paling layak dijadikan reseptor pada
tahap molecular docking.
"""

from pathlib import Path

import networkx as nx
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
EDGES = DATA / "17_string_interactions_protein target.tsv"

SEEDS = {"DPP4", "MMP9", "MMP2", "MIF"}


def main():
    edges = pd.read_csv(EDGES, sep="\t")
    edges = edges.rename(columns={"#node1": "node1"})

    g = nx.Graph()
    for row in edges.itertuples():
        g.add_edge(row.node1, row.node2, weight=float(row.combined_score))

    print(f"simpul: {g.number_of_nodes()}   sisi: {g.number_of_edges()}")
    print(f"seed  : {sorted(SEEDS)}")
    print(f"tetangga tambahan: {sorted(set(g.nodes) - SEEDS)}\n")

    metrics = pd.DataFrame(
        {
            "degree": dict(g.degree()),
            "weighted_degree": dict(g.degree(weight="weight")),
            "betweenness": nx.betweenness_centrality(g, weight=None),
            "closeness": nx.closeness_centrality(g),
            "clustering": nx.clustering(g, weight="weight"),
        }
    )
    metrics["is_seed"] = metrics.index.isin(SEEDS)
    metrics = metrics.sort_values(["weighted_degree", "degree"], ascending=False)
    metrics.index.name = "protein"
    metrics = metrics.round(4)
    metrics.to_csv(HERE / "18.5_ppi_hub_ranking.csv")

    print(metrics.to_string())
    print("\nDensitas jaringan:", round(nx.density(g), 4))

    hubs = metrics[metrics.is_seed].index.tolist()
    print("\nPeringkat seed berdasarkan weighted degree:", " > ".join(hubs))
    print(f"-> {HERE / '18.5_ppi_hub_ranking.csv'}")


if __name__ == "__main__":
    main()
