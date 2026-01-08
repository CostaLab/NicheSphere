import itertools
import numpy as np
import igraph as ig
import networkx as nx
from igraph.clustering import VertexClustering


def _adjacency_from_graph(g: nx.Graph) -> tuple[np.ndarray, np.ndarray]:
    """
    Return positive and negative adjacency matrices (A_pos, A_neg)
    """
    nodes = list(g.nodes())
    n = len(nodes)
    index = {node: i for i, node in enumerate(nodes)}

    A_pos = np.zeros((n, n), dtype=float)
    A_neg = np.zeros((n, n), dtype=float)

    for u, v, data in g.edges(data=True):
        w = data.get("weight", 1.0)

        i, j = index[u], index[v]

        if w > 0:
            A_pos[i, j] = w
            A_pos[j, i] = w
        elif w < 0:
            A_neg[i, j] = -w
            A_neg[j, i] = -w

    return A_pos, A_neg, nodes


def signed_modularity_gomez(g: nx.Graph, partition: dict):
    # positive and negative adjacency matrices
    A_pos, A_neg, nodes = _adjacency_from_graph(g)
    n = len(nodes)
    node_index = {nodes[i]: i for i in range(n)}

    # positive and negative strengths
    k_pos = A_pos.sum(axis=1)
    k_neg = A_neg.sum(axis=1)

    # positive and negative total strengths
    w_pos = k_pos.sum() / 2
    w_neg = k_neg.sum() / 2

    # iterate over all pairs in the same community
    Q_pos = 0.0
    Q_neg = 0.0

    for i_node, j_node in itertools.product(nodes, repeat=2):
        if partition[i_node] != partition[j_node]:
            continue

        i = node_index[i_node]
        j = node_index[j_node]

        # compute deviation against null case random network
        if w_pos > 0:
            expected_pos = (k_pos[i] * k_pos[j]) / (2 * w_pos)
            Q_pos += A_pos[i, j] - expected_pos

        if w_neg > 0:
            expected_neg = (k_neg[i] * k_neg[j]) / (2 * w_neg)
            Q_neg += A_neg[i, j] - expected_neg

    if w_pos > 0:
        Q_pos /= (2 * w_pos)
    if w_neg > 0:
        Q_neg /= (2 * w_neg)

    # calculate final modularity Q
    denom = 2 * w_pos + 2 * w_neg
    Q = (2 * w_pos / denom) * Q_pos - (2 * w_neg / denom) * Q_neg

    return Q


def evaluate_signed_clustering(g: nx.Graph, partition: dict):

    return {"signed_modularity": signed_modularity_gomez(g, partition)}




