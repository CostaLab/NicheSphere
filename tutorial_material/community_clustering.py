import pandas as pd
import numpy as np
import os
from typing import Dict, Any
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.lines

from collections import defaultdict
from community_layout.layout_class import CommunityLayout
from cdlib import algorithms, evaluation, NodeClustering
import igraph as ig
import networkx as nx
import leidenalg as la
import sklearn

import sys
sys.path.insert(1, os.path.join(os.getcwd()+"/signed_louvain")) # change path to signed louvain folder location
import util
import random
import community_detection as signed_louvain
from signet.cluster import Cluster

import nichesphere


# Available community detection methods
UNSIGNED_METHODS = ['louvain_unsigned', 'cdlib_louvain', 'cdlib_leiden']
SIGNED_METHODS = ['leidenalg_mod', 'louvain_signed', 'leidenalg_CPM', 
                      'spinglass_weights', 'sponge_no_weights']

def _get_unsigned_methods():
    return UNSIGNED_METHODS

def _get_signed_methods():
    return SIGNED_METHODS

def _choose_method(n: int) -> str:
    '''
    returns the method name corresponding to the index position in the list of all available methods
    '''
    method_list = UNSIGNED_METHODS + SIGNED_METHODS

    if n > len(method_list)-1:
        raise ValueError(f"Invalid method index '{n}'. Choose an index in the range of 0 to {len(method_list)-1}.")

    return method_list[n]
    

# Default parameters for each community detection method 
DEFAULT_PARAMS = {
    # Unsigned methods
    "louvain_unsigned": {"seed": 12, "resolution": 1.1},
    "cdlib_louvain": {"seed": 12, "resolution": 1.1},
    "cdlib_leiden": {"seed": 12, "resolution": 1.1},

    # Signed methods
    "leidenalg_mod": {"seed": 12},
    "louvain_signed": {"seed": 12, "k": 2, "pass_max": 40},
    "leidenalg_CPM": {"seed": 12, "resolution_CPM": 0.017},
    "spinglass_weights": {"seed": 12, "spins": 25, "gamma": 1.0, "lambda": 0.8},
    "sponge_no_weights": {"seed": 12, "n_clusters": 3},
}

def _get_default_params(name: str) -> Dict[str, Any]:
    if name not in UNSIGNED_METHODS + SIGNED_METHODS:
        raise ValueError(f"Invalid method name '{name}'.")

    return DEFAULT_PARAMS[name]


def niche_coloc_plot(
    HM: pd.DataFrame,
    HMsimm: pd.DataFrame,
    method_name: str,
    community_dict: dict,
    cl: CommunityLayout,
    clist: list,
    save_plot: bool = False
) -> None:
    """
    Generates and optionally saves a colocalization network plot using Nichesphere.

    Args:
        HM (pd.DataFrame): Differential co-localization matrix.
        HMsimm (pd.DataFrame): Adjacency matrix.
        method_name (str): Name of the clustering method used (for title/filename).
        community_dict (Dict[str, List[str]]): Dictionary mapping community labels to node lists.
        cl (CommunityLayout): Computed layout object containing node positions.
        clist (List[str]): List of hex color strings for the communities.
        save_plot (bool, optional): Whether to save the plot to disk. Defaults to False.
    """
    
    # Use context manager to prevent global style pollution
    with plt.rc_context({'axes.facecolor': 'None'}):
        gCol=nichesphere.coloc.colocNW(x_diff=HM, 
                                       adj=HMsimm,
                                       cell_group=community_dict, 
                                       clist=clist, 
                                       nodeSize='betweeness', 
                                       layout=None,                         #layout needs to be set to None if we provide node positions
                                       lab_spacing=0.05, 
                                       thr=1, 
                                       alpha=0.4, 
                                       fsize=(10,10), 
                                       pos=cl.full_positions,                             #node positions (from the CommunityLayout function)
                                       edge_scale=1,                        #edge width
                                       legend_ax=[0.7, 0.05, 0.15, 0.2])    #legend position
        # Custom Legend construction
        community_labels = list(community_dict.keys())
        legend_elements = [
            plt.Line2D([0], [0], marker="o", color='w', markerfacecolor=clist[i],
                       lw=4, label=community_labels[i], ms=10)
            for i in range(len(community_labels))
        ]
        
        plt.gca().add_artist(
            plt.legend(handles=legend_elements, loc='lower left', fontsize=13,
                       title='Niches', alignment='left')
        )
    
        plt.title(f"Community Detection ({method_name})", fontsize=16, weight='bold')
    
        if save_plot:
            save_dir = Path("figures")
            save_dir.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_dir / f"community_detection_{method_name}.jpg")


def niche_plot_community_clustering(
    HM: pd.DataFrame,
    HMsimm: pd.DataFrame,
    method_name: str,
    gCol: nx.Graph,
    vc_map: dict,
    cell_types: list,
    plot: bool = False,
    save_plot: bool = False,
    community_names: list = None
) -> pd.DataFrame:
    """
    Coordinates the layout calculation, coloring, and plotting of community clusters.

    Args:
        HM (pd.DataFrame): Differential co-localization matrix.
        HMsimm (pd.DataFrame): Adjacency matrix.
        method_name (str): Name of the clustering method.
        gCol (nx.Graph): The network graph.
        vc_map (Dict[str, int]): Mapping of vertex names to cluster IDs.
        cell_types (List[str]): List of cell types involved.
        plot (bool, optional): If True, generates the plot. Defaults to False.
        save_plot (bool, optional): If True, saves the plot to disk. Defaults to False.
        community_names (List[str], optional): Custom names for communities. Defaults to None.

    Returns:
        pd.DataFrame: A dataframe containing cell types and their assigned niche colors.
    """

    # Invert vc_map to group nodes by community ID
    communities = defaultdict(list)
    for vertex, community_id in vc_map.items():
        communities[community_id].append(vertex)

    # Ensure communities are processed in sorted order of their IDs
    sorted_cluster_ids = sorted(communities.keys())
    community_list = [communities[i] for i in sorted_cluster_ids]

    # Calculate community layout
    cl = CommunityLayout(gCol,
        community_compression = 0.4,
        layout_algorithm = nx.spring_layout,
        layout_kwargs = {"k":75, "iterations":1000},
        community_algorithm = community_list)

    # Retrieve final communities from layout 
    d_cluster = {index: list(value) for index, value in enumerate(cl.communities())}
    num_communities = len(d_cluster)

    # Replace community labels if community_names is passed with the correct length
    if community_names is None or len(community_names) == 0:
        labels = [f"Community {i + 1}" for i in range(num_communities)]
    elif len(community_names) == num_communities:
        labels = community_names
    else:
        raise ValueError(
            f"Length of community_names ({len(community_names)}) does not match "
            f"number of detected communities ({num_communities})."
        )

    # Map labels to members
    community_dict = {labels[i]: d_cluster[i] for i in range(num_communities)}
    
    # Sssign colors to communities
    cmap = plt.get_cmap("rainbow")   
    clist = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, num_communities)]
    
    community_cols = pd.Series(clist, index=list(community_dict.keys()))
    
    # Generate result DataFrame
    community_df = nichesphere.tl.cells_niche_colors(
        CTs=cell_types,
        niche_colors=community_cols,
        niche_dict=community_dict
    )
    
    if plot:
        niche_coloc_plot(HM, HMsimm, method_name, community_dict, cl, clist, save_plot)
        
        # Debug output
        print(f"Communities detected ({method_name}):")
        print(pd.DataFrame.from_dict(community_dict, orient='index').T.to_string(index=False))

    return community_df
    

def _check_params(method_name: str, params: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """
    Returns a complete parameter dictionary for a given method.
    Only includes relevant parameters for that method.
    Input argument `params` overrides defaults.
    """

    if method_name not in DEFAULT_PARAMS:
        raise ValueError(f"Unknown method '{method_name}'.")

    defaults = DEFAULT_PARAMS[method_name]
    
    # Keep only input keys that are in the defaults for the method
    filtered_params = {k: v for k, v in (params or {}).items() if k in defaults}

    return {**defaults, **filtered_params}
    

def community_detection(
    *,
    gCol: nx.Graph,
    HM : pd.DataFrame,
    HMsimm: pd.DataFrame,
    method_name: str,
    cell_types: list,
    community_names: list[str] | None = None,
    plot: bool = False,
    save_plot: bool = False,
    params: Dict[str, Any] | None = None,
) -> tuple:
    """
    Performs community detection using various signed or unsigned algorithms and visualizes the results.

    Args:
        gCol (nx.Graph): NetworkX graph object.
            Note: For unsigned methods, the graph is often rebuilt from `hm_simm`.
        HM (pd.DataFrame): Differential co-localization matrix (used for weights/plotting).
        HMsimm (pd.DataFrame): Similarity/Adjacency matrix.
        method_name (str): The algorithm to use. Options include:
            Unsigned: ['louvain_unsigned', 'cdlib_louvain', 'cdlib_leiden']
            Signed: ['leidenalg_mod', 'louvain_signed', 'leidenalg_CPM', 'spinglass_weights', 'sponge_no_weights']
        plot (bool): Whether to generate plots.
        save_plot (bool): Whether to save plots to file.
        cell_types (List[str]): List of cell type labels.
        **params (dict): Fine tuning parameters for the community detection.
        community_names (List[str], optional): Custom labels for the resulting communities.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary where keys are method names and values are DataFrames of cell niche colors.
    """
    
    if method_name not in UNSIGNED_METHODS and method_name not in SIGNED_METHODS:
        raise ValueError(f"Method '{method_name}' is not supported.")

    # Check fine tuning parameters and use default
    all_params = _check_params(method_name, params)
   
    
    vc_map: dict[str, int] = {}

    ## Unsigned methods
    if method_name in UNSIGNED_METHODS:
        # Reconstruct unsigned graph from adjacency matrix
        gCol_unsigned = nx.from_pandas_adjacency(HMsimm, create_using=nx.Graph)
        g = gCol_unsigned      

        if method_name == 'louvain_unsigned':
            comms = nx.community.louvain_communities(gCol_unsigned, resolution = all_params['resolution'], seed = all_params['seed'])
            vc_map = {node: i for i, members in enumerate(comms) for node in members}
            
        elif method_name == 'cdlib_louvain':
            clustering = algorithms.louvain(gCol_unsigned, weight='weight', resolution =  all_params['resolution'], randomize = False)
            vc_map = {node: i for i, members in enumerate(clustering.communities) for node in members}
            
        elif method_name == 'cdlib_leiden':
            clustering = algorithms.leiden(gCol_unsigned, weights='weight', resolution_parameter =  all_params['resolution'], seed = all_params['seed'])
            vc_map = {node: i for i, members in enumerate(clustering.communities) for node in members}

    ## Signed methods 
    else:
        g = gCol
        gCol_ig = ig.Graph.from_networkx(gCol) # iGraoh
    
        if method_name == 'leidenalg_mod':
            partition = la.find_partition(gCol_ig, la.ModularityVertexPartition, seed=all_params['seed'])
            vc_map = {v['_nx_name']: partition.membership[v.index] for v in gCol_ig.vs}
            
        elif method_name == 'louvain_signed':
            # Sub-graph with no negative edges 
            G_pos=gCol.copy()
            to_remove=[(a,b) for a, b, attrs in G_pos.edges(data=True) if attrs["weight"] <= 0]
            G_pos.remove_edges_from(to_remove)
            
            # Sub-graph with no positive edges 
            G_neg=gCol.copy()
            to_remove=[(a,b) for a, b, attrs in G_neg.edges(data=True) if attrs["weight"] >= 0]
            G_neg.remove_edges_from(to_remove)

            comms = signed_louvain.best_partition(layers=[G_pos, G_neg],
                            layer_weights=[1., -1.],
                            resolutions=[1., 1.],
                            masks=[False, True],
                            k=all_params['k'],
                            initial_membership=None,
                            weight='weight',
                            random_state=all_params['seed'],
                            pass_max=all_params['pass_max'],
                            return_dendogram=False,
                            silent=True)
            
            vc_map = comms
            
        elif method_name == 'leidenalg_CPM':
            # Vectorized weight scaling
            weights = np.array(gCol_ig.es['weight'])
            min_w, max_w = weights.min(), weights.max()

            # Shift positive to handle CPM requirements (weights in [0, 1])
            if max_w > min_w:
                weights_scaled = (weights - min_w) / (max_w - min_w)
            else:
                weights_scaled = np.zeros_like(weights)
                
            gCol_ig.es['weight_pos'] = weights_scaled.tolist()

            partition = la.find_partition(
                g_ig, la.CPMVertexPartition, 
                weights = 'weight_pos', 
                resolution_parameter = all_params['resolution_CPM'], 
                seed = all_params['seed']
            )

            vc_map = {v['_nx_name']: partition.membership[v.index] for v in gCol_ig.vs}
            
        elif method_name == 'spinglass_weights':
            random.seed(all_params['seed'])
            partition = gCol_ig.community_spinglass(weights="weight", spins=all_params['spins'], gamma=all_params['gamma'], lambda_=all_params['lambda'], implementation='neg')
            vc_map = {v['_nx_name']: partition.membership[v.index] for v in gCol_ig.vs}
            
        elif method_name == 'sponge_no_weights':
            # Extract sparse adjacency (only signs, not weights)
            A = gCol_ig.get_adjacency_sparse(attribute='weight').sign()
            
            cluster_algo = Cluster((A.multiply(A>0), -A.multiply(A<0)))
            # Update method name to reflect K (affecting output key)
            method_name = (f'sponge_no_weights_k{all_params["n_clusters"]}')
            np.random.seed(all_params['seed'])
            sponge_labels = cluster_algo.SPONGE(all_params['n_clusters'])

            vc_map = {v['_nx_name']: int(sponge_labels[v.index]) for v in gCol_ig.vs}


    communities_df =  niche_plot_community_clustering(
        HM, HMsimm, method_name, g, vc_map,
        cell_types, plot, save_plot, community_names
    )

    return communities_df, vc_map, all_params

    