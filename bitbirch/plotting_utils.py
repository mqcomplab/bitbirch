import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
from bitbirch.bitbirch import jt_isim, BitBirch
from bitbirch.cluster_analysis import ClusterAnalysis
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from rdkit.Chem import Draw
from rdkit import Chem

def tsne_plot(brc, fps, method='', title=None):
    clusters = brc.get_cluster_mol_ids()
    sort_clusters = sorted(clusters, key=lambda x: len(x), reverse=True)[:20]
    # #     # Map each mol ID to its cluster ID
    n_molecules = sum([len(x) for x in sort_clusters])


    cluster_labels = [0] * n_molecules

    fps_clusters = []#[0] * n_molecules
    cluster_labels = []#[0] * n_molecules

    for cluster_id, cluster in enumerate(sort_clusters):
        for idx in cluster:
            fps_clusters.append(fps[idx])
            cluster_labels.append(int(cluster_id))


    scaler = StandardScaler()
    fps_scaled = scaler.fit_transform(fps_clusters)

    tsne = TSNE(n_components=2, perplexity=30, random_state=42)
    fps_tsne = tsne.fit_transform(fps_scaled)

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(fps_tsne[:, 0], fps_tsne[:, 1], c=cluster_labels, cmap='tab20', alpha=0.9)
    cbar = plt.colorbar(scatter, label="Cluster ID")
    cbar.set_ticks(np.arange(20))  # Set ticks at the center of each color
    cbar.set_ticklabels(np.arange(1, 21))  # Set labels from 1 to 20
    plt.xlabel("t-SNE Component 1")
    plt.ylabel("t-SNE Component 2")
    if title:
        plt.title(f"t-SNE of Top 20 Largest Clusters for {title}")
    else:
        plt.title(f"t-SNE of Top 20 Largest Clusters for {method[0].upper() + method[1:]}")
    plt.show()

def mol_relocation_plots(brc_1, brc_2, top_n = 2, title=None):
    
    clusters_1 = brc_1.get_cluster_mol_ids()
    clusters_2 = brc_2.get_cluster_mol_ids()

    list_1 = sorted(clusters_1, key=lambda x: len(x), reverse=True)[:top_n]
    list_2 = sorted(clusters_2, key=lambda x: len(x), reverse=True)[:top_n]

    # Map elements to clusters
    cluster_map_1 = {elem: f"Cluster 1-{i}" for i, cluster in enumerate(list_1) for elem in cluster}
    cluster_map_2 = {elem: f"Cluster 2-{i}" for i, cluster in enumerate(list_2) for elem in cluster}

    # Track flow of elements between clusters
    flows = {}
    for elem in cluster_map_1:
        source = cluster_map_1[elem]
        target = cluster_map_2.get(elem, "Other smaller cluster")
        flow_key = (source, target)
        flows[flow_key] = flows.get(flow_key, 0) + 1

    # Prepare data for Sankey diagram
    labels = list(set([src for src, _ in flows] + [tgt for _, tgt in flows]))
    source_indices = [labels.index(src) for src, _ in flows]
    target_indices = [labels.index(tgt) for _, tgt in flows]
    values = list(flows.values())

    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=source_indices,
            target=target_indices,
            value=values
        )
    )])

    # Update layout
    fig.update_layout(
        title_text=title,
        font_size=10
    )
    fig.show()

def get_cluster_cai(brc, data, smiles):
        cai = ClusterAnalysis(brc, data, smiles)
        clusters = cai.get_clusters()
        return clusters, cai


def get_plot_metrics(clusters, cai, fps, smiles, top=20):
    sort_clusters = sort_clusters = sorted(clusters, key=lambda x: len(x), reverse=True)
    biggest_clusters = sort_clusters[:top]
    fps_biggest_clusters = [[fps[i] for i in c] for c in biggest_clusters]

    n_mol_bcs = [len(x) for x in biggest_clusters] # #molecules


    n_scaff = [cai.scaffold_analysis(c)[0] for c in biggest_clusters]
    isim_scaff = [cai.scaffold_analysis(c)[1] for c in biggest_clusters]


    isim_clusters = [jt_isim(np.sum(fpsc, axis=0), n) for fpsc,n in zip(fps_biggest_clusters, n_mol_bcs)]

    return n_mol_bcs, n_scaff, isim_scaff, isim_clusters


def save_metrics_to_csv(brc, fps, smiles, top=20, method=''):
    clusters, cai = get_cluster_cai(brc, fps, smiles)

    # Extract metrics
    n_mol_bcs, n_scaff, isim_scaff, isim_clusters  = get_plot_metrics(clusters, cai, fps, smiles, top=top)

    # Create DataFrame
    df = pd.DataFrame({
        "Cluster": list(range(1, top + 1)),
        "Number of Molecules": n_mol_bcs,
        "Number of Unique Scaffolds": n_scaff,
        "iSIM of Unique Scaffolds": isim_scaff,
        "iSIM": isim_clusters,
    })

    # Ensure output directory exists
    output_dir = "data/"
    os.makedirs(output_dir, exist_ok=True)

    # Save to CSV
    output_file = os.path.join(output_dir, f"{method}_cluster_metrics.csv")
    df.to_csv(output_file, index=False)
    print(f"Saved cluster metrics for {method} to {output_file}")

def init_plot(method, title=None):
    csv_file = f"data/{method}_cluster_metrics.csv"

    if not os.path.exists(csv_file):
        print(f"Error: CSV file not found at {csv_file}")
    else:
        # Read CSV
        df = pd.read_csv(csv_file)

        # Ensure Cluster is treated as a categorical variable
        df["Cluster"] = df["Cluster"].astype(str)

        # Create figure and axes
        plt.figure(figsize=(10, 5))

        # Plot the number of molecules
        plt.bar(df["Cluster"], df["Number of Molecules"], color="blue", alpha=1, label="Molecules")
        
        # Annotate the number of molecules
        for i, mol in enumerate(df["Number of Molecules"]):
            plt.text(i, mol, f"{mol}", ha="center", va="bottom", fontsize=10, color="black")

        # Plot the number of unique scaffolds
        plt.bar(df["Cluster"], df["Number of Unique Scaffolds"], color="orange", alpha=1, label="Scaffolds")

        # Annotate the number of unique scaffolds
        for i, scaf in enumerate(df["Number of Unique Scaffolds"]):
            plt.text(i, scaf, f"{scaf}", ha="center", va="bottom", fontsize=8, color="white")

        # Labels
        plt.xlabel("Cluster", fontsize=12)
        plt.ylabel("Count", fontsize=12)
        plt.xticks(range(len(df["Cluster"])))
        plt.legend(loc="upper right")
     
     
        # Plot iSIM
        plt.twinx()
        plt.plot(df["Cluster"], df["iSIM"], color="green", marker="o", label="iSIM", linewidth=2)
        plt.ylabel("iSIM", fontsize=12)
        plt.yticks(np.arange(0, 1.1, 0.1))
        plt.ylim(0, 1)

        plt.tight_layout()

        plt.legend(loc="upper right", bbox_to_anchor=(0.80, 1))

        # Title
        if title:
            plt.title(f'Top 20 Cluster Metrics for {title}', fontsize=14)
        else:
            plt.title(f'Top 20 Cluster Metrics for {method[0].upper() + method[1:]}', fontsize=14)


        plt.show() 

def mol_visualization(smiles, brc, cluster = 0):
    clustered_ids = brc.get_cluster_mol_ids()

    mols = [Chem.MolFromSmiles(smiles[i]) for i in clustered_ids[cluster-1]]
    print('Number of molecules:', len(mols))

    for i in range(0, len(mols), 30):

        img_1 = Draw.MolsToGridImage(mols[i:i+30], molsPerRow=5)

        # Save the images from the .data attribute
        with open(f'cluster_{cluster}_{i}.png', 'wb') as f:
            f.write(img_1.data)