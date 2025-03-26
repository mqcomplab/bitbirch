from bitbirch.cluster_control import chi, dbi, dunn, birch_analysis, radius, jt_isim
import numpy as np
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
from iSIM.sampling import outlier_sampling, medoid_sampling
from iSIM.utils import binary_fps

class ClusterAnalysis:
    def __init__(self, brc, fingerprints, smiles, 
            double=False):
        self.brc = brc
        self.fps = fingerprints
        self.smiles = smiles
        self.double = double

    def clustering_metrics(self, min_size=10):
        cluster_analysis = birch_analysis(self.brc, self.fps, min_size)
        return chi(cluster_analysis), dbi(cluster_analysis), dunn(cluster_analysis)

    def get_biggest_cluster(self):
        clusters = sorted(self.brc.get_cluster_mol_ids(), key=lambda x: len(x), reverse=True)
        biggest_cluster = np.array(clusters[0])
        fps_biggest_cluster = self.fps[biggest_cluster]
        return biggest_cluster, fps_biggest_cluster

    def biggest_cluster_metrics(self, biggest_cluster, fps_biggest_cluster, out_percentage=5, med_percentage=5):
        n_biggest_cluster = len(biggest_cluster)

        #Get iSIM of the biggest cluster
        isim_biggest_cluster = jt_isim(np.sum(fps_biggest_cluster, axis=0), n_biggest_cluster)

        #Get iSIM of the outliers, percentage defaulted to 5%
        outliers = outlier_sampling(fps_biggest_cluster, percentage=out_percentage, n_ary="JT")
        outliers_fps = fps_biggest_cluster[outliers]
        isim_outliers = jt_isim(np.sum(outliers_fps, axis=0), len(outliers))

        #Get iSIM of the medoids, percentage defaulted to 5%
        medoids = medoid_sampling(fps_biggest_cluster, percentage=med_percentage, n_ary="JT")
        medoids_fps = fps_biggest_cluster[medoids]
        isim_medoids = jt_isim(np.sum(medoids_fps, axis=0), len(medoids))

        #Get radius of biggest_cluster
        radius_value = radius(fps_biggest_cluster)

        return n_biggest_cluster, isim_biggest_cluster, isim_outliers, isim_medoids, radius_value

    def scaffold_analysis(self, biggest_cluster):
        smiles_biggest_cluster = self.smiles[biggest_cluster]
        scaffolds = [MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smile) for smile in smiles_biggest_cluster]
        unique_scaffolds = set(scaffolds)
        n_unique = len(unique_scaffolds)
        scaffolds_fps = binary_fps(unique_scaffolds, fp_type = 'RDKIT')
        scaffolds_iSIM = jt_isim(np.sum(scaffolds_fps, axis = 0), n_unique)
        return n_unique, scaffolds_iSIM #Returning number of unique scaffolds and the iSIM

    def get_clusters(self):
        clusters = self.brc.get_cluster_mol_ids()
        return clusters

    def forward(self):
        print('Running Cluster Analysis: \n')
        print('Computing clustering metrics...')
        self.chi_val, self.dbi_val, self.dunn_val = self.clustering_metrics()
        print(f'Chi:{self.chi_val}')
        print(f'DBI:{self.dbi_val}')
        print(f'Dunn:{self.dunn_val}\n')

        print('Analyzing the biggest cluster...')
        biggest_cluster, fps_biggest_cluster = self.get_biggest_cluster()
        (self.n_biggest_cluster, self.isim_biggest_cluster, 
        self.isim_outliers, self.isim_medoids, self.radius_value) = self.biggest_cluster_metrics(biggest_cluster, fps_biggest_cluster)
        print(f'Size of biggest cluster: {self.n_biggest_cluster}')
        print(f'iSIM: {self.isim_biggest_cluster}')
        print(f'iSIM of the outliers: {self.isim_outliers}')
        print(f'iSIM of the medoids: {self.isim_medoids}')
        print(f'Radius of biggest_cluster: {self.radius_value} \n')

        print('Looking for unique scaffolds...')
        self.unique_scaffolds, self.scaffolds_iSIM = self.scaffold_analysis(biggest_cluster)
        print(f'Number of unique scaffolds: {self.unique_scaffolds}')
        print(f'iSIM of unique scaffolds: {self.scaffolds_iSIM} \n')
        
        print('Getting total clusters...')
        clusters = self.get_clusters()
        self.total_clusters = len(clusters)
        self.larger_clusters = [x for x in [len(cluster) > 10 for cluster in clusters] if x]
        print(f'Total number of clusters: {self.total_clusters}')
        print(f'Total number of clusters with more than 10 molecules: {len(self.larger_clusters)} \n')
        
        print('Analysis Complete!')
        
        return clusters
        #return chi_val, dbi_val, dunn_val, n_biggest_cluster, isim_biggest_cluster, isim_outliers, isim_medoids, radius_value


def create_pandas_dataframe(cluster_analysis_instance, method):
    data = {
        "Chi": cluster_analysis_instance.chi_val,
        "DBI": cluster_analysis_instance.dbi_val,
        "Dunn": cluster_analysis_instance.dunn_val,
        "Biggest Cluster Size": cluster_analysis_instance.n_biggest_cluster,
        "Biggest Cluster iSIM": cluster_analysis_instance.isim_biggest_cluster,
        "Outliers iSIM": cluster_analysis_instance.isim_outliers,
        "Medoids iSIM": cluster_analysis_instance.isim_medoids,
        "Cluster Radius": cluster_analysis_instance.radius_value,
        "Number of Unique Scaffolds": cluster_analysis_instance.unique_scaffolds,
        "Scaffolds iSIM": cluster_analysis_instance.scaffolds_iSIM,
        "Total Clusters": cluster_analysis_instance.total_clusters, 
        "Total Clusters > 10 mols": len(cluster_analysis_instance.larger_clusters)
    }

    # Convert to a DataFrame
    df = pd.DataFrame([data], index=[method])

    return df


def check_cluster_structure(brc):
    total_mol = 0
    mol_j_vals = []
    for index, node in enumerate(brc._get_leaves()):
        for c in node.subclusters_:
            mol_j = c.mol_indices
            total_mol+=c.n_samples_
            mol_j_vals.append(mol_j)
            assert c.n_samples_ == len(mol_j), f"Fatal! N_j:{c.n_samples_} and len(mol_j):{len(mol_j)} are incompatible."
    print('Check to ensure all N_j values are identical to all len(mol_j) values passed.')


    print('Total number of molecules:', total_mol)

    mol_j_vals = [i for cl in mol_j_vals for i in cl]

    unique = []
    for m in mol_j_vals:
        if m in unique:
            print(f'Fatal! Index {m} is repeating!')
        else:
            msg = 'There are no repeated indices.'
        unique.append(m)
    if msg:
        print(msg)

    print(f'The range of molecular indices is {min(mol_j_vals)} thru {max(mol_j_vals)}')
