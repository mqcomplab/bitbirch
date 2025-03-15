import argparse
from bitbirch import BitBirch
import pickle as pkl
import time
import numpy as np

# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--threshold', help='The threshold for the BitBirch model', required=True)

# Parse the arguments
args = parser.parse_args()

centroids_concat = []
centroids_dictionary = {}
centroids_indexes_fps = {}
centroids_indexes_centroids = {}
tracker = 0

# Read the fps_split indexes from the .txt file
with open('fps_splits.txt', 'r') as f:
    splits = f.read().splitlines()

start = time.time()
brc = BitBirch(threshold=float(args.threshold), branching_factor=50)  #threshold could be changed by the user
for split in splits:
    f = open(f'centroids/centroids_{split}.pkl', 'rb')
    data = pkl.load(f)
    brc.fit(np.array(data))
    for i, centroid in enumerate(data):
        centroids_indexes_fps[tracker + i] = split
        centroids_indexes_centroids[tracker + i] = i       
    tracker += len(data)

f.close()
del data

clustered_centroids_ids = brc.get_cluster_mol_ids()

with open('final_centroids_parallel.pkl', 'wb') as f:
    pkl.dump(brc.get_centroids(), f)

del brc

# Obtain the final clustering
final_clusters = []
for centroid_cluster in clustered_centroids_ids:
    cluster = []
    for centroid_id in centroid_cluster:
        # Open the file corresponding to the centroid
        file_mol_ids = open(f'mol_ind/mol_ind_{centroids_indexes_fps[centroid_id]}.pkl', 'rb')
        mol_ids = pkl.load(file_mol_ids)
        file_mol_ids.close()
        cluster.extend(mol_ids[centroids_indexes_centroids[centroid_id]])
    final_clusters.append(cluster)

del mol_ids

# Save the clustered_ids
with open('clustered_ids_parallel.pkl', 'wb') as f:
    pkl.dump(final_clusters, f)

# Save the time in a .txt file
with open('time_centroid_clustering.txt', 'w') as f:
    f.write(str(time.time() - start))
    f.write('\n')
