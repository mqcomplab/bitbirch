import argparse
from bitbirch import BitBirch
import numpy
import time
import pickle as pkl
import os

# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='The file containing the fingerprints', required=True)

# Parse the arguments
args = parser.parse_args()

# Check that directory centroids and mol_ind exists
if not os.path.exists('centroids'):
    os.makedirs('centroids')
if not os.path.exists('mol_ind'):
    os.makedirs('mol_ind')

# Initiate the BIRCH model
birch = BitBirch(threshold=0.65, branching_factor=50)

# Read the fingerprints
fps_npy = numpy.load(args.file, mmap_mode='r')

# Get the initial index and the final index
name_without_extension = args.file.split('.')[0]
name_without_path = name_without_extension.split('/')[-1]

initial_index = int(name_without_path.split('_')[-2])
final_index = int(name_without_path.split('_')[-1])

# Start the timer
start = time.time()

# Fit the BIRCH model
birch.fit_reinsert(fps_npy, list(range(initial_index, final_index+1)))

# End the timer
end = time.time()

# Get the centroids
centroids = birch.get_centroids()

# Get the indexes
indexes = birch.get_cluster_mol_ids()

# Save the centroids and indexes
pkl.dump(centroids, open(f'centroids/centroids_{initial_index}_{final_index}.pkl', 'wb'))
pkl.dump(indexes, open(f'mol_ind/mol_ind_{initial_index}_{final_index}.pkl', 'wb'))

if not os.path.exists('times'):
    os.makedirs('times')
    
with open(f'times/time_{initial_index}_{final_index}.txt', 'w') as f:
    f.write(str(end - start))
    f.write('\n')

