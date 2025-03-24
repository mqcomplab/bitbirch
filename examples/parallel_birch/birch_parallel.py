import argparse
import bitbirch.bitbirch as bb
import numpy
import time
import pickle as pkl
import os

# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='The file containing the fingerprints', required=True)
parser.add_argument('-t', '--threshold', help='The threshold for the BIRCH algorithm', required=True)

# Parse the arguments
args = parser.parse_args()

# Check that directory centroids and mol_ind exists
if not os.path.exists('BFs'):
    os.makedirs('BFs')

# Set the merging criterion
bb.set_merge('diameter')

# Initiate the BIRCH model
birch = bb.BitBirch(threshold=float(args.threshold), branching_factor=50)

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
birch.fit_reinsert(fps_npy, list(range(initial_index, final_index+1)), singly=True)

# End the timer
end = time.time()

# Get the BFs of the leave nodes
BFs, BFs_problematic = birch.prepare_data_BFs(fps_npy, initial_mol=initial_index)

# Delete the fingerprints
del fps_npy

# Save the centroids and indexes
pkl.dump(BFs, open(f'BFs/BFs_{initial_index}_{final_index}.pkl', 'wb'))
pkl.dump(BFs_problematic, open(f'BFs/BFs_problematic_{initial_index}_{final_index}.pkl', 'wb'))

if not os.path.exists('times'):
    os.makedirs('times')
    
with open(f'times/time_{initial_index}_{final_index}.txt', 'w') as f:
    f.write(str(end - start))
    f.write('\n')
