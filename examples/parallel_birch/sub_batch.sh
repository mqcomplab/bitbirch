#!/bin/bash
file_path="fps_splits.txt"
while IFS= read -r line; do
    filename="fps/fps_${line}.npy"
    job_script="job_${line}.sh"
    cat > $job_script << EOL
#!/bin/bash
#SBATCH --job-name=${line}
#SBATCH --mail-user=klopezperez@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=7-00:00:00
#SBATCH --output=${line}.log

module load python
python birch_parallel.py -f $filename
EOL
    # Submit the job
    bash $job_script

    # Remove the job script
    rm $job_script
done < "$file_path"
