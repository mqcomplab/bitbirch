#!/bin/bash
file_path="fps_splits.txt"
while IFS= read -r line; do
    filename="fps/fps_${line}.npy"
    job_script="job_${line}.sh"
    cat > $job_script << EOL
#!/bin/bash
python3 birch_parallel.py -f $filename -t 0.65
EOL
    # Submit the job
    bash $job_script

    # Remove the job script
    rm $job_script
done < "$file_path"