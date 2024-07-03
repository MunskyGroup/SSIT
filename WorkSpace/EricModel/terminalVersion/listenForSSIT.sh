#!/bin/bash

# This script is used to listen for files with the extension _ready.txt in the
# directory 'ready_files'. When such a file is found, it is read. The first line of the
# file is used to determine the string for a variable called 'data_set'. 
# The second line is used to determine the string for a variable called 'ssit_method'.
# The third line is used to determine the string for a variable called 'save_name'.

# The script then renames the file to remove  '_ready.txt' from the filename and renames
# the file to have the extension '_run.txt'. The script then moves the file to the directory
# 'run_files'. 

# The script then runs the matlab function SSIT with the arguments 'data_set' and 'ssit_method'.
# The script does not wait for the matlab function to finish before moving on to the next file.

# The script will continue to listen for files until the user kills the script.

# The script is intended to be run in the background.

pathtoSSIT="../../../src"


while true
do
    # Find all files in the directory 'ready_files'
    files=$(ls ready_files/*)

    # Truncate the list to only include files with the extension _ready.txt
    files=$(echo $files | tr ' ' '\n' | grep '_ready.txt')

    for file in $files
    # The loop will run for each file with the extension _ready.txt
    do
        # if there are no files with the extension _ready.txt, the script will sleep for 10 seconds
        if [ -z "$file" ] 
        then
            # If there are no files with the extension _ready.txt, the script will sleep for 10 seconds
            sleep 10
            continue
        fi

        # print the file name
        echo $file

        # read the first line of the file and assign it to the variable 'model_template'. 
        # This string should only go to the first return carriage.
        model_template=$(head -n 1 $file)
        model_name=$(head -n 2 $file | tail -n 1)
        ssit_method=$(head -n 3 $file | tail -n 1)
        data_set=$(head -n 4 $file | tail -n 1)
        save_name=$(head -n 5 $file | tail -n 1)

        # Move the file to the directory 'run_files'
        mv $file run_files
        
        # Print the matlab function call command
        echo "/Applications/MATLAB_R2023a.app/bin/matlab -nojvm -nodisplay -nosplash -nodesktop -r \"addpath(genpath('$pathtoSSIT')); addpath(genpath('../tmpPropensityFunctions')); SSIT('$model_template', '$model_name', '$ssit_method', '$data_set', '$save_name'); exit &\""

        # Run the matlab function SSIT with the arguments model_template, model_name, data_set, ssit_method, and save_name
        /Applications/MATLAB_R2023a.app/bin/matlab -nojvm -nodisplay -nosplash -nodesktop -r "addpath(genpath('$pathtoSSIT')); addpath(genpath('../tmpPropensityFunctions')); SSIT('$model_template', '$model_name', '$ssit_method', '$data_set', '$save_name'); exit" > log.txt 2>&1 &
    done
done