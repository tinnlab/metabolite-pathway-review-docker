#!/bin/sh
set -e

# Enable core dumps
ulimit -c unlimited

# Path to the license file
LIC_FILE=/code/license/license.lic

# Check if the license file exists and run
if test -f "$LIC_FILE"; then
    export MLM_LICENSE_FILE=$LIC_FILE
else
    echo "Please place your license.lic file in the following path: /license/license.lic"
    exit 1
fi

# Set Java path
export MATLAB_JAVA=/opt/matlab/sys/java/jre/glnxa64/jre

# Directory for temporary files
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd "$tmp_dir"

# Your Matlab script to run the Matlab function
cat <<EOF > run_ogpls.m
% Load your dataset from an XLSX file
dataset = readtable('/example_data/urine.xlsx');  % Use the fixed dataset path

% Set parameters
xData = dataset{:, 1:end-1};  % Adjust this based on your dataset structure
yData = dataset{:, end};  % Adjust this based on your dataset structure

% Create a dummy Meta_Path and debiasing_Coef for demonstration purposes
% Replace these with actual values as needed
Meta_Path = randi([0, 1], size(xData, 2), 10);  % Example Meta_Path matrix
lambda = 0.1;  % Example lambda value
debiasing_Coef = ones(1, size(Meta_Path, 2));  % Example debiasing_Coef vector

% Call the function
[U, T, subU, regressionCoef] = OGPLS(xData, yData, Meta_Path, lambda, debiasing_Coef);

% Save the results
save('/code/output/results.mat', 'U', 'T', 'subU', 'regressionCoef');
EOF

# ensure the script can find the ogPLS.m file
cd /code

# Run Matlab
matlab -nodisplay -nodesktop -batch "run('ogPLS.m'); exit"

# Cleanup
cd ..
rm -rf "$tmp_dir"