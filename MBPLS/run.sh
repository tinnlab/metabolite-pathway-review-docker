#!/bin/sh
set -e

# Path to the license file
LIC_FILE=/code/license/license.lic

# Check if the license file exists and run
if test -f "$LIC_FILE"; then
    export MLM_LICENSE_FILE=$LIC_FILE
else
    echo "Please place your license.lic file in the following path: /license/license.lic"
    exit 1
fi

# Directory for temporary files
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd "$tmp_dir"

# Your Matlab script to run the Matlab function
cat <<EOF > run_mbpls.m
% Load your dataset from an XLSX file
dataset = readtable('/example_data/urine.xlsx');  % Replace our example data with your actual dataset file name

% Assuming xData and yLabel are parts of your dataset
xData = dataset(:, 1:end-1);  % Adjust this based on your dataset structure
yLabel = dataset(:, end);  % Adjust this based on your dataset structure

% Convert tables to arrays
xData = table2cell(xData);
yLabel = table2array(yLabel);

compNum = 3;  % Set the number of components, replace with your desired number

% Call the function
[PIP,tb,t,wb,w] = MBPLS_PIP(xData, yLabel, compNum);

% Save the results
save('/output/results.mat', 'PIP', 'tb', 't', 'wb', 'w');
EOF

# ensure the script can find the MBPLS_PIP.m file
cd /code

# Run Matlab
matlab -nodisplay -nodesktop -batch "run('MBPLS_PIP.m'); exit"

# Cleanup
cd ..
rm -rf "$tmp_dir"