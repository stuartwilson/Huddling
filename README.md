# Huddling

Agent based and other models of thermoregulatory huddling

Requires morphologica to be installed, to include tools for visualization, json, hdf5, etc.

To use the template... 

cd template

define a new eperiment by creating a new .json file in the config directory, e.g., expt2.json

Then run using e.g., 

./model hello config/expt2.json logs

the results should be stored in logs directory, which will be created if it does not exist

then you can run python analyse.py to display the results

