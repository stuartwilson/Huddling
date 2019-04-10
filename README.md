# Huddling

Agent based and other models of thermoregulatory huddling

Requires morphologica to be installed, to include tools for visualization, json, hdf5, etc.

to run an eperiment create a new .json file in the config directory, e.g., expt2.json

Then run using e.g., 

./test hello config/expt2.json logs

the results should be stored in logs directory, which will be created if it does not exist

then you can run python analyse.py to display the results

