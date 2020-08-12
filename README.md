# Huddling

Agent based and other models of thermoregulatory huddling

Requires morphologica to be installed, to include tools for visualization, json, hdf5, etc.

Make in the usual cmake way:

```
cd template
mkdir build
cd build
cmake ..
make
cd ..
```
 
Create a model in a subdirectory, e.g., test, by populating it with test/config.json (there is an example in test already). 

Then run with 

```
./build/agents test 100000 1
```

to run for 100000 steps with seed 1.

the results should be stored in test/agents.h5, which you can visualize by running
 
```
python analyse.py
```

