### COMPILE CPP VERSION TERNIFY

# 1. Create an env for ternify
```
conda create -n ternify python=3.12
conda install mamba -c conda-forge

mamba install rdkit=2025.03.2 -c conda-forge
mamba install librdkit-dev=2025.03.2 -c conda-forge

mamba install libboost-devel=1.86.0 -c conda-forge
mamba install libboost-headers=1.86.0 -c conda-forge
mamba install libboost-python-devel=1.86.0 -c conda-forge # may be not usable

mamba install eigen=3.4.0 -c conda-forge
conda install cmake -c conda-forge # make sure cmake > 3.0. The default cmake version is 4.0.2
```

# 2. Modify the rdkit and rdkit lib path in CMakeLists.txt

# set rdkit lib !!!!
link_directories(/opt/anaconda3/envs/ternify/lib)
set(RDKIT_DIR "/opt/anaconda3/envs/ternify/")

# 3. Compile ternify
conda activate ternify

mkdir build
cd build 
cmake ..
make -j 8

sudo cp ternify /usr/local/bin

3. Run ternify
```
ternify tcs.inp
```
