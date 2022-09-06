# CanteraSurfReactorOpenFOAM
> This is a demo programe which call Cantera surface reaction module in OpenFOAM. We use a Al surface reaction mechanism to estimate the surface and gas phase mole fractions v.s. time.

## How to install
The installation requires [OpenFOAM-7](https://github.com/OpenFOAM/OpenFOAM-7) and [LibCantera-2.6](https://anaconda.org/conda-forge/libcantera-devel), which is precompiled in the present repository. 

Source your OpenFOAM via the default path below (or your own path for OpenFOAM bashrc)
```
source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc 
```
1. Clone the [CanteraSurfReactorOpenFOAM repository](https://github.com/ZmengXu/CanteraSurfReactorOpenFOAM)
```
git clone https://github.com/ZmengXu/CanteraSurfReactorOpenFOAM.git
```
2. Install precompiled [LibCantera-2.6](https://anaconda.org/conda-forge/libcantera-devel)
```
conda create --name ct-dev -c conda-forge boost-cpp
conda activate ct-dev
git clone https://github.com/ZmengXu/Pre-compiledCantera.git
cd Pre-compiledCantera
git checkout v2.6.0
tar -zxvf cantera_build.tar.gz
export CANTERA_ROOT=$PWD/cantera_build
export CANTERA_DATA=$CANTERA_ROOT/data
export LD_LIBRARY_PATH=$CANTERA_ROOT/lib:$LD_LIBRARY_PATH
cd ../CanteraSurfReactorOpenFOAM
```
3. Copy the mechanism `Al_surf_gas.yaml` file to `$CANTERA_DATA` folder
```
cp Al_surf_gas.yaml $CANTERA_DATA
```
4. Compile the demo solver `Al_surf_gas` codes
```
wmake
```

## How to use the code
Everytime you use the code, you need to source the OpenFOAM, e.g.,
```
source /opt/openfoam7/etc/bashrc
```
Activate the conda environment
```
conda activate ct-dev
```
and then cd to the `Pre-compiledCantera` folder and 
```
export CANTERA_ROOT=$PWD/cantera_build
export CANTERA_DATA=$CANTERA_ROOT/data
export LD_LIBRARY_PATH=$CANTERA_ROOT/lib:$LD_LIBRARY_PATH
```
and then cd to any folder and run
```
Al_surf_gas
```
Here is the results,
![0D surface reaction sulimation](https://github.com/ZmengXu/CanteraSurfReactorOpenFOAM/blob/master/CanteraChemkinOpenFOAM.png)
