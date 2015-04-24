#!/usr/bin/ksh

UNAME=`uname | tr A-Z a-z`

# ---------------------------
# Detect number of processors
# ---------------------------
if [ -e /proc/cpuinfo ]
then
    PROCESSORS=`cat /proc/cpuinfo  | grep processor | wc -l`
    if [[ $PROCESSORS > 1 ]]
    then
        export MAKEFLAGS="-j$PROCESSORS"
    fi
fi

    # -------------------------------------------------
    # Compiler Environment variables for normal gcc/g++
    # -------------------------------------------------
    GCC_PATH=/usr/bin
    LOCAL_PATH=/usr/local/bin
    export PATH=$PATH:$GCC_PATH:$LOCAL_PATH
    export CC="/usr/bin/gcc"
    export CXX="/usr/bin/g++"

CUDA_PATH=/usr/local/cuda/bin
export PATH=$PATH:$CUDA_PATH

CUDA_DIR=${CUDA_DIR:="/usr/local/cuda/lib"}

NPW_CUDA_INSTALL_DIR="$HOME/apps/npwcuda"
NPW_CUDA_LIB_DIR="$HOME/apps/npwcuda/lib"
NPW_CUDA_PLUGIN_LIB_DIR="$HOME/apps/npwcuda/plugins/lib"
NPW_CUDA_DIR="$HOME/Dev/npwcuda"
NPW_CUDA_OUTPUT_DIR="$NPW_CUDA_DIR/build"
export NPW_CUDA_DIR

# Library paths
if [[ `uname` == "Linux" ]];
then
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH":"}${CUDA_DIR}:${NPW_CUDA_LIB_DIR}:${NPW_CUDA_PLUGIN_LIB_DIR}:"
    LIBRARY_PATH=${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH LIBRARY_PATH
fi

# -----------------------------------
# Enable vertical synch'ing for NVidia based cards
# -----------------------------------
export __GL_SYNC_TO_VBLANK=1
# we want to allow muti-threaded GL
unset __GL_SINGLE_THREADED
