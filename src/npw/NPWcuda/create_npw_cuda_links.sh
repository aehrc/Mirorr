#!/bin/bash
#this script looks for non-parametric windowing openGL source files
#and, if found, creates links to them in this folder


#colors!

BLACK='\E[30m'
GREEN='\E[32m'
RED='\E[31m'
YELLOW='\E[33m'
BLUE='\E[34m'
MAGENTA='\E[35m'
CYAN='\E[36m'
WHITE='\E[37m'

cecho ()        # arg 1 : message,    arg 2 : color
{
    local default_msg="No message passed."
    message=${1:-$default_msg}
    color=${2:-$black}
    echo -en "$color"
    echo "$message"
    tput sgr0
    return
}
cechon ()        # arg 1 : message,    arg 2 : color
{
    local default_msg="No message passed."
    message=${1:-$default_msg}
    color=${2:-$black}
    echo -en "$color"
    echo -n "$message"
    tput sgr0
    return
}


#
# set the target directory. If not specified, use the current directory
#
TARGET_DIRECTORY=$(pwd)

if [ $# -ge 1 ] ; then
    TARGET_DIRECTORY=$1
fi


#
# locate the source directory
#
NPW_OPENGL_SOURCE_DIRECTORY="*milxRegistration/NPWopengl"
NPW_CUDA_SOURCE_DIRECTORY="*milxRegistration/NPWcuda/npw"
N_OPENGL_SOURCE_DIRECTORIES=`locate $NPW_OPENGL_SOURCE_DIRECTORY -c`
N_CUDA_SOURCE_DIRECTORIES=`locate $NPW_CUDA_SOURCE_DIRECTORY -c`

OPENGL_SOURCE_DIRECTORY=`locate $NPW_OPENGL_SOURCE_DIRECTORY`
CUDA_SOURCE_DIRECTORY=`locate $NPW_CUDA_SOURCE_DIRECTORY`

if [[ $N_OPENGL_SOURCE_DIRECTORIES == 0 ]] ; then
    cecho "*** ERROR *** : could not find source directory \""$NPW_OPENGL_SOURCE_DIRECTORY"\". Have you checked out milx-view somewhere?" $RED
    exit 0
fi
if [[ $N_CUDA_SOURCE_DIRECTORIES == 0 ]] ; then
    cecho "*** ERROR *** : could not find source directory \""$CUDA_OPENGL_SOURCE_DIRECTORY"\". Have you checked out milx-view somewhere?" $RED
    exit 0
fi

#if there are more, use the one with the most recently modified file in it
if [[ $N_OPENGL_SOURCE_DIRECTORIES -ge 2 ]] ; then
    cecho "*** WARNING *** : multiple possible source directories were found:" $MAGENTA

    FIRST_DIR=`locate $NPW_OPENGL_SOURCE_DIRECTORY | head -n1`

    MOST_RECENT_DIR=$FIRST_DIR
    MOST_RECENT_FILE=`ls $FIRST_DIR -t | head -n1`

    for CUR_DIR in $OPENGL_SOURCE_DIRECTORY ; do

        echo "  "$CUR_DIR

        MOST_RECENT_FILE_HERE=`ls "$CUR_DIR" -t | head -n1`
        if test $CUR_DIR"/"$MOST_RECENT_FILE_HERE -nt $MOST_RECENT_DIR"/"$MOST_RECENT_FILE ; then
            MOST_RECENT_FILE=$MOST_RECENT_FILE_HERE
            MOST_RECENT_DIR=$CUR_DIR
        fi
        
    done

   cechon "using the most recent directory: " $MAGENTA
   cecho $MOST_RECENT_DIR $YELLOW

   OPENGL_SOURCE_DIRECTORY=$MOST_RECENT_DIR
fi

cechon "found npw opengl source directory: " $GREEN
cecho $OPENGL_SOURCE_DIRECTORY $YELLOW

SKIP_OPENGL=0 

if [[ $OPENGL_SOURCE_DIRECTORY == $TARGET_DIRECTORY ]] ; then
    cecho "the target directory \""$TARGET_DIRECTORY"\" is the same as the OpenGL source directory \""$OPENGL_SOURCE_DIRECTORY"\"! Skipping OpenGL files" $RED
    SKIP_OPENGL=1
fi

#if there are more, use the one with the most recently modified file in it
if [[ $N_CUDA_SOURCE_DIRECTORIES -ge 2 ]] ; then
    cecho "*** WARNING *** : multiple possible source directories were found:" $MAGENTA

    FIRST_DIR=`locate $NPW_CUDA_SOURCE_DIRECTORY | head -n1`

    MOST_RECENT_DIR=$FIRST_DIR
    MOST_RECENT_FILE=`ls $FIRST_DIR -t | head -n1`

    for CUR_DIR in $CUDA_SOURCE_DIRECTORY ; do

        echo "  "$CUR_DIR

        MOST_RECENT_FILE_HERE=`ls "$CUR_DIR" -t | head -n1`
        if test $CUR_DIR"/"$MOST_RECENT_FILE_HERE -nt $MOST_RECENT_DIR"/"$MOST_RECENT_FILE ; then
            MOST_RECENT_FILE=$MOST_RECENT_FILE_HERE
            MOST_RECENT_DIR=$CUR_DIR
        fi
        
    done

   cechon "using the most recent directory: " $MAGENTA
   cecho $MOST_RECENT_DIR $YELLOW

   CUDA_SOURCE_DIRECTORY=$MOST_RECENT_DIR
fi

cechon "found cuda opengl source directory: " $GREEN
cecho $CUDA_SOURCE_DIRECTORY $YELLOW

SKIP_CUDA=0
if [[ $CUDA_SOURCE_DIRECTORY == $TARGET_DIRECTORY ]] ; then
    cecho "the target directory \""$TARGET_DIRECTORY"\" is the same as the CUDA source directory\""$CUDA_SOURCE_DIRECTORY"\"! Skipping CUDA files" $RED
    SKIP_CUDA=1
fi



#
# check if all files are present
#
REQUIRED_CUDA_FILES=( NPWCudaExtractGeometryKernel.cu NPWCudaKernelPrototypes.h NPWCudaUtilities.h NPWCudaGeometryFunctions.h NPWCudaOpenGLRenderer.cxx NPWCudaShapeDetermination.h NPWCudaVertexBufferPointer.h NPWCudaConstants.h NPWCudaHistogramCreator.cxx NPWCudaOpenGLRenderer.h NPWCudaShapeExpansion.h NPWCudaWindow.cxx NPWCudaDataPointer.h NPWCudaHistogramCreator.h NPWCudaRenderGeometryDirectlyFromCUDAKernel.cu NPWCudaUtilities.cxx NPWCudaWindow.h )
REQUIRED_OPENGL_FILES=( NPWStopwatch.h MutualInformation.h )


#CUDA
if [[ $SKIP_CUDA -eq 0 ]] ; then

    CUDA_SOURCE_FILES=`ls $CUDA_SOURCE_DIRECTORY`
    MISSING_CUDA_FILES=()

    for ((i=0; i<${#REQUIRED_CUDA_FILES[*]}; i++)); do

        FILEWASFOUND=0

        for FOUND_FILE in $CUDA_SOURCE_FILES ; do
            if [ "$FOUND_FILE" == "${CUDA_REQUIRED_FILES[$i]}" ] ; then
                FILEWASFOUND=1
            fi
        done

        if [ $FILEWASFOUND == 0 ] ; then
            MISSING_FILES[${#MISSING_CUDA_FILES[*]}]=${REQUIRED_CUDA_FILES[$i]}
        fi

    done

    if [[ ${#MISSING_CUDA_FILES[@]} != 0 ]] ; then
        cecho "*** ERROR *** : not all source files are present!" $RED
        echo ${#MISSING_CUDA_FILES[@]}" missing files:"
        for ((i=0; i<${#MISSING_CUDA_FILES[*]}; i++)); do
            cecho ${MISSING_CUDA_FILES[$i]} $RED
        done
    else
        cecho "all "${#REQUIRED_CUDA_FILES[@]}" required CUDA source files were found" $GREEN
    fi
fi

#OPENGL
if [[ $SKIP_OPENGL -eq 0 ]] ; then

    OPENGL_SOURCE_FILES=`ls $OPENGL_SOURCE_DIRECTORY`
    MISSING_OPENGL_FILES=()

    for ((i=0; i<${#REQUIRED_OPENGL_FILES[*]}; i++)); do

        FILEWASFOUND=0

        for FOUND_FILE in $OPENGL_SOURCE_FILES ; do
            if [ "$FOUND_FILE" == "${REQUIRED_OPENGL_FILES[$i]}" ] ; then
                FILEWASFOUND=1
            fi
        done

        if [ $FILEWASFOUND == 0 ] ; then
            MISSING_OPENGL_FILES[${#MISSING_OPENGL_FILES[*]}]=${REQUIRED_OPENGL_FILES[$i]}
        fi

    done

    if [[ ${#MISSING_OPENGL_FILES[@]} != 0 ]] ; then
        cecho "*** ERROR *** : not all source files are present!" $RED
        echo ${#MISSING_OPENGL_FILES[@]}" missing files:"
        for ((i=0; i<${#MISSING_OPENGL_FILES[*]}; i++)); do
            cecho ${MISSING_OPENGL_FILES[$i]} $RED
        done
    else
        cecho "all "${#REQUIRED_OPENGL_FILES[@]}" required OpenGL source files were found" $GREEN
    fi
fi

#
# create the links, it they don't already exist
#
if [[ $SKIP_CUDA -eq 0 ]] ; then
    echo -n "creating symbolic links in "
    cecho $TARGET_DIRECTORY $YELLOW

    for ((i=0; i<${#REQUIRED_CUDA_FILES[*]}; i++)); do

        TARGET_FILE=$TARGET_DIRECTORY"/"${REQUIRED_CUDA_FILES[$i]}

        if [[ -e $TARGET_FILE ]] ; then
        
            echo -n "skipping "
            cechon ${REQUIRED_CUDA_FILES[$i]} $CYAN
            echo -n " because it already exists in "
            cecho $TARGET_DIRECTORY $YELLOW

        else

            COMMAND="ln -s "$CUDA_SOURCE_DIRECTORY"/"${REQUIRED_CUDA_FILES[$i]}" "$TARGET_FILE
            eval $COMMAND
        
            echo -n "created symbolic link "
            cecho $TARGET_FILE $CYAN

        fi

    done
fi

if [[ $SKIP_OPENGL -eq 0 ]] ; then
    for ((i=0; i<${#REQUIRED_OPENGL_FILES[*]}; i++)); do

        TARGET_FILE=$TARGET_DIRECTORY"/"${REQUIRED_OPENGL_FILES[$i]}

        if [ -e $TARGET_FILE ] ; then
        
            echo -n "skipping "
            cechon ${REQUIRED_OPENGL_FILES[$i]} $CYAN
            echo -n " because it already exists in "
            cecho $TARGET_DIRECTORY $YELLOW

        else

            COMMAND="ln -s "$OPENGL_SOURCE_DIRECTORY"/"${REQUIRED_OPENGL_FILES[$i]}" "$TARGET_FILE
            eval $COMMAND
        
            echo -n "created symbolic link "
            cecho $TARGET_FILE $CYAN

        fi

    done
fi

cecho "DONE!" $GREEN
echo ""



