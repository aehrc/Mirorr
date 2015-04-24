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
MILXVIEW_SOURCE_DIRECTORY="*milxRegistration/NPWopengl";
SOURCE_DIRECTORY=`locate $MILXVIEW_SOURCE_DIRECTORY`

if [ -z ${SOURCE_DIRECTORY[0]} ] ; then
    cecho "*** ERROR *** : could not find source directory \""$MILXVIEW_SOURCE_DIRECTORY"\". Have you checked out milx-view somewhere?" $RED
    exit 0
fi

#if there are more, use the one with the most recently modified file in it
if [ ${#SOURCE_DIRECTORY[@]} -ge 2 ] ; then
    cecho "*** WARNING *** : multiple possible source directories were found:" $MAGENTA

    FIRST_DIR=`locate *milxRegistration | head -n1`

    MOST_RECENT_DIR=$FIRST_DIR
    MOST_RECENT_FILE=`ls $FIRST_DIR -t | head -n1`

    for CUR_DIR in $SOURCE_DIRECTORY ; do

        echo "  "$CUR_DIR

        MOST_RECENT_FILE_HERE=`ls "$CUR_DIR" -t | head -n1`
        if test $CUR_DIR"/"$MOST_RECENT_FILE_HERE -nt $MOST_RECENT_DIR"/"$MOST_RECENT_FILE ; then
            MOST_RECENT_FILE=$MOST_RECENT_FILE_HERE
            MOST_RECENT_DIR=$CUR_DIR
        fi
        
    done

   cechon "using the most recent directory: " $MAGENTA
   cecho $MOST_RECENT_DIR $YELLOW

   SOURCE_DIRECTORY=$MOST_RECENT_DIR
fi

cechon "found npw opengl source directory: " $GREEN
cecho $SOURCE_DIRECTORY $YELLOW

if [ $SOURCE_DIRECTORY == $TARGET_DIRECTORY ] ; then
    cecho "*** ERROR *** : the target directory is the same as the source directory!" $RED
fi



#
# check if all files are present
#
REQUIRED_FILES=( MutualInformation.h NPWConstants.h NPWGeometry.cxx  NPWGeometry.h NPWGLImage.cxx NPWGLImage.h NPWJointHistogramCreator.cxx NPWJointHistogramCreator.h NPWOpenGLRenderer.cxx NPWOpenGLRenderer.h NPWStopwatch.h NPWVertices.cxx NPWVertices.h NPWWindow.cxx NPWWindow.h Vector2.h )
SOURCE_FILES=`ls $SOURCE_DIRECTORY`
MISSING_FILES=()

for ((i=0; i<${#REQUIRED_FILES[*]}; i++)); do

    FILEWASFOUND=0

    for FOUND_FILE in $SOURCE_FILES ; do
    #for ((j=0; j<${#SOURCE_FILES[*]}; j++)); do
        if [ "$FOUND_FILE" == "${REQUIRED_FILES[$i]}" ] ; then
            FILEWASFOUND=1
        fi
    done

    if [ $FILEWASFOUND == 0 ] ; then
        MISSING_FILES[${#MISSING_FILES[*]}]=${REQUIRED_FILES[$i]}
    fi

done

if [ ${#MISSING_FILES[@]} != 0 ] ; then
    cecho "*** ERROR *** : not all source files are present!" $RED
    echo ${#MISSING_FILES[@]}" missing files:"
    for ((i=0; i<${#MISSING_FILES[*]}; i++)); do
        cecho ${MISSING_FILES[$i]} $RED
    done
else
    cecho "all "${#REQUIRED_FILES[@]}" required source files were found" $GREEN
fi


#
# create the links, it they don't already exist
#
echo -n "creating symbolic links in "
cecho $TARGET_DIRECTORY $YELLOW

for ((i=0; i<${#REQUIRED_FILES[*]}; i++)); do

    TARGET_FILE=$TARGET_DIRECTORY"/"${REQUIRED_FILES[$i]}

    if [ -e $TARGET_FILE ] ; then
    
        echo -n "skipping "
        cechon ${REQUIRED_FILES[$i]} $CYAN
        echo -n " because it already exists in "
        cecho $TARGET_DIRECTORY $YELLOW

    else

        COMMAND="ln -s "$SOURCE_DIRECTORY"/"${REQUIRED_FILES[$i]}" "$TARGET_FILE
        eval $COMMAND
    
        echo -n "created symbolic link "
        cecho $TARGET_FILE $CYAN

    fi

done

cecho "DONE!" $GREEN
echo ""



