#!/bin/bash

## Flag:
## -b to render coursebook using bookdown.

dobook=false
doslides=false

while getopts ":bs" opt; do
    case $opt in
	b)
	    dobook=true
	    ;;
	\?)
	    echo "invalid option: -${OPTARG}." >&2
	    ;;
    esac
done
folder=${!OPTIND}
## Throwing error for no argument (folder) provided.
if [ -z "$folder" ]; then
    echo "ERROR: No folder given." >&2
    exit 1
fi
## book
if [ "$dobook" = true ]; then
    echo "Rendering book"
    cd "$folder"
    if [ -a index.Rmd ]; then
	rmd_exists=true
    else
	rmd_exists=false
	echo "ERROR: .Rmd doesn't exist." >&2
	exit 2
    fi
    R -e "bookdown::render_book('index.Rmd',output_dir = './')" 2>&1 >/dev/null
    mv _bookdown_files/_main_files/ ./
    rmdir _bookdown_files/
fi
