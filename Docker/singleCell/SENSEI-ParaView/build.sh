#!/bin/bash

# individual users must edit this
# environment variable to match their
# system or docker hub username
export IMAGENAME="singlecell-sensei-paraview:latest"

docker build -t $IMAGENAME .

