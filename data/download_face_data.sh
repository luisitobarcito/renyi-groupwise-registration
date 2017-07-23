#!/bin/bash

## download cifar10 dataset from Krivshensky site

DOWNLOAD_DIR=./

wget -v http://www2.imm.dtu.dk/~aam/datasets/face_data.zip
unzip -d face_data  face_data.zip
