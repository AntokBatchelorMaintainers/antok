#!/bin/bash

../build/bin/generateExampleFile
rm output.root -rf && ../build/bin/treereader example.root output.root ../config/example.yaml
