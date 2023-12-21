#!/bin/bash
cd BP
make clean
make
cd ../cooma
make clean
make A1B
cd ../FT
make clean
make TREE2
