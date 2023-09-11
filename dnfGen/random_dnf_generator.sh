#!/bin/bash

mkdir -p accuracy_exp && cd accuracy_exp
mkdir -p 1 && mkdir -p 2 && mkdir -p 3 && mkdir -p 4 && mkdir -p 5
cd ..
# #python2 random_dnf_generator.py 100 1000 300 0.3 1 0.2 0 0.03 0.75 0.16 3 accuracy_exp/1 0
# #python2 random_dnf_generator.py 1000 10000 3000 0.3 1 0.2 0 0.03 0.37 0.08 3  accuracy_exp/2 0

# #python2 random_dnf_generator.py 100 1000 300 1 9 2 0 0.03 0.75 0.16 3 accuracy_exp/3 0
# #python2 random_dnf_generator.py 1000 10000 3000 1 9 2 0 0.03 0.37 0.08 3  accuracy_exp/4 0

python2 random_dnf_generator.py 200 1000 300 0.3 1 0.2 0 0.03 0.75 0.16 3 accuracy_exp/5 0

