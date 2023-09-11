# Timing experiments benchmarks:
python2 random_dnf_generator.py 100000 100001 2 2 10 2 1 3 50 10 20 /path/to/dir/benchmarks/random/dens_great_1/ 0
python2 random_dnf_generator.py 100000 100001 2 0.1 1 0.2 1 3 50 10 20 /path/to/dir/benchmarks/random/dens_less_1/ 0

# Accuracy experiments benchmarks:

## dens_less_1
python2 random_dnf_generator.py 100 1000 300 0.3 1 0.2 0 0.03 0.75 0.16 3 accurancy_exp/1 0
python2 random_dnf_generator.py 1000 10000 3000 0.3 1 0.2 0 0.03 0.37 0.08 3  accurancy_exp/2 0

## dens_great_1
python2 random_dnf_generator.py 100 1000 300 1 9 2 0 0.03 0.75 0.16 3 accurancy_exp/3 0
python2 random_dnf_generator.py 1000 10000 3000 1 9 2 0 0.03 0.37 0.08 3  accurancy_exp/4 0
