#!/bin/bash

TIMES=5

# parameters
n=0
p=0
a=0
g=0
t=0

function run () {
    for ((t = 0; t < $TIMES; t++ )); do
        echo exp $t:
        ./target/BoleyEtAl "./data/data-$n-$p-$a-$g-$t.ptn" "./data/data-$n-$p-$a-$g-$t.grh" "./data/data-$n-$p-$a-$g-$t.id_phenotype" "./data/data-$n-$p-$a-$g-$t.id_population" 1 $a -outname "results/ebg-$n-$p-$a-$g-$t.csv" -common-outname "results/results.csv"
        echo
    done

}

function run_exp() {
    for gtmp in ${gs[@]}; do
        g=$gtmp
        run
    done
}

make clean
make ./target/BoleyEtAl

mkdir results

# Speed 1
echo -e "\e[31mSpeed1:\e[m"; echo

# configure
n=100
p=0.5
a=0.05

gs=(100 200 500 1000 2000 5000 10000 20000)
# gs=(100 200)

run_exp

# Speed 2
echo; echo -e "\e[31mSpeed2:\e[m"; echo

# configure
n=50
p=0.5
a=0.05

gs=(100 200 500 1000 2000 5000 10000 20000)
# gs=(100 200)

run_exp

# Speed 3
echo; echo -e "\e[31mSpeed3:\e[m"; echo

# configure
n=100
p=0.2
a=0.05

gs=(100 200 500 1000 2000 5000 10000 20000)
# gs=(100 200)

run_exp

# Speed 4
echo; echo -e "\e[31mSpeed4:\e[m"; echo

# configure
n=100
p=0.2
a=0.0001

gs=(100 200 500 1000 2000 5000 10000 20000)
# gs=(100 200)

run_exp