#!/bin/bash

TIMES=1

function convert () {
    python3 analysis/convert_CALDERAinstance_to_ours/convert_all_rows_to_ptn.py "./data/tmp$1.binary" ./data/data.ptn
    python3 analysis/convert_CALDERAinstance_to_ours/convert_dbg.py "./data/tmp$1.edges.dbg" ./data/data.grh
}

function run () {
    for ((i = 0; i < $TIMES; i++ )); do
        echo exp $i:
        convert $i
        ./target/BoleyEtAl './data/data.ptn' './data/data.grh' "./data/tmp$i.id_phenotype" "./data/tmp$i.id_population" 1 "$1"
        # ./target/BoleyEtAl './data/data.ptn' './data/data.grh' "./data/tmp$i.id_population" "./data/tmp$i.id_phenotype"  1 "$1"
        echo
    done

}

make clean
make ./target/BoleyEtAl

# Speed 1
echo -e "\e[31mSpeed1:\e[m"; echo

python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 100 -t $TIMES
run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 200 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 500 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 1000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 2000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 5000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 10000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.5 -a 0.05 -g 20000 -t $TIMES
# run '0.05' $TIMES

# Speed 2
echo; echo -e "\e[31mSpeed2:\e[m"; echo
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 100 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 200 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 500 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 1000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 2000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 5000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 10000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 50 -p 0.5 -a 0.05 -g 20000 -t $TIMES
# run '0.05' $TIMES

# Speed 3
echo; echo -e "\e[31mSpeed3:\e[m"; echo
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 100 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 200 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 500 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 1000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 2000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 5000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 10000 -t $TIMES
# run '0.05' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.05 -g 20000 -t $TIMES
# run '0.05' $TIMES

# Speed 4
echo; echo -e "\e[31mSpeed4:\e[m"; echo
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 100 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 200 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 500 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 1000 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 2000 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 5000 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 10000 -t $TIMES
# run '0.0001' $TIMES
# python3 analysis/exp-caldera/main.py -n 100 -p 0.2 -a 0.0001 -g 20000 -t $TIMES
# run '0.0001' $TIMES