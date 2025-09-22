# Enum-Support-Closed-CIS
A fast algorithm for enumerating statistically significant support-closed subgraphs

The programs are used in:
Daiki Watanabe, Takumi Tada, Kazuya Haraguchi: Comparison of Algorithms and Implementations for Enumerating Support-closed Connected Induced Subgraphs. Journal of Information Processing, Vol.32, 2024, pp. 894-902. https://doi.org/10.2197/ipsjjip.32.894

To cite this repository, refer to the above paper. 

## Repository hierarchy

- RandomExp: enumerating solutions of random graph
- BacterialGWAS: directory for experiment of enumerating significant subgraphs
- commands: commands for computational experiments

## Experiment results

Random Graph Exp: `RandomExp/BA/result.zip` and `RandomExp/ER/result.zip`

Bacterial GWAS Exp: `BacterialGWASExp/result.zip`

## execution

### for Random Graph Exp

### Compile:

```
./commands/build_RG.sh
```

### for Bacterial-GWAS

#### Run Exp:
```
$ ./commands/run_BGWAS.sh
```

#### Compile and Run binary:
```
$ ./commads/build_BGWAS.sh
$ ./target/BeleyEtAl # with options
```

options:
```
list of params:
  -seed <INT> ... random seed (1)
  -sigma <INT> ... minimum size of connector (1)
  -delta <DOUBLE> ... minimum density of a connector (0)
  -tlim <DOUBLE> ... time limit in sec; if <0, it is inf. (-1)
  -reduce <BOOL> ... whether reduction is invoked (1)
  -vtable <STR> ... file that describes relation between vertex id & name; itable should be also specified
  -itable <STR> ... file that describes relation between item id & name; vtable should be also specified
  -outname <STR> ... filename of output from vtable & itable (out.txt)
  -distname <STR> ... filename for distribution
  -ramub <INT> ... upper bound on used amount of RAM (-1)
  -common-outname <STR> ... filename of common output from vtable & itable (results.csv)
```

CALDERA: https://github.com/HectorRDB/Caldera_ISMB
