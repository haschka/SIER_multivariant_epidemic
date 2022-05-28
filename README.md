# A multivariant, multi-class (i.e. ages) SEIR epidemic simulator

This model simulates a multivariant SEIR epidemic with cross immunity between variants, a
preimmune population and an ongoing vaccination campain, the later is untested. 

In order to use this program you need a c compilier and lapack libraries, i.e. http://www.netlib.org/lapack/

### Compilation i.e. with gcc
```
gcc -O3 -march=native execute_model.c epidemic.c rkf.c parser.c -lm -llapack -o execute-model
```

### Epidemic simulation
Performing an epidemic simulation requires a prepared input file in keyword-value format.
An example can be found in example.input. 
The simulator generates trajectories written to standard out in tab seperated format.
The first column referes to the simulation time,
follwing by the S I E R compartments.
Depending on how many age classes, variants you simulate you will encounter multiple columns.
```
./execute-model sample.input
```
should yield an epidemic trajectory. An a description of the model implemented as well as
of the input files parameters is found in [documentation/doc.pdf](https://raw.githubusercontent.com/haschka/SIER_multivariant_epidemic/main/documentation/doc.pdf)

### Fitting a model to data
This is experimental and was only developed for in house in modelling the SARS-CoV-2 
epedemic.

First compile the fit program:
```
gcc -O3 -march=native fit-model.c epidemic.c rkf.c parser.c -lm -llapack -o fit-model
```
You can then fit to the example data:
```
fit-model fit-sample.input fit-sample.table fit-sample.out
```
