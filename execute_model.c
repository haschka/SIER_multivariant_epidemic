#include<stdio.h>
#include<stdlib.h>
#include"runge-kutta-fehlberg.h"
#include"epidemic.h"

int main(int argc, char** argv) {

  size_t i,j;

  FILE* input_file;

  model_parameters mp;

  trajectory trj;

  size_t* print_indices;
  
  if ( NULL == (input_file = fopen(argv[1],"r"))) {
    printf("Could not open input file %s\n",argv[1]);
    return(1);
  }

  parse_model_properties(input_file,&mp);

  trj = epidemic_execute(mp);
    
  print_indices = (size_t*)malloc(sizeof(size_t)*trj.n_states);

  for(i=0;i<trj.n_states;i++) {
    print_indices[i] = i;
  }
  
  print_trajectory(stdout,trj,print_indices,trj.n_states); 

}
