#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"runge-kutta-fehlberg.h"

void print_trajectory(FILE* out, trajectory trj,
		      size_t* indices, size_t n_indices) {

  size_t i,j,k;

  for(i=0;i<trj.n_values;i++) {
    fprintf(out,"%12.6E\t",trj.t[i]);
    for(j=0;j<n_indices-1;j++) {
      fprintf(out,"%12.6E\t",trj.states_at_t[indices[j]][i]);
    }
    fprintf(out,"%12.6E\n",trj.states_at_t[indices[n_indices-1]][i]);
  }
}

size_t get_interval_for_t(double t, trajectory a) {
  size_t L = 0;
  size_t R = a.n_values-2;

  size_t midpoint;
  
  while(L<=R) {
    midpoint = (L+R)/2;
    if (a.t[midpoint] < t && a.t[midpoint+1] < t) {
      L = midpoint+1;
    } else if (a.t[midpoint] > t && a.t[midpoint+1] > t) {
      R = midpoint-1;
    } else {
      return (midpoint);
    }
  }
  return(-1);
}

void get_interpolition_coefficients(size_t position_index,
				    size_t trajectory_component_index,
				    trajectory a,
				    double *coeff) {

  size_t i,j,k;

  double m[16];
  double x_v[4];
  double x;

  double* b = coeff;

  double c;
  
  size_t start_pos = position_index;

  if (position_index == 0) start_pos = 1;
  if (position_index > a.n_values-3) start_pos = a.n_values-3;

  for(j=0;j<4;j++) {

    x = a.t[start_pos-1+j]+1;
    x_v[0] = 1;
    b[j] = a.states_at_t[trajectory_component_index][start_pos-1+j];
    for(i=1;i<4;i++) {
      x_v[i] = x_v[i-1]*x;
    }
    
    for(i=0;i<4;i++) {
      m[j*4+i] = x_v[i];
    }

  }

  /* solve m*r = b : results are stored in b ( gauss method ) */
  
  for(j=0;j<4;j++) {
    c = m[j*4+j];
    for(i=0;i<4;i++) {
      m[j*4+i] /= c;
    }
    b[j] /= c;
    for(i=0;i<4;i++) {
      c = m[i*4+j];
      for(k=0;k<4;k++) {
	if (i != j) {
	  m[i*4+k] -= c*m[j*4+k];
	}
      }
      if(i != j) {
	b[i] -= c*b[j];
      }
    }
  }
}

void trajectory_to_state_at_t(double t, trajectory a, double* state) {

  size_t i;

  for(i = 0; i< a.n_states;i++) {
    state[i] = get_interpolated_value_at_t(t,a,i);
  }
}

void trajectory_to_state_at_end(trajectory a, double* state) {

  size_t i;
  for(i = 0; i< a.n_states;i++) {
    state[i] = a.states_at_t[i][a.n_values-1];
  }
}
    
double get_interpolated_value_at_t(double t, trajectory a,
				   size_t trajectory_component_index) {

  size_t index = get_interval_for_t(t,a);

  double coeffs[4];

  double s = t+1;
  
  get_interpolition_coefficients(index,
				 trajectory_component_index,
				 a,
				 coeffs);

  return(coeffs[0]+s*coeffs[1]+s*s*coeffs[2]+s*s*s*coeffs[3]);
}

trajectory combine_trajectories_overwrite(trajectory a, trajectory b) {

  size_t i;

  size_t combined_size = a.n_values+b.n_values;

  a.t = (double*)realloc(a.t,sizeof(double)*(a.n_values+b.n_values));

  for(i=0;i<a.n_states;i++) {
  
    a.states_at_t[i] = (double*)realloc(a.states_at_t[i],sizeof(double)
					*combined_size);

    memcpy(a.states_at_t[i]+a.n_values,
	   b.states_at_t[i],sizeof(double)*b.n_values);
    free(b.states_at_t[i]);
  }
  free(b.states_at_t);

  memcpy(a.t+a.n_values,b.t,sizeof(double)*b.n_values);
  free(b.t);
 
  a.n_values = combined_size;

  return(a);
}
	 
static inline double run_state_eqn(size_t eqn_idx, state_eqn* state_equations,
				 double t, double* states) {

  state_eqn se = state_equations[eqn_idx];
  
  return(se.eqn(t, states, se.params));
}
  

trajectory rkf (double* states,
		state_eqn * state_equations,
		size_t n_states,
		double accepted_error,
		double initial_delta,
		double initial_t,
		double final_t) {

  size_t i,j;

  double delta = initial_delta;
  double t = initial_t;

  double e_tmp;
  double error;

  double a_one = 0.;
  double a_two = 1./4.;
  double a_three = 3./8.;
  double a_four = 12./13.;
  double a_five = 1.;
  double a_six = 0.5;

  double b_two_one = 1./4.;
  double b_three_one = 3./32.;
  double b_four_one = 1932./2197.;
  double b_five_one = 439./216.;
  double b_six_one = -8./27.;

  double b_three_two = 9./32. ;
  double b_four_two = -7200./2197.;
  double b_five_two = -8.;
  double b_six_two = 2.;

  double b_four_three = 7296./2197.;
  double b_five_three = 3680./513. ;
  double b_six_three = -3544./2565. ;

  double b_five_four = -845./4104. ;
  double b_six_four = 1859./4104.;

  double b_six_five = -11./40.;

  double c_one = 25./216.;
  double c_two = 0.;
  double c_three = 1408./2565.;
  double c_four = 2197./4104.;
  double c_five = -1./5.;

  double c_h_one = 16./135.;
  double c_h_two = 0.;
  double c_h_three = 6656./12825. ;
  double c_h_four = 28561./56430. ;
  double c_h_five = -9./50.;
  double c_h_six = 2./55.;

  double c_t_one = 1./360.;
  double c_t_two = 0.;
  double c_t_three = -128./4275.;
  double c_t_four = -2197./75240. ;
  double c_t_five = 1./50. ;
  double c_t_six = 2./55.;

  double * k_one = (double*)malloc(sizeof(double)*n_states);
  double * k_two = (double*)malloc(sizeof(double)*n_states);
  double * k_three = (double*)malloc(sizeof(double)*n_states);
  double * k_four = (double*)malloc(sizeof(double)*n_states);  
  double * k_five = (double*)malloc(sizeof(double)*n_states);
  double * k_six = (double*)malloc(sizeof(double)*n_states);
  
  double * st_buffer = (double*)malloc(sizeof(double)*n_states);

  double * probable_states = (double*)malloc(sizeof(double)*n_states);

  trajectory trj;
  
  trj.n_values = 1;
  trj.states_at_t = (double**)malloc(sizeof(double*)*n_states);
  trj.t = (double*)malloc(sizeof(double*)*n_states);

  for(j=0;j<n_states;j++) {
    trj.states_at_t[j] = (double*)malloc(sizeof(double)*10000);
    trj.states_at_t[j][0] = states[j];
  }
  trj.t = (double*)malloc(sizeof(double)*10000);
  trj.t[0] = initial_t;
  
  while(t<final_t) {
    error = 0;
    
    for(j=0;j<n_states;j++) {
      k_one[j] = delta*run_state_eqn(j,state_equations,
				     t+a_one*delta,states);
      st_buffer[j] = states[j]
	+b_two_one*k_one[j];
    }
    for(j=0;j<n_states;j++) {
      k_two[j] = delta*run_state_eqn(j,state_equations,
				     t+a_two*delta,st_buffer);
      st_buffer[j] = states[j]
	+b_three_one*k_one[j]
	+b_three_two*k_two[j];
    }
    for(j=0;j<n_states;j++) {
      k_three[j] = delta*run_state_eqn(j,state_equations,
				       t+a_three*delta,st_buffer);
      st_buffer[j] = states[j]
	+b_four_one*k_one[j]
	+b_four_two*k_two[j]
	+b_four_three*k_three[j];
    }
    for(j=0;j<n_states;j++) {
      k_four[j] = delta*run_state_eqn(j,state_equations,
				      t+a_four*delta,st_buffer);
      st_buffer[j] = states[j]
	+b_five_one*k_one[j]
	+b_five_two*k_two[j]
	+b_five_three*k_three[j]
	+b_five_four*k_four[j];
    }
    for(j=0;j<n_states;j++) { 
      k_five[j] = delta*run_state_eqn(j,state_equations,
				      t+a_five*delta,st_buffer);

      st_buffer[j] = states[j]
	+b_six_one*k_one[j]
	+b_six_two*k_two[j]
	+b_six_three*k_three[j]
	+b_six_four*k_four[j]
	+b_six_five*k_five[j];
    }
    for(j=0;j<n_states;j++) {
      k_six[j] = delta*run_state_eqn(j,state_equations,
				     t+a_six*delta,st_buffer);
    }

    for(j=0;j<n_states;j++) {
      e_tmp =
	states[j]
	+c_h_one*k_one[j]
	+c_h_two*k_two[j]
	+c_h_three*k_three[j]
	+c_h_four*k_four[j]
	+c_h_five*k_five[j]
	+c_h_six*k_six[j];
      
      probable_states[j] =
	states[j]
	+c_one*k_one[j]
	+c_two*k_two[j]
	+c_three*k_three[j]
	+c_four*k_four[j]
	+c_five*k_five[j];

      e_tmp -= probable_states[j];
      e_tmp = fabs(e_tmp);
      
      if(error < e_tmp) {
	error = e_tmp;
      }
    }
    
    if( accepted_error >= error ) {

      if(trj.n_values%10000 == 0) {
	 for(j=0;j<n_states;j++) {
	   trj.states_at_t[j] =
	     (double*)realloc(trj.states_at_t[j],
			      sizeof(double)*(trj.n_values+10000));
	 }
	 trj.t = (double*)realloc(trj.t,sizeof(double)*(trj.n_values+10000));
      }

      for(j=0;j<n_states;j++) {
	trj.states_at_t[j][trj.n_values] = states[j] = probable_states[j];
      }
      
      t += delta;
      trj.t[trj.n_values] = t;
      trj.n_values++;
    } 
    delta = delta*(double)pow(accepted_error/(2*error),0.25);
  }
  for(j=0;j<n_states;j++) {
    trj.states_at_t[j] = (double*)realloc(trj.states_at_t[j],
					  sizeof(double)*(trj.n_values));
  }

  trj.t = (double*)realloc(trj.t,
			   sizeof(double)*(trj.n_values));

  trj.n_states = n_states;
  
  free(k_one);
  free(k_two);
  free(k_three);
  free(k_four);
  free(k_five);
  free(k_six);
  free(st_buffer);
  free(probable_states);

  return(trj);
}
    
