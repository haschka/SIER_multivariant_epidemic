#define _GNU_SOURCE

#include<unistd.h>
#include<stdio.h>
#include<string.h>
#include<float.h>
#include<limits.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"runge-kutta-fehlberg.h"
#include"epidemic.h"

typedef struct {
  double** points;
  double* t;
  size_t n_values;
  size_t n_variants;
  double preinfected;
  double age_classes[3];
  double total_population;
  double r_first_question[3];
  double r_first;
  double Inc_first_question[3];
  double Inc_first;
  int windowsize;
  double start_date_offset;
} track;

typedef struct {
  track* tracks;
  size_t n_tracks;
} track_dataset;

typedef struct {
  trajectory* trjs;
  size_t n_trjs;
} trajectory_dataset;

double measure_error(track_dataset tra, trajectory_dataset trjs,
		     size_t* trajectory_compartment_indices,
		     double* buffer_n_variant_size, size_t age_groups,
		     double* conversion_rate) {

  size_t i,j,k,l;

  double error = 0.;

  double** track_points; 
  double* track_times;

  double total_track, percentage_track;
  double total_trj, percentage_trj;

  for(i=0;i<tra.n_tracks;i++) {
    track_points = tra.tracks[i].points;
    track_times = tra.tracks[i].t;
    for(k=0;k<tra.tracks[i].n_values;k++) {
      total_track = 0;
      total_trj = 0;
      for(j=0;j<tra.tracks[i].n_variants;j++) {
	total_track +=  track_points[j][k];
        buffer_n_variant_size[j] = 0;
	for(l=0;l<age_groups;l++) {
	  
	  buffer_n_variant_size[j] +=
	    get_interpolated_value_at_t(track_times[k],
					trjs.trjs[i],
					trajectory_compartment_indices[j*
								    age_groups
								       +l])
	    *conversion_rate[j];
	}
	total_trj += buffer_n_variant_size[j];
      }
      for(j=0;j<tra.tracks[i].n_variants;j++) {
	if (total_track != 0 && total_trj != 0) {
	  error +=
	    fabs(track_points[j][k]/total_track
		 - buffer_n_variant_size[j]/total_trj);
	} else {
	  error += 1.;
	}
      }
    }
  }
  //printf("ERROR: %lf \n",error);
  return(error);
}

int determine_Inc_R_at_departure(track* t) {

  int i;
  struct tm times[3];
  struct tm initial_time;

  time_t times_t[3];
  time_t initial_time_t;

  double delta_time[3];

  double interval_time_delta;
  double delta_R;
  double delta_Inc;
  double delta_t;
  
  const char* dates[] = { "2021-12-06", "2021-12-13", "2021-12-20" };

  memset(&initial_time,0,sizeof(struct tm));
  strptime("2021-12-01","%Y-%m-%d", &initial_time);
  initial_time_t = mktime(&initial_time);
  
  for(i=0;i<3;i++) {
    memset(times+i,0,sizeof(struct tm));
    strptime(dates[i],"%Y-%m-%d", times+i);
    times_t[i] = mktime(times+i);

    delta_time[i] = ((double)times_t[i]- (double)initial_time_t)/86400;
  }

  /* Handle edge cases */
  
  if(t->start_date_offset <= delta_time[0]) {
    t->r_first = t->r_first_question[0];
    t->Inc_first = t->Inc_first_question[0];
    goto Inc_R_determinated;
  }

  if(t->start_date_offset == delta_time[1]) {
    t->r_first = t->r_first_question[1];
    t->Inc_first = t->Inc_first_question[1];
    goto Inc_R_determinated;
  }
    
  if(t->start_date_offset >= delta_time[2]) {
    t->r_first = t->r_first_question[2];
    t->Inc_first = t->Inc_first_question[2];
    goto Inc_R_determinated;
  }

  /* Handle inbetween cases */
  
  for(i=0;i<2;i++) {
  
    if (delta_time[i] < t->start_date_offset &&
	t->start_date_offset < delta_time[i+1]) {
      
      interval_time_delta = delta_time[i] - t->start_date_offset;
      
      delta_R = t->r_first_question[i] - t->r_first_question[i+1];
      delta_Inc = t->Inc_first_question[i] - t->Inc_first_question[i+1];
      
      delta_t = delta_time[i+1]-delta_time[i];
      
      t->r_first =
	t->r_first_question[i] + delta_R*(interval_time_delta/delta_t);
      t->Inc_first =
	t->Inc_first_question[i] + delta_Inc*(interval_time_delta/delta_t);

      goto Inc_R_determinated;
    }
    
  }
  
 Inc_R_determinated:
  return(0);
  
}
  
void read_track_point_file(char* filename, track* tracks,
			   size_t n_variants, int windowsize) {

  FILE* trackfile = fopen(filename, "r");

  size_t i;
  
  char * line = NULL;
  size_t line_size = 0;

  size_t n_values = 0;

  char* timebuffer = (char*)malloc(sizeof(char)*11);

  struct tm initial_time;
  time_t initial_time_t;

  size_t check_variant_scan;
  
  struct tm current_time;
  time_t current_time_t;

  int line_scan_offset;
  int offset_adder;

  size_t midpoint;

  track replacement_track;
  double initial_replacement_time;
  
  double* variant_percentage = (double*)malloc(sizeof(double)*n_variants);

  memset(&initial_time,0,sizeof(struct tm));
  memset(&current_time,0,sizeof(struct tm));
      
  rewind(trackfile);
  strptime("2021-12-01","%Y-%m-%d", &initial_time);
  /*  strptime("2021-02-15","%Y-%m-%d", &initial_time); */
  initial_time_t = mktime(&initial_time);
  
  tracks->points = (double**)malloc(sizeof(double*)*n_variants);
  for(i=0;i<n_variants;i++) {
    tracks->points[i] = (double*)malloc(sizeof(double)*25);
  }
  tracks->t = (double*)malloc(sizeof(double)*25);

  getline(&line,&line_size,trackfile);
  
  while(-1 != getline(&line,&line_size,trackfile)) {
    check_variant_scan = 0;
    if( 1 == sscanf(line,"%*s\t%s%n",timebuffer,&line_scan_offset)) {
      variant_percentage[0] = 1;
      for(i=1;i<n_variants;i++) {
	check_variant_scan += sscanf(line+line_scan_offset,"%lf\t%n",
				     variant_percentage+i,&offset_adder);
	variant_percentage[i]/=100.;
	line_scan_offset+=offset_adder;
	variant_percentage[0] -= variant_percentage[i];
      }
      if (check_variant_scan == n_variants-1) {
	
	if((n_values+1)%25 == 0) {
	
	  for(i=0;i<n_variants;i++) {
	  
	    tracks->points[i] =
	      (double*)realloc(tracks->points[i], sizeof(double)*(n_values+25));
	  }
	  tracks->t =
	    (double*)realloc(tracks->t, sizeof(double)*(n_values+25));
	}
      
	for(i=0;i<n_variants;i++) {
	  tracks->points[i][n_values] = variant_percentage[i];
	}
      
	strptime(timebuffer,"%Y-%m-%d", &current_time);
	current_time_t = mktime(&current_time);
      
	tracks->t[n_values] =
	  ((double)current_time_t - (double)initial_time_t)/86400.;
      
	n_values++;
      }
    }	
  }

  for(i=0;i<n_values;i++) {
    if(tracks->points[0][i] < 0.5) {
      midpoint = i;
      goto found_midpoint;
    }
  }

  /* if midpoint is not found */
  midpoint = 32;
  
  fprintf(stderr,"Warning! No midpoint found in dataset set midpoint to 0+32");
  
 found_midpoint:

  if( (long)((long)midpoint-(long)windowsize) >= 0 &&
      midpoint+windowsize <=n_values) {
    replacement_track.t = (double*)malloc(sizeof(double)*2*windowsize);
    replacement_track.points = (double**)malloc(sizeof(double*)*n_variants);
    memcpy(replacement_track.t,
	   tracks->t+(midpoint-windowsize),
	   sizeof(double)*2*windowsize);

    free(tracks->t);
    initial_replacement_time = replacement_track.t[0];
    printf("Inital_time: %lf for track %s", initial_replacement_time, filename);

    for(i=0;i<2*windowsize;i++) {
      replacement_track.t[i] -= initial_replacement_time;
    }
    tracks->t = replacement_track.t;
    
    for(i=0;i<n_variants;i++) {
      replacement_track.points[i] =
	(double*)malloc(sizeof(double)*2*windowsize);
      memcpy(replacement_track.points[i],
	     tracks->points[i]+(midpoint-windowsize),
	     sizeof(double)*2*windowsize);
      free(tracks->points[i]);
      tracks->points[i] = replacement_track.points[i];
    }
  } else {
    fprintf(stderr,"Warning! Midpoint found out of range"
	    "resetting to generic point\n");
    midpoint = 32;
    goto found_midpoint;
  }
  tracks->start_date_offset = initial_replacement_time;
  fclose(trackfile);
  free(line);
  free(timebuffer);
  tracks->n_values = 2*windowsize;
}

int initialize_gradient_state(size_t index,
			      model_parameters* mp,
			      double* state) {

  size_t n_variants = mp->n_variants;
  size_t n_classes = mp->n_classes;

  size_t n_variants_x_n_classes = n_variants*n_classes;

  size_t i,j;

  size_t epsilon_index;

  double sum_I_classes, sum_E_classes, sum_I_variants, sum_E_variants; 
  
  /* Rzero fit variable */
  
  if(index >= 0 && index < n_variants-1) {
    state[index] = mp->Rzero[index+1]; 
    return(0);
  }

  /* varaint distribution */
  
  if(index >= (n_variants-1) && index < (n_variants+(n_variants-1)))  {
    /*
    for(j=0;j<n_variants;j++) {
      sum_I_classes = 0;
      sum_E_classes = 0;
      for(i=0;i<n_classes;i++) {
	sum_I += mp->initial_infected_fraction[i*n_variants+j];
	sum_E += mp->initial_intermediate_fraction[i*n_variants+j];
      }
      
    for(j=0;j<n_variants;j++) {
	state[n_variants-1+j] =
	  mp->initial_infected_fraction[i*n_variants+j]/sum_I;
        state[n_variants-1+j] =
	  mp->initial_intermediate_fraction[i*n_variants+j]/sum_E;
      } 


    */
    state[index] = 1.;
    
    return(0);
  }

  /* lagrange modifier to keep sum of variant distribution to 1 */
  
  if(index == 2*n_variants-1) {
    state[index] = 1.;
    return(0);
  }

  /* epsilon indices */
  
  if(index >= 2*n_variants &&
     index < 2*n_variants+n_variants*(n_variants-1)/2) {
    
    epsilon_index = 0;
    
    for(i=1;i<n_variants;i++) {
      for(j=0;j<i;j++) {
	if (epsilon_index == (index-2*n_variants)) {
	  state[index] = mp->epsilon[i*n_variants+j];
	  state[index] = mp->epsilon[j*n_variants+i];
	}
	epsilon_index++;
      } 
    } 
    return(0);
  }

  /* langrage modifiers to constrain : 2* epsilon indices to 0 and 1 */
  
  if(index >= 2*n_variants+n_variants*(n_variants-1)/2 &&
     index < 2*n_variants+n_variants*(n_variants-1)/2
     +n_variants*(n_variants-1)) {

    state[index] = 1.;
    return(0);
  } 
  return(1);
}

int gradient_state_to_mp(size_t index, model_parameters* mp,
			 double** initial_intermediate_fraction,
			 double** initial_infected_fraction,
			 double value,
			 double* state) {

  size_t n_variants = mp->n_variants;
  size_t n_classes = mp->n_classes;

  size_t n_variants_x_n_classes = n_variants*n_classes;

  size_t i,j;

  size_t local_index;
  
  double sum_I;
  double sum_E;

  double *buffer;
  
  size_t epsilon_index;
  
  if(index >= 0 && index < n_variants-1) {
    mp->Rzero[index+1] = fabs(value);
    return(0);
  }

  if(index >= (n_variants-1) && index < (2*(n_variants-1))) {
    local_index = (index-(n_variants-1))+1;
    for(i=0;i<n_classes;i++) {
      mp->initial_infected_fraction[i*n_variants+local_index] =
	initial_infected_fraction[0][i*n_variants+local_index]*value;
      mp->initial_intermediate_fraction[i*n_variants+local_index] =
	initial_intermediate_fraction[0][i*n_variants+local_index]*value;
    }
    return(0);
  }
  /*
  if(index >= 2*n_variants &&
     index < 2*n_variants+n_variants*(n_variants-1)/2) {
        
    epsilon_index = 0;

    for(i=1;i<n_variants;i++) {
      for(j=0;j<i;j++) {
	if (epsilon_index == (index-2*n_variants)) {
	  mp->epsilon[i*n_variants+j] = fabs(value);
	  mp->epsilon[j*n_variants+i] = fabs(value);
	}
	epsilon_index++;
      }
    } 

    return(0); 
  }
  */
  return(1);
}

void free_trajectory_dataset(trajectory_dataset trjds) {

  size_t i,j;

  for(i=0;i<trjds.n_trjs;i++) {
    free(trjds.trjs[i].t);
    for(j=0;j<trjds.trjs[i].n_states;j++) {
      free(trjds.trjs[i].states_at_t[j]);
    }
    free(trjds.trjs[i].states_at_t);
  }
}



track_dataset read_tracks(FILE* in, size_t n_variants,
			  size_t n_classes) {

  char * line = NULL;
  size_t line_size = 0;

  int i;

  track* tracks=NULL;

  double preinf_buffer;

  double Inc_first_buffer[3];
  double r_first_buffer[3];
  
  double age_class_fraction_buffer[3];
  
  double total_population_buffer;

  char track_filename[PATH_MAX];
  
  size_t n_tracks =0;

  int windowsize_buffer;
  
  track_dataset tds;
  
  rewind(in);
  getline(&line,&line_size,in);

  while(-1 != getline(&line,&line_size,in)) {
    if (10 == sscanf(line,
		     "%*i\t%s\t%lf\t%lf\t%lf\%lf\t%lf\%lf\t%lf\t%i\t%lf",
		     track_filename,
		     &preinf_buffer,
		     r_first_buffer+0,
		     r_first_buffer+1,
		     r_first_buffer+2,
		     Inc_first_buffer+0,
		     Inc_first_buffer+1,
		     Inc_first_buffer+2,		
		     &windowsize_buffer,
		     &total_population_buffer)) {

      tracks = (track*)realloc(tracks,sizeof(track)*(n_tracks+1));
      tracks[n_tracks].preinfected = preinf_buffer;
      
      if(n_classes == 1) {
	tracks[n_tracks].age_classes[0] = total_population_buffer;
	for(i=0;i<3;i++) {
	  tracks[n_tracks].Inc_first_question[i] = Inc_first_buffer[i];
	  tracks[n_tracks].r_first_question[i] = r_first_buffer[i];
	}
	tracks[n_tracks].windowsize = windowsize_buffer;
      } else {
	fprintf(stderr,
		"Fatal: This version only supports a single age class!\n");
	_exit(1);
	/*
	for(i=0;i<n_classes;i++) {
	  tracks[n_tracks].age_classes[i] = age_class_fraction_buffer[i];
	}
	*/
      }
      tracks[n_tracks].total_population = total_population_buffer;
      tracks[n_tracks].n_variants = n_variants;
      read_track_point_file(track_filename, tracks+n_tracks,
			    n_variants, windowsize_buffer);
      n_tracks++;

    } else {
      fprintf(stderr,"Fatal! Malformed Track File.");
      _exit(1);
    }
  }

  tds.tracks = tracks;
  tds.n_tracks = n_tracks;
  free(line);
  return(tds);
}

void execute_epidemic_with_track_dataset_conditions(trajectory_dataset* trjds,
						    model_parameters* mp,
						    track_dataset* tds) {

  size_t i,j;
  
  for(i=0;i<tds->n_tracks;i++) {
    mp->n_variants = tds->tracks[i].n_variants;
    for(j=0;j<mp->n_classes;j++) {
      mp->pre_immune_fraction[j] = tds->tracks[i].preinfected/100.;

    }

    if (mp->n_classes > 1) {
   /*      mp->population_class_distribution = tds->tracks[i].age_classes; */
      fprintf(stderr,
	      "Fatal! This code does not run with multiple age classes\n");
      _exit(1);
    } else {
      mp->population_class_distribution[0] = 1.;
    }
    mp->total_population = tds->tracks[i].total_population;
    mp->Rzero[0] = tds->tracks[i].r_first;
    

    if(trjds->trjs[i].n_values != 0) {
      for(j=0;j<trjds->trjs[i].n_states;j++) {
	free(trjds->trjs[i].states_at_t[j]);
      }
      free(trjds->trjs[i].states_at_t);
      free(trjds->trjs[i].t);
    }
    
    trjds->trjs[i] = epidemic_execute(mp[0]);
  }
}

double lagrange_error_addition(double* state, size_t n_variants) {

  double sum = 0;
  size_t i, start_eps_constraint_one, start_eps_constraint_two;

  double lagrange = 0.;
  
  for(i=0;i<n_variants-1;i++) {
    sum += state[n_variants-1+i];
  }
  
  lagrange -= state[2*n_variants-1]*(1-sum);

  start_eps_constraint_one = 2*n_variants+n_variants*(n_variants-1)/2;
  start_eps_constraint_two = 2*n_variants+n_variants*(n_variants-1);
  for(i=0;i< n_variants*(n_variants-1)/2;i++) {
    lagrange -= state[start_eps_constraint_one+i]*state[2*n_variants+i];
    lagrange += state[start_eps_constraint_two+i]*(state[2*n_variants+i]-1);
  }
  return(0);
} 

void test_line_search(size_t n_gradient_variables, double* gradient,
		      double** initial_intermediate_fraction,
		      double** initial_infected_fraction,
		      double* gradient_state,
		      double* possible_new_gradient_state,
		      double* current_step_size,
		      double* current_error,
		      double previous_step_error,
		      model_parameters* mp,
		      trajectory_dataset* trjds,
		      track_dataset* tds,
		      size_t* compartment_indices,
		      double* buffer_n_variant_size) {
  
  size_t j, failed_error_count;
  
  for(j = 0;j<n_gradient_variables;j++) {
    possible_new_gradient_state[j] =
      gradient_state[j]-current_step_size[0]*gradient[j];
  }
  for(j = 0;j<n_gradient_variables;j++) {
    gradient_state_to_mp(j,mp,initial_intermediate_fraction,
			 initial_infected_fraction,
			 possible_new_gradient_state[j],
			 possible_new_gradient_state);
  }
  
  execute_epidemic_with_track_dataset_conditions(trjds,
						 mp,
						 tds);
  
  
  current_error[0] = measure_error(tds[0],trjds[0],compartment_indices,
				   buffer_n_variant_size,mp->n_classes,
				   mp->conversion_rate);

    current_error[0] +=
     lagrange_error_addition(possible_new_gradient_state,mp->n_variants);

    //      current_error[0] *= current_error[0];
  
  if(current_error[0] < previous_step_error) { 
    current_step_size[0]*=1.1;
    printf("current step size : %12.6E\n", current_step_size[0]);

    printf("New Gradient ");
    for(j = 0;j<n_gradient_variables;j++) {
      gradient_state[j] = possible_new_gradient_state[j];
      printf(" %lf ",gradient_state[j]);
    }
    printf("\n");
  
    return;
  } else {
    failed_error_count = 0;
    while (current_error[0] >= previous_step_error) {
      
      current_step_size[0]*=0.1;
      printf("current step size : %12.6E\n", current_step_size[0]);

      for(j = 0;j<n_gradient_variables;j++) {
	possible_new_gradient_state[j] =
	  gradient_state[j]-current_step_size[0]*gradient[j];
      }
      for(j = 0;j<n_gradient_variables;j++) {
	gradient_state_to_mp(j,mp,
			     initial_intermediate_fraction,
			     initial_infected_fraction,
			     possible_new_gradient_state[j],
			     possible_new_gradient_state);
      }
      
      execute_epidemic_with_track_dataset_conditions(trjds,
						     mp,
						     tds);
  

      printf("Too far, decr step size: "
	     "Error - current %12.6E - previous %12.6E :\n",
	     current_error[0], previous_step_error);
      
      current_error[0] = measure_error(tds[0],trjds[0],compartment_indices,
				       buffer_n_variant_size,mp->n_classes,
				       mp->conversion_rate);

      current_error[0] +=
	lagrange_error_addition(possible_new_gradient_state,mp->n_variants);
      
      //           current_error[0] *= current_error[0];

      failed_error_count++;
      if(failed_error_count == 4) {
	break;
      }	
    }
    printf("New Gradient ");
    for(j = 0;j<n_gradient_variables;j++) {
      gradient_state[j] = possible_new_gradient_state[j];
      printf(" %lf ",gradient_state[j]);
    }
    printf("\n");
    return;
  }
}


						    
int main(int argc, char** argv) {

  size_t i,j,k;

  FILE* input_file;
  FILE* track_file;
  FILE* fitted_model;
  
  model_parameters mp;

  trajectory_dataset trjds;

  size_t* print_indices;

  size_t* compartment_indices;
  
  track_dataset tds;

  double initialization_sum;
  
  double machine_epsilon_sqrt = 1.4832396974191326e-08 ;
  double* gradient_state;
  double* evaluative_gradient_state;
  double* possible_new_gradient_state;
  double* gradient;
  size_t n_gradient_variables;
  double derivative_delta, forward_error, backward_error;

  double epsilon=10e-10;
  double error=1.;
  double this_step_error, previous_step_error = DBL_MAX;

  double gradient_sum;
  
  double errors[30];
  size_t gradient_iterations = 0;
  
  double current_step_size = 0.0000001;

  double* buffer_n_variant_size;

  double** initial_infected_fraction;
  double** initial_intermediate_fraction;

  if ( NULL == (input_file = fopen(argv[1],"r"))) {
    printf("Could not open input file %s\n",argv[1]);
    return(1);
  }
  if ( NULL == (track_file = fopen(argv[2],"r"))) {
    printf("Could not open input file %s\n",argv[2]);
    return(1);
  }
  if ( NULL == (fitted_model = fopen(argv[3],"w"))) {
    printf("Could not open output file %s\n",argv[3]);
    return(1);
  }

  parse_model_properties(input_file,&mp);

  n_gradient_variables =
    2*(mp.n_variants-1);

  gradient = (double*)malloc(sizeof(double)*n_gradient_variables);
  gradient_state = (double*)malloc(sizeof(double)*n_gradient_variables);
  possible_new_gradient_state =
    (double*)malloc(sizeof(double)*n_gradient_variables);
  evaluative_gradient_state =
    (double*)malloc(sizeof(double)*n_gradient_variables);
  
  tds = read_tracks(track_file, mp.n_variants, mp.n_classes);

  for(i=0;i<tds.n_tracks;i++) {
    determine_Inc_R_at_departure(tds.tracks+i);
  }
  
  initial_intermediate_fraction =
    (double**)malloc(sizeof(double*)*tds.n_tracks);
  initial_infected_fraction =
    (double**)malloc(sizeof(double*)*tds.n_tracks);
  
  for(i=0;i<tds.n_tracks;i++) {
    initial_intermediate_fraction[i]=
      (double*)malloc(sizeof(double)*mp.n_variants*mp.n_classes);
    initial_infected_fraction[i]=
      (double*)malloc(sizeof(double)*mp.n_variants*mp.n_classes);
  }
    
  /* needs to be reworked for multiple tracks */
  for(i=0;i<tds.n_tracks;i++) {
    initialization_sum = 0;
    for(j=0;j<mp.n_variants;j++) {
      initialization_sum += tds.tracks[i].points[j][0];
    }
    for(k=0;k<mp.n_classes;k++) {
      for(j=0;j<mp.n_variants;j++) {
	initial_intermediate_fraction[i][k*mp.n_variants+j] =
	  (tds.tracks[i].Inc_first/mp.conversion_rate[j])
	  *(tds.tracks[i].points[j][0]/initialization_sum);
	initial_infected_fraction[i][k*mp.n_variants+j] =
	  (tds.tracks[i].Inc_first*mp.recovery_rate[j])
	  *(tds.tracks[i].points[j][0]/initialization_sum);
      }
    }
  }

  /* warning special startup conditions - not really portable */
  mp.initial_intermediate_fraction[0] = initial_intermediate_fraction[0][0];
  mp.initial_infected_fraction[0] = initial_infected_fraction[0][0];
  
  trjds.n_trjs = tds.n_tracks;
  trjds.trjs = (trajectory*)malloc(sizeof(trajectory)*trjds.n_trjs);
  for(i=0;i<trjds.n_trjs;i++) {
    trjds.trjs[i].n_values=0;
  }
  
  buffer_n_variant_size = (double*)malloc(sizeof(double)*mp.n_variants);
  
  compartment_indices =
    (size_t*)malloc(sizeof(size_t)*mp.n_variants*mp.n_classes);

  mp.S_offset = 0;
  mp.E_offset = mp.n_variants*mp.n_classes;
  mp.bigI_offset = 2*mp.n_variants*mp.n_classes;
  mp.Inc_offset = 3*mp.n_variants*mp.n_classes;
  mp.V_offset = 4*mp.n_variants*mp.n_classes;

  for(i=0;i<mp.n_classes*mp.n_variants;i++) {
    compartment_indices[i] = mp.E_offset+i;
  }
  
  for(i=0;i<n_gradient_variables;i++) {
    initialize_gradient_state(i,
			      &mp,
			      gradient_state);
  }
  
  for(i=0;i<n_gradient_variables;i++) {
    gradient_state_to_mp(i,&mp,
			 initial_intermediate_fraction,
			 initial_infected_fraction,
			 gradient_state[i]+derivative_delta,
			 evaluative_gradient_state);
  }
  
  execute_epidemic_with_track_dataset_conditions(&trjds,
						 &mp,
						 &tds);

  
  while(error > epsilon) {

    for(j=0;j<n_gradient_variables;j++) {

      derivative_delta = machine_epsilon_sqrt*gradient_state[j];

      memcpy(evaluative_gradient_state,gradient_state,
	     sizeof(double)*n_gradient_variables);
      
      evaluative_gradient_state[j]+=derivative_delta;

      gradient_state_to_mp(j,&mp,
			   initial_intermediate_fraction,
			   initial_infected_fraction,
			   gradient_state[j]+derivative_delta,
			   evaluative_gradient_state);

      
      execute_epidemic_with_track_dataset_conditions(&trjds,
						     &mp,
						     &tds);
      
      forward_error = measure_error(tds,trjds,
				    compartment_indices,buffer_n_variant_size,
				    mp.n_classes,mp.conversion_rate);

      //forward_error += lagrange_error_addition(evaluative_gradient_state,
      //				       mp.n_variants);
      //      forward_error *= forward_error;

      evaluative_gradient_state[j]-= 2*derivative_delta;
      
      gradient_state_to_mp(j,&mp,
			   initial_intermediate_fraction,
			   initial_infected_fraction,
			   gradient_state[j]-derivative_delta,
			   evaluative_gradient_state);

      execute_epidemic_with_track_dataset_conditions(&trjds,
						     &mp,
						     &tds);
      
      backward_error = measure_error(tds,trjds,
				     compartment_indices,buffer_n_variant_size,
				     mp.n_classes,mp.conversion_rate);

      //backward_error += lagrange_error_addition(evaluative_gradient_state,
      //						mp.n_variants);
      //      backward_error *= backward_error;

      
      gradient[j] = (forward_error - backward_error)/(2*derivative_delta);
    }
    gradient_sum = 0;
    for(j=0;j<n_gradient_variables;j++) {
      gradient_sum += gradient[j]*gradient[j];
    }
    for(j=0;j<n_gradient_variables;j++) {
      gradient[j]/=sqrt(gradient_sum);
    }
    
    test_line_search(n_gradient_variables, gradient,
		     initial_infected_fraction,
		     initial_intermediate_fraction,
		     gradient_state,
		     possible_new_gradient_state,
		     &current_step_size,
		     &this_step_error,
		     previous_step_error, &mp, &trjds, &tds,
		     compartment_indices,
		     buffer_n_variant_size);
    printf("current error : %12.6E\n", this_step_error);
    errors[gradient_iterations%30] =
      fabs(previous_step_error - this_step_error);

    if(gradient_iterations > 30) {
      for(i=0;i<30;i++) {
	error += errors[i];
	error *= (1./30.);
      }
      printf("Delta-E %E \n", error);
    }
    previous_step_error = this_step_error;
    gradient_iterations++;
  }
  print_model_parameters(mp);
  for(i=0;i<trjds.n_trjs;i++) {
    print_trajectory(fitted_model,trjds.trjs[i],
		     compartment_indices,mp.n_variants*mp.n_classes);
  }
  fclose(fitted_model);
}
