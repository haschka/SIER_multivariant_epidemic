#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"runge-kutta-fehlberg.h"
#include"epidemic.h"

void dgeev_(char*, char*, int*, double*, int*, double*, double*, double*,
	    int*, double*, int*, double*, int*, int*);


void print_model_parameters(model_parameters p) {

  size_t i,j;
  
  printf("n_variants        %12li\n",p.n_variants);
  printf("n_classes         %12li\n",p.n_classes);
  printf("runtime           %12.6E\n",p.runtime);
  printf("error_tolerance   %12.6E\n",p.error_tolerance);

  printf("Rzero             ");
  for(i=0;i<p.n_variants;i++) {
    printf("%12.6E ",p.Rzero[i]);
  }
  printf("\n");
  
  printf("epsilon:\n");
  for(i=0;i<p.n_variants;i++) {
    for(j=0;j<p.n_variants-1;j++) {
      printf("%12.6E ",p.epsilon[i*p.n_variants+j]);
    }
    printf("%12.6E\n",p.epsilon[i*p.n_variants+p.n_variants-1]);
  }
  printf("\n");

  printf("Initial E Compartment:\n");
  for(i=0;i<p.n_classes;i++) {
    for(j=0;j<p.n_variants-1;j++) {
      printf("%12.6E ",p.initial_intermediate_fraction[i*p.n_variants+j]);
    }
    printf("%12.6E\n",
	   p.initial_intermediate_fraction[i*p.n_variants+p.n_variants-1]);
  }

  printf("Initial I Compartment:\n");
  for(i=0;i<p.n_classes;i++) {
    for(j=0;j<p.n_variants-1;j++) {
      printf("%12.6E ",p.initial_infected_fraction[i*p.n_variants+j]);
    }
    printf("%12.6E\n",
	   p.initial_infected_fraction[i*p.n_variants+p.n_variants-1]);
  }

  printf("Pre Immune Fraction:\n");
  for(i=0;i<p.n_classes;i++) {
    printf("%12.6E ",p.pre_immune_fraction[i]);
  }
  printf("\n");
}

static inline double vaccine_rate(int variant, int class, double t,
				  model_parameters p) {
  
  int i = variant;
  int j = class;
  size_t location = j*p.n_variants+i;

  double rate = p.vaccine_logistic_maximum[location]*1.
	  /(1+exp(-p.vaccine_logistic_growth[location]
		  *(t-p.vaccine_logistic_midpoint[location])));

  rate /= p.total_population;
  
  return (rate);
  
}

static inline double* build_population_per_class(model_parameters p) {

  int i;

  double* pop = (double*)malloc(sizeof(double)*p.n_classes);

  for(i=0;i<p.n_classes;i++) {
    pop[i] = p.total_population*p.population_class_distribution[i];
  }

  return(pop);
}

static inline double* build_rho(model_parameters p) {

  int i,j,k;

  char JOBVL = 'N';
  char JOBVR = 'N';
  int N = (int)p.n_classes;
  double*A = (double*)malloc(sizeof(double)*p.n_classes*p.n_classes);
  int LDA = (int)p.n_classes;
  double*WR = (double*)malloc(sizeof(double)*p.n_classes);
  double*WI = (double*)malloc(sizeof(double)*p.n_classes);
  double*VL = NULL;
  int LDVL = 1;
  double* VR  = NULL;
  int LDVR = 1;
  double*WORK = (double*)malloc(sizeof(double));
  int LWORK = -1;
  int INFO;

  double MAX;
  
  double optimal;

  double *rho = (double*)malloc(sizeof(double)*p.n_variants);

  dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR ,
	 WORK, &LWORK, &INFO);
  
  optimal = WORK[0];

  LWORK = (int)optimal;
  
  WORK = (double*)realloc(WORK,sizeof(double)*optimal);
  
  for(k=0;k<p.n_variants;k++) {
  
    for(i=0;i<p.n_classes;i++) {
      for(j=0;j<p.n_classes;j++) {
	A[i*p.n_classes+j] =
	  (1-p.pre_immune_fraction[i*p.n_variants+k])*p.sigma[i*p.n_variants+k]
	  *p.contacts_between_classes[i*p.n_classes+j]
	  *p.pop_per_class[i]/p.pop_per_class[j];
      }
    }

    dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR ,
	   WORK, &LWORK, &INFO);

    MAX = WR[0];
    for(i=0;i<p.n_classes;i++) {
      if ( MAX < WR[i] ) MAX = WR[i];
    }
    
    rho[k] = MAX;
  }
  
  free(A);
  free(WR);
  free(WI);
  free(WORK);
  return(rho);
}
    
static inline double* build_lambda(model_parameters p) {

  int i;

  double* lambda = (double*)malloc(sizeof(double)*p.n_variants);
  
  for(i = 0;i < p.n_variants; i++) {
    lambda[i] = p.Rzero[i]*p.recovery_rate[i]/p.rho[i];
  }
  return (lambda);
}

static inline double* build_beta(model_parameters p) {

  int j,l,k;

  double * beta =
    (double*)malloc(sizeof(double)*p.n_classes*p.n_classes*p.n_variants);

  double total_over_class;

  for(j=0;j<p.n_classes;j++) {
    for(l=0;l<p.n_classes;l++) {
      for(k=0;k<p.n_variants;k++) {
	
	total_over_class = p.total_population/p.pop_per_class[l];
	
	beta[k*p.n_classes*p.n_classes+l*p.n_classes+j] =
	  p.lambda[k]*p.sigma[j*p.n_variants+k]*total_over_class
	  *p.contacts_between_classes[j*p.n_classes+l];
      }
    }
  }
  
  return(beta);
}

static inline double dS_dt (double t, double* states, void*parameters) {

  state_function_parameters* p = (state_function_parameters*)parameters;
  
  int i=p->variant;
  int j=p->class;
  int k,l;

  model_parameters* mp = p->mp;
  
  int n_variants = mp->n_variants;
  int n_classes = mp->n_classes;

  double*S = states+(mp->S_offset);
  double*E = states+(mp->E_offset);
  double*bigI = states+(mp->bigI_offset);
  double*V = states+(mp->V_offset);

  double*beta = mp->beta;
  double*epsilon = mp->epsilon;

  double* vaccine_efficiency = mp->vaccine_efficiency;
  
  double S_i_j = S[i*n_classes+j];
  double E_i_j = E[i*n_classes+j];
  double bigI_i_j = bigI[i*n_classes+j];
  double V_i_j = V[i*n_classes+j];
  
  double asymptomatic_i_j = mp->asymptomatics_fraction[i*n_classes+j];
  double reporting_rate = mp->reporting_rate;
  
  double R_i_j =
    mp->population_class_distribution[j]- S_i_j - E_i_j - bigI_i_j - V_i_j;
  
  double D_i_j =
    S_i_j + E_i_j + asymptomatic_i_j*bigI_i_j
    +(1 - reporting_rate)*R_i_j+V_i_j;

  double vaccination;
    
  double sum = 0;
  
  for(k=0;k<n_variants;k++) {
    for(l=0;l<n_classes;l++) {
      sum += epsilon[i*n_variants+k]
	*beta[k*n_classes*n_classes+l*n_classes+j]
	*bigI[k*n_classes+l];
    }
  }

  vaccination = (vaccine_rate(i,j,t+1,mp[0])
		 -vaccine_rate(i,j,t,mp[0]))*(S_i_j/D_i_j)
    *vaccine_efficiency[j*n_variants+i]*S_i_j;
  
  return(mp->immune_drain[j*n_variants+i]*R_i_j
	 +mp->vaccine_drain[j*n_variants+i]*V_i_j-sum*S_i_j-vaccination);
}

static inline double dE_dt (double t, double* states, void*parameters) {

  state_function_parameters* p = (state_function_parameters*)parameters;
  model_parameters* mp = p->mp;
  
  int i=p->variant, j=p->class,l;

  double * beta = mp->beta;
  double * S = states+(mp->S_offset);
  double * E = states+(mp->E_offset);
  double * bigI = states+(mp->bigI_offset);

  int n_classes = mp->n_classes;
  int n_variants = mp->n_variants;
  
  double sum = 0;
  for(l=0;l<n_classes;l++) {
    sum += beta[i*n_classes*n_classes+l*n_classes+j]
      *bigI[i*n_classes+l];
  }
  return(sum*S[i*n_classes+j]-mp->conversion_rate[i]*E[i*n_classes+j]);
}

static inline double dI_dt (double t, double* states, void*parameters) {

  state_function_parameters* p = (state_function_parameters*)parameters;
  model_parameters* mp = p->mp;
  
  int i=p->variant, j=p->class;

  int n_classes = mp->n_classes;
  int n_variants = mp->n_variants;
  
  double* E = states+(mp->E_offset);
  double* bigI = states+(mp->bigI_offset);

  double* conversion_rate = mp->conversion_rate;
  double* recovery_rate = mp->recovery_rate;
  
  return(conversion_rate[i]*E[i*n_classes+j]
	 -recovery_rate[i]*bigI[i*n_classes+j]);
}

static inline double dInc_dt (double t, double* states, void* parameters) {

  state_function_parameters* p = (state_function_parameters*)parameters;
  model_parameters* mp = p->mp;
  
  int i=p->variant, j=p->class;

  int n_classes = mp->n_classes;
  int n_variants = mp->n_variants;
  
  double* E = states+(mp->E_offset);
  
  return((mp->conversion_rate[i])*E[i*n_classes+j]);
}

static inline double dV_dt(double t, double* states, void* parameters) {
  state_function_parameters* p = (state_function_parameters*)parameters;
  model_parameters* mp = p->mp;
  
  int i=p->variant, j=p->class;

  int n_classes = mp->n_classes;
  int n_variants = mp->n_variants;

  double* vaccine_efficiency = mp->vaccine_efficiency;

  double * S = states+(mp->S_offset);
  double * E = states+(mp->E_offset);
  double * bigI = states+(mp->bigI_offset);
  double* V = states+(mp->V_offset);
  
  double S_i_j = S[i*n_classes+j];
  double E_i_j = E[i*n_classes+j];
  double bigI_i_j = bigI[i*n_classes+j];
  double V_i_j = V[i*n_classes+j];
  
  double asymptomatic_i_j = mp->asymptomatics_fraction[i*n_classes+j];
  double reporting_rate = mp->reporting_rate;
  
  double R_i_j =
    mp->population_class_distribution[j]- S_i_j - E_i_j - bigI_i_j - V_i_j;
  
  double D_i_j =
    S_i_j + E_i_j + asymptomatic_i_j*bigI_i_j
    +(1 - reporting_rate)*R_i_j+V_i_j;

  double vaccination;

  vaccination =
    (vaccine_rate(i,j,t+1,mp[0])
     -vaccine_rate(i,j,t,mp[0]))*(S_i_j/D_i_j)
    *vaccine_efficiency[j*n_variants+i]*S_i_j;

  return(vaccination - mp->vaccine_drain[j*n_variants+i]*V_i_j);
}

trajectory epidemic_execute(model_parameters mp) {

  size_t i,j;
  
  int n_variants, n_classes, n_variants_times_n_classes;
  int n_eqns;

  double addition_to_E, addition_to_I;
  
  state_eqn* state_eqns;
  state_function_parameters* ps;

  size_t idx_S,idx_E,idx_I,idx_Inc,idx_V;

  trajectory* partial_trajectory;
  
  trajectory trj;

  double* initial_state;

  n_variants = mp.n_variants;
  n_classes = mp.n_classes;
  
  n_eqns = n_variants*n_classes*5;

  partial_trajectory = (trajectory*)malloc(sizeof(trajectory)*n_variants);
  
  mp.pop_per_class = build_population_per_class(mp);
  mp.rho = build_rho(mp);
  mp.lambda = build_lambda(mp);
  mp.beta = build_beta(mp);

  state_eqns = (state_eqn*)malloc(sizeof(state_eqn)*n_eqns);

  n_variants_times_n_classes = n_variants*n_classes;

  mp.S_offset = 0;
  mp.E_offset = n_variants_times_n_classes;
  mp.bigI_offset = 2*n_variants_times_n_classes;
  mp.Inc_offset = 3*n_variants_times_n_classes;
  mp.V_offset = 4*n_variants_times_n_classes;
  
  ps =
    (state_function_parameters*)malloc(sizeof(state_function_parameters)
 				       *n_eqns);
  
  initial_state = (double*)malloc(sizeof(double)*n_eqns);
  
  for(i=0;i<n_variants;i++) {
    for(j=0;j<n_classes;j++) {
      idx_S = i*n_classes+j;
      idx_E = n_variants_times_n_classes+idx_S;
      idx_I = 2*n_variants_times_n_classes+idx_S;
      idx_Inc = 3*n_variants_times_n_classes+idx_S;
      idx_V = 4*n_variants_times_n_classes+idx_S;
      
      initial_state[idx_S] =
	mp.total_population*mp.population_class_distribution[j];

      state_eqns[idx_S].eqn = &dS_dt;
      ps[idx_S].variant = i;
      ps[idx_S].class = j;
      ps[idx_S].mp = &mp;
      state_eqns[idx_S].params =
	(void*)((state_function_parameters*)(ps+idx_S));

      initial_state[idx_E] = 0.;

      state_eqns[idx_E].eqn = &dE_dt;
      ps[idx_E].variant = i;
      ps[idx_E].class = j;
      ps[idx_E].mp = &mp;
      state_eqns[idx_E].params =
	(void*)((state_function_parameters*)(ps+idx_E));

      initial_state[idx_I] = 0.;
      
      state_eqns[idx_I].eqn = &dI_dt;
      ps[idx_I].variant = i;
      ps[idx_I].class = j;
      ps[idx_I].mp = &mp;
      state_eqns[idx_I].params =
	(void*)((state_function_parameters*)(ps+idx_I));

      initial_state[idx_Inc] = 0.;
      
      state_eqns[idx_Inc].eqn = &dInc_dt;
      ps[idx_Inc].variant = i;
      ps[idx_Inc].class = j;
      ps[idx_Inc].mp = &mp;
      state_eqns[idx_Inc].params =
	(void*)((state_function_parameters*)(ps+idx_Inc));

      initial_state[idx_V] = 0.;

      state_eqns[idx_V].eqn = &dV_dt;
      ps[idx_V].variant = i;
      ps[idx_V].class = j;
      ps[idx_V].mp = &mp;
      state_eqns[idx_V].params =
	(void*)((state_function_parameters*)(ps+idx_V));

      
      initial_state[idx_S] -=
	initial_state[idx_S]*mp.pre_immune_fraction[j*mp.n_variants+i];
      
      initial_state[idx_S] /= mp.total_population;
      initial_state[idx_E] /= mp.total_population;
      initial_state[idx_I] /= mp.total_population;
      
    }
  }
  
  for(i=0;i<n_variants;i++) {
    for(j=0;j<n_classes;j++) {

      idx_S = i*n_classes+j;
      idx_E = n_variants_times_n_classes+idx_S;
      idx_I = 2*n_variants_times_n_classes+idx_S;
      idx_Inc = 3*n_variants_times_n_classes+idx_S;
      
      addition_to_E = mp.population_class_distribution[j]
	*mp.initial_intermediate_fraction[j*n_variants+i];
      
      addition_to_I = mp.population_class_distribution[j]
	*mp.initial_infected_fraction[j*n_variants+i];
      
      initial_state[idx_E] += addition_to_E;
      initial_state[idx_I] += addition_to_I;
      initial_state[idx_Inc] += addition_to_I;
       
      initial_state[idx_S] -= addition_to_E;
      initial_state[idx_S] -= addition_to_I;
      
    }
    
    if (i < n_variants-1) {
      partial_trajectory[i] = rkf(initial_state,
				  state_eqns,
				  n_eqns,
				  mp.error_tolerance,
				  1.,
				  mp.variant_introduction_date[i],
				  mp.variant_introduction_date[i+1]);
      trajectory_to_state_at_end(partial_trajectory[i],initial_state);
      partial_trajectory[i].n_values -=1;
    } else {
      partial_trajectory[i] = rkf(initial_state,
				  state_eqns,
				  n_eqns,
				  mp.error_tolerance,
				  1.,
				  mp.variant_introduction_date[i],
				  mp.runtime);
    }
  }
  
  free(state_eqns);
  free(ps);
  free(mp.beta);
  free(mp.lambda);
  free(mp.pop_per_class);
  free(mp.rho);
  free(initial_state);
  
  for(i=1;i<n_variants;i++) {
    partial_trajectory[0] =
      combine_trajectories_overwrite(partial_trajectory[0],
				     partial_trajectory[i]);
  }
  trj = partial_trajectory[0];
  free(partial_trajectory);
  return(trj);
}
