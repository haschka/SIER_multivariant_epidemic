typedef struct {
  size_t n_variants;
  size_t n_classes;

  double runtime;
  double error_tolerance;
  
  double total_population;
  double reporting_rate;

  double* epsilon;
  
  size_t S_offset;
  size_t E_offset;
  size_t bigI_offset;
  size_t Inc_offset;
  size_t V_offset;
  
  double* initial_intermediate_fraction;
  double* initial_infected_fraction;
  
  double* vaccine_logistic_maximum;
  double* vaccine_logistic_growth;
  double* vaccine_logistic_midpoint;
  double* vaccine_efficiency;

  double* immune_drain;
  double* vaccine_drain;
  
  double* population_class_distribution;
  double* contacts_between_classes;
  double* pre_immune_fraction;
  double* sigma;

  double* Rzero;

  double* variant_introduction_date;
  
  double* asymptomatics_fraction;

  double* conversion_rate;
  double* recovery_rate;

  double* pop_per_class;
  double* rho;
  double* lambda;
  double* beta;
  
} model_parameters;

void print_model_parameters(model_parameters p);

typedef struct {
  int variant;
  int class;
  model_parameters* mp;
} state_function_parameters;

void parse_model_properties(FILE* in, model_parameters* mp);
trajectory epidemic_execute(model_parameters mp);
