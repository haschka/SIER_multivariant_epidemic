typedef struct {
  double ** states_at_t;
  double * t;
  size_t n_states;
  size_t n_values;
} trajectory;

typedef struct {
  double (*eqn) (double t, double *states, void * parameters);
  void * params;
} state_eqn;

double get_interpolated_value_at_t(double t, trajectory a,
				   size_t trajectory_component_index);

void trajectory_to_state_at_end(trajectory a, double* state);
void trajectory_to_state_at_t(double t, trajectory a, double* state);

trajectory combine_trajectories_overwrite(trajectory a, trajectory b);

trajectory rkf (double* states,
		state_eqn* state_equations,
		size_t n_states,
		double accepted_error,
		double initial_delta,
		double initial_t,
		double final_t);

void print_trajectory(FILE* out,trajectory trj,
		      size_t* indices, size_t n_indices);
