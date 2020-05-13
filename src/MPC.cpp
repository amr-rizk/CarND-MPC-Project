#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
//#include <math.h> 

using CppAD::AD;
using Eigen::VectorXd;
using namespace std;

/**
 * TODO: Set the timestep length and duration
 */
size_t N = 15;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
//   simulator around in a circle with a constant steering angle and velocity on
//   a flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
//   presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

//define input limits
const double a_max = 1;
const double a_min = -1;
const double del_max = 0.43633;
const double del_min = -0.43633;

//define refrences
//const double cte_ref = 0.0;
//const double epsi_ref = 0.0;
const double v_ref = 50.0;

//define start indices for the vars vector
const int x_start = 0;
const int y_start = x_start + N;
const int psi_start = y_start + N;
const int v_start = psi_start + N;
const int cte_start = v_start + N;
const int epsi_start = cte_start + N;
const int delta_start = epsi_start + N;
const int a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  VectorXd coeffs;
  FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    /**
     * TODO: implement MPC
     * `fg` is a vector of the cost constraints, `vars` is a vector of variable 
     *   values (state & actuators)
     * NOTE: You'll probably go back and forth between this function and
     *   the Solver function below.
     */
     fg[0] = 0;
     
     for (unsigned int i = 0; i < N; i++){ 
     	fg[0] += 1000*CppAD::pow(vars[cte_start+i], 2) + 1000*CppAD::pow(vars[epsi_start + i], 2) + CppAD::pow(vars[v_start + i] - v_ref, 2);
     	
     	if (i < N - 1){
     		fg[0] += 50*CppAD::pow(vars[a_start + i], 2) + 50*CppAD::pow(vars[delta_start + i], 2);
     	}
     	if (i < N - 2){
     		fg[0] += 250000*CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2) + 5000*CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
     	}
     }
     
    for (unsigned int n = 0; n < N; n++){
      if (n == 0){
      	fg[1 + x_start] = vars[x_start];
    	fg[1 + y_start] = vars[y_start];
    	fg[1 + psi_start] = vars[psi_start];
    	fg[1 + v_start] = vars[v_start];
    	fg[1 + cte_start] = vars[cte_start];
    	fg[1 + epsi_start] = vars[epsi_start];
    	continue;
      }
      
      // The state at time t+1
      AD<double> x_n1 = vars[x_start + n];
      AD<double> y_n1 = vars[y_start + n];
      AD<double> psi_n1 = vars[psi_start + n];
      AD<double> v_n1 = vars[v_start + n];
      AD<double> cte_n1 = vars[cte_start + n];
      AD<double> epsi_n1 = vars[epsi_start + n];

      // The state at time t
      AD<double> x_n = vars[x_start + n - 1];
      AD<double> y_n = vars[y_start + n - 1];
      AD<double> psi_n = vars[psi_start + n - 1];
      AD<double> v_n = vars[v_start + n - 1];
      AD<double> cte_n = vars[cte_start + n - 1];
      AD<double> epsi_n = vars[epsi_start + n - 1];

      // Only consider the actuation at time t.
      AD<double> delta_n = vars[delta_start + n - 1];
      AD<double> a_n = vars[a_start + n - 1];

      AD<double> fx = coeffs[0] + coeffs[1]*x_n + coeffs[2]*CppAD::pow(x_n, 2) + coeffs[3]*CppAD::pow(x_n, 3);
      AD<double> psi_des_n = CppAD::atan(coeffs[1] + 2*coeffs[2]*x_n + 3*coeffs[3]*CppAD::pow(x_n, 2));

      fg[1 + x_start + n] = x_n1 - (x_n + v_n*CppAD::cos(psi_n)*dt);
      fg[1 + y_start + n] = y_n1 - (y_n + v_n*CppAD::sin(psi_n)*dt);
      fg[1 + psi_start + n] = psi_n1 - (psi_n - (v_n/Lf)*delta_n*dt);
      fg[1 + v_start + n] = v_n1 - (v_n + a_n*dt);
      fg[1 + cte_start + n] = cte_n1 - (fx - y_n + v_n*CppAD::sin(epsi_n)*dt);
      fg[1 + epsi_start + n] = epsi_n1 - (psi_n - psi_des_n - (v_n/Lf)*delta_n*dt);
    } 
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

std::vector<double> MPC::Solve(const VectorXd &state, const VectorXd &coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  /**
   * TODO: Set the number of model variables (includes both states and inputs).
   * For example: If the state is a 4 element vector, the actuators is a 2
   *   element vector and there are 10 timesteps. The number of variables is:
   *   4 * 10 + 2 * 9
   */
  size_t n_vars = N*6 + (N - 1)*2;
  /**
   * TODO: Set the number of constraints
   */
  size_t n_constraints = N*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (unsigned int i = 0; i < n_vars; ++i) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  /**
   * TODO: Set lower and upper limits for variables.
   */
  
  //no limits on the states
  for (int i = 0; i < delta_start; i++){
    vars_lowerbound[i] = -std::numeric_limits<int>::max();
    vars_upperbound[i] = std::numeric_limits<int>::max();
  }
  
  
  //limits on actuators
  for (int i = delta_start; i < a_start; i++ ) {
    vars_lowerbound[i] = del_min;
    vars_upperbound[i] = del_max;
  }
  for (unsigned int i = a_start; i < n_vars; i++ ) {
    vars_lowerbound[i] = a_min;
    vars_upperbound[i] = a_max;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (unsigned int i = 0; i < n_constraints; ++i) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  //equality constraints
  constraints_lowerbound[x_start] = state[0];
  constraints_lowerbound[y_start] = state[1];
  constraints_lowerbound[psi_start] = state[2];
  constraints_lowerbound[v_start] = state[3];
  constraints_lowerbound[cte_start] = state[4];
  constraints_lowerbound[epsi_start] = state[5];

  constraints_upperbound[x_start] = state[0];
  constraints_upperbound[y_start] = state[1];
  constraints_upperbound[psi_start] = state[2];
  constraints_upperbound[v_start] = state[3];
  constraints_upperbound[cte_start] = state[4];
  constraints_upperbound[epsi_start] = state[5];

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // NOTE: You don't have to worry about these options
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  //   of sparse routines, this makes the computation MUCH FASTER. If you can
  //   uncomment 1 of these and see if it makes a difference or not but if you
  //   uncomment both the computation time should go up in orders of magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  /**
   * TODO: Return the first actuator values. The variables can be accessed with
   *   `solution.x[i]`.
   *
   * {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
   *   creates a 2 element double vector.
   */
   
  vector<double> result;

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (unsigned int i = 0; i < N - 2; i++ ) {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }
  return result;
}
