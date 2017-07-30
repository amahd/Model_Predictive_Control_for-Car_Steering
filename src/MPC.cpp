#include "MPC.h"

#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

using CppAD::AD;



// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + HORIZON;
size_t psi_start = y_start + HORIZON;
size_t v_start = psi_start + HORIZON;
size_t cte_start = v_start + HORIZON;
size_t epsi_start = cte_start + HORIZON;
size_t delta_start = epsi_start + HORIZON;
size_t a_start = delta_start + HORIZON - 1;

class FG_eval {
public:
// Fitted polynomial coefficients
Eigen::VectorXd coeffs;
FG_eval(Eigen::VectorXd coeffs) {
this->coeffs = coeffs;
}

typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
void operator()(ADvector& fg, const ADvector& vars) {

	// `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
	// NOTE: You'll probably go back and forth between this function and
	// the Solver function below.

	// The cost is stored is the first element of `fg`.
	// Any additions to the cost should be added to `fg[0]`.
	fg[0] = 0;

	// The part of the cost based on the reference state.
	for (unsigned int t = 0; t < HORIZON; t++) {
		fg[0] += W_CTE * CppAD::pow(vars[cte_start + t], 2);
		fg[0] += W_ORE * CppAD::pow(vars[epsi_start + t], 2);
		fg[0] += W_VEL * CppAD::pow(vars[v_start + t] - REFV, 2);
	}

	// Minimize the use of actuators.
	for (unsigned int t = 0; t < HORIZON - 1; t++) {
		fg[0] += W_DELUSE * CppAD::pow(vars[delta_start + t], 2);
		fg[0] += W_THRUSE * CppAD::pow(vars[a_start + t], 2);
	}

	// Minimize the value gap between sequential actuations.
	for (unsigned int t = 0; t < HORIZON - 2; t++) {
		fg[0] += W_DELDIFF * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
		fg[0] += W_THRDIFF * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
	}

	//
	// Setup Constraints
	//
	// NOTE: In this section you'll setup the model constraints.

	// Initial constraints
	//
	// We add 1 to each of the starting indices due to cost being located at
	// index 0 of `fg`.
	// This bumps up the position of all the other values.
	fg[1 + x_start] = vars[x_start];
	fg[1 + y_start] = vars[y_start];
	fg[1 + psi_start] = vars[psi_start];
	fg[1 + v_start] = vars[v_start];
	fg[1 + cte_start] = vars[cte_start];
	fg[1 + epsi_start] = vars[epsi_start];

	// The rest of the constraints
	for (int t = 0; t < HORIZON - 1; t++) {
		// The state at time t+1 .
		AD<double> x1 = vars[x_start + t + 1];
		AD<double> y1 = vars[y_start + t + 1];
		AD<double> psi1 = vars[psi_start + t + 1];
		AD<double> v1 = vars[v_start + t + 1];
		AD<double> cte1 = vars[cte_start + t + 1];
		AD<double> epsi1 = vars[epsi_start + t + 1];

		// The state at time t.
		AD<double> x0 = vars[x_start + t];
		AD<double> y0 = vars[y_start + t];
		AD<double> psi0 = vars[psi_start + t];
		AD<double> v0 = vars[v_start + t];
		AD<double> cte0 = vars[cte_start + t];
		AD<double> epsi0 = vars[epsi_start + t];

		// Only consider the actuation at time t.
		AD<double> delta0 = vars[delta_start + t];
		AD<double> a0 = vars[a_start + t];

		// AD<double> f0 = coeffs[0] + coeffs[1] * x0;
		// AD<double> psides0 = CppAD::atan(coeffs[1]);

		AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
		AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0);

		// Here's `x` to get you started.
		// The idea here is to constraint this value to be 0.
		//
		// Recall the equations for the model:
		// x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
		// y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
		// psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
		// v_[t+1] = v[t] + a[t] * dt
		// cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
		// epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
		fg[2 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * STEP);
		fg[2 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * STEP);
		fg[2 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / LF * STEP);
		fg[2 + v_start + t] = v1 - (v0 + a0 * STEP);
		fg[2 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * STEP));
		fg[2 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / LF * STEP);
		}
	  }
};

//
// MPC class definition implementation.
//
//MPC::MPC() {
//}
MPC::~MPC() {
}

vector<double>  MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
	bool ok = true;
	typedef CPPAD_TESTVECTOR(double) Dvector;

	double x = state[0];
	double y = state[1];
	double psi = state[2];
	double v = state[3];
	double cte = state[4];
	double epsi = state[5];


	// For example: If the state is a 4 element vector, the actuators is a 2
	// element vector and there are 10 timesteps. The number of variables is:
	//
	// 4 * 10 + 2 * 9
	// number of independent variables
	// N timesteps == N - 1 actuations
	size_t n_vars = HORIZON * 6 + (HORIZON - 1) * 2;

	// Number of constraints
	size_t n_constraints = HORIZON * 6;

	// Initial value of the independent variables.
	// SHOULD BE 0 besides initial state.
	Dvector vars(n_vars);
	for (unsigned int i = 0; i < n_vars; i++)
		vars[i] = 0;

	Dvector vars_lowerbound(n_vars);
	Dvector vars_upperbound(n_vars);


	// Set all non-actuators upper and lowerlimits
	// to the max negative and positive values.
	for (unsigned int i = 0; i < delta_start; i++) {
		vars_lowerbound[i] = -1.0e19;
		vars_upperbound[i] = 1.0e19;
	}

	// The upper and lower limits of delta are set to -25 and 25
	// degrees (values in radians).
	// NOTE: Feel free to change this to something else.
	for (unsigned int i = delta_start; i < a_start; i++) {
		vars_lowerbound[i] = -0.436332;
		vars_upperbound[i] = 0.436332;
	}

	// Acceleration/decceleration upper and lower limits.
	// NOTE: Feel free to change this to something else.
	for (unsigned int i = a_start; i < n_vars; i++) {
		vars_lowerbound[i] = -1.0;
		vars_upperbound[i] = 1.0;
	}

	// Lower and upper limits for the constraints
	// Should be 0 besides initial state.
	Dvector constraints_lowerbound(n_constraints);
	Dvector constraints_upperbound(n_constraints);
	for (unsigned int i = 0; i < n_constraints; i++) {
		constraints_lowerbound[i] = 0;
		constraints_upperbound[i] = 0;
	}

	constraints_lowerbound[x_start] = x;
	constraints_lowerbound[y_start] = y;
	constraints_lowerbound[psi_start] = psi;
	constraints_lowerbound[v_start] = v;
	constraints_lowerbound[cte_start] = cte;
	constraints_lowerbound[epsi_start] = epsi;

	constraints_upperbound[x_start] = x;
	constraints_upperbound[y_start] = y;
	constraints_upperbound[psi_start] = psi;
	constraints_upperbound[v_start] = v;
	constraints_upperbound[cte_start] = cte;
	constraints_upperbound[epsi_start] = epsi;

	// object that computes objective and constraints
	FG_eval fg_eval(coeffs);

	//
	// NOTE: You don't have to worry about these options
	//
	// options for IPOPT solver
	std::string options;
	// Uncomment this if you'd like more print information
	options += "Integer print_level  0\n";
	// NOTE: Setting sparse to true allows the solver to take advantage
	// of sparse routines, this makes the computation MUCH FASTER. If you
	// can uncomment 1 of these and see if it makes a difference or not but
	// if you uncomment both the computation time should go up in orders of
	// magnitude.
	options += "Sparse  true        forward\n";
	options += "Sparse  true        reverse\n";
	// NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
	// Change this as you see fit.
	options += "Numeric max_cpu_time          0.5\n";

	// place to return solution
	CppAD::ipopt::solve_result < Dvector > solution;

	// solve the problem
	CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
			constraints_upperbound, fg_eval, solution);

	// Check some of the solution values
	ok &= solution.status == CppAD::ipopt::solve_result < Dvector > ::success;

	// Cost
	auto cost = solution.obj_value;
	std::cout << "Cost " << cost << std::endl;



	vector<double> sol(2*HORIZON + 2);   // Return the optimized actuator inputs plus the location of x and y coordinates
	  	  	  	  	  	  	  	   // over the horizon
	sol[0] = solution.x[delta_start];
	sol[1] = solution.x[a_start];
    for (unsigned int k = 0; k < HORIZON; ++k){
		  sol[k + 2] =   solution.x[x_start + k];
		  sol[k + 2 + HORIZON] =   solution.x[y_start + k];
		}

   return sol;
}

