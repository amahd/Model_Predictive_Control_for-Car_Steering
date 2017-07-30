#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "defines.h"

using namespace std;


class MPC {
 public:

 size_t N ;
 double dt ;
 const double Lf;
  MPC(): N(HORIZON), dt(STEP), Lf(LF){} ;

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
