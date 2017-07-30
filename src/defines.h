/*
 * defines.h
 *
 *  Created on: Jul 30, 2017
 *      Author: s
 */

#ifndef DEFINES_H_
#define DEFINES_H_

// the timestep length and duration
# define HORIZON 10
# define STEP  0.2

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
#define  LF  2.67

//Converting mph to m/s, required when using simulator latency
#define MILE2MPS 0.447


// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 30 mph, but it can be set upto 40mph
#define  REFV  30*MILE2MPS   // converted already in m/sec

// Latency in the actuation system
#define LATENCY 0.1


#define  W_CTE  200   	// weight for CTE
#define  W_ORE  200		//weight for orientation error
#define  W_VEL	20 		//weight for difference between current and desired velocity


#define  W_DELUSE	8000	//weight for minimising steering
#define  W_THRUSE	1		//weight for throttle use minimization


#define  W_DELDIFF	16000	//weight for stopping sudden steering changes
#define  W_THRDIFF	1000  	//weight for stopping sudden throttle changes





#endif /* DEFINES_H_ */
