# IEKF
This repository contains the code for state estimation for a unicycle system. The states are estimated using Unscented Kalman Filter (UKF), Extended Kalman Filter (EKF) and Invariant Extended Kalman Filter (IEKF). The code was developed in an effort to regenerate the results of [1].

# System equations
The system considered in this problem is a unicycle system with the following equations
$$\dot{x} = v\cos(\theta) + C_xd$$
$$\dot{y} = v\sin(\theta) + C_yd$$
$$\dot{\theta} = \omega,$$
where $(x,y)$ are the position of the robot, $\theta$ is the heading, $v$ is the linear velocity, and $\omega$ is the turning rate. $d$ is the disturbance acting on the system through the matrices $C_x$ and $C_y$. The disturbance dynamics is given by
$$\dot{d} = Ad.$$
For more information see section 3 in [1].

## Undisturbed system
In the first example, there is no disturbance acting on the system, i.e., $d=0$. Try running `runMe.m` in the `Undisturbed` folder to get the following results.

<img src="/images/1.png" width="70%" height="70%">
