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
In the first example, there is no disturbance acting on the system, i.e., $d=0$. Try running `runMe.m` in the `undisturbed` folder to get the following results.

<img src="/images/x.png" width="50%" height="50%">
<img src="/images/y.png" width="50%" height="50%">
<img src="/images/theta.png" width="50%" height="50%">

As it could be seen from these plots, the IEKF outperforms the EKF and UKF in terms of RMSE of estimation. 

## Disturbed system
In the second example, a disturbance $d \in \mathbb{R}^4$ acts on the system. There are two different IEKF designs provided in [1]. Run `runMe.m` in the `disturbed` folder to get the results. Some results are shown here.

<img src="/images/theta.png" width="50%" height="50%">
<img src="/images/d2.png" width="50%" height="50%">
<img src="/images/d4.png" width="50%" height="50%">

# To do
The result from the second example does not exactly match the result from the paper.

# References
1- Coleman, K., Bai, H. and Taylor, C.N., 2021. Extended invariant-EKF designs for state and additive disturbance estimation. Automatica, 125, p.109464.
