# IEKF
This repository contains the code for state estimation for a unicycle system. The states are estimated using Unscented Kalman Filter (UKF), Extended Kalman Filter (EKF) and Invariant Extended Kalman Filter (IEKF). The code was developed in an effort to regenerate the results of [1].

# System equations
The system considered in this problem is a unicycle system with the following equations
$$\begin{aligned}
\dot{x} & = v\cos(\theta) + C_xd \\
\dot{y} & = v\sin(\theta) + C_yd\\
\dot{\theta} & = \omega,$$

<img src="/images/1.png" width="70%" height="70%">
