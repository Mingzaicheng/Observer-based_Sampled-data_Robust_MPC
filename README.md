# Observer-based Sampled-data RMPC

This simulation validates the proposed observer-based sampled-data robust MPC approach for LPV systems by the continuous-time continuous stirred tank reactor (CSTR) model and Numerical model.

## Requirement

1) MATLAB R2024a
2) LMI toolbox
3) ellipsoidal toolbox

## Short introduction 

This simulation includes two folders, corresponding to Simulation Examples 1 and 2. Each folder contains simulation files for Case (a) (considering saturated control inputs) and Case (b) (neglecting input saturation) as described in Algorithm 1.

By running the `main.m` file, response plots of the system states, estimated states, and estimated error sets can be obtained, along with the system state trajectory and control input response plots. The `main.m` file first designs the CSTR/Numerical model (with experimental parameters detailed in the `system_model.m` file), then performs offline optimization of the state observer gain by solving the LMI problem in the `observer_off.m` file. The state observer system used for state estimation is defined in the `observer_model.m` file. Subsequently, the feedback controller gain is optimized online by solving the `optimization_on.m` file, while the estimated error set is updated in real time by calling the `error_refresh.m` file.

