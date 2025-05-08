#!/bin/bash

export OMP_NUM_THREADS=1

# # -------------------------------- # 
# # Automatische Differentiezug (AD) #
# # -------------------------------- # 
# echo -e "\nExecuting: ad_odeint"
# ./../build/bin/ad_odeint_example

# echo -e "\nExecuting: ad_square"
# ./../build/bin/ad_square_example

# # ------------------------------------ # 
# # Ordinary Differential Equation (ODE) #
# # ------------------------------------ # 
# echo -e "\nExecuting: ode_seair"
# ./../build/bin/ode_seair_example

# echo -e "\nExecuting: ode_secirvvs"
# ./../build/bin/ode_secirvvs_example

# ------------ # 
# Optimization #
# ------------ # 
echo -e "\nExecuting: ode_seair_optimization"
./../build/bin/ode_seair_optimization

# echo -e "\nExecuting: ode_secirvvs_optimzation"
# ./../build/bin/ode_secirvvs_optimzation

# echo -e "\nExecuting: ode_secirvvs_optimzation2"
# ./../build/bin/ode_secirvvs_optimzation2


# echo -e "\nExecuting: IPOPT_1"
# ./../build/bin/ipopt_1

# echo -e "\nExecuting: IPOPT_2"
# ./../build/bin/ipopt_2

# echo -e "\nExecuting: IPOPT_3"
# ./../build/bin/ipopt_3