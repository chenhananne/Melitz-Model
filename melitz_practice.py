import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve
import os

# Given parameters value

theta=3.8  # elasticity of substitution
k=3.4  # firms' productivity homogeneity/Pareto shape
c_bar=2  # fixed labor cost for entry
c=0.1 # fixed labor cost for operation
b=0.2  # minimum value of firms' productivity
beta=0.975 # firm survive probability every period

# Calculate the parameters
mu = k * theta / (theta - 1)
Theta = (theta - 1)**(theta - 1) / theta**theta
B = c_bar * (1 - beta) * ((c / Theta)**(mu /theta)) / b**k
A = (1 / c) * (1 / theta - 1 / mu)

# print the results
print(f"mu = {mu}")
print(f"Theta = {Theta}")
print(f"B = {B}")
print(f"A = {A}")

# Given variables in billion euros
Y_h=17070 # GDP in EU
Y_j=16540 # GDP in China
Y_R=59890 # GDP rest of world
X_hj=514 # export from EU to China
X_jh=514 # export from China to EU
X_hR=2337 # export from EU to rest of world
X_Rh=2337 # exports from rest of world to EU
X_jR=2626 # exports from China to rest of world
X_Rj=2626 # exports from rest of world to China

#calculate domestic sales

X_hh =Y_h - X_hj - X_hR # domestic sales in China
X_jj =Y_j - X_jh - X_jR # domestic sales in EU
X_RR =Y_R - X_Rh - X_Rj # domestic sales in rest of world


# print the results
print(f"X_hh = {X_hh}")
print(f"X_jj = {X_jj}")
print(f"X_RR = {X_RR} ")

# Define the variables
D_h, D_j, D_R, P_h, P_j, P_R, delta_hj, delta_hR, delta_jR, W_h,  W_j, W_R, S_h,  S_j, S_R, M_h, M_j, M_R = symbols (
    'D_h D_j D_R P_h P_j P_R delta_hj delta_hR delta_jR W_h  W_j W_R S_h  S_j S_R M_h M_j M_R ')

# Equations
eq1 = Eq(D_h, (Y_h**2) / X_hh)
eq2 = Eq(D_j, (Y_j**2) / X_jj)
eq3 = Eq(D_R, (Y_R**2) / X_RR)

eq4 = Eq(P_h, 1)
eq5 = Eq(P_j, (D_j / D_h)**(1 / (2 * (1 - mu))))
eq6 = Eq(P_R, (D_R / D_h)**(1 / (2 * (1 - mu))))

eq7 = Eq(delta_hj,  ((X_hj * X_jh) / (X_jj * X_hh))** (1 / (2 * (1 - mu))))
eq8 = Eq(delta_hR,  ((X_hR * X_Rh) / (X_RR * X_hh))** (1 / (2 * (1 - mu))))
eq9 = Eq(delta_jR,  ((X_jR * X_Rj) / (X_RR * X_jj))** (1 / (2 * (1 - mu))))

eq10 = Eq(W_h, (B*A*mu / D_h) **(- 1/mu))
eq11 = Eq(W_j, (B*A*mu*(P_j**(1 - mu))/D_j )** (- 1 / mu))
eq12 = Eq(W_R, (B*A*mu*(P_R**(1 - mu))/D_R )** (- 1 / mu))

eq13 = Eq(S_h, Y_h / W_h)
eq14 = Eq(S_j, Y_j / W_j)
eq15 = Eq(S_R, Y_R / W_R)

eq16 = Eq(M_h, A * S_h * P_h / D_h)
eq17 = Eq(M_j, A * S_j * P_j / D_j)
eq18 = Eq(M_R, A * S_R * P_R / D_R)

#solve the equations
solutions = solve ([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18],
                  (D_h, D_j, D_R, P_h, P_j, P_R, delta_hj, delta_hR, delta_jR, W_h,  W_j, W_R, S_h,  S_j, S_R, M_h, M_j, M_R))

for sol in solutions:
    print(sol)
#initialize an empty list to store the evaluated solutions
solution_list=[]

# Iterate through each solution, check if it is a tuple, and evaluate each element if necessary
for sol in solutions:
    if isinstance(sol, tuple):  # If the solution is a tuple, we need to handle each element
        solution_list.extend([s.evalf() for s in sol])  # Apply evalf to each element
    else:
        solution_list.append(sol.evalf())  # If it's a single solution, just apply evalf
#list of variable names
variables =['D_h', 'D_j', 'D_R', 'P_h', 'P_j', 'P_R', 'delta_hj', 'delta_hR', 'delta_jR', 'W_h', 'W_j', 'W_R', 'S_h', 'S_j', 'S_R', 'M_h', 'M_j', 'M_R']

results = {}
for i in range(len(variables)):
    globals()[variables[i]] = solution_list[i]
    results[variables[i]] = solution_list[i]


print(results)

# Check
LHS_eq1= D_h
RHS_eq1= Y_h + (delta_hj/P_j)**(1-mu)*Y_j + (delta_hR/P_R)**(1-mu)*Y_R
LHS_eq2= D_j
RHS_eq2= Y_j + (delta_hj*P_j)**(1-mu)*Y_h + (delta_jR*P_j/P_R)**(1-mu)*Y_R
LHS_eq3= D_R
RHS_eq3= Y_R + (delta_hR*P_R)**(1-mu)*Y_h + (delta_jR*P_R/P_j)**(1-mu)*Y_j


#PRINT THE RESULTS AND COMPARE
print(f"Equation 1: LHS = {LHS_eq1}, RHS = {RHS_eq1}, Difference = {abs(LHS_eq1 - RHS_eq1)}")
print(f"Equation 2: LHS = {LHS_eq2}, RHS = {RHS_eq2}, Difference = {abs(LHS_eq2 - RHS_eq2)}")
print(f"Equation 3: LHS = {LHS_eq3}, RHS = {RHS_eq3}, Difference = {abs(LHS_eq3 - RHS_eq3)}")



# Check if the differences are negligible
tolerance = 1e-5  # Define an acceptable tolerance level
if abs(LHS_eq1 - RHS_eq1) < tolerance and abs(LHS_eq2 - RHS_eq2) < tolerance and abs(LHS_eq3 - RHS_eq3) < tolerance:
    print("Model is calibrated successfully!")
else:
    print("Model calibration check failed!")

#endogenous variables
 #average avenue, profit, labor
r_h=A**(-1)*W_h*D_h/P_h
pi_h=mu**(-1)*r_h
L_h=r_h/W_h*(1-1/mu)

r_j=A**(-1)*W_j*D_j/P_j
pi_j=mu**(-1)*r_j
L_j=r_j/W_j*(1-1/mu)

r_R=A**(-1)*W_R*D_R/P_R
pi_R=mu**(-1)*r_R
L_R=r_R/W_R*(1-1/mu)
print (r_h,pi_h,L_h,r_j,pi_j,L_j,r_R,pi_R,L_R)

#Minimum productivity of a firm located in h for producing into market h and j and R
Theta_hh=(c/Theta)**(1/(theta-1))*(W_h/P_h)**(theta/(theta-1))
Theta_hj=(c/Theta)**(1/(theta-1))*(delta_hj*W_h/P_j)**(theta/(theta-1))
Theta_hR=(c/Theta)**(1/(theta-1))*(delta_hR*W_h/P_R)**(theta/(theta-1))
print(Theta_hh,Theta_hj,Theta_hR)

#P(exist_h) value
P_exist_h= c_bar*(1-beta) / (S_h/M_h - L_h )

print(P_exist_h)

#average productivity
varphi_bar=abs(k*b/(1-k))

print (varphi_bar)