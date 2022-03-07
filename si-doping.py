import numpy as np
import decimal
import matplotlib.pyplot as plt
from math import e
from decimal import Decimal,localcontext
import math
from scipy.optimize import fsolve
import math
from math import exp
from scipy.optimize import root
from scipy.interpolate import make_interp_spline, BSpline

def fun(x):
    # x[0]= V = x , y=x[1]= Oi , z=D =x[2]
    
    kt = kb*T

    return [Oinit - ((exp(2.21/kt))*x[0]*x[1]) - (2 *exp(4.05/kt) * x[0] *(x[1])*x[1]) - (exp(2.4/kt) * x[0]*x[1]*x[2]) -x[1],
            Vinit - (exp(2.21/kt)*x[0]*x[1]) - (exp(4.05/kt) * x[0] *(x[1])*x[1])  - (2 *exp(2.52/kt)*(x[0])*x[0] ) - (exp(0.4/kt)* x[0] *x[2]) - (exp(2.4/kt) * x[0]*x[1]*x[2] )- x[0],   
            Dinit - (exp(0.4/(kt))* x[0] *x[2]) - ((exp(2.4/kt)) * x[0]*x[1]*x[2]) - x[2] ]


#Boltzmann constant
decimal.getcontext().prec = 10000
#kb = 8.617333262e-5
kb = 8.6e-5

#array initialization
Oi_arr = []
V_arr = []
VV_arr = []
VOi_arr = []
DV_arr = []
DVOi_arr = []
VO2_arr = []
T_arr = []
D_arr= []

#dedomena
Vinit = 10**(18)
Oinit = 10**(18)
Dinit = 10**(19)


for T in range(800,1500,5):
    T_arr.append((T))
    kt = kb * T
    #LYSI SYSTIMATOS
    x0 = [10**10,10**(2),5*10**11]
    sol = root(fun,x0,method='hybr',tol=10e-15)

    V= sol['x'][0]
    Oi = sol['x'][1]
    D = sol['x'][2]
  
    D_arr.append(np.log10(D))
    V_arr.append(np.log10(V))
    Oi_arr.append(np.log10(Oi))
    
    # [VV]/[V][V]
   
    VV = exp(2.52/kt) * V * V
    VV = np.log10(VV)
    VV_arr.append(VV)

    #[VOi]/[V][Oi]

    VOi = exp(2.21/kt) * V * Oi
    VOi = np.log10(VOi)
    VOi_arr.append(VOi)

    #[DV]/[D][V]
    DV = exp(0.4/kt) * D * V
    DV = np.log10(DV)
    DV_arr.append(DV)

    #[DVOi]/[D][V][Oi]
    
    DVOi = exp(2.4/kt) * D * V * Oi
    DVOi = np.log10(DVOi)
    DVOi_arr.append(DVOi)
    
    #[VO2]/[V][Oi][Oi]

    VO2 = exp(4.05/kt) * V * Oi *Oi
    VO2 = np.log10(VO2)
    VO2_arr.append(VO2)
   
    print("T=" + str(T))
    print("[DVOi] = " + str(DVOi))
    print("[VV] = " + str(VV))
    print("[VOi]= " +str(VOi))
    print("[VO2]= " + str(VO2))
    print("[DV] = " + str(DV))
    print("[DVOi] = " + str(DVOi))
    print("Vacancy(V): " + str(V))
    print("Dopant(D): "+ str(D))
    print("Oxygen(Oi): " + str(Oi))
    print("\n")



plt.plot(T_arr,D_arr,label = '[D]')
plt.plot(T_arr,V_arr, label = '[V]')
plt.plot(T_arr,Oi_arr, label = '[Oi]]')
plt.plot(T_arr,VV_arr, label = '[VV]')
plt.plot(T_arr,VOi_arr,label = '[VOi]')
plt.plot(T_arr,DV_arr,label = '[DV]')
plt.plot(T_arr,DVOi_arr,label = '[DVOi]')
plt.plot(T_arr,VO2_arr,label= '[VO2]')
plt.plot(T_arr,V_arr,label= '[V]')
plt.xlabel("Temperature (K)")
plt.ylabel("[V] concetration (lοg10 scale cm^3)")
plt.legend()
plt.show()
