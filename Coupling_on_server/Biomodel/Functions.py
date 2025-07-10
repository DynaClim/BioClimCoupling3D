import numpy as np
import pandas as pd
import pickle
from math import *

#Loads the parameters of the biological model infered from Corkrey et al.  
file = open('Bio_model_parameter', 'rb')
Par = pickle.load(file)

#Returns trait values according to the parameters loaded above
def Run_def_FT_random(FT,a_Qc,b_Qc,a_qmax,b_qmax,c_qmax,a_mg,b_mg,c_mg,a_ks,a_kd,a_mort,a_gmax,T,rc=1):
    Topt       = T                                       # Optimal temperature 
    Vc         = (4/3)*pi*rc**3                          # Cell Volume (µm3)
    Qc         = (a_Qc*Vc**(b_Qc))                       # Cell C content (molX.Cell-1)
    qmax       = a_qmax*np.exp(c_qmax*Topt)*Vc**b_qmax   # Maximum metabolic rate (d-1) (Gillooly 2001)
    mg         = 24*a_mg*np.exp(c_mg*Topt)*Vc**(b_mg)    # Maintenance (J.Cell-1d-1) (Tijhuis 93 + Litchman 2007)
    ks         = a_ks*1e-8                               # Metabolic half-saturation constant (mol.L-1)
    kd         = a_kd*0.5                                # Maximum decay rate (d-1)
    #mort       = a_mort*0.1                              # Basal mortality rate (d-1) INITIAL VALUE
    mort       = a_mort*1                                # Basal mortality rate (d-1)
    thresh     = Qc                                      # Inflexion point of the division function (molX.Cell-1)
    slope      = 10                                      # Slope of the division function
    gmax       = a_gmax*1                                       # Maximum division rate (-d)
    
    #print(Topt,Vc,Qc,qmax,mg,ks,kd,'mort =',mort,thresh,slope,'gmax =',gmax)

    return ([Topt,rc,Vc,Qc,ks,qmax,mg,kd,mort,thresh,slope,gmax])

def ReturnTraits(T,Par,rc):
    [a_qmax_fact,c_qmax,a_mg_fact,c_mg] = Par
    
    a_Qc        = 18E-15
    b_Qc        = 0.94
    #c_mg        = 0.08
    a_mg        = np.exp(-46.72+0.08*298 - c_mg*298)*a_mg_fact
    b_qmax      = 0.82
    a_qmax      = np.exp(-55.76+0.1*298 - c_qmax*298)*a_qmax_fact
    b_mg        = 0.67    
    a_ks        = 1
    a_kd        = 1
    a_mort      = 0.1    
    a_gmax      = 20
    Topt,rc,Vc,Qc,ks,qmax,mg,kd,mort,thresh,slope,gmax = Run_def_FT_random('meth',a_Qc,b_Qc,a_qmax,b_qmax,c_qmax,a_mg,b_mg,c_mg,a_ks,a_kd,a_mort,a_gmax,T,rc=1)
    return([Topt,rc,Vc,Qc,ks,qmax,mg,kd,mort,thresh,slope,gmax])

