import numpy as np
import time as libtime
from math import *
from Constants_GCM import *


def step_Bio(vec,T,FT,NC,X0,traits,p,growth_lim=False,debug=False):

    [H,C,G,N2]  = vec
    alphaC     = np.exp(9345.17/T-167.8108+23.3585*np.log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35.0)

    if   FT == 'meth'      : import Methanogens as ft        ; dgcat = ft.DeltaGcat(T,H/alphaH,C/alphaC,G/alphaG)
#     if   FT == 'meth'      : import Methanogens as ft        ; dgcat = ft.DeltaGcat(T,H,C,G)
    elif FT == 'NO3metht'  : import NO3Methanotrophs as ft   ; dgcat = ft.DeltaGcat(T,C,G,NO3,NO2)
    elif FT == 'NO2metht'  : import NO2Methanotrophs as ft   ; dgcat = ft.DeltaGcat(T,C,G,NO2,N2)
    elif FT == 'H2SO4metht': import H2SO4Methanotrophs as ft ; dgcat = ft.DeltaGcat(T,C,G,H2SO4,H2S)
    elif FT == 'acet'      : import Acetogens1 as ft         ; dgcat = ft.DeltaGcat(T,CO,C,CH3COOH)
    elif FT == 'acet2'     : import Acetogens2 as ft         ; dgcat = ft.DeltaGcat(T,C,H,CH3COOH)
    elif FT == 'acett'     : import Acetotrophs as ft        ; dgcat = ft.DeltaGcat(T,C,G,CH3COOH)
    elif FT == 'ferm'      : import Fermentors as ft         ; dgcat = ft.DeltaGcat(T,X0D,N,H,CH3COOH,C)
    elif FT == 'photoH2'   : import PhotoH2 as ft

    Vc     = traits[2]
    Qc     = traits[3]
    ks     = traits[4]
    mg     = traits[6]
    kd     = traits[7]
    mort   = traits[8]
    thresh = traits[9]
    qmax   = traits[5]
#    qmax   = traits[5]*(10*thresh-X0)/(8*thresh)
    slope  = traits[10]
    gmax   = traits[11]

    if dgcat < 0:
        if   FT == 'meth'      : qcat = ft.QCat(dgcat,H,C,qmax,ks)
        elif FT == 'acet'      : qcat = ft.QCat(dgcat,CO,qmax,ks)
        elif FT == 'acett'     : qcat = ft.QCat(dgcat,CH3COOH,qmax,ks)
        mreq     = ft.Mreq(mg,dgcat)                                      # Rate of cell maintenance
        decay    = ft.Decay(mreq,qcat,dgcat,kd)                           # Rate of decay
        if(debug == True): print(T,H,C,N2,X0/Vc)   
        dgana    = ft.DeltaGana(T,H,C,N2,X0/Vc)     # Energy required for the anabolic reaction
        Lam      = max(-((dgana+ft.dgdiss)/dgcat),0)                      # Metabolic coupling
        Y        = ft.Yl(Lam)                                             # Metabolic stochiometry
        if(debug == True): print(dgcat,dgana)   
        slim     = ft.Slim([H,C,N2],[Y[i] for i in [0,1,9]])              # Limiting substrate
        QMet_t   = ft.QMet(dgcat,qmax,ks,slim)                            # Metabolic rate
        qana     = ft.QAna(dgcat,dgana,Lam,qcat,QMet_t,mreq,qmax,ks,slim) # Anabolic rate
        qcat     = qcat                                                   # Catabolic rates
        if growth_lim == True: new_cell = ft.Gamma(thresh,slope,gmax,X0,Qc)
        else: new_cell = qana/(Qc)
        dNC      = (new_cell - decay - mort)*NC
        dX0      = qana - new_cell*X0

    else:
        mreq     = ft.Mreq(mg,dgcat)                                      # Rate of cell maintenance        
        qcat     = 0
        qana     = 0
        dgana    = ft.DeltaGana(T,H,C,N2,X0/Vc)     # Energy required for the anabolic reaction
        decay    = kd
        new_cell = ft.Gamma(thresh,slope,gmax,X0,Qc)
        dNC      = (new_cell - decay - mort)*NC
        dX0      = qana - new_cell*X0

    dX0D         = (decay+mort)*(X0*NC)

    #print(' T=',T,' dgc=',dgcat,' dga=',dgana,' qc=',qcat,' qa=',qana,'\n')
    return(dgcat,dNC,dX0,qana,qcat,mreq,new_cell,ft.Catabolism,ft.Anabolism)

