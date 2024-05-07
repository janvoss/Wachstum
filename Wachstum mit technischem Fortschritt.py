import plotnine as gg
import numpy as np
import pandas as pd

# Parameter
alpha = 0.4
delta = 0.1
s = 0.2
L = 1 # einfachster Fall
n = 0 # Bevölkerungswachstum

#neu
g=0.012

# Zeitschritte
t = np.linspace(0, 100, 101)

# Kapitalakkumulation
##Vorbeteitung

K_t = np.zeros(len(t))
Y_t = np.zeros(len(t))
I_t = np.zeros(len(t))
L_t = np.zeros(len(t))
k_t = np.zeros(len(t))
y_t = np.zeros(len(t))


#neu
A_t =np.zeros(len(t))


## Ausgangswerte

K_t[0] = 2

#neu
A_t[0] = 1
Y_t[0] = A_t[0]*(K_t[0]**alpha)*(L**(1-alpha)) #CD-Funktion
I_t[0] = s*Y_t[0]
L_t[0] = L
y_t[0] = Y_t[0]/L_t[0]
k_t[0] = K_t[0]/L_t[0]

## Dynamik


for i in range(1, len(t)):
    
    A_t[i] = (1+g)A_t[i-1]
    
    K_t[i] = (1-delta)*K_t[i-1] + I_t[i-1]
    Y_t[i] = A_t[i]*(K_t[i]**alpha) * (L**(1-alpha))
    I_t[i] = s*Y_t[i]
    L_t[i] = (1+n)*L_t[i-1]
    k_t[i] = K_t[i]/L_t[i]
    y_t[i] = Y_t[i]/L_t[i]



## In Datenframe speichern
df = pd.DataFrame({'t': t, 'Kapital': K_t, 'Y': Y_t,
                    'Investitionen': I_t, 'Abschreibungen': delta*K_t,
                    'Bevölkerung': L_t, 'Kapital pro Kopf': k_t,
                  'y': y_t})
