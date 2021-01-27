
import numpy as np
from spectrum import *

N=1000 #q4
S = SpectrumPiano(N=N) #detuning of 100cents occurs at cents[18] ish want self.nu[17]
S.mk_figs(7, nu_max=20000)

nu, cents = S.anharmonicity(nu_max=20000, fname="-")
print("N={}, node threshold: nu={:.3f}, {:.3f} cents".format(N, nu, cents))
     
    
