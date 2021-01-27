
import numpy as np
from spectrum import *

N=1000
S = SpectrumGuitar(N=N)
S.mk_figs(7, nu_max=20000)

nu, cents = S.anharmonicity(nu_max=20000, fname="-")
print("N={}, node threshold: nu={:.3f}, {:.3f} cents".format(N, nu, cents))
     
    
