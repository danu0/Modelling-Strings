import numpy as np
import matplotlib.pyplot as plt
from spectrumbase import SpectrumBase

class Spectrum(SpectrumBase):
        
    def index_of_nu(self, nu):
        """
        Function that returns the index of nu in the list self.nu
        nu: the frequency to consider
        """ #add comments innit
        eps = 1e-12
        for i in range(1,len(self.nu)):
            if nu-self.nu[i]<eps:
                return i-1
        return len(self.nu)-1

            
    def anharmonicity(self, nu_max, fname="-"):
        """
        Calculate and plot the detuning of nu_max in cents.
        nu_max: maximum frequency considered
        fname: name of the figure/savefile
        """
        imax = Spectrum.index_of_nu(self, nu_max)
        new_nu = self.nu[:imax+1]
        cents = []
        print(new_nu[18])
        for i in range(0,len(new_nu)):
            cents.append(1200*(np.log(new_nu[i]/((i+1)*new_nu[0]))/np.log(2)))
        if fname != "-":
            plt.plot(range(0,len(new_nu)),cents)
            plt.xlabel('i', fontsize=24) #change fontsize
            plt.ylabel('cents',fontsize=24) #change fontsize
        if fname != ("-" or " "):
            plt.savefig(fname+".png")
        plt.show()
        print(cents)
        return self.nu[imax], cents[-1]

    
class SpectrumGuitar(Spectrum):
    """ Class which simulates a guitar string."""

    def __init__(self, N):
        """ Initialiser, with N the number of beads."""
        self.N=N
        m = 1
        c = 329.6
        L = 0.63
        dx = L/(N+1)
        k = c**2/dx**2
        
        A = -2*np.eye(N)               # diagonal is -2
        A[1:N,   0:N-1] += np.eye(N-1) # a[i,j-1]=1
        A[0:N-1, 1:N]   += np.eye(N-1) # a[i-1,j]=1
        C = k*A/m

        super().__init__(C, "guitar")

class SpectrumPiano(Spectrum):
    """Class which simulates a piano string."""

    def __init__(self,N):
        """ Initialiser, with N the number of beads."""
        self.N = N
        m = 1
        c_1 = 329.6
        c_2 = 1.25
        L = 0.63
        dx = L/(N+3)
        
        A_2 = -2*np.eye(N)
        A_2[1:N,   0:N-1] += np.eye(N-1)
        A_2[0:N-1, 1:N]   += np.eye(N-1)
        
        A_4 = 6*np.eye(N)
        A_4[1:N, 0:N-1] += -4*np.eye(N-1)
        A_4[0:N-1, 1:N] += -4*np.eye(N-1)
        A_4[2:N, 0:N-2] += np.eye(N-2)
        A_4[0:N-2, 2:N] += np.eye(N-2)
        
        B_2 = np.zeros((N,N))
        B_2[0,0] = 0.5
        B_2[N-1,N-1] = 0.5 
        
        B_4 = np.zeros((N,N))
        B_4[0,0] = -2
        B_4[N-1,N-1] = -2
        B_4[0,1] = 2
        B_4[N-1, N-2] = 2
        B_4[1,1] = -1
        B_4[N-2,N-2] = -1
        
        C = m*((c_1/dx)**2 * (A_2 + B_2))-((c_2)**2 /(dx**4) * (A_4 + B_4))

        super().__init__(C, "piano")
