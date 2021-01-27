import numpy as np
import math
import matplotlib.pyplot as plt

class SpectrumBase:
    def __init__(self, C, name=""):
        """ Solve the eigenvalue problem C v = lambda v.
            C:    Input matrix.
            name: Name for this type of string (used for figures).
        """
        self.Mat= C
        (Ev , V) = np.linalg.eig(C)
        idx = Ev.argsort()[::-1]   
        self.Ev = Ev[idx]
        self.V = V[:,idx]
        self.nu = np.sqrt(-self.Ev)/(2*math.pi)
        self.pltformats = ["r-", "g-", "b-", "c-", "m-", "y-", "k-"]
        self.name = name

    def index_of_nu(self, nu):
        pass

    def anharmonicity(self, nu_max, fname="-"):
        pass

    def anharmonicity_summary(self, nu_max):
        pass

    def plot_modes(self, Nm, fname=""):
        """ Plot the lowest Nm mode profiles of the string.
            Nm:    Number of modes to display.
            fname: Filename for figure output (if non-empty).
        """
        for i in range(0, Nm):
            # We must add the (fixed) end nodes by hand.
            W = np.pad(self.V[:,i], (1,1), 'constant', constant_values=(0, 0))
            plt.plot(W, self.pltformats[i%len(self.pltformats)],
                    label=r'$\omega={1:.3g}$'.format(i,math.sqrt(-self.Ev[i])))

        plt.legend(loc='lower right', bbox_to_anchor=(1.45, 0.0), prop={'size':10})

        ax = plt.gca()
        ax.set_aspect(1.0/ax.get_data_ratio()*0.5)
        plt.tight_layout(rect=[0, 0, 1.1, 1])
        plt.xlabel("i", fontsize=24)
        plt.ylabel("y", fontsize=24)
        plt.tight_layout(pad=0.1)
        if(fname != ""):
            plt.savefig(fname)
        plt.show()
        
    def plot_spectrum(self, nu_max, fname=""):   
        """ Plot the spectrum of the string in blue.
            Plot multiples of the fundamental frequency in red.
            nu_max: Largest frequency mode to display.
            fname:  Filename for figure output (if non-empty).
        """
        Nmax = self.index_of_nu(nu_max)
        nu_index =  self.nu/self.nu[0]

        # Plot extected index in blue
        plt.plot(nu_index[0:Nmax], "b-")

        # Plot actual index in red
        plt.plot(range(1, Nmax+1), "r-")
        plt.xlabel("i", fontsize=24)
        plt.ylabel(r'$\nu/\nu_0$', fontsize=24)
        plt.tight_layout(pad=0.4)
        if(fname != ""):
           plt.savefig(fname)
        plt.show()

    def mk_figs(self, Nmodes, nu_max):
        """ Plot normal modes, spectrum and anharmonicity.
            Nmodes: Number of normal mode profiles to plot.
            nu_max: Max frequency for frequency plot and anharmonicity.
            fname : Filename for figure output (if non-empty).
        """

        self.plot_modes(Nmodes, "lowest_N{}_modes_{}.pdf".format(Nmodes, self.name))
        self.plot_spectrum(nu_max, "spectrum_numax{}_{}.pdf".format(nu_max, self.name))
        self.anharmonicity(nu_max, "anharmonicity_numax{}_{}.pdf".format(nu_max, self.name))


