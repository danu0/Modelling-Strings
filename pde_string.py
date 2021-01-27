import numpy as np
import matplotlib.pyplot as plt
from pde_base import PDE_Base

class PDE_String(PDE_Base):
    def F(self, t, v):
        """ Equation to solve:
        dy[i]/dt = doty[i]
        d doty[i]/dt = k*(y[i+1]+y[i-1]-2*y[i])
        t: current time
        v: current function as a vector of dim N
        """
        F = np.zeros(2*self.N)
        F[0:self.N] = v[self.N:2*self.N]
        F[self.N+1:2*self.N] = self.c**2/self.dx**2*(v[2:self.N+1] + v[:self.N-1] - 2*v[1:self.N])
        return(F)
        

    def plot_snapshots(self, Nsnapshot, T, fig_dt, fname=""):
        """
        Generates the plots of y(x,t_i) from t_i=0 to t_Ns = T/2 on the same plot
        Nsnapshot: How many times to plot
        T:
        fig_dt: A small change in time
        fname: Name of the figure
        """
        t_int = len(self.t_list)/Nsnapshot
        x = np.linspace(0,L,N)
        for i in range(Nsnapshot):
            t = int(i*t_int)
            y = self.V_list[t][:N]
            t_label = self.t_list[t]
            plt.plot(x,y, self.formats[i], label=t_label)
        plt.legend()
        plt.xlabel('x', fontsize = 24)
        plt.ylabel('y', fontsize = 24)
        plt.show()
        if fname != "":
            plt.savefig(fname)

    def FFT_mod(self, x):
        """
        Returns the co-efficients of the fourier transform as an array.
        x: Position of the particle
        """

        eps = self.L/(2*self.N)
        for i in range(self.N):
            if abs(i*self.L/self.N -x) <eps:
                index =  i
        z = np.array(self.V_list)
        v = z[:,[index]]
        v = np.fft.fft(v)

        modfft = np. abs (v [0: len(self.V_list)//2+1]) *2/ self.N
        modfft [0] /= 2 # constant term
        modfft [self.N //2] /= 2 # higher frequency
        print(len(modfft))
        return modfft

if __name__ == "__main__":  
    N = 200    # number of points

    # STRING PARAMETERS
    L = 0.63      # length of string: 63cm
    c = 329.6     # sound speed on the string
    T = 2*L/c     # period
    tmax = 20*T   # integration time 
    fig_dt=T/200  # interval between figure capture
    spr = PDE_String(c, L, N)
    
    # 1) lowest mode
    x = np.array(range(0,N))*spr.dx # x coordinates
    y = np.sin(x*np.pi/L)
    doty = np.zeros(N)
    spr = PDE_String(c, L, N)
    spr.reset(0,y,doty)
    spr.iterate(tmax,fig_dt)
    print("N T=",spr.t_list[-1]/T)
    spr.mk_figs(7, T, L/3, 20000, 5, "Mode1")

    # 2) plucked string
    y = []
    for i in x:
        if 0 <= i <= L/3:
            y.append(3*i/L)
        else:
            y.append(1.5*(1-i/L))
    doty = np.zeros(N)
    spr = PDE_String(c, L, N)
    spr.reset(0,y,doty)
    spr.iterate(tmax, fig_dt)
    print("N T=",spr.t_list[-1]/T)
    spr.mk_figs(7, T, L/3, 20000, 5, "Mode2")


    # 3) hit string
    # TO COMPLETE (coding task 3)
    x = np.array(range(0,N))*spr.dx # x coordinates
    y = np.zeros(N)
    doty = np.exp(-((x-L/3)/0.05)**2)
    spr = PDE_String(c, L, N)
    spr.reset(0,y,doty)
    spr.iterate(tmax,fig_dt)
    print("N T=",spr.t_list[-1]/T)
    spr.mk_figs(7, T, L/3, 20000, 5, "Mode13")
