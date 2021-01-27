import sys
import ode_rk4
import numpy as np
import matplotlib.pyplot as plt

class PDE_Base(ode_rk4.ODE_RK4):
  """ A class to solve partial differential equations by
      rewriting them as a system of ordinary differential
      equations.
  """
  
  def __init__(self, c, L, N):
      """Initialiser for our class. 
         c : string constant
         L : string length
         N : number of nodes
      """  
      self.c = c
      self.N = N # number of points
      self.L = L
      self.dx = L/(N-1)
      self.dt = self.dx/(self.c*4)
      self.formats=["r-","g-","b-","c-","m-","y-","k-", "r-.","g-.","b-.","c-.","m-.","y-.","k-."]
      # we use an empty array for now
      super().reset(np.empty(2*self.N, dtype='float64'), self.dt, t0=0)

      
  def reset(self, t, y, doty):
      """ Initialise the paramaters (t and initial conditions)
         y : spring extension : array of size N
         doty : dy/dt         : array of size N 
      """
      self.t = t
     
      if (len(y) != self.N) or (len(doty) != self.N):
         print("spring.__init__ ERROR: invalid size for y or doty")
         sys.exit(5)

      # Create array to hold the functions y and dy/dt
      V = np.empty(2*self.N, dtype='float64')
      # v[0:N]   : y
      # v[N:2*N] : dy/dt 
      V[0:self.N] = y
      V[self.N:2*self.N] = doty
      
      super().reset(V, self.dt, t0=t)

  def F(self, t, v):
      """ equation to solve: 
          dy[i]/dt = doty[i]  
          d doty[i]/dt = k*(y[i+1]+ y[i-1]-2*y[i])
          t: current time 
          v: current function as a vector
      """
      F = np.zeros(2*self.N)
      for i in range(self.N):
        # dv[i]/dt = d doty[i]
        F[i] = v[self.N+i]
      
        # d doty[i]/dt = c^2/dx^2(y[i+1]+y[i-1]-2y[i])
        if(i == 0) or (i == self.N-1):
          F[self.N+i] = 0
        else:
          F[self.N+i] = self.c**2/self.dx**2*np.array(v[i+1]+v[i-1]-2*v[i])
      return(F)

  def plot_follow_node(self, x, fname=""):
      """ Plot the oscillation of one node, at position x, during the 
          integration.
          x: position of the node
          fname: filename for figure outout
      """
      i = int((self.N-1)*x/self.L)+1 # add 1 because plot() substracts it
      self.plot(i,0) # add 1 because plot() substracts it
      plt.xlabel("t",fontsize=24)
      plt.ylabel("y["+str(x)+"]",fontsize=24)
      plt.tight_layout(pad=0.4)
      if (fname != ""):
        plt.savefig(fname)
      plt.show()
    

  def plot_modfft(self, freq_max, NFFT, fname=""):
      """ Plot modFFT for up to freq_max
          freq_max: upper frequency for spectrum graph
          NFFT: number of nodes to follow
          fname: filename for figure outout
      """
      modFFT = 0
      for frac in np.linspace(0.5/NFFT,0.5,NFFT):
          modFFT += self.FFT_mod(self.L*frac)
      modFFT /= NFFT
      
      samplesize = len(self.V_list)
      freq = np.fft.fftfreq(samplesize,self.t_list[-1]/samplesize)
      # select frequencies of interest
      freq = np.abs(freq)[0:samplesize//2+1]
      sel_m = modFFT[freq < freq_max]
      sel_f = freq[freq < freq_max]
      plt.xlabel(r'$\nu$ (Hz)',fontsize=24)
      plt.ylabel("|FT|",fontsize=24)
      plt.loglog(sel_f,sel_m,"b.")
      plt.tight_layout(pad=0.4)
      if (fname != ""):
        plt.savefig(fname)
      plt.show()

  def mk_figs(self, Nsnapshot, T, node_x, freq_max, NFFT, name=""):
      """
      """
      #  snapshots of string profile
      fig_dt = self.t_list[1]-self.t_list[0]
      self.plot_snapshots(Nsnapshot , T, fig_dt, "string_profiles_N{}_{}.pdf".format(self.N, name))

      # Plot Profile f(node_x) as a function of time
      self.plot_follow_node(node_x, "node_profiles_N{}_{}.pdf".format(self.N, name))

      # Compute average FFT sampled at 5 different places.
      self.plot_modfft(freq_max, NFFT, "string_spectrum_N{}_{}.pdf".format(self.N, name))
      
  def plot_snapshots(self, Ns, T, fig_dt, fname=""):
      pass

  def FFT_mod(self, x):
      pass

############################################################
        
if __name__ == "__main__":
  N = 200    # number of points

  # STRING PARAMETERS
  L = 0.63      # length of string: 63cm
  c = 329.6     # sound speed on the string
  T = 2*L/c     # period
  tmax = 20*T   # integration time 
  fig_dt=T/200  # interval between figure capture
  
  # Needed to initialise dx
  spr = PDE_String(c, L, N)

  # 1) lowest mode
  x = np.array(range(0,N))*spr.dx # x coordinates
  y = np.sin(x*np.pi/L)
  doty = np.zeros(N)

  spr = PDE_String(c, L, N)
  spr.reset(0, y, doty)
  spr.iterate(tmax, fig_dt)
  print("N T=", spr.t_list[-1]/T)
  spr.mk_figs(7, T, L/3, 20000, 5, "Mode1")
      
  # 2) plucked string 
  # TO COMPLETE

  # 3) hit string 
  # TO COMPLETE
