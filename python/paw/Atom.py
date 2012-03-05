
from __future__ import division
from numpy import array, linspace, zeros, sign
from pylab import plot, ion, ioff

DEBUG = 1
REALLY_SMALL = 1.0e-6

def numerov(F, h, grid, x0, x1, direction="forward"):
  
  x = zeros(grid.shape)
  x[0] = x0
  x[1] = x1
  if direction == "backward":
    F = F[::-1]
  
  for i, t in enumerate(grid[2:]):   
    #x[i+2] = ((2.0 + 5.0 * (h ** 2) * F[i+1] / 6.0) * x[i+1] - (1.0 - (h ** 2) * F[i] / 12.0) * x[i]) / (1.0 + (h ** 2) * F[i+2] / 12.0)
    #if DEBUG: print "Index %d: x[i+1] = %f, x[i] = %f, F[i+2] = %f, F[i+1] = %f, F[i] = %f." % (i, x[i+1], x[i], F[i+2], F[i+1], F[i])
    x[i+2] = (x[i+1] * (2.0 + (5.0 / 6.0) * (h ** 2) * F[i+1]) - x[i] * (1.0 - (h ** 2) * F[i] / 12.0)) / (1 - (h ** 2) * F[i+2] / 12.0)
  
  if direction == "forward":
    return x
  else:
    return x[::-1]

def nodecount(y):
  
  # We use a dirty hack here - if the array y happens to actually contain a zero,
  # we will step over it twice, so just add a 0.5 each time. Yuck!
  
  count = 0
  for i in range(len(y[1:])):
    if sign(y[i]) == 0 or sign(y[i-1]) == 0:
      count += 0.5 # Dirty hack!
    elif sign(y[i]) != sign(y[i-1]):
      count += 1.0
  
  return count
      
def eigensearch(guess, n, l, f, h, grid):

  # Initially just use half-interval method and only shoot from the centre out.

  if DEBUG: print "Searching for eigenvalue correponding to n=%d, l=%d." % (n, l)
  if DEBUG: print "Initial energy guess: %f Ha." % guess
  
  # From n and l we can decide how many nodes we expect.
  nodes = n - l - 1
  if DEBUG: print "Expecting %d nodes" % nodes
  
  # Search window is from 0 down to twice the guess.
  Emin = 2.0 * guess
  Emax = 0.0
  Eupper = Emax
  Elower = Emin
  Ecentre = (Emax + Emin) / 2.0
  iteration = 0
  
  not_found = True
  if DEBUG: ion()
  
  while not_found:
    if DEBUG: print "Entering iteration %d:" % iteration
    Fupper = f(grid, l, Eupper)
    Fcentre = f(grid, l, Ecentre)
    Flower = f(grid, l, Elower)
    Xupper = numerov(Fupper, h, grid, 0, 1)
    Xcentre = numerov(Fcentre, h, grid, 0, 1)
    Xlower = numerov(Flower, h, grid, 0, 1)
    
    if DEBUG: plot(grid, Xlower, 'b-')
    if DEBUG: plot(grid, Xcentre, 'g-')
    if DEBUG: plot(grid, Xupper, 'r-')
    
    if abs(Xcentre[-1]) < REALLY_SMALL:
      # Converged to solution - check it's the correct solution.
      if DEBUG: print "Found solution at E=%f Ha..." % Ecentre
      cur_nodes = nodecount(Xcentre)
      if cur_nodes == nodes:
        if DEBUG: print "Solution has the correct number of nodes."
        if DEBUG: ioff()
        return Ecentre
      elif cur_nodes > nodes:
        # Too many nodes - set Eupper to Ecentre and restart.
        if DEBUG: print "Solution has too many nodes (%d)." % cur_nodes
        Eupper = Ecentre
        Elower = Emin
      else:
        # Not enough nodes: set Elower to Ecentre and restart.
        if DEBUG: print "Solution has too few nodes (%d)." % cur_nodes
        Elower = Ecentre
        Eupper = Emax
    elif sign(Xcentre[-1]) == sign(Xupper[-1]):
      # Solution is not between the centre and the top.
      if DEBUG: print "Solution is below the range %f to %f." % (Ecentre, Eupper)
      Eupper = Ecentre
    elif sign(Xcentre[-1]) == sign(Xlower[-1]):
      if DEBUG: print "Solution is above the range %f to %f." % (Elower, Ecentre)
      Elower = Ecentre
      
    Ecentre = (Eupper + Elower) / 2
    iteration += 1
    
  return None
    
class Atom:
  
  nuclear_charge = 1.0
  n_electrons = 1
  n_states = 2
  occupancies = array([1.0, 0.0])
  quantum_numbers = array([[1,0,0], [2,0,0]])
  num_gridpoints = 2000
  grid = []
  step = 0.01
  eigenvalues = array([-1.0, -0.1])
  wavefunctions = []
  
  def __init__(self):
    
    self.init_states()
    
  def init_states(self):
    # Task: take the number of electrons and the occupancies, find the
    # first-pass eigenvalues and workfunctions. Note this is not SCF at all.
    
    # We work on a linear grid
    self.grid = linspace(0.01,10.01,self.num_gridpoints)
    self.step = self.grid[1] - self.grid[0]
    # For each state, generate a numerov F and run an eigenvalue search
    for i in range(self.n_states):
      self.eigenvalues[i] = eigensearch(self.eigenvalues[i], self.quantum_numbers[i,0], self.quantum_numbers[i,1], self.numerov_f, self.step, self.grid)
      
    print self.eigenvalues  
   
  def numerov_f(self, r, l, E):
  
    if type(r) == type(array([])):
      return array([self.numerov_f(x, l, E) for x in r])
    else:
      return 2.0 * self.potential(r) + l * (l + 1) / (r ** 2) - 2.0 * E
      
  def potential(self, r):
    
    if type(r) == type(array([])):
      return array([self.potential(x) for x in r])
    else:
      # Nuclear-electron potential
      #return -1.0 * self.nuclear_charge / r
      # Free particle
      return 0.0
      # Harmonic Oscillator
      #return 0.5 * r ** 2
  
  
