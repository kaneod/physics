
import numpy as np

tol = 1e-6

def preedge_calculate(x,y):
  """ P = pes_utils.preedge_calculate(x,y)
  
  Calculates the best-fit linear pre-edge for a dataset (x,y). Finds the biggest peak,
  then finds the pre-edge region using a sequence of linear fits starting from the end
  point.
  
  """

  # First ensure the energy values are *decreasing* in the array,
  # if not, reverse them.
  
  if x[0] < x[-1]:
    print "Found x,y arrays badly ordered: Reversing order!"
    E = x[::-1]
    J = y[::-1]
    
  else:
    E = x
    J = y
    
  # Locate the biggest peak.
  maxidx = np.abs(y - np.amax(y)).argmin()
  
  # Find the gradient of every possible linear fit between the lowest binding energy
  # and the biggest peak.
  grads = []
  for i in range(2, len(x) - maxidx):
    
    # Best linear fit to the last i values
    xs = x[-i:]
    ys = y[-i:]
    p = np.polyfit(xs,ys,1)
    grads.append(p[0])
    
  # Differentiate the gradient array.
  dgrads = []
  for i in range(len(grads)-1):
    dgrads.append(grads[i+1] - grads[i])
  
  dgrads = np.array(dgrads)
  
  # Find the minimum index of the absolute of the gradient of gradients.
  mingrad = np.abs(dgrads).argmin()
  
  # Make a best linear fit from this number of pre-edge points, generate linear
  # pre-edge.
  p = np.polyfit(x[-mingrad:], y[-mingrad:], 1)
  return np.polyval(p,x)
def shirley_calculate(x,y):
  """ S = pes_utils.shirley_calculate(x,y)
  
  Calculate the best auto-Shirley background S for a dataset (x,y). Finds the biggest peak 
  and then uses the minimum value either side of this peak as the terminal points of the
  Shirley background.
  """
  
  # First ensure the energy values are *decreasing* in the array,
  # if not, reverse them.
  
  if x[0] < x[-1]:
    print "Found x,y arrays badly ordered: Reversing order!"
    E = x[::-1]
    J = y[::-1]
    
  else:
    E = x
    J = y
    
  # Locate the biggest peak.
  maxidx = np.abs(y - np.amax(y)).argmin()
  
  # Locate the minima either side of maxidx.
  lmidx = np.abs(y[0:maxidx] - np.amin(y[0:maxidx])).argmin()
  rmidx = np.abs(y[maxidx:] - np.amin(y[maxidx:])).argmin() + maxidx
  # Get Eleft, Eright, Ileft, Iright
  Eleft = E[lmidx]
  Ileft = J[lmidx]
  Eright = E[rmidx]
  Iright = J[rmidx]
  
  # Max integration index
  imax = rmidx - 1

  # Initial value of B. The total background S = Iright + B,
  # and B is equal to (Ileft - Iright) below lmidx and initially
  # zero above.
  B = np.zeros(E.shape)
  B[:lmidx] = Ileft - Iright
  Bnew = B.copy()
  
  while True:
    
    # Calculate new kn = (Ileft - Iright) / (int_(Eleft)^(Eright) J(E') - Iright - Bn-1(E') dE')
    ksum = 0.0
    for i in range(lmidx,imax):
      ksum += (E[i] - E[i+1]) * 0.5 * (J[i] + J[i+1] - 2 * Iright - B[i] - B[i+1])
    k = (Ileft - Iright) / ksum
    
    # Calculate new B
    for i in range(lmidx,rmidx):
      Isum = 0.0
      for j in range(i,imax):
        Isum += (E[j] - E[j+1]) * 0.5 * (J[j] + J[j+1] - 2 * Iright - B[j] - B[j+1])
      Bnew[i] = k * Isum
      
      
    # If Bnew is close to B, exit.
    if np.linalg.norm(Bnew - B) < tol:
      B = Bnew.copy()
      break
    else:
      B = Bnew.copy()
      
  return Iright + B

    
