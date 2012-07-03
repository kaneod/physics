
import numpy as np

tol = 1.0e-6
HUGE = 1e10

def shirley_integral(E,F,eidx):
  
  # I(E) = int_{E}^{Emax} F(E') dE'
  
  # Here "Emax" is the smallest binding energy,
  # so we want the energy scale to be *decreasing*. 
  
  Isum = 0.0
  
  for i in range(eidx,len(E)-1):
    Isum += F[i]
    

def shirley_calculate(x,y,emin):

  # First ensure the energy values are *decreasing* in the array,
  # if not, reverse them.
  
  if x[0] < x[-1]:
    E = x[::-1]
    J = y[::-1]
  else:
    E = x
    J = y
    
  # Figure out the index closest to emin
  emidx = 0
  edist = HUGE
  for i,energy in enumerate(E):
    if abs(energy - emin) < edist:
      edist = abs(energy - emin)
      emidx = i
  
  # Get Eleft, Eright, Ileft, Iright
  Eleft = E[emidx]
  Ileft = J[emidx]
  Eright = E[-1]
  Iright = J[-1]
  
  # Max integration index
  imax = len(E) - 1

  # Initial value of B. The total background S = Iright + B,
  # and B is equal to Ileft - Iright below emin and initially
  # zero above.
  B = zeros(E.shape)
  B[:emidx] = Ileft - Iright
  Bnew = B
  
  Bs = []
  
  while True:
    
    # Calculate new kn = (Ileft - Iright) / (int_(Eleft)^(Eright) J(E') - Iright - Bn-1(E') dE')
    ksum = 0.0
    for i in range(emidx,imax):
      ksum += (E[i] - E[i+1]) * 0.5 * (J[i] + J[i+1] - 2 * Iright - B[i] - B[i+1])
    k = (Ileft - Iright) / ksum
    
    print "Found new k: ", k
    
    # Calculate new B
    for i in range(emidx,len(E)):
      Isum = 0.0
      for j in range(i,imax):
        Isum += (E[i] - E[i+1]) * 0.5 * (J[i] + J[i+1] - 2 * Iright - B[i] - B[i+1])
      print "Calculated isum for E: ", E[i], Isum
      Bnew[i] = k * Isum
      
      
    # If Bnew is close to B, exit.
    print "Bnew difference: ", np.linalg.norm(Bnew - B)
    if np.linalg.norm(Bnew - B) < tol:
      Bs.append(Bnew)
      B = Bnew
      break
    else:
      Bs.append(Bnew)
      B = Bnew
      
  return Iright, Bs

    
