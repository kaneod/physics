
from __future__ import division
from numpy import array, linspace, zeros, sign, max, polyfit, roots, exp,arange
from pylab import plot, ion, semilogy, figure

DEBUG = 1
REALLY_SMALL = 1.0e-10
REALLY_BIG  = 1.0e6


def numerov(F, grid, x0, x1, direction="forward"):
  
  x = zeros(grid.shape)
  x[0] = x0
  x[1] = x1
  
  # Create repeated items
  A = (2.0 + (5.0 / 6.0) * F)
  B = (1.0 - F / 12.0)
  C = B ** -1.0
  if direction == "backward":
    F = F[::-1]
  
  for i in range(len(grid[2:])):   
    x[i+2] = (x[i+1] * A[i+1] - x[i] * B[i]) * C[i+2]
    #if DEBUG: print "%10f" % x[i+2]
    if abs(x[i+2]) > REALLY_BIG:
      return x,x[i+2], i+2
  
  if direction == "forward":
    return x, x[-1], len(x) - 1
  else:
    return x[::-1], x[-1], len(x) - 1

def numerov_f(l, r, E):
  
  return 2.0 * V(r) + l * (l + 1) / (r ** 2) - 2.0 * E
  
def V(r):
  
  return -1.0 / r

rmax = 100
jmax = 10000
delta = 0.001
rp = rmax / (exp(jmax * delta) - 1)
grid = array([rp * (exp(j * delta) - 1) for j in range(jmax+1)])
grid = grid[1:]

if DEBUG: print len(grid)

l = 0
Elower = -0.75
Eupper = -0.25
Ecentre = (Elower + Eupper) / 2.0

iteration = 0
not_done = True
while not_done:
  Fl = array([((rp ** 2) * (delta ** 2) * exp(2 * j * delta)) * numerov_f(l,x,Elower) + 0.25 * (delta ** 2) for j, x in enumerate(grid)])

  xl, xlm, ilm = numerov(Fl, grid, 0.0, 0.0001, "forward")
  
  Fu = array([((rp ** 2) * (delta ** 2) * exp(2 * j * delta)) * numerov_f(l,x,Eupper) + 0.25 * (delta ** 2) for j, x in enumerate(grid)])
  xu, xum, ium = numerov(Fu, grid, 0.0, 0.0001, "forward")
  
  Fc = array([((rp ** 2) * (delta ** 2) * exp(2 * j * delta)) * numerov_f(l,x,Ecentre) + 0.25 * (delta ** 2) for j, x in enumerate(grid)])
  xc, xcm, icm = numerov(Fc, grid, 0.0, 0.0001, "forward")
  
  
  # Convert x (really v from page 116) to u.
  if DEBUG: print len(xl), len(arange(jmax+1))
  ul = xl * exp(arange(jmax+1)[1:] * delta)
  uu = xu * exp(arange(jmax+1)[1:] * delta)
  uc = xc * exp(arange(jmax+1)[1:] * delta)
  ucm = xcm * exp(icm * delta)
  uc = uc / max(abs(uc))
  
  if abs(ucm) < REALLY_SMALL:
    print "Finished with energy %f, residual %f." % (Ecentre, ucm)
    break
  elif iteration > 50:
    print "Exceeded 200 iterations."
    break
  elif sign(xcm) == sign(xlm):
    Elower = Ecentre
  elif sign(xcm) == sign(xum):
    Eupper = Ecentre
  print "Iteration %d: Energy %f, residual %f." %(iteration, Ecentre, ucm)
  
  
  Ecentre = (Elower + Eupper) / 2.0 
  iteration += 1

ion()
plot(grid, uc)
figure()
plot(grid, xc)
figure()
plot(arange(len(grid)), grid)
#semilogy()
