import numpy as np
import time
from scipy import sparse as sp
import scipy.sparse.linalg as spla

theta = np.array([ -10.843917078696988026 + 19.277446167181652284j,
                   -5.2649713434426468895 + 16.220221473167927305j,
                   5.9481522689511774808 + 3.5874573620183222829j,
                   3.5091036084149180974 + 8.4361989858843750826j,
                   6.4161776990994341923 + 1.1941223933701386874j,
                   1.4193758971856659786 + 10.925363484496722585j,
                   4.9931747377179963991 + 5.9968817136039422260j,
                   -1.4139284624888862114 + 13.497725698892745389j])

alpha = np.array([ -.0000005090152186522491565 - .00002422001765285228797j,
                      .00021151742182466030907 + .0043892969647380673918j,
                      113.39775178483930527 + 101.9472170421585645j,
                      15.059585270023467528 - 5.7514052776421819979j,
                      -64.500878025539646595 - 224.59440762652096056j,
                      -1.4793007113557999718 + 1.7686588323782937906j,
                      -62.518392463207918892 - 11.19039109428322848j,
                      .041023136835410021273 - .15743466173455468191j])

alpha_0 = np.array([2.1248537104952237488e-16])

def solveSystem(A, t, n_0):
    s = len(theta)
    A = A*t
    n = 0*n_0
    #A = sp.csc_matrix(A)
    ident = sp.identity(np.shape(A)[0],format="csc")
    #ident = np.identity(np.shape(A)[0])

    for j in xrange(s):
        n = n + spla.spsolve(A - theta[j]*ident, alpha[j]*n_0)
        #n = n + np.linalg.solve(A - theta[j]*ident, alpha[j]*n_0)

    n = 2.*n.real
    n = n + alpha_0*n_0
    return n

if __name__ == "__main__":
    for x in xrange(1):
        n = 1000000
        A = sp.random(n,n,density=0.01, format="csc")
        #A = sp.random(n,n)
        #A = np.random.rand(n,n)
        t = 0.1
        n0 = np.ones((n,1))
        
        start = time.time()
        sol = solveSystem(A, t, n0)
        end = time.time()
        
        runTime = end - start
        print runTime
