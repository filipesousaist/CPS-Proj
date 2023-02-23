import numpy as np
import sympy as sp
from numpy.linalg.linalg import LinAlgError
import numpy.random as rand

j : int = 0

def keygen(q: int, *V : int):
    u: int = len(V)
    n: int = V[u - 1]
    O: tuple = ()

    rand.seed(1)

    # Generate oil variables
    for i in range(u - 1):
        O += (V[i + 1] - V[i],)

    F = ()
    for i in range(u - 1):
        Layer : tuple = ()
        
        for _ in range(O[i]):
            Alpha : np.ndarray  = rand.randint(0, q, size=(O[i], V[i]), dtype=int)
            Beta :  np.ndarray  = rand.randint(0, q, (V[i], V[i]), dtype=int)
            Gamma : np.ndarray  = rand.randint(0, q, (V[i + 1],), dtype=int)
            Eta :   int         = rand.randint(0, q)
            Layer += ((Alpha, Beta, Gamma, Eta),)

        F += (Layer,)
    
    # L1 and L2 are affine transformations (Y = A*X + B)
    L : tuple = ()
    for size in (n - V[0], n):
        while True:
            #A : np.ndarray = rand.randint(0, q, (size, size))
            A : np.ndarray = np.identity(size, dtype=int)
            try:
                print ("Before inverse")
                inv_mod(A, q)
                print("After inverse")
            except ValueError:
                print("Not invertible") 
                print(A)
                continue
            print("Invertible")
            break
        B : np.ndarray = rand.randint(0, q, (size,), dtype=int)
        L += ((A, B),)
    

    return (L[0], L[1], F)


def sign(q : int, Y : np.ndarray, s_k : tuple):
    L1, L2, F = s_k
    u = len(F) + 1

    Y_ : np.ndarray = (inv_mod(L1[0], q) @ (Y - L1[1])) % q

    while True:
        # Generate random x1, ..., x(v1)
        v : int = F[0][0][0].shape[1] 
        X : np.ndarray = rand.randint(0, q, (v,), dtype=int)
        o : int = 0
        print("X")
        success : bool = True
        for i in range(u - 1):
            o_i = len(F[i])
            try:
                X_ = plug_coefficients(q, F[i], X, Y_[o : o + o_i])
            except ValueError or sp.matrices.common.NonInvertibleMatrixError:
                print("Fail", i)
                success = False
                break
            print("Success", i)
            X = np.concatenate((X.reshape((X.shape[0], 1)), X_))
            o += o_i
        
        if success:
            return inv_mod(L2[0], q) @ (X - L2[1]) % q


def verify(q: int, F: tuple, X: np.ndarray, Y: np.ndarray):
    nLayers = len(F)
    Y2 = np.array([])
    for i in range(nLayers):
        Y2 = np.concatenate(Y2, computeF(q, F[i], X))

    for i in range(len(Y2)):
        if (Y[i] != Y2[i]):
            return False
    return True


def computeF(q, Layer, X):
    Y = np.array([0]*len(Layer))
    for i in range(len(Layer)):
        (Alpha, Beta, Gamma, Eta) = Layer[i]
        Xo = X[Beta.shape[0]:Gamma.shape[0]]
        Xv = X[:Beta.shape[0]]
        Xv2 = X[:Gamma.shape[0]]
        Y[i] = (Xo.T @ Alpha @ Xv + Xv.T @ Beta @ Xv + Gamma.dot(Xv2) + Eta) % q
    return Y

def plug_coefficients(q: int, Layer: tuple, X: np.ndarray, Y: np.ndarray):
    o : int = len(Layer)
    v : int = X.shape[0]
    A : np.ndarray = np.zeros((o, o), dtype=int) # Coefficients of the oil variables
    B : np.ndarray = np.zeros((o,), dtype=int) # Free parameters
    
    for l in range(o):
        print("plug_coefficients", l)
        Alpha, Beta, Gamma, Eta = Layer[l]
        for i in range(o):
            A[l, i] = (Alpha[i].dot(X) + Gamma[v + i]) % q

        B[l] = (X.T @ Beta @ X + Gamma[:v].dot(X) + Eta) % q
    
    return (inv_mod(A, q) @ (Y - B)) % q


def inv_mod3(M : np.ndarray, n : int):
    Adj = np.matrix(M).getH()
    det = determinantMod(M, n)
    print(det)
    r = modinv(det, n)
    return (-r * Adj) % n

def inv_mod(M : np.ndarray, n : int):
    return np.array(sp.Matrix(M).inv_mod(n))


def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, _ = egcd(a, m)
    if g != 1:
        print (a, "and", m, "are not coprime")
        raise ValueError
    else:
        return x % m

def determinantMod(M: np.ndarray, n: int):
    size = M.shape[0]
    if (size == 1):
        return M[0,0]
    total : int = 0
    mult : int = 1
    for i in range(size):  
        newM : np.ndarray = \
            np.concatenate((M[1:,:i], M[1:,i+1:]), axis = 1)
        total += (mult * M[0,i] * determinantMod(newM, n)) % n
        mult *= -1
    return total % n