


from utils import *
#Operations of Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) B'*y^2=C*x^3+A*x^2+C*x.

# Doubling in Montgomery curve with cost (4M + 2S + 4a)
def montgomery_double(X1, Z1, A24, C24, p, counter):
    """Doubling of a Montgomery point in projective coordinates (X:Z)


       Input: projective Montgomery x-coordinates P = (X1:Z1), 
              Montgomery curve constants A24=A+2C and C24=4C, 
              a prime number p and an initial counter

       Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2) 
               and updates the total number of operations in counter.
       
       """

    t0 = (X1 - Z1) % p
    t1 = (X1 + Z1) % p
    t0 = pow(t0, 2, p)
    t1 = pow(t1, 2, p)
    Z2 = (C24 * t0) % p
    X2 = (Z2 * t1) % p
    t1 = (t1 - t0) % p
    t0 = (A24 * t1) % p
    Z2 = (Z2 + t0) % p
    Z2 = (Z2 * t1) % p
    counter[0] += 4; counter[1] += 2; counter[2] += 4  
    return X2, Z2, counter


# Addition in Montgomery curve with cost (4M + 2S + 6a)
def montgomery_add(X1, Z1, X2, Z2, X0, Z0, A, C, p, counter):
    """ Addition of two Montgomery points in projective coordinates (X:Z)


       Input: projective Montgomery x-coordinates x(P) = (X1:Z1),
              x(Q)=(X2:Z2), x(P-Q)=(X0:Z0), Montgomery curve constants
              A and C, a prime number p and an initial counter

       Output: projective Montgomery x-coordinates of x(P+Q)=(X_plus:Z_plus) 
               and updates the total number of operations in counter.
       
       """
    V0 = (X1 + Z1) % p
    V1 = (X2 - Z2) % p
    V1 = (V1 * V0) % p
    V0 = (X1 - Z1) % p
    V2 = (X2 + Z2) % p
    V2 = (V2 * V0) % p
    V3 = (V1 + V2) % p
    V3 = pow(V3, 2, p)
    V4 = (V1 - V2) % p
    V4 = pow(V4, 2, p)
    X_plus = (Z0 * V3) % p
    Z_plus = (X0 * V4) % p
    counter[0] += 4; counter[1] += 2; counter[2] += 6 
    return X_plus, Z_plus, counter

# Optimized doubling when Z=1 with cost (4M + 1S + 5a)
def montgomery_double_z_1(X, A24, C24, p, counter):
    """
        Input: projective Montgomery x-coordinates P = (X1:1), 
              Montgomery curve constants A24=A+2C and C24=4C, 
              a prime number p and an initial counter
        
            Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2) 
                    and updates the total number of operations in counter.
    """
    t0 = (X - 1) % p
    X4 = (2 * X) % p
    X4 = (2 * X4) % p
    t0 = pow(t0, 2, p)
    t1 = (t0 + X4) % p
    Z2 = (C24 * t0) % p
    X2 = (Z2 * t1) % p
    t0 = (A24 * X4) % p
    Z2 = (Z2 + t0) % p
    Z2 = (Z2 * X4) % p
    counter[0] += 4; counter[1] += 1; counter[2] += 5 
    return X2, Z2, counter

#  Combined Add and Double operation with cost (8M + 4S + 8a)
def montgomery_add_double(X1, Z1, X2, Z2, X0, Z0, A24, C24, p, counter):
    """Simultaneous doubling and addition of two Montgomery points 
       in projective coordinates (X:Z)


       Input: projective Montgomery x-coordinates x(P) = (X1:Z1),
              x(Q)=(X2:Z2), x(P-Q)=(X0:Z0), Montgomery curve constants
              A24=A+2C and C24=4C, a prime number p and an initial counter

        Output: projective Montgomery x-coordinates of x(P+Q)=(X3:Z3),  
               x(2P)=(X4:Z4) and updates the total number of operations in counter.
              
    """
    V0 = (X1 + Z1) % p
    t1 = pow(V0, 2, p)
    V1 = (X2 - Z2) % p
    V1 = (V1 * V0) % p
    V0 = (X1 - Z1) % p
    t0 = pow(V0, 2, p)
    V2 = (X2 + Z2) % p
    V2 = (V2 * V0) % p
    V3 = (V1 + V2) % p
    V3 = pow(V3, 2, p)
    V4 = (V1 - V2) % p
    V4 = pow(V4, 2, p)
    X3 = (Z0 * V3) % p
    Z3 = (X0 * V4) % p
    Z4 = (C24 * t0) % p
    X4 = (Z4 * t1) % p
    t1 = (t1 - t0) % p
    t0 = (A24 * t1) % p
    Z4 = (Z4 + t0) % p
    Z4 = (Z4 * t1) % p
    counter[0] += 8; counter[1] += 4; counter[2] += 8 
    return X4, Z4, X3, Z3, counter

# Optimized combined Add and Double when Z0=1 with cost (7M + 4S + 8a)
def montgomery_add_double_z_1(X1, Z1, X2, Z2, X0, A24, C24, p, counter):
    """optimized simultaneous doubling and addition of two Montgomery points 
       in projective coordinates (X:Z)


       Input: projective Montgomery x-coordinates x(P) = (X1:Z1),
              x(Q)=(X2:Z2), x(P-Q)=(X0:1), Montgomery curve constants
              A24=A+2C and C24=4C, a prime number p and an initial counter

        Output: projective Montgomery x-coordinates of x(P+Q)=(X3:Z3),  
               x(2P)=(X4:Z4) and updates the total number of operations in counter.
              
    """
    V0 = (X1 + Z1) % p
    t1 = pow(V0, 2, p)
    V1 = (X2 - Z2) % p
    V1 = (V1 * V0) % p
    V0 = (X1 - Z1) % p
    t0 = pow(V0, 2, p)
    V2 = (X2 + Z2) % p
    V2 = (V2 * V0) % p
    V3 = (V1 + V2) % p
    V3 = pow(V3, 2, p)
    V4 = (V1 - V2) % p
    V4 = pow(V4, 2, p)
    X3 = V3
    Z3 = (X0 * V4) % p
    Z4 = (C24 * t0) % p
    X4 = (Z4 * t1) % p
    t1 = (t1 - t0) % p
    t0 = (A24 * t1) % p
    Z4 = (Z4 + t0) % p
    Z4 = (Z4 * t1) % p
    counter[0] += 7; counter[1] += 4; counter[2] += 8 
    return X4, Z4, X3, Z3, counter


# Scalar multiplication using the Montgomery ladder
def montgomery_ladder(X, Z, n, A24, C24, p, counter):
    """Input: projective Montgomery x-coordinates P = (X:Z), 
              an integer n, Montgomery curve constants A24=A+2C and C24=4C, 
              a prime number p and an initial counter

       Output: projective Montgomery x-coordinates Q = n*P = (X1:Z1) 
               and updates the total number of operations in counter.
       """
    X1, Z1 = 1, 0
    X2, Z2 = X, Z
    LIST = bin_expansion(n)
    for i in range(len(LIST) - 1, -1, -1):
        if LIST[i] == 0:
            X1, Z1, X2, Z2, counter = montgomery_add_double(X1, Z1, X2, Z2, X, Z, A24, C24, p, counter)
        else:
            X2, Z2, X1, Z1, counter = montgomery_add_double(X2, Z2, X1, Z1, X, Z, A24, C24,  p, counter)
    return X1, Z1, counter

# Optimized calar multiplication using the Montgomery ladder when initial Z=1.
def montgomery_ladder_z_1(X, n, A24, C24, p, counter):
    """Input: projective Montgomery x-coordinates P = (X:1), 
              an integer n, Montgomery curve constants A24=A+2C and C24=4C, 
              a prime number p and an initial counter
              
        Output: projective Montgomery x-coordinates Q = n*P = (X1:Z1) 
               and updates the total number of operations in counter.
       """
    X1, Z1 = X, 1
    X2, Z2, counter = montgomery_double_z_1(X, A24, C24,  p, counter)

    LIST = [int(x) for x in bin(n)[2:][::-1]]
    for i in range(len(LIST) - 2, -1, -1):
        if LIST[i] == 0:
            X1, Z1, X2, Z2, counter = montgomery_add_double_z_1(X1, Z1, X2, Z2, X, A24, C24,  p, counter)
        else:
            X2, Z2, X1, Z1, counter = montgomery_add_double_z_1(X2, Z2, X1, Z1, X, A24, C24,  p, counter)
    return X1, Z1, counter


# Pohlig-Hellman algorithm
def pohlig_hellman(r, P1, P2, Q1, Q2, A, C, p, counter):
    """Solves the discrete log problem."""
    A24 = (A + 2 * C) % p
    C24 = (4 * C) % p
    counter[0] += 2; counter[2] += 1
    m = 1
    P_List = [[P1, P2]]
    Q_List = [[Q1, Q2]]
    for l in range(r - 2):
        pp1, pp2, counter = montgomery_double(P_List[l][0], P_List[l][1], A24, C24, p, counter)
        P_List.append([pp1, pp2])
        qq1, qq2, counter = montgomery_double(Q_List[l][0], Q_List[l][1], A24, C24, p, counter)
        Q_List.append([qq1, qq2])

    for l in range(2, r):
        QQ1, QQ2 = Q_List[r - l - 1]
        QQ1, QQ2, counter = montgomery_ladder(QQ1, QQ2, m, A24, C24, p, counter)
        PP1, PP2 = P_List[r - l - 1]
        if (PP1 * QQ2) % p != (QQ1 * PP2) % p:
            m += 2**l
    return m, counter

# Computing the points in the kernel
def kernel_points(X1, Z1, A24, C24, d, p, counter):
    """Computes points in the isogeny kernel."""
    KernelList = [[X1, Z1]]
    if d == 2:
        X2, Z2, counter = montgomery_double(X1, Z1, A24, C24, p, counter)
        KernelList.append([X2, Z2])
    elif d > 2:
        X2, Z2, counter = montgomery_double(X1, Z1, A24, C24,  p, counter)
        KernelList.append([X2, Z2])
        x1, z1, x2, z2 = X1, Z1, X2, Z2
        for i in range(3, d + 1):
            x3, z3, counter = montgomery_add(x2, z2, X1, Z1, x1, z1, A24, C24, p, counter)
            KernelList.append([x3, z3])
            x1, z1, x2, z2 = x2, z2, x3, z3
    return KernelList, counter

# Computing an odd degree isogeny
def odd_isogeny_montgomery(KernelList, X, Z, ell, A, C2, p, counter):
    """Computes an odd degree isogeny."""
    phi_X, phi_Z, phi_a, phi_d = 1, 1, 1, 1
    t1, t2 = (X + Z) % p, (X - Z) % p
    a, d_val = (A + C2) % p, (A - C2) % p
    s = len(KernelList)
    for i in range(s):
        t3 = (KernelList[i][1] + KernelList[i][0]) % p
        t4 = (KernelList[i][0] - KernelList[i][1]) % p
        F = (t2 * t3) % p
        G = (t1 * t4) % p
        phi_X = (phi_X * (F + G)) % p
        phi_Z = (phi_Z * (F - G)) % p
        phi_a = (phi_a * t3) % p
        phi_d = (phi_d * t4) % p

    phi_X = (pow(phi_X, 2, p) * X) % p
    phi_Z = (pow(phi_Z, 2, p) * Z) % p

    phi_a = pow(phi_a, 4, p)
    as_val, counter = ppower(a, s, p, counter)
    phi_a = (pow(phi_a * as_val, 2, p) * a) % p

    phi_d = pow(phi_d, 4, p)
    ds_val, counter = ppower(d_val, s, p, counter)
    phi_d = (pow(phi_d * ds_val, 2, p) * d_val) % p

    phi_A = (2 * (phi_a + phi_d)) % p
    phi_C = (phi_a - phi_d) % p

    counter[0] += (6 * s + 6)
    counter[1] += 8
    counter[2] += (4 * s + 7)
    return phi_X, phi_Z, phi_A, phi_C, counter

# Generating the distinguished point
def generating_distinguished_point(A, C, p, PrimeList, counter):
    """Generates a distinguished point on the curve."""
    A24 = (A + 2 * C) % p
    C24 = (4 * C) % p
    counter[2] += 4
    i = 1
    while True:
        i += 1
        P1 = (-i) % p
        # R = P1^3 * C^2 + P1^2 * A * C + P1 * C^2
        R = (pow(P1, 3, p) * pow(C, 2, p) + pow(P1, 2, p) * A * C + P1 * pow(C, 2, p)) % p
        counter[0] += 5; counter[1] += 4; counter[2] += 4
        if jacobi_symbol(R, p) == 1:
            break
    P2 = 1
    for l in PrimeList:
        P1, P2, counter = montgomery_ladder(P1, P2, l, A24, C24, p, counter)
    return P1, P2, counter