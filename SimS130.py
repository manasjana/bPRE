##########################################################################################################
#                                                                                                        
# This  code was obtained by modifying that of SiGamal available online at :                             
# https://github.com/BorisFouotsa/SimS/                                                                  
#                                                                                                        
#                                                                                                        
################################################################################
# Some functions in this file are based on those in SIDH Library.
# 
# SIDH Library: https://www.microsoft.com/en-us/research/project/sidh-library/
# SIDH Library
# 
# Copyright (c) Microsoft Corporation
# All rights reserved. 
# 
# MIT License
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
# associated documentation files (the ""Software""), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
# 
# THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
#############################################################################################################




import random

# Helper for Jacobi Symbol 
def jacobi_symbol(a, n):
    """Computes the Jacobi symbol (a/n)."""
    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0

# Binary expansion
def two_expansion(n):
    """Returns the binary expansion of n as a list."""
    LIST = []
    while n != 0:
        m = n % 2
        n = (n - m) // 2
        LIST.append(m)
    return LIST

# Compute D^s
def ppower(D, s, counter, p):
    """Computes D^s mod p and updates the operation counter."""
    LIST = two_expansion(s)
    Ds = D
    for i in range(len(LIST) - 2, -1, -1):
        if LIST[i] == 0:
            Ds = pow(Ds, 2, p)
            counter[1] += 1  # Squaring (S)
        else:
            Ds = pow(Ds, 2, p)
            Ds = (Ds * D) % p
            counter[0] += 1  # Multiplication (M)
            counter[1] += 1  # Squaring (S)
    return Ds, counter

# Doubling in Montgomery curves with cost (4M + 2S + 4a)
def montgomery_double(X1, Z1, A24, C24, counter, p):
    """Doubling of a Montgomery point in projective coordinates (X:Z)


       Input: projective Montgomery x-coordinates P = (X1:Z1), 
       where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.

       Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
       
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

def montgomery_double_z_1(X, A24, C24, counter, p):
    """Optimized doubling when Z=1."""
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
    counter[0] += 4; counter[1] += 1; counter[2] += 5 # 4M + 1S + 5a
    return X2, Z2, counter

# Addition in Montgomery curves
def montgomery_add(X1, Z1, X2, Z2, X0, Z0, A, C, counter, p):
    """Performs point addition on a Montgomery curve."""
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
    counter[0] += 4; counter[1] += 2; counter[2] += 6 # 4M + 2S + 6a
    return X_plus, Z_plus, counter

def montgomery_add_double(X1, Z1, X2, Z2, X0, Z0, A24, C24, counter, p):
    """Combined Add and Double operation."""
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
    counter[0] += 8; counter[1] += 4; counter[2] += 8 # 8M + 4S + 8a
    return X4, Z4, X3, Z3, counter

def montgomery_add_double_z_1(X1, Z1, X2, Z2, X0, A24, C24, counter, p):
    """Optimized combined Add and Double when Z=1."""
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
    counter[0] += 7; counter[1] += 4; counter[2] += 8 # 7M + 4S + 8a
    return X4, Z4, X3, Z3, counter

# Montgomery ladder
def montgomery_ladder(X, Z, n, A24, C24, counter, p):
    """Scalar multiplication using the Montgomery ladder."""
    X1, Z1 = 1, 0
    X2, Z2 = X, Z
    LIST = two_expansion(n)
    for i in range(len(LIST) - 1, -1, -1):
        if LIST[i] == 0:
            X1, Z1, X2, Z2, counter = montgomery_add_double(X1, Z1, X2, Z2, X, Z, A24, C24, counter, p)
        else:
            X2, Z2, X1, Z1, counter = montgomery_add_double(X2, Z2, X1, Z1, X, Z, A24, C24, counter, p)
    return X1, Z1, counter

def montgomery_ladder_z_1(X, n, A24, C24, counter, p):
    """Optimized ladder when initial Z=1."""
    X1, Z1 = X, 1
    X2, Z2, counter = montgomery_double_z_1(X, A24, C24, counter, p)
    # IntegerToSequence(n, 2) in Magma is equivalent to bin(n) bits in reverse order
    LIST = [int(x) for x in bin(n)[2:][::-1]]
    for i in range(len(LIST) - 2, -1, -1):
        if LIST[i] == 0:
            X1, Z1, X2, Z2, counter = montgomery_add_double_z_1(X1, Z1, X2, Z2, X, A24, C24, counter, p)
        else:
            X2, Z2, X1, Z1, counter = montgomery_add_double_z_1(X2, Z2, X1, Z1, X, A24, C24, counter, p)
    return X1, Z1, counter

# Pohlig-Hellman algorithm
def pohlig_hellman(r, P1, P2, Q1, Q2, A, C, count, p):
    """Solves the discrete log problem."""
    A24 = (A + 2 * C) % p
    C24 = (4 * C) % p
    count[0] += 2; count[2] += 1
    m = 1
    P_List = [[P1, P2]]
    Q_List = [[Q1, Q2]]
    for l in range(r - 2):
        pp1, pp2, count = montgomery_double(P_List[l][0], P_List[l][1], A24, C24, count, p)
        P_List.append([pp1, pp2])
        qq1, qq2, count = montgomery_double(Q_List[l][0], Q_List[l][1], A24, C24, count, p)
        Q_List.append([qq1, qq2])

    for l in range(2, r):
        QQ1, QQ2 = Q_List[r - l - 1]
        QQ1, QQ2, count = montgomery_ladder(QQ1, QQ2, m, A24, C24, count, p)
        PP1, PP2 = P_List[r - l - 1]
        if (PP1 * QQ2) % p != (QQ1 * PP2) % p:
            m += 2**l
    return m, count

# Computing the points in the kernel
def kernel_points(X1, Z1, A24, C24, d, counter, p):
    """Computes points in the isogeny kernel."""
    KernelList = [[X1, Z1]]
    if d == 2:
        X2, Z2, counter = montgomery_double(X1, Z1, A24, C24, counter, p)
        KernelList.append([X2, Z2])
    elif d > 2:
        X2, Z2, counter = montgomery_double(X1, Z1, A24, C24, counter, p)
        KernelList.append([X2, Z2])
        x1, z1, x2, z2 = X1, Z1, X2, Z2
        for i in range(3, d + 1):
            x3, z3, counter = montgomery_add(x2, z2, X1, Z1, x1, z1, A24, C24, counter, p)
            KernelList.append([x3, z3])
            x1, z1, x2, z2 = x2, z2, x3, z3
    return KernelList, counter

# Computing an odd degree isogeny
def odd_isogeny_montgomery(KernelList, X, Z, ell, A, C2, counter, p):
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
    as_val, counter = ppower(a, s, counter, p)
    phi_a = (pow(phi_a * as_val, 2, p) * a) % p

    phi_d = pow(phi_d, 4, p)
    ds_val, counter = ppower(d_val, s, counter, p)
    phi_d = (pow(phi_d * ds_val, 2, p) * d_val) % p

    phi_A = (2 * (phi_a + phi_d)) % p
    phi_C = (phi_a - phi_d) % p

    counter[0] += (6 * s + 6)
    counter[1] += 8
    counter[2] += (4 * s + 7)
    return phi_X, phi_Z, phi_A, phi_C, counter

# Generating the distinguished point
def generating_distinguished_point(A, C, counter, p, PrimeList):
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
        P1, P2, counter = montgomery_ladder(P1, P2, l, A24, C24, counter, p)
    return P1, P2, counter

# Class group action evaluation
def evaluating_the_class_group_action_montgomery(A, C, IntegersList, PrimeList, mm, nn, counter, p):
    """Evaluates the class group action using SIMBA."""
    loopcount1 = loopcount2 = 0
    kList = [0] * len(PrimeList)
    SIMBA = [[0] * len(PrimeList) for _ in range(mm)]
    P_val = p + 1

    for i in range(len(PrimeList)):
        j = (i + 1) % mm
        if j == 0: j = mm
        SIMBA[j-1][i] = IntegersList[i]

    # SIMBA1
    for _ in range(nn):
        for nu in range(mm):
            if any(x != 0 for x in SIMBA[nu]):
                S = []
                s = 0
                X = 0
                while not S:
                    loopcount1 += 1
                    X = random.randint(0, p - 1)
                    # symb = C*X + A*X^2 + C*X^3
                    symb = (C * X + A * pow(X, 2, p) + C * pow(X, 3, p)) % p
                    counter[0] += 4; counter[1] += 1; counter[2] += 3
                    s = jacobi_symbol(symb, p)
                    if s == 1:
                        S = [t for t in range(len(SIMBA[nu])) if SIMBA[nu][t] > 0]
                    elif s == -1:
                        S = [t for t in range(len(SIMBA[nu])) if SIMBA[nu][t] < 0]

                k = 1
                for j in S:
                    kList[j] = k
                    k *= PrimeList[j]
                
                S.reverse()
                C2, A24, C24 = (2 * C) % p, (A + 2 * C) % p, (4 * C) % p
                counter[2] += 3
                QX, QZ, counter = montgomery_ladder_z_1(X, P_val // k, A24, C24, counter, p)

                for l in S:
                    RX, RZ, counter = montgomery_ladder(QX, QZ, kList[l], A24, C24, counter, p)
                    if RX != 0 and RZ != 0:
                        loopcount2 += 1
                        d = (PrimeList[l] - 1) // 2
                        KernelList, counter = kernel_points(RX, RZ, A24, C24, d, counter, p)
                        QX, QZ, A, C, counter = odd_isogeny_montgomery(KernelList, QX, QZ, PrimeList[l], A, C2, counter, p)
                        C2, A24, C24 = (2 * C) % p, (A + 2 * C) % p, (4 * C) % p
                        counter[2] += 3
                        SIMBA[nu][l] -= s

    # SIMBA2
    SIMBA2 = [0] * len(PrimeList)
    for i in range(len(PrimeList)):
        j = (i + 1) % mm
        if j == 0: j = mm
        SIMBA2[i] = SIMBA[j-1][i]

    for i in range(len(IntegersList)):
        while SIMBA2[i] != 0:
            S = []
            s = 0
            X = 0
            while not S:
                loopcount1 += 1
                X = random.randint(0, p - 1)
                symb = (C * X + A * pow(X, 2, p) + C * pow(X, 3, p)) % p
                counter[0] += 4; counter[1] += 1; counter[2] += 2
                s = jacobi_symbol(symb, p)
                if s == 1:
                    S = [t for t in range(len(SIMBA2)) if SIMBA2[t] > 0]
                elif s == -1:
                    S = [t for t in range(len(SIMBA2)) if SIMBA2[t] < 0]
            
            k = 1
            for j in S:
                kList[j] = k
                k *= PrimeList[j]
            
            S.reverse()
            C2, A24, C24 = (2 * C) % p, (A + 2 * C) % p, (4 * C) % p
            counter[2] += 3
            QX, QZ, counter = montgomery_ladder_z_1(X, P_val // k, A24, C24, counter, p)

            for l in S:
                RX, RZ, counter = montgomery_ladder(QX, QZ, kList[l], A24, C24, counter, p)
                if RX != 0 and RZ != 0:
                    loopcount2 += 1
                    d = (PrimeList[l] - 1) // 2
                    KernelList, counter = kernel_points(RX, RZ, A24, C24, d, counter, p)
                    QX, QZ, A, C, counter = odd_isogeny_montgomery(KernelList, QX, QZ, PrimeList[l], A, C2, counter, p)
                    C2, A24, C24 = (2 * C) % p, (A + 2 * C) % p, (4 * C) % p
                    counter[2] += 3
                    SIMBA2[l] -= s
    
    return A, C, counter, loopcount1, loopcount2

# Main Simulation Test
def SimS130_test():
    """Runs the SimS130 cryptographic test"""
    # Public parameters 
    PrimeList = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 569]
    p = pow(2, 130)
    for l in PrimeList:
        p *= l
    p -= 1
    
    rr, n, m_sec = 130, len(PrimeList), 10
    A_init, C_init = 0, 1

    # Alice Key Gen 
    Alice_IntegersList = [random.randint(-m_sec, m_sec) for _ in range(n)]
    counter_kgen = [0, 0, 0]
    Alice_A, Alice_C, counter_kgen, _, _ = evaluating_the_class_group_action_montgomery(A_init, C_init, Alice_IntegersList, PrimeList, 1, 0, counter_kgen, p)
    cost_kgen = counter_kgen[0] + 0.8 * counter_kgen[1] + 0.05 * counter_kgen[2]

    # Bob Encryption 
    counter_enc = [0, 0, 0]
    Bob_IntegersList = [random.randint(-m_sec, m_sec) for _ in range(n)]
    plaintext = 0
    for _ in range(128):
        plaintext = plaintext * 2 + random.randint(0, 1)
    plaintext = plaintext * 2 + 1
    
    Bob_A, Bob_C, counter_enc, _, _ = evaluating_the_class_group_action_montgomery(A_init, C_init, Bob_IntegersList, PrimeList, 1, 0, counter_enc, p)
    Bob_A2, Bob_C2, counter_enc, _, _ = evaluating_the_class_group_action_montgomery(Alice_A, Alice_C, Bob_IntegersList, PrimeList, 1, 0, counter_enc, p)
    
    P1, P2, counter_enc = generating_distinguished_point(Bob_A2, Bob_C2, counter_enc, p, PrimeList)
    P11, P22, counter_enc = montgomery_ladder(P1, P2, plaintext, Bob_A2 + Bob_C2 * 2, Bob_C2 * 4, counter_enc, p)
    
    Bob_mont2, counter_enc = ppower(Bob_C2, p - 2, counter_enc, p)
    Bob_mont2 = (Bob_A2 * Bob_mont2) % p
    counter_enc[0] += 1
    
    xP4, counter_enc = ppower(P22, p - 2, counter_enc, p)
    xP4 = (P11 * xP4) % p
    counter_enc[0] += 1
    
    randomized_xP4 = int(Bob_mont2) ^ int(xP4)
    cost_enc = counter_enc[0] + 0.8 * counter_enc[1] + 0.05 * counter_enc[2]

    # Alice Decryption 
    counter_dec = [0, 0, 0]
    Alice_A2, Alice_C2, counter_dec, _, _ = evaluating_the_class_group_action_montgomery(Bob_A, Bob_C, Alice_IntegersList, PrimeList, 1, 0, counter_dec, p)
    P1_dec, P2_dec, counter_dec = generating_distinguished_point(Alice_A2, Alice_C2, counter_dec, p, PrimeList)
    
    Alice_mont2, counter_dec = ppower(Alice_C2, p - 2, counter_dec, p)
    Alice_mont2 = (Alice_A2 * Alice_mont2) % p
    counter_dec[0] += 1
    
    xP4_dec = int(Alice_mont2) ^ int(randomized_xP4)
    P11_dec, P22_dec = xP4_dec % p, 1
    
    mm, counter_dec = pohlig_hellman(rr, P11_dec, P22_dec, P1_dec, P2_dec, Alice_A2, Alice_C2, counter_dec, p)
    cost_dec = counter_dec[0] + 0.8 * counter_dec[1] + 0.05 * counter_dec[2]

    if mm == plaintext or (pow(2, rr) - mm) == plaintext:
        print("Successful, we have the following costs in terms of field multiplications in Fp")
        return f"KGen. cost: {round(cost_kgen)}, Enc. cost: {round(cost_enc)}, Dec. cost: {round(cost_dec)}"
    else:
        print("Unsuccessful")
        return 0

if __name__ == "__main__":
    result = SimS130_test()
    if result:
        print(result)


def MSimS130_test():
    # Public parameters 
    PrimeList = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 569]
    p = pow(2, 130)
    for l in PrimeList:
        p *= l
    p -= 1
    
    rr, n, m_sec = 130, len(PrimeList), 10
    A_init, C_init = 0, 1

    # Alice Key Gen 
    Alice_IntegersList = [random.randint(-m_sec, m_sec) for _ in range(n)]
    counter_kgen = [0, 0, 0]
    Alice_A, Alice_C, counter_kgen, _, _ = evaluating_the_class_group_action_montgomery(A_init, C_init, Alice_IntegersList, PrimeList, 1, 0, counter_kgen, p)
    cost_kgen = counter_kgen[0] + 0.8 * counter_kgen[1] + 0.05 * counter_kgen[2]
    
    # Bob Encryption 
    counter_enc = [0, 0, 0]
    Bob_IntegersList = [random.randint(-m_sec, m_sec) for _ in range(n)]
    plaintext = 0
    for _ in range(128):
        plaintext = plaintext * 2 + random.randint(0, 1)
    plaintext = plaintext * 2 + 1
    
    Bob_A, Bob_C, counter_enc, _, _ = evaluating_the_class_group_action_montgomery(A_init, C_init, Bob_IntegersList, PrimeList, 1, 0, counter_enc, p)
    Bob_A2, Bob_C2, counter_enc, _, _ = evaluating_the_class_group_action_montgomery(Alice_A, Alice_C, Bob_IntegersList, PrimeList, 1, 0, counter_enc, p)
    
    P1, P2, counter_enc = generating_distinguished_point(Bob_A, Bob_C, counter_enc, p, PrimeList)
    P11, P22, counter_enc=montgomery_ladder(P1, P2, plaintext, Bob_A+Bob_C *2, Bob_C*4, counter_enc, p)

    Bob_mont1, counter_enc= ppower(Bob_C, p-2, counter_enc, p)
    Bob_mont1=(Bob_A * Bob_mont1) %p
    counter_enc[0]+=1

    xP3, counter_enc= ppower(P22,p-2, counter_enc, p)
    xP3=(P11*xP3) % p

    randomized_xP3=int(Bob_mont1) ^ int(xP3)
    cost_enc=counter_enc[0]+0.8 * counter_enc [1] + 0.05 * counter_enc[2]

    # Alice Decryption
    counter_dec =[0,0,0]

    update_Alice_IntegerList=[]
    for i in Alice_IntegersList:
        update_Alice_IntegerList.append(-i)
        counter_dec[2]+=1

    Alice_A, Alice_C, counter_dec, _, _ =evaluating_the_class_group_action_montgomery(Bob_A2, Bob_C2, update_Alice_IntegerList, PrimeList, 1, 0, counter_dec, p)
    P1_dec, P2_dec, counter_dec= generating_distinguished_point(Alice_A, Alice_C,counter_dec, p, PrimeList)

    Alice_mont1, counter_dec=ppower(Alice_C, p-2, counter_dec, p)
    Alice_mont1= (Alice_A*Alice_mont1) % p
    counter_dec[0]+=1

    xP3_dec=int(Alice_mont1) ^ int(randomized_xP3)
    P11_dec, P22_dec= xP3_dec % p, 1

    mm, counter_dec= pohlig_hellman(rr, P11_dec, P22_dec, P1_dec, P2_dec, Alice_A, Alice_C, counter_dec, p)
    cost_dec=counter_dec[0] + 0.8 * counter_dec[1] + 0.05 * counter_dec[2]


    # Re-Key Generation

    # Alice's new key generation
    re_key_count=[0,0,0]
    re_key_IntegerList=[]
    Alice_new_IntegerList=[random.randint(-m_sec, m_sec) for _ in range(n)]
    random_List=[random.randint(-m_sec, m_sec) for _ in range(n)]
    for i in range(len(random_List)):
        re_key_IntegerList.append(random_List[i]+ Alice_new_IntegerList[i]-Alice_IntegersList[i])
        re_key_count[2]+=2
    
    print(f"{re_key_IntegerList}")

    
    re_keygen_cost=re_key_count[0] + 0.8 * re_key_count[1] + 0.05 * re_key_count[2]


    # Re-Encryption
    counter_ReEnc=[0,0,0]
    Bob_new_A2, Bob_new_C2, counter_ReEnc,_,_= evaluating_the_class_group_action_montgomery(Bob_A2,Bob_C2,re_key_IntegerList,PrimeList,1,0,counter_ReEnc,p)
    cost_ReEnc=counter_ReEnc[0]+ 0.8 * counter_ReEnc[1]+0.05 * counter_ReEnc[2]

    if mm == plaintext or (pow(2, rr) - mm) == plaintext:
        print("Successful, we have the following costs in terms of field multiplications in Fp")
        return f"KGen. cost: {round(cost_kgen)}, Enc. cost: {round(cost_enc)}, Dec. cost: {round(cost_dec)}, reKeyGen. cost:{round(re_keygen_cost)}, ReEnc. cost:{round(cost_ReEnc)}"
    else:
        print("Unsuccessful")
        return 0
    


    
if __name__ == "__main__":
    result = MSimS130_test()
    if result:
        print(result)