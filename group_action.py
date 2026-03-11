import random
from utils import *
from montgomery import *



# Evaluation of Class group action
def evaluating_the_class_group_action_montgomery(A, C, IntegersList, PrimeList, mm, nn, p, counter):
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
                QX, QZ, counter = montgomery_ladder_z_1(X, P_val // k, A24, C24, p, counter)

                for l in S:
                    RX, RZ, counter = montgomery_ladder(QX, QZ, kList[l], A24, C24, p, counter)
                    if RX != 0 and RZ != 0:
                        loopcount2 += 1
                        d = (PrimeList[l] - 1) // 2
                        KernelList, counter = kernel_points(RX, RZ, A24, C24, d, p, counter)
                        QX, QZ, A, C, counter = odd_isogeny_montgomery(KernelList, QX, QZ, PrimeList[l], A, C2, p, counter)
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
            QX, QZ, counter = montgomery_ladder_z_1(X, P_val // k, A24, C24, p, counter)

            for l in S:
                RX, RZ, counter = montgomery_ladder(QX, QZ, kList[l], A24, C24, p, counter)
                if RX != 0 and RZ != 0:
                    loopcount2 += 1
                    d = (PrimeList[l] - 1) // 2
                    KernelList, counter = kernel_points(RX, RZ, A24, C24, d, p, counter)
                    QX, QZ, A, C, counter = odd_isogeny_montgomery(KernelList, QX, QZ, PrimeList[l], A, C2, p, counter)
                    C2, A24, C24 = (2 * C) % p, (A + 2 * C) % p, (4 * C) % p
                    counter[2] += 3
                    SIMBA2[l] -= s
    
    return A, C, counter, loopcount1, loopcount2