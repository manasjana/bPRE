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



# Binary expansion of an Integer n
def bin_expansion(n):
    """
       Input: An integer n

       Output: Binary expansion of n as a list. 
    
    """
    Binary_List = []
    while n != 0:
        m = n % 2
        n = (n - m) // 2
        Binary_List.append(m)
    return Binary_List


# Jacobi Symbol 
def jacobi_symbol(a, n):
    """
      Input: An integer a and a positive odd integer n
    
      Output : The Jacobi symbol (a/n).
    """
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



# Compute D^s (mod p)
def ppower(D, s,p, counter):
    """
       Input: An integer D, a positive intger s, a prime number p
             and a initial list counter=[0,0,0]
       
       Output: Returns D^s (mod p) and countes the total number of operations
               in counter where counter[0] represents total number of field multiplication,
               counter[1] represents the total number field squaring and counter[2] represents
               total number of field addition, substraction or doubling.
    """
    Binary_List = bin_expansion(s)
    Ds = D
    for i in range(len(Binary_List) - 2, -1, -1):
        if Binary_List[i] == 0:
            Ds = pow(Ds, 2, p)
            counter[1] += 1  
        else:
            Ds = pow(Ds, 2, p)
            Ds = (Ds * D) % p
            counter[0] += 1  
            counter[1] += 1  
    return Ds, counter
