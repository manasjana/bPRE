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
############################################################################################################





import random
from utils import *
from montgomery import *
from group_action import *
import time
import sys


def MSimS130_test():
    # Generation of Public parameters 
    PrimeList = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 569]
    r = 130
    p = pow(2, r)
    for l in PrimeList:
        p *= l
    p -= 1
    
    r, n, s = 130, len(PrimeList), 10
    A_init, C_init = 0, 1
    

     # Key generation of Receiver (Alice) 
    kgen_time=time.time()
    Alice_IntegersList = [random.randint(-s, s) for _ in range(n)]
    
    counter_kgen = [0, 0, 0]
    Alice_A, Alice_C, counter_kgen, _, _ = evaluating_the_class_group_action_montgomery(A_init, C_init, Alice_IntegersList, PrimeList, 1, 0, p, counter_kgen)
    cost_kgen = counter_kgen[0] + 0.8 * counter_kgen[1] + 0.05 * counter_kgen[2]
    print(f"KeyGen Time:{(time.time()-kgen_time)*1e3} ms")
    pk_size=(Alice_A.bit_length()+7)//8 + (Alice_C.bit_length()+7)//8
    print(f"pk size: {pk_size}")
    
    
    # Encryption by Sender (Bob)
    enc_time=time.time()
    counter_enc = [0, 0, 0]
    Bob_IntegersList = [random.randint(-s, s) for _ in range(n)]
    plaintext = 0
    for _ in range(r-2):
        plaintext = plaintext * 2 + random.randint(0, 1)
    plaintext = plaintext * 2 + 1
    
    Bob_A, Bob_C, counter_enc, _, _ = evaluating_the_class_group_action_montgomery(A_init, C_init, Bob_IntegersList, PrimeList, 1, 0, p, counter_enc)
    Bob_A2, Bob_C2, counter_enc, _, _ = evaluating_the_class_group_action_montgomery(Alice_A, Alice_C, Bob_IntegersList, PrimeList, 1, 0, p, counter_enc)
    
    P1, P2, counter_enc = generating_distinguished_point(Bob_A, Bob_C,  p, PrimeList, counter_enc)
    P11, P22, counter_enc=montgomery_ladder(P1, P2, plaintext, Bob_A+Bob_C *2, Bob_C*4, p, counter_enc)

    Bob_mont1, counter_enc= ppower(Bob_C, p-2, p, counter_enc)
    Bob_mont1=(Bob_A * Bob_mont1) %p
    counter_enc[0]+=1

    xP3, counter_enc= ppower(P22,p-2, p, counter_enc)
    xP3=(P11*xP3) % p

    randomized_xP3=int(Bob_mont1) ^ int(xP3)
    cost_enc=counter_enc[0]+0.8 * counter_enc [1] + 0.05 * counter_enc[2]
    print(f"Enc. Time: {(time.time()-enc_time)*1e3} ms")
    
    #print(f"Ciphertext size:{(Bob_mont1.bit_length()+7)//8 + (randomized_xP3.bit_length()+7)//8}")
    
    # Decryption by Receiver (Alice)
    dec_time=time.time()
    counter_dec =[0,0,0]

    update_Alice_IntegerList=[]
    for i in Alice_IntegersList:
        update_Alice_IntegerList.append(-i)
        counter_dec[2]+=1

    Alice_A, Alice_C, counter_dec, _, _ =evaluating_the_class_group_action_montgomery(Bob_A2, Bob_C2, update_Alice_IntegerList, PrimeList, 1, 0, p, counter_dec)
    P1_dec, P2_dec, counter_dec= generating_distinguished_point(Alice_A, Alice_C, p, PrimeList, counter_dec)

    Alice_mont1, counter_dec=ppower(Alice_C, p-2, p, counter_dec)
    Alice_mont1= (Alice_A*Alice_mont1) % p
    counter_dec[0]+=1

    xP3_dec=int(Alice_mont1) ^ int(randomized_xP3)
    P11_dec, P22_dec= xP3_dec % p, 1

    mm, counter_dec= pohlig_hellman(r, P11_dec, P22_dec, P1_dec, P2_dec, Alice_A, Alice_C, p, counter_dec)
    cost_dec=counter_dec[0] + 0.8 * counter_dec[1] + 0.05 * counter_dec[2]
    print(f"Dec. Time: {(time.time()-dec_time)*1e3} ms")

    # # Re-Key Generation

    # # Alice's new key generation
    # re_key_count=[0,0,0]
    # re_key_IntegerList=[]
    # Alice_new_IntegerList=[random.randint(-s, s) for _ in range(n)]
    # random_List=[random.randint(-s, s) for _ in range(n)]
    # for i in range(len(random_List)):
    #     re_key_IntegerList.append(random_List[i]+ Alice_new_IntegerList[i]-Alice_IntegersList[i])
    #     re_key_count[2]+=2
    
    

    
    # re_keygen_cost=re_key_count[0] + 0.8 * re_key_count[1] + 0.05 * re_key_count[2]


    # # Re-Encryption
    # counter_ReEnc=[0,0,0]
    # Bob_new_A2, Bob_new_C2, counter_ReEnc,_,_= evaluating_the_class_group_action_montgomery(Bob_A2,Bob_C2,re_key_IntegerList,PrimeList,1,0,p, counter_ReEnc)
    # cost_ReEnc=counter_ReEnc[0]+ 0.8 * counter_ReEnc[1]+0.05 * counter_ReEnc[2]

    if mm == plaintext or (pow(2, r) - mm) == plaintext:
        print("Decryption is Successful, we have the following costs in terms of field multiplications in Fp")
        return f"KGen. cost: {round(cost_kgen)}, Enc. cost: {round(cost_enc)}, Dec. cost: {round(cost_dec)}"
    else:
        print("Decryption is Unsuccessful")
        return 0
    


    
if __name__ == "__main__":
    
    result = MSimS130_test()
    
    if result:
        print(result)