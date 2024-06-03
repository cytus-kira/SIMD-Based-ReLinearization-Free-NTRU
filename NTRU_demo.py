import time

import numpy as np

from BFV import *
from Ntru import Ntru
from ckks_encode import CKKSEncoder
from helper import *

from random import randint
from math import log, ceil

# This implementation follows the description at https://eprint.iacr.org/2012/144.pdf
# Brakerski/Fan-Vercauteren (BFV) somewhat homomorphic encryption scheme
#
# Polynomial arithmetic on ciphertext domain is performed in Z[x]_q/x^n+1
# Polynomial arithmetic on plaintext domain is performed in Z[x]_t/x^n+1
# * n: ring size
# * q: ciphertext coefficient modulus
# * t: plaintext coefficient modulus (if t is equal to 2, no negative values is accepted)
# * psi,psiv,w,wv: polynomial arithmetic parameters
#
# Note that n,q,t parameters together determine the multiplicative depth.

# Parameter generation (pre-defined or generate parameters)
PD = 1  # 0: generate -- 1: pre-defined

if PD == 0:
    # Select one of the parameter sets below
    t = 16
    n, q, psi = 1024, 132120577, 73993  # log(q) = 27
    # t = 256;  n, q, psi = 2048 , 137438691329      , 22157790             # log(q) = 37
    # t = 1024; n, q, psi = 4096 , 288230376135196673, 60193018759093       # log(q) = 58

    # other necessary parameters
    psiv = modinv(psi, q)
    w = pow(psi, 2, q)
    wv = modinv(w, q)
else:
    # Enter proper parameters below
    # t, n, logq = 16, 1024, 27
    t, n, logq = 16, 2048, 37
    # t, n, logq = 1024, 4096, 58

    # other necessary parameters (based on n and log(q) determine other parameter)
    q, psi, psiv, w, wv = ParamGen(n, logq)

# Determine mu, sigma (for discrete gaussian distribution)

mu = 0
sigma = 0.5 * 3.2

# Generate polynomial arithmetic tables
w_table = [1] * n
wv_table = [1] * n
psi_table = [1] * n
psiv_table = [1] * n
for i in range(1, n):
    w_table[i] = ((w_table[i - 1] * w) % q)
    wv_table[i] = ((wv_table[i - 1] * wv) % q)
    psi_table[i] = ((psi_table[i - 1] * psi) % q)
    psiv_table[i] = ((psiv_table[i - 1] * psiv) % q)

qnp = [w_table, wv_table, psi_table, psiv_table]

print("--- Starting BFV Demo")

# Generate BFV evaluator
Evaluator = Ntru(n, q, t, mu, sigma, qnp)

# Generate Keys
start_time = time.time()
Evaluator.SecretKeyGen(True)
Evaluator.PublicKeyGen(True)
#
start_enc_time = time.time()
print(start_enc_time - start_time)

# print system parameters
print(Evaluator)

# Generate random message
# n1, n2 = 15, -5
# n1, n2 = randint(-(2 ** 15), 2 ** 15 - 1), randint(-(2 ** 15), 2 ** 15 - 1)

# print("--- Random integers n1 and n2 are generated.")
# print("* n1: {}".format(n1))
# print("* n2: {}".format(n2))
# print("* n1+n2: {}".format(n1 + n2))
# print("* n1-n2: {}".format(n1 - n2))
# print("* n1*n2: {}".format(n1 * n2))
# print("")

# Encode random messages into plaintext polynomials
# print("--- n1 and n2 are encoded as polynomials m1(x) and m2(x).")
# m1 = Evaluator.IntEncode(n1)
# m2 = Evaluator.IntEncode(n2)
n1 = [1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1,
      3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7,
      1, 3, 4, 2, 1, 1, 3, 4, 5, 1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 2, 1, 1, 3, 4, 5, 1]
# n2 = np.random.rand(n)
print(n1)
# print(n2)

# encode = CKKSEncoder(len(n1) * 4, 2 ** 16)

# e1 = encode.encode(n1)
# print(e1)
# d_1 = encode.decode(e1)
# print(d_1.tolist())
# m1 = encode.encode(n1).coef
m1 = n1
# m2 = encode.encode(n2).coef
# print("* m1(x): {}".format(m1))
# print("* m2(x): {}".format(m2))
# print("")

# Encrypt message

ct1 = Evaluator.Encryption(m1)
print(ct1)
# ct2 = Evaluator.Encryption(m2)
ct1_2 = ct1 + ct1
start_dec_time = time.time()
print('enc_time', start_dec_time - start_enc_time)
dec_1 = Evaluator.Decryption(ct1_2)
print(dec_1.F)
end_dec_time = time.time()
print('enc_time', end_dec_time - start_dec_time)
# dec_2 = Evaluator.Decryption(ct2)
# dec_1 = np.poly1d(dec_1.F)
# dec_2 = np.poly1d(dec_2.F)
# decode_1 = encode.decode(dec_1)
# decode_2 = encode.decode(dec_2)
#
# data1 = Evaluator.IntDecode(dec1)
# data2 = Evaluator.IntDecode(dec2)
# print(decode_1)
# print(decode_2)
#
# print("--- m1 and m2 are encrypted as ct1 and ct2.")
# print("* ct1[0]: {}".format(ct1[0]))
# print("* ct1[1]: {}".format(ct1[1]))
# print("* ct2[0]: {}".format(ct2[0]))
# print("* ct2[1]: {}".format(ct2[1]))
# print("")
#
# # Homomorphic Addition
# ct = Evaluator.HomomorphicAddition(ct1,ct2)
# mt = Evaluator.Decryption(ct)
#
# nr = Evaluator.IntDecode(mt)
# ne = (n1+n2)
#
# print("--- Performing ct_add = Enc(m1) + Enc(m2)")
# print("* ct_add[0] :{}".format(ct[0]))
# print("* ct_add[1] :{}".format(ct[1]))
# print("--- Performing ct_dec = Dec(ct_add)")
# print("* ct_dec    :{}".format(mt))
# print("--- Performing ct_dcd = Decode(ct_dec)")
# print("* ct_dcd    :{}".format(nr))
#
# if nr == ne:
#     print("* Homomorphic addition works.")
# else:
#     print("* Homomorphic addition does not work.")
# print("")
#
# # Homomorphic Subtraction
# ct = Evaluator.HomomorphicSubtraction(ct1,ct2)
# mt = Evaluator.Decryption(ct)
#
# nr = Evaluator.IntDecode(mt)
# ne = (n1-n2)
#
# print("--- Performing ct_sub = Enc(m1) - Enc(m2)")
# print("* ct_sub[0] :{}".format(ct[0]))
# print("* ct_sub[1] :{}".format(ct[1]))
# print("--- Performing ct_dec = Dec(ct_sub)")
# print("* ct_dec    :{}".format(mt))
# print("--- Performing ct_dcd = Decode(ct_dec)")
# print("* ct_dcd    :{}".format(nr))
#
# if nr == ne:
#     print("* Homomorphic subtraction works.")
# else:
#     print("* Homomorphic subtraction does not work.")
# print("")

# Multiply two message (no relinearization)
# ct = Evaluator.HomomorphicMultiplication(ct1, ct2)
# mt = Evaluator.DecryptionV2(ct)
#
# nr = Evaluator.IntDecode(mt)
# ne = (n1 * n2)
#
# print("--- Performing ct_mul = Enc(m1) * Enc(m2) (no relinearization)")
# print("* ct_mul[0] :{}".format(ct[0]))
# print("* ct_mul[1] :{}".format(ct[1]))
# print("--- Performing ct_dec = Dec(ct_sub)")
# print("* ct_dec    :{}".format(mt))
# print("--- Performing ct_dcd = Decode(ct_dec)")
# print("* ct_dcd    :{}".format(nr))
#
# if nr == ne:
#     print("* Homomorphic multiplication works.")
# else:
#     print("* Homomorphic multiplication does not work.")
# print("")
#
# # Multiply two message (relinearization v1)
# ct = Evaluator.HomomorphicMultiplication(ct1, ct2)
# ct = Evaluator.RelinearizationV1(ct)
# mt = Evaluator.Decryption(ct)
#
# nr = Evaluator.IntDecode(mt)
# ne = (n1 * n2)
#
# print("--- Performing ct_mul = Enc(m1) * Enc(m2) (with relinearization v1)")
# print("* ct_mul[0] :{}".format(ct[0]))
# print("* ct_mul[1] :{}".format(ct[1]))
# print("--- Performing ct_dec = Dec(ct_sub)")
# print("* ct_dec    :{}".format(mt))
# print("--- Performing ct_dcd = Decode(ct_dec)")
# print("* ct_dcd    :{}".format(nr))

# if nr == ne:
#     print("* Homomorphic multiplication works.")
# else:
#     print("* Homomorphic multiplication does not work.")
# print("")

"""
# Multiply two message (relinearization v2)
ct = Evaluator.HomomorphicMultiplication(ct1,ct2)
ct = Evaluator.RelinearizationV2(ct)
mt = Evaluator.Decryption(ct)

nr = Evaluator.IntDecode(mt)
ne = (n1*n2)

print("--- Performing ct_mul = Enc(m1) * Enc(m2) (with relinearization v2)")
print("* ct_mul[0] :{}".format(ct[0]))
print("* ct_mul[1] :{}".format(ct[1]))
print("--- Performing ct_dec = Dec(ct_sub)")
print("* ct_dec    :{}".format(mt))
print("--- Performing ct_dcd = Decode(ct_dec)")
print("* ct_dcd    :{}".format(nr))

if nr == ne:
    print("* Homomorphic multiplication works.")
else:
    print("* Homomorphic multiplication does not work.")
"""
#
