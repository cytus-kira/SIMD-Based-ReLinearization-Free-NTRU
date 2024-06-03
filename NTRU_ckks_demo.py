import pickle
import time
import numpy as np
from Ntru_ckks_encode import Ntru_Ckks
from helper import *
from math import log

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
use = 0
if use:
    with open('NTRU.pkl', 'rb') as file:
        Evaluator = pickle.load(file)

else:
    PD = 1  # 0: generate -- 1: pre-defined
    if PD == 0:
        # Select one of the parameter sets below
        # scale = 20
        # n, q, psi = 1024, 132120577, 73993  # log(q) = 27
        # t = 256;  n, q, psi = 2048 , 137438691329      , 22157790             # log(q) = 37
        scale = 40
        n, q, psi = 4096, 288230376151130113, 205582747241994398  # log(q) = 58

        # other necessary parameters
        psiv = modinv(psi, q)
        w = pow(psi, 2, q)
        wv = modinv(w, q)
    else:
        # Enter proper parameters below
        # scale, n, logq = 40, 16, 27
        # scale, n, logq = 20, 2048, 37
        scale, n, logq = 40, 8, 200
        level = 2
        q_ = log(n * (2 ** scale * 2 ** 7 * n) ** level * (6.5 * n ** 3.5 * log(n, 2)) ** level, 2)
        # print(q_)
        if q_ > logq:
            q = q_
        # other necessary parameters (based on n and log(q) determine other parameter)
        q, psi, psiv, w, wv = ParamGen(n, logq)

    # Determine mu, sigma (for discrete gaussian distribution)

    mu = 0
    sigma = int(6 * n ** 3.5 * log(n, 2))

    qnp1, q1 = get_qnp(n, logq)
    qnp2, q2 = get_qnp(n, logq - scale)
    qnp3, q3 = get_qnp(n, logq - scale - scale)
    # n2 = np.random.rand(n)
    print("--- Starting NTRU Demo")
    start_init_time = time.time()
    # Generate BFV evaluator
    Evaluator = Ntru_Ckks(n, q1, q2, q3, mu, sigma, qnp1, qnp2, qnp3, scale)
    print(Evaluator)

    # Generate Keys
    start_time = time.time()
    print('init_time', start_time - start_init_time)
    Evaluator.SecretKeyGen(read=False)
    Evaluator.PublicKeyGen(read=False)
    key_time = time.time()
    print(key_time - start_time)
    #
    # with open('NTRU.pkl', 'wb') as file:
    #     pickle.dump(Evaluator, file)

n1 = np.random.randint(-7, 7, int(Evaluator.n / 2), int).tolist()
# n1 = [1, 2, 3, 4]
# a = Evaluator.encode.encode(n1)
# b = Evaluator.encode.decode(a)
# print(n1)
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

# Encrypt message
m1 = Evaluator.encode.encode(n1)
decode_1 = Evaluator.encode.decode(np.polynomial.Polynomial(m1))
print(np.array(n1))
print(np.real(decode_1))
ct1 = Evaluator.Encryption(m1)
dec_1 = Evaluator.decryption(ct1)
# end_dec_time = time.time()
# print('dec_time', end_dec_time - start_dec_time)
#
# #
decode_1 = Evaluator.encode.decode(np.polynomial.Polynomial(dec_1.F))
f_v = Evaluator.encode.no_scale_decode(np.polynomial.Polynomial(Evaluator.sk.F))
v = decode_1 / f_v
decode_time = time.time()
# print('decode_time', decode_time - end_dec_time)
print('dec', np.real(v))
# start_init_time = time.time()
#
a = time.time()
mul_data = ct1 * ct1
rescale_data = Evaluator.rescale(mul_data)
print('mul_time', time.time() - a)
dec_2 = Evaluator.decryption_mul(rescale_data)

decode_mul = Evaluator.encode.decode(np.polynomial.Polynomial(dec_2.F))
#
f_v_dou = Evaluator.encode.no_scale_decode(np.polynomial.Polynomial(Evaluator.der_key.F))
v = decode_mul / f_v_dou
decode_time = time.time()
# print('decode_time', decode_time - end_dec_time)
print('dec', np.real(v))
# data = data * ct1
#
# # Generate Keys
# start_time = time.time()
# print('ct_mul', start_time - start_init_time)
#
# m_data = Poly(n, q, qnp)
# m_data.F = m1
# start_init_time = time.time()
# # Generate BFV evaluator
#
# data2 = ct1 * m_data
#
# # Generate Keys
# start_time = time.time()
# print('pc_mul', start_time - start_init_time)
