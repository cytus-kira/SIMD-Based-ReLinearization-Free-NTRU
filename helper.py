import gmpy2

from generate_prime import *
from random import randint


# Modular inverse of an integer
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return x % m


# GCD of two integers
def gcd(n1, n2):
    a = n1
    b = n2
    while b != 0:
        a, b = b, a % b
    return a


# Bit-Reverse integer
def intReverse(a, n):
    b = ('{:0' + str(n) + 'b}').format(a)
    return int(b[::-1], 2)


# Bit-Reversed index
def indexReverse(a, r):
    n = len(a)
    b = [0] * n
    for i in range(n):
        rev_idx = intReverse(i, r)
        b[rev_idx] = a[i]
    return b


# Reference Polynomial Multiplication
# with f(x) = x^n + 1
def RefPolMul(A, B, M):
    C = [0] * (2 * len(A))
    D = [0] * (len(A))
    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % M

    for i in range(len(A)):
        D[i] = (C[i] - C[i + len(A)]) % M
    return D


def get_qnp(n, logq):
    q, psi, psiv, w, wv = ParamGen(n, logq)
    w_table = [1] * n
    wv_table = [1] * n
    psi_table = [1] * n
    psiv_table = [1] * n
    for i in range(1, n):
        w_table[i] = mod_q((w_table[i - 1] * w), q)
        wv_table[i] = mod_q((wv_table[i - 1] * wv), q)
        psi_table[i] = mod_q((psi_table[i - 1] * psi), q)
        psiv_table[i] = mod_q((psiv_table[i - 1] * psiv), q)
    qnp = [w_table, wv_table, psi_table, psiv_table]
    return qnp, q


# Reference Polynomial Multiplication (w/ modulus)
# with f(x) = x^n + 1
def RefPolMulv2(A, B):
    C = [0] * (2 * len(A))
    D = [0] * (len(A))
    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB)

    for i in range(len(A)):
        D[i] = (C[i] - C[i + len(A)])
    return D


# Check if input is m-th (could be n or 2n) primitive root of unity of q
def isrootofunity(w, m, q):
    if w == 0:
        return False
    elif pow(w, m // 2, q) == (q - 1):
        return True
    else:
        return False


# Returns a proper NTT-friendly prime
def GetProperPrime(n, logq):
    factor = 2 * n
    value = (1 << logq) - factor + 1
    lbound = (1 << (logq - 1))
    while (value > lbound):
        if is_prime(value) == True:
            return value
        else:
            value = value - factor
    raise Exception("Failed to find a proper prime.")


# Returns a primitive root
def FindPrimitiveRoot(m, q):
    g = (q - 1) // m

    if (q - 1) != g * m:
        return False

    attempt_ctr = 0
    attempt_max = 100

    while (attempt_ctr < attempt_max):
        a = randint(2, q - 1)
        b = pow(a, g, q)
        # check 
        if isrootofunity(b, m, q):
            return True, b
        else:
            attempt_ctr = attempt_ctr + 1

    return True, 0


# Generate necessary BFV parameters given n and log(q)
def ParamGen(n, logq):
    pfound = False
    while (not (pfound)):
        q = GetProperPrime(n, logq)
        pfound, psi = FindPrimitiveRoot(2 * n, q)
    psiv = modinv(psi, q)
    w = pow(psi, 2, q)
    wv = modinv(w, q)
    return q, psi, psiv, w, wv


def is_prime(n, k=5):
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0:
        return False

    # Find r and s such that n-1 = 2^r * s
    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2

    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, s, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False

    return True


def is_2_power(n):
    return n != 0 and (n & (n - 1) == 0)


def mod_q(x, q):
    """
    计算 value 对 q 取模，结果在 -q/2 到 q/2 之间。
    """
    x = gmpy2.mpz(x)
    q = gmpy2.mpz(q)
    adjusted_x = gmpy2.f_mod(x + q // 2, q) - q // 2
    return adjusted_x
