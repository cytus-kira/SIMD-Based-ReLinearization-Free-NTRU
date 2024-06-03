# BFV
import json

import numpy as np
from sympy import symbols

from new_test import invert_poly
from poly import *


class Ntru:
    # Definitions
    # Z_q[x]/f(x) = x^n + 1 where n=power-of-two

    # Operations
    # -- SecretKeyGen
    # -- PublicKeyGen
    # -- Encryption
    # -- Decryption
    # -- EvaluationKeyGenV1
    # -- EvaluationKeyGenV2 (need to be fixed)
    # -- HomAdd
    # -- HomMult
    # -- RelinV1
    # -- RelinV2 (need to be fixed)

    # Parameters
    # (From outside)
    # -- n (ring size)
    # -- q (ciphertext modulus)
    # -- t (plaintext modulus)
    # -- mu (distribution mean)
    # -- sigma (distribution std. dev.)
    # -- qnp (NTT parameters: [w,w_inv,psi,psi_inv])
    # (Generated with parameters)
    # -- sk
    # -- pk
    # -- rlk1, rlk2

    def __init__(self, n, q, t, mu, sigma, qnp):
        self.n = n
        self.q = q
        self.t = t
        self.mu = mu
        self.sigma = sigma
        self.qnp = qnp  # array NTT parameters: [w,w_inv,psi,psi_inv]
        #
        self.sk = Poly(self.n, self.q, self.qnp)
        self.pk = Poly(self.n, self.q, self.qnp)
        self.rlk1 = Poly(self.n, self.q, self.qnp)
        self.rlk2 = Poly(self.n, self.q, self.qnp)

    #
    def __str__(self):
        str = "\n--- Parameters:\n"
        str = str + "n    : {}\n".format(self.n)
        str = str + "q    : {}\n".format(self.q)
        str = str + "t    : {}\n".format(self.t)
        str = str + "mu   : {}\n".format(self.mu)
        str = str + "sigma: {}\n".format(self.sigma)
        return str

    #
    def SecretKeyGen(self, read=False):
        """
        sk <- R_2
        """
        if not read:
            s = Poly(self.n, self.q, self.qnp)
            s.randomize(0, domain=False, type=1, mu=self.mu, sigma=self.sigma)
            data = [data * self.t for data in s.F]
            data[0] += 1
            with open('sk.txt', 'w') as f:
                f.write(json.dumps(data))
            s.F = data
            self.sk = s
        else:
            with open('sk.txt', 'r') as f:
                self.sk.F = json.loads(f.read())

    #
    def PublicKeyGen(self, read=False):
        """
        a <- R_q
        e <- X
        pk[0] <- (-(a*sk)+e) mod q
        pk[1] <- a
        """
        if not read:
            g = Poly(self.n, self.q, self.qnp)
            g.randomize(0, domain=False, type=1, mu=self.mu, sigma=self.sigma)
            x = symbols('x')
            from sympy import ZZ
            from sympy import Poly as spPoly
            self.R_poly = spPoly(x ** self.n + 1, x).set_domain(ZZ)
            f = spPoly(self.sk.F[::-1], x).set_domain(ZZ)
            f_q = [int(data) for data in invert_poly(f, self.R_poly, self.q).all_coeffs()][::-1]
            f_q_poly = Poly(self.n, self.q, self.qnp)
            f_q_poly.F = f_q
            gf_q = f_q_poly * g
            self.pk.F = [(self.t * data) % self.q for data in gf_q.F]
            with open('pk.txt', 'w') as f:
                f.write(json.dumps(self.pk.F))
        else:
            with open('pk.txt', 'r') as f:
                self.pk.F = json.loads(f.read())

    #
    def Encryption(self, m):
        """
        """
        m_data = Poly(self.n, self.q, self.qnp)
        m_data.F = m
        s = Poly(self.n, self.q, self.qnp)
        s.randomize(1)
        e = Poly(self.n, self.q, self.qnp)
        e.randomize(1)
        e.F = [data * self.t for data in e.F]
        ct = self.pk * s + m_data + e
        return ct

    #
    def decryption(self, ct):
        """
        """
        m = ct * self.sk
        m = m % self.q % self.t

        return m