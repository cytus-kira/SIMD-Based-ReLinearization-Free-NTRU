import json
import pickle

import gmpy2
import numpy as np
from sympy import symbols

from ckks_encode import CKKSEncoder
from new_test import invert_poly
from poly import *


class Ntru_Ckks:

    def __init__(self, n, q1, q2, q3, mu, sigma, qnp1, qnp2, qnp3, scale):
        self.n = n
        self.q = q1
        self.q2 = q2
        self.q3 = q3
        self.mu = mu
        self.sigma = sigma
        self.scale = scale
        self.scale_data = 2 ** self.scale
        self.qnp = qnp1  # array NTT parameters: [w,w_inv,psi,psi_inv]
        self.qnp2 = qnp2  # array NTT parameters: [w,w_inv,psi,psi_inv]
        self.qnp3 = qnp3  # array NTT parameters: [w,w_inv,psi,psi_inv]
        #

        self.sk = Poly(self.n, self.q, self.qnp)
        self.pk = Poly(self.n, self.q, self.qnp)
        self.f_q_poly = Poly(self.n, self.q, self.qnp)
        self.encode = CKKSEncoder(self.n * 2, 2 ** scale)

    #

    def __str__(self):
        str = "\n--- Parameters:\n"
        str = str + "n    : {}\n".format(self.n)
        str = str + "q    : {}\n".format(self.q)
        # str = str + "t    : {}\n".format(self.t)
        str = str + "mu   : {}\n".format(self.mu)
        str = str + "sigma: {}\n".format(self.sigma)
        return str

    #
    def SecretKeyGen(self, read=False):
        """
        sk <- R_2
        """
        if not read:
            # self.sk.F = [8, 7, 6, 5, 4, 3, 2, 1]
            self.sk.randomize(0, domain=False, type=0, mu=self.mu, sigma=self.sigma)
            name = 'sk_q' + str(self.q) + 'n' + str(self.n) + '.pkl'
            with open(name, 'wb') as f:
                pickle.dump(self.sk.F, f)
        else:
            name = 'sk_q' + str(self.q) + 'n' + str(self.n) + '.pkl'
            with open(name, 'rb') as f:
                self.sk.F = pickle.load(f)

    #
    def PublicKeyGen(self, read=False):
        """
        """
        if not read:
            g = Poly(self.n, self.q, self.qnp)
            g.randomize(0, domain=False, type=1, mu=self.mu, sigma=self.sigma)
            # g.F = [1, 2, 3, 4, 5, 6, 7, 8]
            x = symbols('x')
            from sympy import ZZ
            from sympy import Poly as spPoly
            self.R_poly = spPoly(x ** self.n + 1, x).set_domain(ZZ)
            f = spPoly(self.sk.F[::-1], x).set_domain(ZZ)
            f_q = [mod_q(data, self.q) for data in invert_poly(f, self.R_poly, self.q).all_coeffs()][::-1]
            self.f_q_poly.F = f_q
            self.pk = self.f_q_poly * g
            name = 'pk_q' + str(self.q) + 'n' + str(self.n) + '.pkl'
            with open(name, 'wb') as f:
                pickle.dump(self.pk.F, f)
            # with open('pk_8192_73q.txt', 'w') as f:
            #     f.write(json.dumps(self.sk.F))
        else:
            name = 'pk_q' + str(self.q) + 'n' + str(self.n) + '.pkl'
            with open(name, 'rb') as f:
                self.pk.F = pickle.load(f)
            # with open('ntru_pk.pkl', 'ab') as f:
            #     pickle.dump(self.pk.F, f)
            # with open('pk_8192_73q.txt', 'r') as f:
            #     self.pk.F = json.loads(f.read())

    #
    def Encryption(self, m):
        """
        """
        m_data = Poly(self.n, self.q, self.qnp)
        m_data.F = m

        s = Poly(self.n, self.q, self.qnp)
        s.randomize(2, type=1)
        # s.F = [1,1,1,1,2,2,2,2]
        e = Poly(self.n, self.q, self.qnp)
        e.randomize(2, domain=False, type=1, mu=self.mu, sigma=3.2)
        # e.F = [2,2,2,2,1,1,1,1]
        ct = (self.pk * s + m_data)
        return ct

        #

    def decryption(self, ct):
        """
        """
        m = ct * self.sk
        return m

    def decryption_mul(self, ct):
        data = (self.sk * self.sk).F
        self.der_key = Poly(self.n, self.q2, self.qnp2)
        self.der_key.F = data
        m = ct * self.der_key
        return m

    def rescale(self, data):
        data.F = [x // self.scale_data for x in data.F]
        data.q = self.q2
        data.np = self.qnp2

        return data

    def encode(self, data):
        return self.encode.encode(data)

    def decode(self, data):
        return self.encode.decode(data)

    def no_scale_decode(self, data):
        return self.encode.no_scale_decode(data)
