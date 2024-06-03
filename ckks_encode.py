import json

import gmpy2
import numpy as np
from joblib import Parallel, delayed
from numpy.polynomial import Polynomial
from scipy.linalg import lu_factor, lu_solve


def round_coordinates(coordinates):
    """Gives the integral rest."""
    coordinates = coordinates - np.floor(coordinates)
    return coordinates


def coordinate_wise_random_rounding(coordinates):
    """Rounds coordinates randonmly."""
    r = round_coordinates(coordinates)
    f = np.array([np.random.choice([c, c - 1], 1, p=[1 - c, c]) for c in r]).reshape(-1)

    rounded_coordinates = coordinates - f
    rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
    return rounded_coordinates


class CKKSEncoder:
    def __init__(self, poly_degree: int, scale: float):
        """Initialization of the encoder for poly_degree a power of 2.

        xi, which is an poly_degree-th root of unity will, be used as a basis for our computations.
        """

        self.xi = np.exp(2 * np.pi * 1j / poly_degree)
        self.poly_degree = poly_degree

        self.create_sigma_R_basis()
        self.scale = scale

        self.roots = [self.xi ** (2 * i + 1) for i in range(poly_degree // 2)]

        self.A_inv = lu_factor(CKKSEncoder.vandermonde(self.xi, self.poly_degree))


    def pi(self, z: np.array) -> np.array:
        """Projects a vector of H into C^{N/2}."""

        N = self.poly_degree // 4
        return z[:N]

    def pi_inverse(self, z: np.array) -> np.array:
        """Expands a vector of C^{N/2} by expanding it with its
        complex conjugate."""

        z_conjugate = z[::-1]
        z_conjugate = [np.conjugate(x) for x in z_conjugate]
        return np.concatenate([z, z_conjugate])

    def create_sigma_R_basis(self):
        """Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))."""

        self.sigma_R_basis = np.array(self.vandermonde(self.xi, self.poly_degree)).T

    def compute_basis_coordinates(self, z):
        """Computes the coordinates of a vector with respect to the orthogonal lattice basis."""
        output = np.array([np.real(np.vdot(z, b) / np.vdot(b, b)) for b in self.sigma_R_basis])
        return output

    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding."""
        coordinates = self.compute_basis_coordinates(z)

        rounded_coordinates = coordinate_wise_random_rounding(coordinates)
        y = np.matmul(self.sigma_R_basis.T, rounded_coordinates)
        return y

    def encode(self, z: np.array) -> Polynomial:
        """Encodes a vector by expanding it first to H,
        scale it, project it on the lattice of sigma(R), and performs
        sigma inverse.
        """
        pi_z = self.pi_inverse(z)
        scaled_pi_z = self.scale * pi_z
        rounded_scale_pi_zi = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scale_pi_zi)

        # We round it afterwards due to numerical imprecision
        coef = np.round(np.real(p))

        return [gmpy2.mpz(data) for data in coef]

    def decode(self, p: Polynomial) -> np.array:
        """Decodes a polynomial by removing the scale,
        evaluating on the roots, and project it on C^(N/2)"""
        # if type(p) != np.polynomial:
        #     p = Polynomial(p)
        rescaled_p = p / self.scale
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        return pi_z

    def decode_mul(self, p: Polynomial,degree) -> np.array:
        """Decodes a polynomial by removing the scale,
        evaluating on the roots, and project it on C^(N/2)"""
        # if type(p) != np.polynomial:
        #     p = Polynomial(p)
        for i in range(degree):
            p = p / self.scale
        z = self.sigma(p)
        pi_z = self.pi(z)
        return pi_z

    def no_scale_decode(self, p: Polynomial) -> np.array:
        """Decodes a polynomial by removing the scale,
        evaluating on the roots, and project it on C^(N/2)"""
        # rescaled_p = p / self.scale
        z = self.sigma(p)
        pi_z = self.pi(z)
        return pi_z

    @staticmethod
    def vandermonde(xi: np.complex128, M: int) -> np.array:
        """Computes the Vandermonde matrix from a m-th root of unity."""

        N = M // 2
        matrix = []
        # We will generate each row of the matrix
        for i in range(N):
            # For each row we select a different root
            root = xi ** (2 * i + 1)
            row = []

            # Then we store its powers
            for j in range(N):
                row.append(root ** j)
            matrix.append(row)
        return matrix

    def sigma_inverse(self, b: np.array) -> Polynomial:
        """Encodes the vector b in a polynomial using an M-th root of unity."""

        # First we create the Vandermonde matrix

        # Then we solve the system
        coeffs = lu_solve(self.A_inv, b)


        return coeffs

    def sigma(self, p: Polynomial) -> np.array:
        """Decodes a polynomial by applying it to the M-th roots of unity."""
        # outputs = []
        # outputs = Parallel(n_jobs=-1)(delayed(p)(root) for root in self.roots)
        # outputs = np.polyval(p.coef, self.roots)
        # sorted_roots = np.sort(self.roots)
        outputs = np.polyval(p.coef[::-1], self.roots)
        return np.array(outputs)


if __name__ == '__main__':
    a = 64
    p1 = np.random.random(a)
    print(p1)
    encode = CKKSEncoder(a * 4, 2 ** 40)
    e1 = encode.encode(p1)
    # print(e1)
    d_1 = encode.decode(e1)
    print(np.real(d_1))
