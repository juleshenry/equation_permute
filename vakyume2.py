from scipy import sqrt


class LinearEqn:
    @staticmethod
    def eqn___c(a: float, b: float):
        # a + b = c
        result = []
        c = a - b
        result.append(c)
        return c

    @staticmethod
    def eqn___b(a: float, c: float):
        # a + b = c
        result = []
        b = a - c
        result.append(b)
        return b

    @staticmethod
    def eqn___a(b: float, c: float):
        # a + b = c
        result = []
        a = b + c
        result.append(a)
        return a


class QuadEqn:
    @staticmethod
    def eqn_1_3__v(T: float, m: float):
        # .5 * m * v**2 = 1.5 * (1.380649e-23) * T
        result = []
        v = -6.4357959880655e-12 * sqrt(T / m)
        result.append(v)
        v = 6.4357959880655e-12 * sqrt(T / m)
        result.append(v)
        return v

    @staticmethod
    def eqn_1_3__m(T: float, v: float):
        # .5 * m * v**2 = 1.5 * (1.380649e-23) * T
        result = []
        m = 4.141947e-23 * T / v**2
        result.append(m)
        return m

    @staticmethod
    def eqn_1_3__T(m: float, v: float):
        # .5 * m * v**2 = 1.5 * (1.380649e-23) * T
        result = []
        T = 2.41432350534664e22 * m * v**2
        result.append(T)
        return T


class UniversalGasLaw:
    @staticmethod
    def eqn___V(R: float, T: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        V = R * T * n / p
        result.append(V)
        return V

    @staticmethod
    def eqn___n(R: float, T: float, V: float, p: float):
        # p * V = n * R * T
        result = []
        n = V * p / (R * T)
        result.append(n)
        return n

    @staticmethod
    def eqn___R(T: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        R = V * p / (T * n)
        result.append(R)
        return R

    @staticmethod
    def eqn___p(R: float, T: float, V: float, n: float):
        # p * V = n * R * T
        result = []
        p = R * T * n / V
        result.append(p)
        return p

    @staticmethod
    def eqn___T(R: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        T = V * p / (R * n)
        result.append(T)
        return T
