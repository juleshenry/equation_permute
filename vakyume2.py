class Eqn:
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
    def eqn_1_3__T(m: float, v: float):
        # .5 * m * v**2 = 1.5 * (1.380649e-23) * T
        result = []
        T = 2.41432350534664e22 * m * v**2
        result.append(T)
        return T

    @staticmethod
    def eqn_1_3__m(T: float, v: float):
        # .5 * m * v**2 = 1.5 * (1.380649e-23) * T
        result = []
        m = 4.141947e-23 * T / v**2
        result.append(m)
        return m
