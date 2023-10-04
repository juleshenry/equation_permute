class GasLaw:
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

    @staticmethod
    def eqn___R(T: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        R = V * p / (T * n)
        result.append(R)
        return R


class GasLaw:
    @staticmethod
    def eqn___p(R: float, T: float, V: float, n: float):
        # p * V = n * R * T
        result = []
        p = R * T * n / V
        result.append(p)
        return p

    @staticmethod
    def eqn___V(R: float, T: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        V = R * T * n / p
        result.append(V)
        return V

    @staticmethod
    def eqn___R(T: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        R = V * p / (T * n)
        result.append(R)
        return R
        pass  # NotImplementedError

    @staticmethod
    def eqn___T(R: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        T = V * p / (R * n)
        result.append(T)
        return T

    @staticmethod
    def eqn___n(R: float, T: float, V: float, p: float):
        # p * V = n * R * T
        result = []
        n = V * p / (R * T)
        result.append(n)
        return n


class GasLaw:
    @staticmethod
    def eqn___T(R: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        T = V * p / (R * n)
        result.append(T)
        return T

    @staticmethod
    def eqn___R(T: float, V: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        R = V * p / (T * n)
        result.append(R)
        return R

    @staticmethod
    def eqn___V(R: float, T: float, n: float, p: float):
        # p * V = n * R * T
        result = []
        V = R * T * n / p
        result.append(V)
        return V

    @staticmethod
    def eqn___p(R: float, T: float, V: float, n: float):
        # p * V = n * R * T
        result = []
        p = R * T * n / V
        result.append(p)
        return p

    @staticmethod
    def eqn___n(R: float, T: float, V: float, p: float):
        # p * V = n * R * T
        result = []
        n = V * p / (R * T)
        result.append(n)
        return n
