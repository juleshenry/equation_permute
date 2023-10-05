from scipy import sqrt
from sympy import Piecewise, Eqn


class pressMgmt:
    @staticmethod
    def eqn_3_1__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return Vacuum

    @staticmethod
    def eqn_3_1__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return Abs_Pressure

    @staticmethod
    def eqn_3_1__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return BarometricPressure

    @staticmethod
    def eqn_3_2__G_C(G: float, H: float, P: float, rho: float):
        # P = G / (G_C * rho * H)
        result = []
        G_C = G / (H * P * rho)
        result.append(G_C)
        return G_C

    @staticmethod
    def eqn_3_2__H(G: float, G_C: float, P: float, rho: float):
        # P = G / (G_C * rho * H)
        result = []
        H = G / (G_C * P * rho)
        result.append(H)
        return H

    @staticmethod
    def eqn_3_2__G(G_C: float, H: float, P: float, rho: float):
        # P = G / (G_C * rho * H)
        result = []
        G = G_C * H * P * rho
        result.append(G)
        return G

    @staticmethod
    def eqn_3_2__P(G: float, G_C: float, H: float, rho: float):
        # P = G / (G_C * rho * H)
        result = []
        P = G / (G_C * H * rho)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_2__rho(G: float, G_C: float, H: float, P: float):
        # P = G / (G_C * rho * H)
        result = []
        rho = G / (G_C * H * P)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_3_3__H_2(H_1: float, P: float, P_P: float):
        # P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 + P + P_P
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_3__P_P(H_1: float, H_2: float, P: float):
        # P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 - P
        result.append(P_P)
        return P_P

    @staticmethod
    def eqn_3_3__P(H_1: float, H_2: float, P_P: float):
        # P_P - P = H_2 - H_1
        result = []
        P = -H_1 + H_2 - P_P
        result.append(P)
        return P

    @staticmethod
    def eqn_3_3__H_1(H_2: float, P: float, P_P: float):
        # P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 - P - P_P
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_4__V(KAPPA: float, P: float):
        # P * V = KAPPA
        result = []
        V = KAPPA / P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_4__P(KAPPA: float, V: float):
        # P * V = KAPPA
        result = []
        P = KAPPA / V
        result.append(P)
        return P

    @staticmethod
    def eqn_3_4__KAPPA(P: float, V: float):
        # P * V = KAPPA
        result = []
        KAPPA = P * V
        result.append(KAPPA)
        return KAPPA

    @staticmethod
    def eqn_3_5__P_P(P: float, V: float, V_P: float):
        # P_P = P * (V / V_P)
        result = []
        P_P = P * V / V_P
        result.append(P_P)
        return P_P

    @staticmethod
    def eqn_3_5__V(P: float, P_P: float, V_P: float):
        # P_P = P * (V / V_P)
        result = []
        V = P_P * V_P / P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_5__P(P_P: float, V: float, V_P: float):
        # P_P = P * (V / V_P)
        result = []
        P = P_P * V_P / V
        result.append(P)
        return P

    @staticmethod
    def eqn_3_5__V_P(P: float, P_P: float, V: float):
        # P_P = P * (V / V_P)
        result = []
        V_P = P * V / P_P
        result.append(V_P)
        return V_P

    @staticmethod
    def eqn_3_6__V(H_1: float, H_2: float, P: float, V_P: float):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P * (-H_1 + H_2 + P) / P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_6__V_P(H_1: float, H_2: float, P: float, V: float):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P * V / (-H_1 + H_2 + P)
        result.append(V_P)
        return V_P

    @staticmethod
    def eqn_3_6__H_2(H_1: float, P: float, V: float, V_P: float):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P * V / V_P - P
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_6__P(H_1: float, H_2: float, V: float, V_P: float):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P * (-H_1 + H_2) / (V - V_P)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_6__H_1(H_2: float, P: float, V: float, V_P: float):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P * V / V_P + P
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_8__H_2(A_C: float, V_P: float):
        # V_P = A_C * H_2
        result = []
        H_2 = V_P / A_C
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_8__A_C(H_2: float, V_P: float):
        # V_P = A_C * H_2
        result = []
        A_C = V_P / H_2
        result.append(A_C)
        return A_C

    @staticmethod
    def eqn_3_8__V_P(A_C: float, H_2: float):
        # V_P = A_C * H_2
        result = []
        V_P = A_C * H_2
        result.append(V_P)
        return V_P

    @staticmethod
    def eqn_3_9__A_C(H_1: float, H_2: float, P: float, V: float):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        A_C = P * V / (H_2 * (-H_1 + H_2 + P))
        result.append(A_C)
        return A_C

    @staticmethod
    def eqn_3_9__V(A_C: float, H_1: float, H_2: float, P: float):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        V = A_C * H_2 * (-H_1 + H_2 + P) / P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_9__H_2(A_C: float, H_1: float, P: float, V: float):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_2 = (
            A_C * (H_1 - P)
            - sqrt(
                A_C * (A_C * H_1**2 - 2 * A_C * H_1 * P + A_C * P**2 + 4 * P * V)
            )
        ) / (2 * A_C)
        result.append(H_2)
        H_2 = (
            A_C * (H_1 - P)
            + sqrt(
                A_C * (A_C * H_1**2 - 2 * A_C * H_1 * P + A_C * P**2 + 4 * P * V)
            )
        ) / (2 * A_C)
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_9__P(A_C: float, H_1: float, H_2: float, V: float):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        P = A_C * H_2 * (H_1 - H_2) / (A_C * H_2 - V)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_9__H_1(A_C: float, H_2: float, P: float, V: float):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_1 = H_2 + P - P * V / (A_C * H_2)
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_11__H_2(A_C: float, P: float, V: float):
        # P = A_C / V * (H_2) ** 2
        result = []
        H_2 = -sqrt(P * V / A_C)
        result.append(H_2)
        H_2 = sqrt(P * V / A_C)
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_11__A_C(H_2: float, P: float, V: float):
        # P = A_C / V * (H_2) ** 2
        result = []
        A_C = P * V / H_2**2
        result.append(A_C)
        return A_C

    @staticmethod
    def eqn_3_11__V(A_C: float, H_2: float, P: float):
        # P = A_C / V * (H_2) ** 2
        result = []
        V = A_C * H_2**2 / P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_11__P(A_C: float, H_2: float, V: float):
        # P = A_C / V * (H_2) ** 2
        result = []
        P = A_C * H_2**2 / V
        result.append(P)
        return P

    @staticmethod
    def eqn_3_12__H_2(KAPPA_1: float, P: float):
        # P = KAPPA_1 * H_2 ** 2
        result = []
        H_2 = -sqrt(P / KAPPA_1)
        result.append(H_2)
        H_2 = sqrt(P / KAPPA_1)
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_12__KAPPA_1(H_2: float, P: float):
        # P = KAPPA_1 * H_2 ** 2
        result = []
        KAPPA_1 = P / H_2**2
        result.append(KAPPA_1)
        return KAPPA_1

    @staticmethod
    def eqn_3_12__P(H_2: float, KAPPA_1: float):
        # P = KAPPA_1 * H_2 ** 2
        result = []
        P = H_2**2 * KAPPA_1
        result.append(P)
        return P

    @staticmethod
    def eqn_3_13__H_2(H_1: float, KAPPA_2: float, P: float):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_2 = H_1 + P / KAPPA_2
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_13__P(H_1: float, H_2: float, KAPPA_2: float):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        P = KAPPA_2 * (-H_1 + H_2)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_13__H_1(H_2: float, KAPPA_2: float, P: float):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_1 = H_2 - P / KAPPA_2
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_13__KAPPA_2(H_1: float, H_2: float, P: float):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        KAPPA_2 = -P / (H_1 - H_2)
        result.append(KAPPA_2)
        return KAPPA_2

    @staticmethod
    def eqn_3_15__V_PMIN():
        # V_PMIN = 3.141592653589793 / 4
        result = []
        V_PMIN = 0.785398163397448
        result.append(V_PMIN)
        return V_PMIN

    @staticmethod
    def eqn_3_16__V_div_V_P_MAX():
        # V_div_V_P_MAX = 200000 / (3.141592653589793 / 4)
        result = []
        V_div_V_P_MAX = 254647.908947033
        result.append(V_div_V_P_MAX)
        return V_div_V_P_MAX

    @staticmethod
    def eqn_3_17__P_MIN():
        # P_MIN = (3.141592653589793 / 4) / (200000)
        result = []
        P_MIN = 0.00000392699081698724
        result.append(P_MIN)
        return P_MIN
