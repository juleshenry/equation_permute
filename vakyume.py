

class VacuumTheory:
    @staticmethod
    def eqn_1_3__k(T: float, m: float, v ** 2: float):
    # .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return k
    @staticmethod
    def eqn_1_3__m(T: float, k: float, v ** 2: float):
    # .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return m
    @staticmethod
    def eqn_1_3__v ** 2(T: float, k: float, m: float):
    # .5 * m * v**2 = 1.5 * k * T
        pass # unable  to solve
    @staticmethod
    def eqn_1_3__T(k: float, m: float, v ** 2: float):
    # .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return T
    @staticmethod
    def eqn_1_7__n(R: float, T: float, V: float, p: float):
    # p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return n
    @staticmethod
    def eqn_1_7__p(R: float, T: float, V: float, n: float):
    # p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return p
    @staticmethod
    def eqn_1_7__V(R: float, T: float, n: float, p: float):
    # p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return V
    @staticmethod
    def eqn_1_7__R(T: float, V: float, n: float, p: float):
    # p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return R
    @staticmethod
    def eqn_1_7__T(R: float, V: float, n: float, p: float):
    # p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return T
    @staticmethod
    def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float):
    # P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return m
    @staticmethod
    def eqn_1_8__M(P: float, R: float, T: float, V: float, m: float):
    # P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return M
    @staticmethod
    def eqn_1_8__P(M: float, R: float, T: float, V: float, m: float):
    # P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return P
    @staticmethod
    def eqn_1_8__R(M: float, P: float, T: float, V: float, m: float):
    # P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return R
    @staticmethod
    def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float):
    # P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return V
    @staticmethod
    def eqn_1_8__T(M: float, P: float, R: float, V: float, m: float):
    # P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return T
    @staticmethod
    def eqn_1_9__M(P: float, R: float, T: float, rho: float):
    # rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return M
    @staticmethod
    def eqn_1_9__P(M: float, R: float, T: float, rho: float):
    # rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return P
    @staticmethod
    def eqn_1_9__R(M: float, P: float, T: float, rho: float):
    # rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return R
    @staticmethod
    def eqn_1_9__rho(M: float, P: float, R: float, T: float):
    # rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_1_9__T(M: float, P: float, R: float, rho: float):
    # rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return T
    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
    # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return P_1
    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
    # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return T_1
    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
    # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return P_2
    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
    # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return T_2
    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
    # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return V_2
    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
    # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return V_1
    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
    # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        pass # unable  to solve
    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
    # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        W = 0
        result.append(W)
        return W
    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
    # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        pass # unable  to solve
    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
    # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        pass # unable  to solve
    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
    # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        T = 0
        result.append(T)
        return T
    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
    # Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return Total_P
    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
    # Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return sum_partial_pressures
    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
    # y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return y_a
    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
    # y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return n_a
    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
    # y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return n
    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
    # y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return y_a
    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
    # y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return p_a
    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
    # y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return P


class FluidFlowVacuumLines:
    @staticmethod
    def eqn_2_1__Re(D: float, mu: float, rho: float, v: float):
    # Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return Re
    @staticmethod
    def eqn_2_1__v(D: float, Re: float, mu: float, rho: float):
    # Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return v
    @staticmethod
    def eqn_2_1__mu(D: float, Re: float, rho: float, v: float):
    # Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_1__rho(D: float, Re: float, mu: float, v: float):
    # Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_2_1__D(Re: float, mu: float, rho: float, v: float):
    # Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return D
    @staticmethod
    def eqn_2_2__lambd(delta: float, psi: float):
    # lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
        return lambd
    @staticmethod
    def eqn_2_2__delta(lambd: float, psi: float):
    # lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        delta = -0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        delta = 0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        return delta
    @staticmethod
    def eqn_2_2__psi(delta: float, lambd: float):
    # lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return psi
    @staticmethod
    def eqn_2_3__lambd(D: float, kn: float):
    # kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return lambd
    @staticmethod
    def eqn_2_3__kn(D: float, lambd: float):
    # kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return kn
    @staticmethod
    def eqn_2_3__D(kn: float, lambd: float):
    # kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return D
    @staticmethod
    def eqn_2_4__vel_grad(_beta: float, mu: float):
    # _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return vel_grad
    @staticmethod
    def eqn_2_4___beta(mu: float, vel_grad: float):
    # _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return _beta
    @staticmethod
    def eqn_2_4__mu(_beta: float, vel_grad: float):
    # _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_5__q(D: float, L: float, delta_P: float, mu: float):
    # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return q
    @staticmethod
    def eqn_2_5__L(D: float, delta_P: float, mu: float, q: float):
    # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_5__delta_P(D: float, L: float, mu: float, q: float):
    # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return delta_P
    @staticmethod
    def eqn_2_5__mu(D: float, L: float, delta_P: float, q: float):
    # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_5__D(L: float, delta_P: float, mu: float, q: float):
    # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        D = -2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = -2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        return D
    @staticmethod
    def eqn_2_6__lambd(mu: float, rho: float, v_a: float):
    # mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return lambd
    @staticmethod
    def eqn_2_6__v_a(lambd: float, mu: float, rho: float):
    # mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return v_a
    @staticmethod
    def eqn_2_6__rho(lambd: float, mu: float, v_a: float):
    # mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_2_6__mu(lambd: float, rho: float, v_a: float):
    # mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_7__v_a(T: float, k: float, m: float):
    # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return v_a
    @staticmethod
    def eqn_2_7__k(T: float, m: float, v_a: float):
    # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return k
    @staticmethod
    def eqn_2_7__m(T: float, k: float, v_a: float):
    # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return m
    @staticmethod
    def eqn_2_7__T(k: float, m: float, v_a: float):
    # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return T
    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return T_c
    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return P_c
    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return M
    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return mu_c
    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return T_c
    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return P_c
    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return M
    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
    # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return mu_c
    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
    # Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return oper_press
    @staticmethod
    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float):
    # Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return Suc_Pres
    @staticmethod
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float):
    # Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return delta_P
    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float):
    # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return h_r
    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
    # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return v
    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
    # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float):
    # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return g_c
    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
    # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return f
    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
    # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return D
    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
        return g
    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return d
    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = -0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        v = 0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        return v
    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155*L*f*rho*v**2/(d*g)
        result.append(delta_P)
        return delta_P
    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838*d*delta_P*g/(L*rho*v**2)
        result.append(f)
        return f
    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
    # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_2_13__d(L: float, delta_P: float, f: float, q: float, rho: float):
    # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        d = 1.16543402167043*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) - 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) + 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) - 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) + 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        return d
    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
    # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return q
    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
    # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, f: float, q: float, rho: float):
    # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        delta_P = 2.15*L*f*q**2*rho/d**5
        result.append(delta_P)
        return delta_P
    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
    # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return f
    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
    # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float):
    # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return v_s
    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float):
    # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return M
    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
    # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return k
    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float):
    # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return R
    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
    # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return g_c
    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
    # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return T
    @staticmethod
    def eqn_2_15__f(Re: float):
    # f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return f
    @staticmethod
    def eqn_2_15__Re(f: float):
    # f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return Re
    @staticmethod
    def eqn_2_16__f(Re: float):
    # f = 64 / Re
        result = []
        f = 64/Re
        result.append(f)
        return f
    @staticmethod
    def eqn_2_16__Re(f: float):
    # f = 64 / Re
        result = []
        Re = 64/f
        result.append(Re)
        return Re
    @staticmethod
    def eqn_2_16__d ** 2(L: float, delta_P: float, mu: float, v: float):
    # delta_P = 0.0345* mu * L * v / d**2
        pass # unable  to solve
    @staticmethod
    def eqn_2_16__v(L: float, d ** 2: float, delta_P: float, mu: float):
    # delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return v
    @staticmethod
    def eqn_2_16__L(d ** 2: float, delta_P: float, mu: float, v: float):
    # delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_16__delta_P(L: float, d ** 2: float, mu: float, v: float):
    # delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return delta_P
    @staticmethod
    def eqn_2_16__mu(L: float, d ** 2: float, delta_P: float, v: float):
    # delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_16__q(L: float, d ** 4: float, delta_P: float, mu: float):
    # delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return q
    @staticmethod
    def eqn_2_16__L(d ** 4: float, delta_P: float, mu: float, q: float):
    # delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952*d**4*delta_P/(mu*q)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_16__delta_P(L: float, d ** 4: float, mu: float, q: float):
    # delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105*L*mu*q/d**4
        result.append(delta_P)
        return delta_P
    @staticmethod
    def eqn_2_16__mu(L: float, d ** 4: float, delta_P: float, q: float):
    # delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952*d**4*delta_P/(L*q)
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_16__d ** 4(L: float, delta_P: float, mu: float, q: float):
    # delta_P = 0.105 * mu * L * q / d**4
        pass # unable  to solve
    @staticmethod
    def eqn_2_18a__D_eq(R_ll: float):
    # D_eq = 4 * R_ll
        result = []
        D_eq = 4*R_ll
        result.append(D_eq)
        return D_eq
    @staticmethod
    def eqn_2_18a__R_ll(D_eq: float):
    # D_eq = 4 * R_ll
        result = []
        R_ll = D_eq/4
        result.append(R_ll)
        return R_ll
    @staticmethod
    def eqn_2_18b__R_ll(h: float, w: float):
    # R_ll = w * h / (2 * (w + h))
        result = []
        R_ll = h*w/(2*(h + w))
        result.append(R_ll)
        return R_ll
    @staticmethod
    def eqn_2_18b__h(R_ll: float, w: float):
    # R_ll = w * h / (2 * (w + h))
        result = []
        h = 2*R_ll*w/(-2*R_ll + w)
        result.append(h)
        return h
    @staticmethod
    def eqn_2_18b__w(R_ll: float, h: float):
    # R_ll = w * h / (2 * (w + h))
        result = []
        w = 2*R_ll*h/(-2*R_ll + h)
        result.append(w)
        return w
    @staticmethod
    def eqn_2_19a__Re(R_ll: float, mu: float, rho: float, v: float):
    # Re = 4 * R_ll * rho * v / mu
        result = []
        Re = 4*R_ll*rho*v/mu
        result.append(Re)
        return Re
    @staticmethod
    def eqn_2_19a__v(R_ll: float, Re: float, mu: float, rho: float):
    # Re = 4 * R_ll * rho * v / mu
        result = []
        v = Re*mu/(4*R_ll*rho)
        result.append(v)
        return v
    @staticmethod
    def eqn_2_19a__mu(R_ll: float, Re: float, rho: float, v: float):
    # Re = 4 * R_ll * rho * v / mu
        result = []
        mu = 4*R_ll*rho*v/Re
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_19a__R_ll(Re: float, mu: float, rho: float, v: float):
    # Re = 4 * R_ll * rho * v / mu
        result = []
        R_ll = Re*mu/(4*rho*v)
        result.append(R_ll)
        return R_ll
    @staticmethod
    def eqn_2_19a__rho(R_ll: float, Re: float, mu: float, v: float):
    # Re = 4 * R_ll * rho * v / mu
        result = []
        rho = Re*mu/(4*R_ll*v)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_2_19b__Re(h: float, mu: float, rho: float, v: float, w: float):
    # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        Re = 2*h*rho*v*w/(mu*(h + w))
        result.append(Re)
        return Re
    @staticmethod
    def eqn_2_19b__w(Re: float, h: float, mu: float, rho: float, v: float):
    # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        w = Re*h*mu/(-Re*mu + 2*h*rho*v)
        result.append(w)
        return w
    @staticmethod
    def eqn_2_19b__v(Re: float, h: float, mu: float, rho: float, w: float):
    # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        v = Re*mu*(h + w)/(2*h*rho*w)
        result.append(v)
        return v
    @staticmethod
    def eqn_2_19b__h(Re: float, mu: float, rho: float, v: float, w: float):
    # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        h = Re*mu*w/(-Re*mu + 2*rho*v*w)
        result.append(h)
        return h
    @staticmethod
    def eqn_2_19b__mu(Re: float, h: float, rho: float, v: float, w: float):
    # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        mu = 2*h*rho*v*w/(Re*(h + w))
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_19b__rho(Re: float, h: float, mu: float, v: float, w: float):
    # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        rho = Re*mu*(h + w)/(2*h*v*w)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_2_20__sum_pipe(L: float, sum_equivalent_length: float):
    # L = sum_pipe + sum_equivalent_length
        result = []
        sum_pipe = L - sum_equivalent_length
        result.append(sum_pipe)
        return sum_pipe
    @staticmethod
    def eqn_2_20__sum_equivalent_length(L: float, sum_pipe: float):
    # L = sum_pipe + sum_equivalent_length
        result = []
        sum_equivalent_length = L - sum_pipe
        result.append(sum_equivalent_length)
        return sum_equivalent_length
    @staticmethod
    def eqn_2_20__L(sum_equivalent_length: float, sum_pipe: float):
    # L = sum_pipe + sum_equivalent_length
        result = []
        L = sum_equivalent_length + sum_pipe
        result.append(L)
        return L
    @staticmethod
    def eqn_2_22__S_p(P_s: float, Q_throughput: float):
    # Q_throughput = S_p * P_s
        result = []
        S_p = Q_throughput/P_s
        result.append(S_p)
        return S_p
    @staticmethod
    def eqn_2_22__P_s(Q_throughput: float, S_p: float):
    # Q_throughput = S_p * P_s
        result = []
        P_s = Q_throughput/S_p
        result.append(P_s)
        return P_s
    @staticmethod
    def eqn_2_22__Q_throughput(P_s: float, S_p: float):
    # Q_throughput = S_p * P_s
        result = []
        Q_throughput = P_s*S_p
        result.append(Q_throughput)
        return Q_throughput
    @staticmethod
    def eqn_2_25__P_1(C: float, P_2: float, Q_throughput: float):
    # C = Q_throughput / (P_1 - P_2)
        result = []
        P_1 = P_2 + Q_throughput/C
        result.append(P_1)
        return P_1
    @staticmethod
    def eqn_2_25__C(P_1: float, P_2: float, Q_throughput: float):
    # C = Q_throughput / (P_1 - P_2)
        result = []
        C = Q_throughput/(P_1 - P_2)
        result.append(C)
        return C
    @staticmethod
    def eqn_2_25__P_2(C: float, P_1: float, Q_throughput: float):
    # C = Q_throughput / (P_1 - P_2)
        result = []
        P_2 = P_1 - Q_throughput/C
        result.append(P_2)
        return P_2
    @staticmethod
    def eqn_2_25__Q_throughput(C: float, P_1: float, P_2: float):
    # C = Q_throughput / (P_1 - P_2)
        result = []
        Q_throughput = C*(P_1 - P_2)
        result.append(Q_throughput)
        return Q_throughput
    @staticmethod
    def eqn_2_26__q(D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, mu: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        q = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*mu)
        result.append(q)
        return q
    @staticmethod
    def eqn_2_26__L(D: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, q: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        L = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(mu*q)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_26__P_p(D: float, L: float, P_downstream: float, P_upstream: float, mu: float, q: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_p = 0.0
        result.append(P_p)
        return P_p
    @staticmethod
    def eqn_2_26__mu(D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, q: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        mu = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*q)
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_26__P_downstream(D: float, L: float, P_p: float, P_upstream: float, mu: float, q: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_downstream = P_upstream - 40.7436654315252*L*mu*q/D**4
        result.append(P_downstream)
        return P_downstream
    @staticmethod
    def eqn_2_26__P_upstream(D: float, L: float, P_downstream: float, P_p: float, mu: float, q: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_upstream = P_downstream + 40.7436654315252*L*mu*q/D**4
        result.append(P_upstream)
        return P_upstream
    @staticmethod
    def eqn_2_26__D(L: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, q: float):
    # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        D = -14953.4878122122*I*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = 14953.4878122122*I*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = -14953.4878122122*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = 14953.4878122122*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        return D
    @staticmethod
    def eqn_2_28__L(C: float, D: float, P_p: float, mu: float):
    # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        L = 0.0245436926061703*D**4*P_p/(C*mu)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_28__P_p(C: float, D: float, L: float, mu: float):
    # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        P_p = 40.7436654315252*C*L*mu/D**4
        result.append(P_p)
        return P_p
    @staticmethod
    def eqn_2_28__mu(C: float, D: float, L: float, P_p: float):
    # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        mu = 0.0245436926061703*D**4*P_p/(C*L)
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_28__C(D: float, L: float, P_p: float, mu: float):
    # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        C = 0.0245436926061703*D**4*P_p/(L*mu)
        result.append(C)
        return C
    @staticmethod
    def eqn_2_28__D(C: float, L: float, P_p: float, mu: float):
    # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        D = -2.52647511098426*I*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = -2.52647511098426*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = 2.52647511098426*(C*L*mu/P_p)**(1/4)
        result.append(D)
        return D
    @staticmethod
    def eqn_2_29__S_1(C: float, S_2: float):
    # S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_1 = C*S_2/(C + S_2)
        result.append(S_1)
        return S_1
    @staticmethod
    def eqn_2_29__C(S_1: float, S_2: float):
    # S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        C = -S_1*S_2/(S_1 - S_2)
        result.append(C)
        return C
    @staticmethod
    def eqn_2_29__S_2(C: float, S_1: float):
    # S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_2 = C*S_1/(C - S_1)
        result.append(S_2)
        return S_2
    @staticmethod
    def eqn_2_31__S_p(C: float, S_pump_speed: float):
    # S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_p = C*S_pump_speed/(C - S_pump_speed)
        result.append(S_p)
        return S_p
    @staticmethod
    def eqn_2_31__C(S_p: float, S_pump_speed: float):
    # S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        C = S_p*S_pump_speed/(S_p - S_pump_speed)
        result.append(C)
        return C
    @staticmethod
    def eqn_2_31__S_pump_speed(C: float, S_p: float):
    # S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_pump_speed = C*S_p/(C + S_p)
        result.append(S_pump_speed)
        return S_pump_speed
    @staticmethod
    def eqn_2_32__geometric_sum_C(C_series: float):
    # 1 / C_series = geometric_sum_C
        result = []
        geometric_sum_C = 1/C_series
        result.append(geometric_sum_C)
        return geometric_sum_C
    @staticmethod
    def eqn_2_32__C_series(geometric_sum_C: float):
    # 1 / C_series = geometric_sum_C
        result = []
        C_series = 1/geometric_sum_C
        result.append(C_series)
        return C_series
    @staticmethod
    def eqn_2_33__C_paralell(arithmetic_sum_C: float):
    # 1 / C_paralell = arithmetic_sum_C
        result = []
        C_paralell = 1/arithmetic_sum_C
        result.append(C_paralell)
        return C_paralell
    @staticmethod
    def eqn_2_33__arithmetic_sum_C(C_paralell: float):
    # 1 / C_paralell = arithmetic_sum_C
        result = []
        arithmetic_sum_C = 1/C_paralell
        result.append(arithmetic_sum_C)
        return arithmetic_sum_C
    @staticmethod
    def eqn_2_34__C_1(C: float, C_2: float, D: float, L: float, P_p: float, mu: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C_1 = mu*(C*L - C_2*D**3)/(D**4*P_p)
        result.append(C_1)
        return C_1
    @staticmethod
    def eqn_2_34__C_2(C: float, C_1: float, D: float, L: float, P_p: float, mu: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C_2 = C*L/D**3 - C_1*D*P_p/mu
        result.append(C_2)
        return C_2
    @staticmethod
    def eqn_2_34__L(C: float, C_1: float, C_2: float, D: float, P_p: float, mu: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        L = D**3*(C_1*D*P_p + C_2*mu)/(C*mu)
        result.append(L)
        return L
    @staticmethod
    def eqn_2_34__P_p(C: float, C_1: float, C_2: float, D: float, L: float, mu: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        P_p = mu*(C*L - C_2*D**3)/(C_1*D**4)
        result.append(P_p)
        return P_p
    @staticmethod
    def eqn_2_34__mu(C: float, C_1: float, C_2: float, D: float, L: float, P_p: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        mu = C_1*D**4*P_p/(C*L - C_2*D**3)
        result.append(mu)
        return mu
    @staticmethod
    def eqn_2_34__C(C_1: float, C_2: float, D: float, L: float, P_p: float, mu: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C = D**3*(C_1*D*P_p + C_2*mu)/(L*mu)
        result.append(C)
        return C
    @staticmethod
    def eqn_2_34__D(C: float, C_1: float, C_2: float, L: float, P_p: float, mu: float):
    # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        D = Piecewise((-sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (-sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((-sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (-sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        return D
    @staticmethod
    def eqn_2_11__C_L(C_T: float, F_p: float):
    # C_T = C_L * F_p
        result = []
        C_L = C_T/F_p
        result.append(C_L)
        return C_L
    @staticmethod
    def eqn_2_11__C_T(C_L: float, F_p: float):
    # C_T = C_L * F_p
        result = []
        C_T = C_L*F_p
        result.append(C_T)
        return C_T
    @staticmethod
    def eqn_2_11__F_p(C_L: float, C_T: float):
    # C_T = C_L * F_p
        result = []
        F_p = C_T/C_L
        result.append(F_p)
        return F_p
    @staticmethod
    def eqn_2_36__F_t(C: float, C_0: float):
    # C = C_0 * F_t
        result = []
        F_t = C/C_0
        result.append(F_t)
        return F_t
    @staticmethod
    def eqn_2_36__C(C_0: float, F_t: float):
    # C = C_0 * F_t
        result = []
        C = C_0*F_t
        result.append(C)
        return C
    @staticmethod
    def eqn_2_36__C_0(C: float, F_t: float):
    # C = C_0 * F_t
        result = []
        C_0 = C/F_t
        result.append(C_0)
        return C_0
    @staticmethod
    def eqn_2_37__F_t(A: float, C: float, M: float, T: float):
    # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        F_t = 0.000681714375311032*C**2*M/(A*T)
        result.append(F_t)
        return F_t
    @staticmethod
    def eqn_2_37__M(A: float, C: float, F_t: float, T: float):
    # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        M = 1466.89*A*F_t*T/C**2
        result.append(M)
        return M
    @staticmethod
    def eqn_2_37__A(C: float, F_t: float, M: float, T: float):
    # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        A = 0.000681714375311032*C**2*M/(F_t*T)
        result.append(A)
        return A
    @staticmethod
    def eqn_2_37__C(A: float, F_t: float, M: float, T: float):
    # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        C = 38.3*sqrt(A*F_t*T/M)
        result.append(C)
        return C
    @staticmethod
    def eqn_2_37__T(A: float, C: float, F_t: float, M: float):
    # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        T = 0.000681714375311032*C**2*M/(A*F_t)
        result.append(T)
        return T


class PressMgmt:
    @staticmethod
    def eqn_3_1__BarometricPressure(Abs_Pressure: float, Vacuum: float):
    # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return BarometricPressure
    @staticmethod
    def eqn_3_1__Abs_Pressure(BarometricPressure: float, Vacuum: float):
    # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return Abs_Pressure
    @staticmethod
    def eqn_3_1__Vacuum(Abs_Pressure: float, BarometricPressure: float):
    # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return Vacuum
    @staticmethod
    def eqn_3_2__G(G_C: float, H: float, P: float, rho: float):
    # P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return G
    @staticmethod
    def eqn_3_2__G_C(G: float, H: float, P: float, rho: float):
    # P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return G_C
    @staticmethod
    def eqn_3_2__P(G: float, G_C: float, H: float, rho: float):
    # P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return P
    @staticmethod
    def eqn_3_2__H(G: float, G_C: float, P: float, rho: float):
    # P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return H
    @staticmethod
    def eqn_3_2__rho(G: float, G_C: float, H: float, P: float):
    # P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return rho
    @staticmethod
    def eqn_3_3__P_P(H_1: float, H_2: float, P: float):
    # P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 - P
        result.append(P_P)
        return P_P
    @staticmethod
    def eqn_3_3__H_1(H_2: float, P: float, P_P: float):
    # P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 - P - P_P
        result.append(H_1)
        return H_1
    @staticmethod
    def eqn_3_3__H_2(H_1: float, P: float, P_P: float):
    # P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 + P + P_P
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_3__P(H_1: float, H_2: float, P_P: float):
    # P_P - P = H_2 - H_1
        result = []
        P = -H_1 + H_2 - P_P
        result.append(P)
        return P
    @staticmethod
    def eqn_3_4__KAPPA(P: float, V: float):
    # P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return KAPPA
    @staticmethod
    def eqn_3_4__V(KAPPA: float, P: float):
    # P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return V
    @staticmethod
    def eqn_3_4__P(KAPPA: float, V: float):
    # P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return P
    @staticmethod
    def eqn_3_5__V_P(P: float, P_P: float, V: float):
    # P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return V_P
    @staticmethod
    def eqn_3_5__P_P(P: float, V: float, V_P: float):
    # P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return P_P
    @staticmethod
    def eqn_3_5__V(P: float, P_P: float, V_P: float):
    # P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return V
    @staticmethod
    def eqn_3_5__P(P_P: float, V: float, V_P: float):
    # P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return P
    @staticmethod
    def eqn_3_6__H_1(H_2: float, P: float, V: float, V_P: float):
    # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return H_1
    @staticmethod
    def eqn_3_6__H_2(H_1: float, P: float, V: float, V_P: float):
    # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_6__P(H_1: float, H_2: float, V: float, V_P: float):
    # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return P
    @staticmethod
    def eqn_3_6__V_P(H_1: float, H_2: float, P: float, V: float):
    # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return V_P
    @staticmethod
    def eqn_3_6__V(H_1: float, H_2: float, P: float, V_P: float):
    # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return V
    @staticmethod
    def eqn_3_8__V_P(A_C: float, H_2: float):
    # V_P = A_C * H_2
        result = []
        V_P = A_C*H_2
        result.append(V_P)
        return V_P
    @staticmethod
    def eqn_3_8__A_C(H_2: float, V_P: float):
    # V_P = A_C * H_2
        result = []
        A_C = V_P/H_2
        result.append(A_C)
        return A_C
    @staticmethod
    def eqn_3_8__H_2(A_C: float, V_P: float):
    # V_P = A_C * H_2
        result = []
        H_2 = V_P/A_C
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_9__A_C(H_1: float, H_2: float, P: float, V: float):
    # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        A_C = P*V/(H_2*(-H_1 + H_2 + P))
        result.append(A_C)
        return A_C
    @staticmethod
    def eqn_3_9__H_1(A_C: float, H_2: float, P: float, V: float):
    # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_1 = H_2 + P - P*V/(A_C*H_2)
        result.append(H_1)
        return H_1
    @staticmethod
    def eqn_3_9__H_2(A_C: float, H_1: float, P: float, V: float):
    # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_2 = (A_C*(H_1 - P) - sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        H_2 = (A_C*(H_1 - P) + sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_9__P(A_C: float, H_1: float, H_2: float, V: float):
    # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        P = A_C*H_2*(H_1 - H_2)/(A_C*H_2 - V)
        result.append(P)
        return P
    @staticmethod
    def eqn_3_9__V(A_C: float, H_1: float, H_2: float, P: float):
    # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        V = A_C*H_2*(-H_1 + H_2 + P)/P
        result.append(V)
        return V
    @staticmethod
    def eqn_3_11__A_C(H_2: float, P: float, V: float):
    # P = A_C / V * (H_2) ** 2
        result = []
        A_C = P*V/H_2**2
        result.append(A_C)
        return A_C
    @staticmethod
    def eqn_3_11__V(A_C: float, H_2: float, P: float):
    # P = A_C / V * (H_2) ** 2
        result = []
        V = A_C*H_2**2/P
        result.append(V)
        return V
    @staticmethod
    def eqn_3_11__H_2(A_C: float, P: float, V: float):
    # P = A_C / V * (H_2) ** 2
        result = []
        H_2 = -sqrt(P*V/A_C)
        result.append(H_2)
        H_2 = sqrt(P*V/A_C)
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_11__P(A_C: float, H_2: float, V: float):
    # P = A_C / V * (H_2) ** 2
        result = []
        P = A_C*H_2**2/V
        result.append(P)
        return P
    @staticmethod
    def eqn_3_12__KAPPA_1(H_2: float, P: float):
    # P = KAPPA_1 * H_2 ** 2
        result = []
        KAPPA_1 = P/H_2**2
        result.append(KAPPA_1)
        return KAPPA_1
    @staticmethod
    def eqn_3_12__H_2(KAPPA_1: float, P: float):
    # P = KAPPA_1 * H_2 ** 2
        result = []
        H_2 = -sqrt(P/KAPPA_1)
        result.append(H_2)
        H_2 = sqrt(P/KAPPA_1)
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_12__P(H_2: float, KAPPA_1: float):
    # P = KAPPA_1 * H_2 ** 2
        result = []
        P = H_2**2*KAPPA_1
        result.append(P)
        return P
    @staticmethod
    def eqn_3_13__KAPPA_2(H_1: float, H_2: float, P: float):
    # P = KAPPA_2 * (H_2 - H_1)
        result = []
        KAPPA_2 = -P/(H_1 - H_2)
        result.append(KAPPA_2)
        return KAPPA_2
    @staticmethod
    def eqn_3_13__H_1(H_2: float, KAPPA_2: float, P: float):
    # P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_1 = H_2 - P/KAPPA_2
        result.append(H_1)
        return H_1
    @staticmethod
    def eqn_3_13__H_2(H_1: float, KAPPA_2: float, P: float):
    # P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_2 = H_1 + P/KAPPA_2
        result.append(H_2)
        return H_2
    @staticmethod
    def eqn_3_13__P(H_1: float, H_2: float, KAPPA_2: float):
    # P = KAPPA_2 * (H_2 - H_1)
        result = []
        P = KAPPA_2*(-H_1 + H_2)
        result.append(P)
        return P
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
