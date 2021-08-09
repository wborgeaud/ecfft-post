# p is the size of the base field of the curve Secp256k1
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
F = GF(p)
# See the post on how to find a and b such that 2^12 divides the order of E.
a, b = 1, 641632526086088189863104799840287255425991771106608433941469413117059106896
E = EllipticCurve(F, [a,b])

log_n = 12
n = 2^log_n
assert E.order() % n == 0
g = E.gens()[0]
G = (g.order()//n) * g
assert G.order() == n


R = E.random_element()
H = [R + i*G for i in range(2^log_n)]
L = [h.xy()[0] for h in H]
S = [L[i] for i in range(0, n, 2)]
S_prime = [L[i] for i in range(1, n, 2)]

def precompute(log_n, S, S_prime, E):
    Ss = {}
    Ss_prime = {}
    matrices = {}
    inverse_matrices = {}
    for i in range(log_n, -1, -1):
        n = 1 << i
        nn = n // 2

        Ss[n] = S
        Ss_prime[n] = S_prime
        matrices[n] = []
        inverse_matrices[n] = []

        for iso in E.isogenies_prime_degree(2):
            psi = iso.x_rational_map()
            if len(set([psi(x) for x in S]))==nn:
                break
        v = psi.denominator()
        q = nn - 1
        for j in range(nn):
            s0, s1 = S[j], S[j + nn]
            assert psi(s0) == psi(s1)
            M = Matrix(F, [[v(s0)^q,s0*v(s0)^q],[v(s1)^q, s1*v(s1)^q]])
            inverse_matrices[n].append(M.inverse())

            s0, s1 = S_prime[j], S_prime[j + nn]
            assert psi(s0) == psi(s1)
            M = Matrix(F, [[v(s0)^q,s0*v(s0)^q],[v(s1)^q, s1*v(s1)^q]])
            matrices[n].append(M)

        S = [psi(x) for x in S[:nn]]
        S_prime = [psi(x) for x in S_prime[:nn]]
        E = iso.codomain()

    return Ss, Ss_prime, matrices, inverse_matrices

# Precompute the data needed to compute EXTEND_S,S'
Ss, Ss_prime, matrices, inverse_matrices = precompute(log_n-1, S, S_prime, E)
        
def extend(P_evals):
    n = len(P_evals)
    nn = n // 2
    if n == 1:
        return P_evals
    S = Ss[n]
    S_prime = Ss_prime[n]
    P0_evals = []
    P1_evals = []
    for j in range(nn):
        s0, s1 = S[j], S[j + nn]
        y0, y1 = P_evals[j], P_evals[j + nn]
        Mi = inverse_matrices[n][j]
        p0, p1 = Mi * vector([y0, y1])
        P0_evals.append(p0)
        P1_evals.append(p1)

    P0_evals_prime = extend(P0_evals)
    P1_evals_prime = extend(P1_evals)

    ansL = []
    ansR = []
    for M, p0, p1 in zip(matrices[n], P0_evals_prime, P1_evals_prime):
        v = M * vector([p0, p1]) 
        ansL.append(v[0])
        ansR.append(v[1])
    return ansL + ansR

# Generate a random polynomial for testing
R.<X> = F[]
P = sum(F.random_element() * X^i  for i in range(1<<(log_n - 1)))

# Evaluate P on S
P_evals = [P(x) for x in S]
# result holds the evaluation of P on S'
result = extend(P_evals)
assert result == [P(x) for x in S_prime]
