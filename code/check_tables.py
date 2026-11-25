#!/usr/bin/env python3
"""Companion code for "Congruences modulo powers of 5 for broken 12-diamond
partitions" (Ramanuja).

Recomputes every table in Appendix A of the paper from the five polynomials
e_1..e_5 of the modular equation (Theorem 5.1) alone, using exact integer
arithmetic and only the Python standard library, and re-runs every finite
check used in the proofs.

  python3 check_tables.py            # tables and all finite checks (seconds)
  python3 check_tables.py --qseries  # also: q-series checks of Lemma 4.1 and
                                     # Theorem 5.1, and numerical tests of
                                     # Theorems 1.1-1.3 (about 5 minutes)

The numerical tests are a sanity check only; they are not used in the paper.
"""
import sys

P6 = 5**6

def v5(n):
    if n == 0:
        return float('inf')
    v = 0
    while n % 5 == 0:
        n //= 5; v += 1
    return v

def w(n): return n // 5

# ---------- polynomials in y: lists of ints, index = degree ----------
def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a
def padd(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def psc(c, a): return trim([c*v for v in a])
def pmul(a, b):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                r[i+j] += u*v
    return trim(r)
def coef(p, i): return p[i] if 0 <= i < len(p) else 0

E = {1: [0, 205, 4300, 34000, 120000, 160000],
     2: [0, -215, -4475, -35000, -122000, -160000],
     3: [0, 85, 1750, 13525, 46500, 60000],
     4: [0, -15, -305, -2325, -7875, -10000],
     5: [0, 1, 20, 150, 500, 625]}

ok = True
def check(cond, msg):
    global ok
    print(('  ok   ' if cond else '  FAIL ') + msg)
    ok = ok and cond

print('Lemma 5.2 (norms):')
n5 = [1]; n4 = [1]
for l in range(1, 6):
    n5 = padd(n5, psc(5**l, E[l])); n4 = padd(n4, psc(4**l, E[l]))
x5 = [1]
for _ in range(5): x5 = pmul(x5, [1, 5])
check(n5 == x5, '1 + sum 5^l e_l = (1+5y)^5')
check(n4 == [1, 4], '1 + sum 4^l e_l = 1+4y')

print('Lemma 6.2 (valuations of e_l):')
check(all(v5(E[l][k]) >= -(-(5*k-l)//6) for l in E for k in range(1, 6)),
      'v5([y^k]e_l) >= ceil((5k-l)/6)')

# ---------- b_j = U5(y^j) by Newton's identities ----------
J = 40
p = [[5]]
for j in range(1, J+1):
    s = []
    for l in range(1, min(j, 5)+1):
        term = pmul(E[l], p[j-l]) if l < j else psc(j, E[l])
        s = padd(s, psc((-1)**(l+1), term))
    p.append(s)
assert all(c % 5 == 0 for q in p for c in q)
b = [[c//5 for c in q] for q in p]

# ---------- bivariate arithmetic mod Q(Y;y) ----------
Qc = [psc(-1, E[5]), E[4], psc(-1, E[3]), E[2], psc(-1, E[1]), [1]]
def bmul(A, B):
    r = [[] for _ in range(len(A)+len(B)-1)]
    for i, u in enumerate(A):
        for j, v in enumerate(B):
            r[i+j] = padd(r[i+j], pmul(u, v))
    return r
def red(A):
    A = [list(t) for t in A]
    while len(A) > 5:
        lc = A.pop(); d = len(A)
        for k in range(5):
            A[d-5+k] = padd(A[d-5+k], pmul(lc, psc(-1, Qc[k])))
    while len(A) < 5: A.append([])
    return A
def div_linear(num, c):
    d = len(num)-1; q = [[] for _ in range(d)]
    rem = [list(t) for t in num]
    for i in range(d, 0, -1):
        assert all(v % c == 0 for v in rem[i])
        q[i-1] = [v//c for v in rem[i]]
        rem[i-1] = padd(rem[i-1], psc(-1, q[i-1]))
    assert trim(rem[0]) == []
    return q
V4 = div_linear([padd(psc(1024, Qc[0]), [1, 4])] + [psc(1024, t) for t in Qc[1:]], 4)
def U5(Pb):
    R = red(Pb); s = []
    for i in range(5): s = padd(s, pmul(R[i], b[i]))
    return s
def Ypow(k): return [[] for _ in range(k)] + [[1]]
a = [U5(bmul(Ypow(j+2), V4)) for j in range(5)]
for j in range(5, J+1):           # recurrence (4.3)
    s = []
    for l in range(1, 6):
        s = padd(s, psc((-1)**(l+1), pmul(E[l], a[j-l])))
    a.append(s)

print('Appendix A.1 and Lemma 6.3 (initial cases):')
check(a[0] == [0, 10, 255, 2560, 12800, 32000, 32000], 'a_0 as printed')
check(all(v5(coef(b[j], i)) >= (5*i-j-1)//6 for j in range(5) for i in range(len(b[j]))),
      'v5([y^i]b_j) >= floor((5i-j-1)/6), j<=4')
check(all(v5(coef(a[j], i)) >= (5*i-j-7)//6 for j in range(5) for i in range(len(a[j]))),
      'v5([y^i]a_j) >= floor((5i-j-7)/6), j<=4')
check(all(v5(coef(b[j], i)) >= (5*i-j-1)//6 for j in range(J+1) for i in range(len(b[j]))),
      '  (consistency) same bound for b_j, j<=%d' % J)
check(all(v5(coef(a[j], i)) >= (5*i-j-7)//6 for j in range(J+1) for i in range(len(a[j]))),
      '  (consistency) same bound for a_j, j<=%d' % J)

# ---------- y-series ----------
K = 25
def series_div_x(poly, m, K):     # poly * (1+5y)^(-m), truncated
    inv = [1] + [0]*K
    for _ in range(m):
        inv = [sum(inv[t-s]*(-5)**s for s in range(t+1)) for t in range(K+1)]
    return [sum(coef(poly, i)*inv[t-i] for i in range(t+1)) for t in range(K+1)]
alpha = {j: series_div_x(a[j], 2, K) for j in range(6)}
beta = lambda k, i: coef(b[i], k)
s_poly = pmul([1, 4], b[1])
check(s_poly == [0, 41, 1024, 10240, 51200, 128000, 128000], 's = (1+4y) b_1 as printed')
sigma = series_div_x(s_poly, 1, K)
i0 = lambda j: -(-(j+2)//5)

def c(k, j):
    return sum(alpha[j][i]*beta(k, i) for i in range(i0(j), 5*k+1))
L1 = {k: sum(sigma[i]*beta(k, i) for i in range(1, 5*k+1)) for k in range(1, 6)}

print('Table 2 (beta_{k,i} mod 5^6):')
for k in range(1, 6): print('  k=%d' % k, [beta(k, i) % P6 for i in range(1, 26)])
print('Table 3 (alpha_{i,j} mod 5^6):')
for j in range(6): print('  j=%d' % j, [alpha[j][i] % P6 for i in range(1, 11)])
print('Table 4 (sigma_i mod 5^6):', [sigma[i] % P6 for i in range(1, 11)])

print('Proposition 6.5 (finite part of the contraction):')
for j in range(1, 5):
    check(v5(c(1, j)) >= 1, 'c_{1,%d} = %d is divisible by 5' % (j, c(1, j)))

print('Lemma 7.1 (L_1):')
check(L1[1] == 656835 and L1[1] % 5 == 0, '[y]L_1 = 656835 = 5*131367')
pairs = [(k, j) for k in range(1, 6) for j in range(1, 6)
         if (k, j) in {(1,1),(1,2),(1,3),(1,4),(1,5),(2,1),(2,2),(2,3),(2,4),(2,5),
                       (3,1),(3,2),(3,3),(3,4),(4,1),(4,2),(5,1),(5,2)}]
print('Table 5 (c_{k,j} mod 5^6):', {pq: c(*pq) % P6 for pq in pairs})
print('Table 6 (L_1 coefficients mod 5^6):', {k: L1[k] % P6 for k in L1})

lam = {k: (L1[k]//5**(1+w(k))) % 5 for k in (1, 2)}
chat = lambda k, j: (c(k, j)//5**(1+w(k)-w(j))) % 5
check(lam == {1: 2, 2: 3}, 'lambda = (2, 3)')
print('Proposition 7.2 (internal congruence):')
for k in (1, 2, 3):
    lhs = sum(chat(k, j)*lam[j] for j in (1, 2)) % 5
    check(lhs == (-lam.get(k, 0)) % 5, '(Cbar lambda)_%d = -lambda_%d' % (k, k))
print('Lemma 8.1 and Proposition 8.2 (three-term congruence):')
Dnum = {k: sum(L1[j]*c(k, j) for j in range(1, 6) if (k, j) in pairs) + 5*L1[k] for k in range(1, 6)}
delta = {}
for k in range(1, 6):
    m = 5**(4+w(k))
    assert Dnum[k] % 5**(3+w(k)) == 0
    delta[k] = (Dnum[k] % m)//5**(3+w(k))
check(delta == {1: 0, 2: 0, 3: 1, 4: 4, 5: 2}, 'delta = (0,0,1,4,2)')
for k in (1, 2, 3):
    lhs = sum(chat(k, j)*delta[j] for j in (3, 4, 5)) % 5
    check(lhs == lam.get(k, 0), '(Cbar delta)_%d = lambda_%d' % (k, k))

if '--qseries' in sys.argv:
    N = 260
    def eta(dexps, N=N):
        f = [0]*(N+1); f[0] = 1
        for d, e in dexps:
            for m in range(1, N//d+1):
                s = d*m
                for _ in range(abs(e)):
                    if e > 0:
                        for n in range(N, s-1, -1): f[n] -= f[n-s]
                    else:
                        for n in range(s, N+1): f[n] += f[n-s]
        return f
    def mul(u, v, N=N):
        r = [0]*(N+1)
        for i, x in enumerate(u):
            if x:
                for j in range(N+1-i): r[i+j] += x*v[j]
        return r
    sh = lambda u, k: ([0]*k+u)[:N+1]
    one = [1]+[0]*N
    y = sh(eta([(2, 1), (10, 3), (1, -3), (5, -1)]), 1)
    x = eta([(2, 5), (5, 1), (1, -5), (10, -1)])
    g = eta([(5, 3), (1, 1), (10, -1), (2, -3)])
    h = sh(eta([(1, -2), (5, -6), (10, 8)]), 2)
    k_ = eta([(5, 2), (1, 6), (2, -8)])
    p4 = [o+4*t for o, t in zip(one, y)]
    print('q-series checks of Lemma 4.1 and Theorem 5.1 to q^%d:' % N)
    check(x == [o+5*t for o, t in zip(one, y)], 'x = 1+5y')
    check(mul(g, x) == p4, 'g x = 1+4y')
    check(mul(h, p4) == mul(y, y), 'h (1+4y) = y^2')
    check(mul(k_, mul(x, x)) == p4, 'k x^2 = 1+4y')
    t = [0]*(N+1)
    for n in range(N//5+1): t[5*n] = y[n]
    tp = [one]
    for _ in range(5): tp.append(mul(tp[-1], t))
    yp = [one]
    for _ in range(5): yp.append(mul(yp[-1], y))
    def ev(poly):
        r = [0]*(N+1)
        for kk, cc in enumerate(poly):
            if cc: r = [u+cc*v for u, v in zip(r, tp[kk])]
        return r
    Qv = list(yp[5])
    for l, sg in ((1, -1), (2, 1), (3, -1), (4, 1), (5, -1)):
        Qv = [u+sg*v for u, v in zip(Qv, mul(ev(E[l]), yp[5-l]))]
    check(all(v == 0 for v in Qv), 'Q(y(tau); y(5 tau)) = 0')
    # Theorems 1.1, 1.2 for small alpha
    NN = 15625*2 + 14324
    print('computing Delta_12 up to n = %d ...' % NN)
    F = eta([(2, 1), (25, 1), (1, -3), (50, -1)], NN)
    r = lambda al: (11*25**al+13)//12
    H = {al: [F[25**al*n+r(al)] for n in range((NN-r(al))//25**al+1)] for al in (1, 2, 3)}
    check(all(v5(H[al][n]) >= al for al in H for n in range(len(H[al]))), '(1.2) for alpha<=3')
    check(all(v5(H[al+1][n]+5*H[al][n]) >= al+2 for al in (1, 2) for n in range(len(H[al+1]))), '(1.3) for alpha<=2')
    check(all(v5(H[3][n]+5*H[2][n]-125*H[1][n]) >= 5 for n in range(len(H[3]))), '(1.4) for alpha=2')
    cser = eta([(2, 3), (5, 2), (10, 1), (1, -4)], 2000)
    dser = eta([(2, 8), (1, -2)], 2000)
    check(all((H[1][n]//5 - 2*cser[n]) % 5 == 0 for n in range(len(H[1]))), 'Theorem 1.3(a), alpha=1')
    G = [(H[2][n]+5*H[1][n])//125 for n in range(len(H[2]))]
    check(all(G[n] % 5 == 0 for n in range(len(G)) if n % 5 != 2) and
          all((G[5*m+2]-dser[m]) % 5 == 0 for m in range((len(G)-3)//5+1)), 'Theorem 1.3(b)')

print('\nALL CHECKS PASSED' if ok else '\nSOME CHECK FAILED')
