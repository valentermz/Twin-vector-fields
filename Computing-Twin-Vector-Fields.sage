
# coding: utf-8

# # Computation of Twins for Quadratic Vector Fields on $\mathbb{C}^2$
#
# ## Introduction
#
# Let $v=P(x,y)\frac{\partial}{\partial x}+Q(x,y)\frac{\partial}{\partial y}$ be a quadratic vector field on $\mathbb{C}^2$. In the generic case $v$ has 4 isolated singularities, we consider only such vector fields.
#
# ### Theorem:
#
# Given a generic quadratic field $v$, there exists a unique quadratic vector field $v^{t}$ which has the same singular set as $v$, at each singular point the same spectrum, yet the vector fields are different. We call these *"twin vector fields"*.
#
# Below we compute $v^{t}$ in terms of the polynomials $P$ and $Q$ defining $v$.
#
# ## Preliminaries
#
# By an affine change of coordinates we may assume $v$ has singularities at $p_0,p_1,p_2$ given by $(0,0)$, $(1,0)$ and $(0,1)$. The position of the fourth singularity $p_3$ is unknown, but may be recovered from the Euler-Jacobi formulae.
#
# We denote by $t_k$ and $d_k$ the trace and determinant of the linear
# part of $v$ at $p_k$.

# ### Variables

# In[1]:

# Coordinates on C2
var('x y')

# Parameters defining P, Q
for j in range(6):
    var('a{}'.format(j))

# Traces
t = [''] * 3

# Determinants
d = [''] * 3

# Auxiliary variables
m0, m1, m2, m3 = var('m0 m1 m2 m3')


# ### Definitions

# In[2]:

P = a0 * x ^ 2 + a1 * x * y + a2 * y ^ 2 - a0 * x - a2 * y
Q = a3 * x ^ 2 + a4 * x * y + a5 * y ^ 2 - a3 * x - a5 * y

v = (P, Q)

Sing = [[x == 0, y == 0], [x == 1, y == 0], [x == 0, y == 1]]

Dv = matrix([
    (diff(P, x), diff(P, y)),
    (diff(Q, x), diff(Q, y))
])

for k in range(3):
    t[k] = Dv.subs(Sing[k]).trace()
    d[k] = Dv.subs(Sing[k]).determinant()


# ### The trick
#
# Suppose $v' = (P', Q')$ has the same singular locus as $v$. By Max Noether's theorem there exist complex numbers $m_0,\ldots,m_3$ such that
#
# $\pmatrix{
#     P' \\
#     Q' \\
# }
# = \pmatrix{
#     m_0 & m_1 \\
#     m_2 & m_3 \\
# }
# \pmatrix{
#     P \\
#     Q \\
# }$.
#
# Thus $v$ and $v'$ will be twins if and only if the following system of equations is satisfied:
#
# $\operatorname{tr}MDv(p_k) = t_k, \text{ for } k=1,2,3, \quad
# \operatorname{det}M = 1$.

# In[3]:

M = matrix([(m0, m1), (m2, m3)])

lin_syst = [''] * 3

for k in range(3):
    lin_syst[k] = (M * Dv).subs(Sing[k]).trace() - t[k]

quad_eq = M.determinant() - 1


# We will first solve the linear system with respect to $m_1,m_2,m_3$, and then we solve the remaining quadratic equation on $m_0$.
#
# ### The linear system

# In[4]:

lin_syst_sol = solve(lin_syst, [m1, m2, m3])[0]

for k in range(3):
    print lin_syst_sol[k], '\n'


# ### The quadratic equation

# In[5]:

quad_eq = quad_eq.subs(lin_syst_sol).simplify_rational().factor()
quad_eq


# Note how the denominator has a factor of $m_0-1$.

# In[6]:

quad_eq_sol = solve(quad_eq.numerator() / (m0 - 1), m0)[0]
quad_eq_sol


# ## The final expressions

# In[7]:

m0_final = m0.subs(quad_eq_sol).simplify_rational().factor()
m1_final = m1.subs(lin_syst_sol).subs(quad_eq_sol).simplify_rational().factor()
m2_final = m2.subs(lin_syst_sol).subs(quad_eq_sol).simplify_rational().factor()
m3_final = m3.subs(lin_syst_sol).subs(quad_eq_sol).simplify_rational().factor()

for k in range(4):
    print 'm{} ='.format(k), eval('m' + str(k) + '_final'), '\n'


# ### Sanity check
#
# We obtain a new vector field $v^t = P^t\frac{\partial}{\partial x} +
# Q^t\frac{\partial}{\partial y}$. Below we verify that what we obtain
# does indeed have the same spectra as $v$.

# In[8]:

Pt = m0_final * P + m1_final * Q
Qt = m2_final * P + m3_final * Q

Dvt = matrix([
    (diff(Pt, x), diff(Pt, y)),
    (diff(Qt, x), diff(Qt, y))
])

# Thse should print only zeroes
for k in range(3):
    print (Dvt.subs(Sing[k]).trace() - t[k]).simplify_rational()
    print (Dvt.subs(Sing[k]).determinant() - d[k]).simplify_rational()


# ## Export final expressions for $P^t$, $Q^t$

# In[9]:

Pt = Pt.simplify_rational().factor()
Qt = Qt.simplify_rational().factor()


# In[10]:

print Pt


# In[11]:

print Qt


# ### A remark
#
# All of the above works fine as long as the denominators don't vanish,
# the linear system has non-zero determinat, and the quadratic equation
# has two different roots. Thus, the explicit *genericity assumptions* can
# be recovered from the above computations.
