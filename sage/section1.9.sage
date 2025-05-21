load("utils1.9.sage")

var("t", domain=RR)

r, theta = literal_function("r"), literal_function("theta")
show(r)

show((r(t), theta(t)), simplify=False)

q_prime = path_function([r, theta])
show(q_prime(t), simplify=False)

show(Gamma(q_prime))

# load("show_expression.sage")
show(Gamma(q_prime)(t), simplify=False)

F = p_to_r
show(compose(F, Gamma(q_prime))(t), simplify=False)

Q = lambda t: compose(p_to_r, Gamma(q_prime))(t)
show(Gamma(Q)(t), simplify=False)

show(f_bar(q_prime)(t))

t_prime = var("tt", domain="positive", latex_name="t'")
q = path_function([literal_function("r"), literal_function("theta")])
local = Gamma(q)(t)
show(osculating_path(local)(t_prime))

show(Gamma_bar(f_bar)(local))

show(F_to_C(p_to_r)(local))

q = path_function([literal_function("x")])
local = Gamma(q, 4)(t)
show(local)

m, k = var("m k", domain="positive")
L = L_harmonic(m, k)
show(Euler_Lagrange_operator(L)(local))
