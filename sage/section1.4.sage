import numpy as np

load(
    "utils.sage",
    "utils1.4.sage",
)

var('m', domain='positive')

q = vector(
    [
        literal_function("q_x"),
        literal_function("q_y"),
        literal_function("q_z"),
    ]
)

show(L_free_particle(m)(Gamma(q)(t)))

var("T", domain="positive")

show(latex(Lagrangian_action(L_free_particle(m), q, 0, T)))

test_path = vector([4 * t + 7, 3 * t + 5, 2 * t + 1])
Lagrangian_action(L_free_particle(mass=3), test_path, 0, 10)

hard_path = vector([4 * t + 7, 3 * t + 5, 2 * exp(-t) + 1])

result = Lagrangian_action(L_free_particle(mass=3), hard_path, 0, 10)
show(result)
show(float(result))

def make_eta(nu, t1, t2):
    return (t - t1) * (t - t2) * nu(t=t)

nu = vector([sin(t), cos(t), t ^ 2])
show(1 / 3 * make_eta(nu, 3, 4)  + test_path)

def varied_free_particle_action(mass, q, nu, t1, t2):
    eta = make_eta(nu, t1, t2)

    def f(eps):
        return Lagrangian_action(L_free_particle(mass), q + eps * eta, t1, t2).n()

    return f

show(varied_free_particle_action(3.0, test_path, nu, 0.0, 10.0)(0.001))

res = find_local_minimum(
    varied_free_particle_action(3.0, test_path, nu, 0, 10), -2, 1
)
show(res)

ts = np.linspace(0, pi / 2, 5)
qs = [cos(t).n() for t in ts]
lp = Lagrangian_polynomial(ts, qs)
ts = np.linspace(0, pi / 2, 20)
Cos = [lp(x=t).n() for t in ts]
Sin = [-lp.derivative(x)(x=t).n() for t in ts]
Zero = [abs(Cos[i] ^ 2 + Sin[i] ^ 2 - 1) for i in range(len(ts))]
show(max(Zero))

def make_path(t0, q0, t1, q1, qs):
    ts = np.linspace(t0, t1, len(qs) + 2)
    qs = np.r_[q0, qs, q1]
    return vector([Lagrangian_polynomial(ts, qs)(t)])

def parametric_path_action(Lagrangian, t0, q0, t1, q1):
    def f(qs):
        path = make_path(t0, q0, t1, q1, qs=qs)
        return Lagrangian_action(Lagrangian, path, t0, t1)

    return f

t0, t1 = 0, pi / 2
q0, q1 = cos(t0), cos(t1)
ts = np.linspace(0, pi / 2, 5)
initial_qs = [cos(t).n() for t in ts][1:-1]
parametric_path_action(L_harmonic(m=1, k=1), t0, q0, t1, q1)(initial_qs)

def find_path(Lagrangian, t0, q0, t1, q1, n):
    ts = np.linspace(t0, t1, n)
    initial_qs = np.linspace(q0, q1, n)[1:-1]
    minimizing_qs = minimize(
        parametric_path_action(Lagrangian, t0, q0, t1, q1),
        initial_qs,
    )
    return make_path(t0, q0, t1, q1, minimizing_qs)

best_path = find_path(L_harmonic(m=1, k=1), t0=0, q0=1, t1=pi / 2, q1=0, n=5)
result = [
    abs(best_path(t=t).n()[0] - cos(t).n()) for t in np.linspace(0, pi / 2, 10)
]
show(max(result))

ts = np.linspace(0, pi / 2, 20)
q = vector([cos(t)])
lvalues = [L_harmonic(m=1, k=1)(Gamma(q)(t)).n() for t in ts]
points = list(zip(ts, lvalues))
plot = list_plot(points, color="black", size=30)
plot.axes_labels(["$t$", "$L$"])
plot.save("figures/Lagrangian.png", figsize=(4, 2))
