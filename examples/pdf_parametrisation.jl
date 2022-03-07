# # Input PDF parametrisation and priors
#
# An important part of PDF fitting is defining a useful parametrisation for the
# PDF shapes, as well as meaningful prior distributions that encode our knowledge
# of the problem.
#
# In this notebook, and in the `PartonDensity` package, we explore two different approaches:
# * *Full Dirichlet*
# * *Valence shape + Dirichlet*
#
# We are still investigating which approach will be most
# effective, and so both options are currently implemented,
# as detailed below.

using Distributions, Plots, SpecialFunctions, Printf
const sf = SpecialFunctions;

# ## "Full Dirichlet" approach
#
# A clean way to ensure the momentum sum rule would be to sample different
# contributions of the momentrum density integral from a Dirichlet distribution,
# then use these weights to set the parameters on the individual Beta distributions.

# 9 components of decreasing importance

dirichlet = Dirichlet([3., 2., 1, 0.5, 0.3, 0.2, 0.1, 0.1, 0.1])
data = rand(dirichlet, 1000);

# Have a look

plot()
for i in 1:9
    histogram!(data[i,:], bins=range(0, stop=1, length=20), alpha=0.7)
end
plot!(xlabel="I_i = A_i B_i")

# This would be great as the sum rule is automatically conserved

sum(data, dims=1)

# But, it could be non-trivial to choose the Dirichlet weights for a
# sensible prior, and connect to the physics of the problem.

I = rand(dirichlet)

# As we are more interested in `K_u` and `K_d` than `λ_u` and `λ_d`
# for our high-x data, we can set the `λ`s according to the `K`s.

K_u = rand(Uniform(2, 10))

# Integral of number density must = 2, and integral of
# momentum density must = `I[1]`. This implies the following:

λ_u = (I[1] * (K_u + 1)) / (2 - I[1])

# Let's check this

A_u = 2 / sf.beta(λ_u, K_u + 1)

I[1] ≈ A_u * sf.beta(λ_u + 1, K_u + 1)

# While this approach might be nice, it could be hard to set priors on the shape
# of the valence distributions, as the `λ_u` and `λ_d` are now derived parameters, dependent
# on the specified Dirichlet weights and `K_u`/`K_d` values. However, sensible prior choices
# could be made using e.g. prior predictive checks.

# ## "Valence shape + Dirichlet" approach
#
# An alternative approach is to specify constraints on the
# valence params through the shape of their Beta distributions, then using a
# Dirichlet to specify the weights of the gluon and sea components. The problem
# here is it isn't clear how to specify that the d contribution must be less than
# the u contribution, but it is possible to do this indirectly through priors on
# the shape parameters. This would however also require some further investigation.

x = range(0, stop=1, length=50)

# High-level priors
# Looks like we maybe want to change lambda and K priors to boost these components

λ_u = 0.7 #rand(Uniform(0, 1))
K_u = 4 #rand(Uniform(2, 10))
λ_d = 0.5 #rand(Uniform(0, 1))
K_d = 6 #rand(Uniform(2, 10))

u_V = Beta(λ_u, K_u+1)
A_u = 2 / sf.beta(λ_u, K_u+1)

d_V = Beta(λ_d, K_d+1)
A_d = 1 / sf.beta(λ_d, K_d+1)

# Integral contributions

I_u = A_u * sf.beta(λ_u+1, K_u+1)
I_d = A_d * sf.beta(λ_d+1, K_d+1)

plot(x, x .* A_u .* x.^λ_u .* (1 .- x).^K_u * 2, alpha=0.7, label="x u(x)", lw=3)
plot!(x, x .* A_d .* x.^λ_d .* (1 .- x).^K_d, alpha=0.7, label="x d(x)", lw=3)
plot!(xlabel="x", legend=:topright)

#

@printf("I_u = %.2f\n", I_u)
@printf("I_d = %.2f\n", I_d)

# The remaining 7 integrals can be dirichlet-sampled with decreasing importance

remaining = 1 - (I_u + I_d)
dirichlet = Dirichlet([3., 2., 1, 0.5, 0.3, 0.2, 0.1])
I = rand(dirichlet) * remaining;
sum(I) ≈ remaining

# Gluon contributions

λ_g1 = rand(Uniform(-1, 0))
λ_g2 = rand(Uniform(0, 1))
K_g = rand(Uniform(2, 10))
A_g2 = I[1] / sf.beta(λ_g2+1, K_g+1)
A_g1 = I[2] / sf.beta(λ_g1+1, 5+1);

# Sea quark contributions

λ_q = rand(Uniform(-1, 0))
A_ubar = I[3] / (2 * sf.beta(λ_q+1, 5+1))
A_dbar = I[4] / (2 * sf.beta(λ_q+1, 5+1))
A_s = I[5] / (2 * sf.beta(λ_q+1, 5+1))
A_c = I[6] / (2 * sf.beta(λ_q+1, 5+1))
A_b = I[7] / (2 * sf.beta(λ_q+1, 5+1));

total = A_u * sf.beta(λ_u+1, K_u+1) + A_d * sf.beta(λ_d+1, K_d+1) 
total += A_g1 * sf.beta(λ_g1+1, 5+1) + A_g2 * sf.beta(λ_g2+1, K_g+1)
total += 2 * (A_ubar + A_dbar + A_s + A_c + A_b) * sf.beta(λ_q+1, 5+1)
total ≈ 1

#

x = 10 .^ range(-2, stop=0, length=500)

# How does it look?

xg2 = A_g2 * x.^λ_g2 .* (1 .- x).^K_g
xg1 = A_g1 * x.^λ_g1 .* (1 .-x).^5
plot(x, x .* A_u .* x.^λ_u .* (1 .- x).^K_u * 2, alpha=0.7, label="x u(x)", lw=3)
plot!(x, x .* A_d .* x.^λ_d .* (1 .- x).^K_d, alpha=0.7, label="x d(x)", lw=3)
plot!(x, xg1 + xg2, alpha=0.7, label="x g(x)", lw=3)
plot!(x, A_ubar * x.^λ_q .* (1.0 .- x).^5, alpha=0.7, label="x ubar(x)", lw=3)
plot!(x, A_dbar * x.^λ_q .* (1.0 .- x).^5, alpha=0.7, label="x dbar(x)", lw=3)
plot!(x, A_s * x.^λ_q .* (1.0 .- x).^5, alpha=0.7, label="x s(x)", lw=3)
plot!(x, A_c * x.^λ_q .* (1.0 .- x).^5, alpha=0.7, label="x c(x)", lw=3)
plot!(x, A_b * x.^λ_q .* (1.0 .-  x).^5, alpha=0.7, label="x b(x)", lw=3)
plot!(xlabel="x", legend=:bottomleft, xscale=:log, ylims=(1e-8, 10), yscale=:log) 

# ### Prior predictive checks
#
# We can start to visualise the type of PDFs that are allowed by the
# combination of the choice of parametrisation and prior distributions
# with some simple prior predictive checks, as done below for the valence
# style parametrisation...

N = 100
alpha = 0.03
total = Array{Float64, 1}(undef, N)
first = true
leg = 0

plot()
for i in 1:N

    λ_u_i = rand(Uniform(0, 1))
    K_u_i = rand(Uniform(2, 10))
    λ_d_i = rand(Uniform(0, 1))
    K_d_i = rand(Uniform(2, 10))
    A_u_i = 2 / sf.beta(λ_u_i, K_u_i+1)
    A_d_i = 1 / sf.beta(λ_d_i, K_d_i+1)
    I_u_i = A_u * sf.beta(λ_u_i+1, K_u_i+1)
    I_d_i = A_d * sf.beta(λ_d_i+1, K_d_i+1)
    u_V_i = Beta(λ_u_i, K_u_i+1)
    d_V_i = Beta(λ_d_i, K_d_i+1)

    remaining_i = 1 - (I_u_i + I_d_i)
    dirichlet_i = Dirichlet([3., 2., 1, 0.5, 0.3, 0.2, 0.1])
    I_i = rand(dirichlet_i) * remaining_i
    
    λ_g1_i = rand(Uniform(-1, 0))
    λ_g2_i = rand(Uniform(0, 1))
    K_g_i = rand(Uniform(2, 10))
    A_g2_i = I_i[1] / sf.beta(λ_g2_i+1, K_g_i+1)
    A_g1_i = I_i[2] / sf.beta(λ_g1_i+1, 5+1)

    λ_q_i = rand(Uniform(-1, 0))
    A_ubar_i = I_i[3] / (2 * sf.beta(λ_q_i+1, 5+1))
    A_dbar_i = I_i[4] / (2 * sf.beta(λ_q_i+1, 5+1))
    A_s_i = I_i[5] / (2 * sf.beta(λ_q_i+1, 5+1))
    A_c_i = I_i[6] / (2 * sf.beta(λ_q_i+1, 5+1))
    A_b_i = I_i[7] / (2 * sf.beta(λ_q_i+1, 5+1))
    
    total[i] = A_u_i * sf.beta(λ_u_i+1, K_u_i+1) + A_d_i * sf.beta(λ_d_i+1, K_d_i+1) 
    total[i] += A_g1_i * sf.beta(λ_g1_i+1, 5+1) + A_g2_i * sf.beta(λ_g2_i+1, K_g_i+1)
    total[i] += 2 * (A_ubar_i + A_dbar_i + A_s_i + A_c_i + A_b_i) * sf.beta(λ_q_i+1, 5+1)
    
    xg2_i = A_g2_i * x.^λ_g2_i .* (1 .- x).^K_g_i
    xg1_i = A_g1_i * x.^λ_g1_i .* (1 .- x).^5
    plot!(x, [x .* A_u_i .* x.^λ_u_i .* (1 .- x).^K_u_i * 2], alpha=alpha, color="blue", lw=3)
    plot!(x, x .* A_d_i .* x.^λ_d_i .* (1 .- x).^K_d_i, alpha=alpha, color="orange", lw=3)
    plot!(x, xg1_i + xg2_i, alpha=alpha, color="green", lw=3)
    plot!(x, A_ubar_i * x.^λ_q_i .* (1 .- x).^5, alpha=alpha, color="red", lw=3)
    plot!(x, A_dbar_i * x.^λ_q_i .* (1 .- x).^5, alpha=alpha, color="purple", lw=3)
    plot!(x, A_s_i * x.^λ_q_i .* (1 .- x).^5, alpha=alpha, color="brown", lw=3)
    plot!(x, A_c_i * x.^λ_q_i .* (1 .- x).^5, alpha=alpha, color="pink", lw=3)
    plot!(x, A_b_i * x.^λ_q_i .* (1 .- x).^5, alpha=alpha, color="grey", lw=3)
end

plot!(xlabel="x", ylabel="x f(x)", xscale=:log, legend=false,
      ylims=(1e-8, 10), yscale=:log)

# Looks like naive priors need some work...

# ## PDF Parametrisation interface
#
# `PartonDensity` provides a handy interface to both the parametrisations through the
# `DirichletPDFParams` and `ValencePDFParams` options.

using PartonDensity

# Let's try the Dirichlet specification...

pdf_params = DirichletPDFParams(K_u=4.0, K_d=6.0, λ_g1=0.7, λ_g2=-0.4,
                         K_g=6.0, λ_q=-0.5, seed=5,
                         weights=[3., 1., 5., 5., 1., 1., 1., 0.5, 0.5]);

plot_input_pdfs(pdf_params)

#

int_xtotx(pdf_params) ≈ 1

# And now the valence specification...

pdf_params = ValencePDFParams(λ_u=0.7, K_u=4.0, λ_d=0.5, K_d=6.0, λ_g1=0.7, λ_g2=-0.4,
                              K_g=6.0, λ_q=-0.5, seed=5, weights=[5., 5., 1., 1., 1., 0.5, 0.5]);

plot_input_pdfs(pdf_params)

#

int_xtotx(pdf_params) ≈ 1
