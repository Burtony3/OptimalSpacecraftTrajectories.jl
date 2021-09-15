using Symbolics, LinearAlgebra, DifferentialEquations

@variables r μ x y z

# SETUP
g = -μ/norm([x, y, z])^3 * [x, y, z]
jac = Array{Any}(undef, 3, 3)
hes = Array{Any}(undef, 3, 3, 3)
sub = Dict(sqrt(abs2(x) + abs2(y) + abs2(z)) => r);
D = [Differential(x), Differential(y), Differential(z)]

# CALCULATING PARTIALS
for i = 1:3, j = 1:3, k = 1:3
	
	# n×n×n HESSIAN
	hes[k, j, i] = substitute(
		expand_derivatives(D[i](-expand_derivatives(D[j](g[k])))), 
		sub)

	# n×n JACOBIAN
	jac[k, j] = substitute(
		expand_derivatives(D[j](g[k])), 
		sub)
end

# CREATING FUNCTIONS
_buildJacobian = build_function(jac, [x, y, z, μ, r], expression=false)[1]
jacobian(r⃗, μ₀) = _buildJacobian([r⃗[1], r⃗[2], r⃗[3], μ₀, norm(r⃗)])
_buildHessian  = build_function(hes, [x, y, z, μ, r], expression=false)[1]
hessian(r⃗, μ₀) = _buildHessian([r⃗[1], r⃗[2], r⃗[3], μ₀, norm(r⃗)])

export jacobian, hessian
