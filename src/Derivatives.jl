using Symbolics, LinearAlgebra, DifferentialEquations

# ================== CUSTOM TYPES ================== #

# FINITE DIFFERENCE
abstract type Difference end
struct ComplexDifference <: Difference end
struct CentralDifference <: Difference 
	points::Int64
	derivative::Int64
	function CentralDifference(points=3, derivative=1)
		return new(points, derivative)
	end
end
struct ForwardDifference <: Difference 
	points::Int64
	derivative::Int64
	function ForwardDifference(points=3, derivative=1)
		return new(points, derivative)
	end
end
struct BackwardDifference <: Difference 
	points::Int64
	derivative::Int64
	function BackwardDifference(points=3, derivative=1)
		return new(points, derivative)
	end
end

# ================== CONSTANTS ================== #
# CENTRAL DIFFERENCE QUERIES
cntr3pt = Dict( 1 => (f, x, δx) -> (-f(x - δx) + f(x + δx))/(2*δx),        # Accuracy: δx²
				2 => (f, x, δx) -> (f(x - δx) -2*f(x) + f(x + δx))/(δx^2)) # Accuracy: δx²

cntr5pt = Dict( 1 => (f, x, δx) -> (f(x - 2*δx) - 8*f(x - δx) + 8*f(x + δx) - f(x + 2*δx))/(12*δx),                # Accuracy: δx⁴
				2 => (f, x, δx) -> (-f(x - 2*δx) + 16*f(x - δx) - 30*f(x) + 16*f(x + δx) - f(x + 2*δx))/(12*δx^2), # Accuracy: δx⁴
				3 => (f, x, δx) -> (-f(x - 2*δx) + 2*f(x - δx) -2*f(x + δx) + f(x + 2*δx))/(2*δx^3))               # Accuracy: δx²

cntr7pt = Dict( 1 => (f, x, δx) -> (-f(x - 3*δx) + 9*f(x - 2*δx) - 45*f(x - δx) + 45*f(x + δx) - 9*f(x + 2*δx) + f(x + 3*δx))/(60*δx),                     # Accuracy: δx⁶
				2 => (f, x, δx) -> (2*f(x - 3*δx) - 27*f(x - 2*δx) + 270*f(x - δx) -490*f(x) + 270*f(x + δx) - 27*f(x + 2*δx) + 2*f(x + 3*δx))/(180*δx^2), # Accuracy: δx⁶
				3 => (f, x, δx) -> (f(x - 3*δx) - 8*f(x - 2*δx) + 13*f(x - δx) - 13*f(x + δx) + 8*f(x + 2*δx) - f(x + 3*δx))/(8*δx^3))                     # Accuracy: δx⁴

cntrDiff = Dict(3 => cntr3pt, 5 => cntr5pt, 7 => cntr7pt)

# FORWARD DIFFERENCE


# BACKWARD DIFFERENCE

# ================== FUNCTIONS ================== #

# FINITE DIFFERENCES
# function gradient(f::Function, x::Vector, h::Union{Float64, Int64}, method::Difference = ComplexDifference())
function gradient(f, x::Vector, h::Union{Float64, Int64}, method = ComplexDifference())
	return _gradient(method, f, x, h)
end

function _gradient(method::ComplexDifference, f, x, h)
	∂f = []
	L = length(x)
	for i = 1:L
		δx = zeros(Complex{Float64}, L)
		δx[i] = h*im
		push!(∂f, imag(f(x + δx))/h)
	end
	return ∂f
end

# TODO: Make work with 2nd+ order derivatives (ie: 1x1 )
function _gradient(method::CentralDifference, f, x, h)
	∂ⁿf = []
	L = length(x)
	g = cntrDiff[method.points][method.derivative]
	for i = 1:L
		δx = zeros(L)
		δx[i] = h
		push!(∂ⁿf, g(f, x, δx))
	end
	return ∂ⁿf
end

# TODO: Make Forward and Backwards Differences

# ------------- OLD CODE ------------- #

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

export jacobian, hessian, cntrDiff, gradient
