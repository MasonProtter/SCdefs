# /ssh:mprotter@cedar.computecanada.ca:/BCS_PH/SCdefs.jl 
module SCdefs

ϵ(k::Number) = @fastmath -2cos(k)
ϵ(k::Array)  = @fastmath -2*∑(cos, k)
ξ(k, μ) = @fastmath ϵ(k) - μ

E(k, μ, Δ) = sqrt(Esq(k, μ, Δ))
Esq(k, μ, Δ) = @fastmath (ξ(k, μ)^2 + Δ^2)

E(k, μ, Δ, Δ̄)  = @fastmath sqrt(Esq(k, μ, Δ, Δ̄))
Esq(k, μ, Δ, Δ̄)  = @fastmath (ξ(k, μ)^2 + Δ*Δ̄)


βξ(k, μ, β) = @fastmath β*ξ(k, μ)
βE(k, μ, Δ, β)  = @fastmath (β*E(k, μ, Δ))

βE(k, μ, Δ, Δ̄, β) = @fastmath (β*E(k, μ, Δ, Δ̄))

f(x) = (       1.0
	/#--------------
	  (exp(x) + 1.0))
# function f(x::Float64)
#     @fastmath (x < -20.0 ? 1.0 :
#                x >  20.0 ? 0.0 :
#                inv(exp(x) + 1.0))
# end

# function b(x::Float64)
#     @fastmath (x < -20.0 ? -1.0 :
#                x >  20.0 ?  0.0 :
#                inv(exp(x) - 1.0))
# end

b(x) = (       1.0
	/#--------------
	   (exp(x) - 1.0))

G₀(k, iωn, μ) = 1/(iωn - ξ(k, μ))
G(k, iωn, μ, Δ) = (iωn + ξ(k, μ))/(iωn^2 - Esq(k, μ, Δ))

G(k, iωn, μ, Δ, Δ̄) = (iωn + ξ(k, μ))/(iωn^2 - Esq(k, μ, Δ, Δ̄))
F(k, iωn, μ, Δ, Δ̄) = -Δ/(iωn^2 - Esq(k, μ, Δ, Δ̄))
F̄(k, iωn, μ, Δ, Δ̄) = -Δ̄/(iωn^2 - Esq(k, μ, Δ, Δ̄))

iω(n, β::T) where {T<:Number} = @fastmath (im*π*(2n+1)/β)::Complex{T}
iΩ(m, β::T) where {T<:Number} = @fastmath (im*π*2m/β)::Complex{T}

∑ = sum

function generate_ks(N::Float64)
    nmax = (N-1.0)/2.0
    kmax = 2.0*π*nmax/N
    kstep = 2*π/N
    -kmax:kstep:kmax
end

generate_ks(N) = generate_ks(Float64(N))

export ∑, ϵ, ξ, E, Esq, βE, βξ, f, b, G₀, G, F, F̄, iω, iΩ, generate_ks

end
