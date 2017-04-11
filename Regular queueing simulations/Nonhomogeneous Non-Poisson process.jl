# Written by Yongkyu Cho.
# Reference: Transforming renewal processes for simulation of nonstationary arrival processes, Gerhardt and Nelson (2009), INFORMS Journal on Computing
# I use the inversion method with Weibull(0.5,0.5) base distribution.

using Distributions, PyPlot, Roots, JuMP, Ipopt

a = 4.0; b = -3.0; c = 1000.0
λ(t) = a + b*sin((π/c)*t)
R(s,x) = a*(x-s)+((b*c)/π)*( cos(((π*s)/c))-cos(((π*x)/c)) ) # Note that R(0.0,x) = R(x)


function inverse_integrated_rate_function(f::Function, s::Float64, val::Float64)  # f: integrated rate function (obtained by hand), s: starting point
    return fzero(x -> f(s,x)-val , s)
end

function inverse_integrated_function(f::Function, val::Float64) # Numerically calculating the inverse of the integrated function
  x = 0.0
  while quadgk(f,0.0,x)[1] < val
    x += 0.00001
  end
  return x
end

function generate_NHWP(R::Function, Θ::Float64, α::Float64, N::Int64) # non-homogeneous Weibull process until N arrivals
  n = 1
  S = inverse_integrated_function(x->2*e^(-sqrt(2x)),rand()) # first sample from the Stationary Excess distribution
  V = [inverse_integrated_rate_function(R,0.0,S)]
  while n < N
    push!(V, inverse_integrated_rate_function(R,V[n],rand(Weibull(Θ,α))) )
    n += 1
  end
  return V
end

function generate_NHWP(R::Function, Θ::Float64, α::Float64, T::Float64) # non-homogeneous Weibull process until T times
  n = 1
  S = inverse_integrated_function( x->2*e^(-sqrt(2*x)) , rand())
  V = [inverse_integrated_rate_function(R,0.0,S)]
  while V[n] < T
    push!(V, inverse_integrated_rate_function(R,V[n],rand(Weibull(Θ,α))) )
    n += 1
  end
  return V
end

function generate_NHPP(λ::Function, T::Float64)
  m = Model(solver = IpoptSolver(print_level = 0))
  #JuMP.registerNLFunction(m, :λ, 1, λ, autodiff=true)
  JuMP.register(m, :λ, 1, λ, autodiff=true)
  @variable(m, t)
  @NLobjective(m, Max, λ(t))
  solve(m)
  λ_max = getobjectivevalue(m)
  x = Float64[]
  t = 0.0
  while t < T
    t -= (1/λ_max)*log(rand())
    if rand() <= λ(t)/λ_max
      push!(x, t)
    end
  end
  return x
end

function generate_NHPP(λ::Function, N::Int64)
  m = Model(solver = IpoptSolver(print_level = 0))
  #JuMP.registerNLFunction(m, :λ, 1, λ, autodiff=true)
  JuMP.register(m, :λ, 1, λ, autodiff=true)
  @variable(m, t)
  @NLobjective(m, Max, λ(t))
  solve(m)
  λ_max = getobjectivevalue(m)
  x = Float64[]
  t = 0.0
  n = 0
  while n < N
    t -= (1/λ_max)*log(rand())
    if rand() <= λ(t)/λ_max
      push!(x, t)
      n += 1
    end
  end
  return x
end

# plotting for verification
plt = PyPlot
range = 5000.0
for k in 1:1
  x = generate_NHWP(R,0.5,0.5,range)
  y = [1]
  for i in 1:length(x)-1
    push!(y, y[i]+1)
  end
  plt.step(x,y)
end
for k in 1:1
  x = generate_NHPP(λ,range)
  y = [1]
  for i in 1:length(x)-1
    push!(y, y[i]+1)
  end
  plt.step(x,y)
end
plt.legend(["NHWP","NHPP"])
plt.savefig("Nonhomogeneous Non-Poisson process.pdf")
