using Distributions
using Plots
pyplot()

function sampleBM(x0::Float64, μ::Float64, σ::Float64, n::Int64)   # n: number of samples
  path = Float64[]
  for i in 1:n
    if i == 1
      push!(path, x0)
    else
      push!(path, path[i-1] + rand(Normal(μ, σ)))
    end
  end
  return path
end

function sup_max(BM_path::Array{Float64}, k::Int64)
  temp_array = Float64[]
  for i in 1:k
    if BM_path[i] < 0.0
      push!(temp_array, -BM_path[i])
    else
      push!(temp_array, 0.0)
    end
  end
  return maximum(temp_array)
end

function sampleRBM(x0::Float64, μ::Float64, σ::Float64, n::Int64)   # n: number of samples
  BM_path = Float64[]
  RBM_path = Float64[]
  for i in 1:n
    if i == 1
      push!(BM_path, x0)
      push!(RBM_path, x0)
    else
      temp_normal = rand(Normal(μ, σ))
      push!(BM_path, BM_path[i-1] + temp_normal)
      push!(RBM_path, BM_path[i-1] + sup_max(BM_path,i))
    end
  end
  return RBM_path
end

function BM_to_RBM(BM_path::Array{Float64})
  RBM_path = Float64[]
  for i in 1:length(BM_path)
    if i == 1
      push!(RBM_path, BM_path[1])
    else
      push!(RBM_path, BM_path[i] + sup_max(BM_path, i))
    end
  end
  return RBM_path
end

BM_path = sampleBM(6.0,-0.01,0.09,1000)
#RBM_path = sampleRBM(6.0,-0.01,0.09,1000)
RBM_path = BM_to_RBM(BM_path)

x_axis = Float64[]
for i in 1:1000
  push!(x_axis, 0.0)
end

plt = PyPlot
#plt.subplot(2,1,1)
plt.xticks(0:100:1000)

plt.plot(BM_path)
plt.plot(RBM_path)
plt.plot(x_axis, color="black")
plt.savefig("(p379) RBM_GG1_approximation.pdf")
show()
