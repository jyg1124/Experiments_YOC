# Written by Yongkyu Cho.
# Reference: Pattern Recognition and Machine Learning, Christopher Bishop.

using Distributions, PyPlot

# Type definitions
type GMM
  X::Matrix{Float64}          # Given data points
  μ::Array{Vector{Float64}}   # Mean vectors
  Σ::Array{Array{Float64,2}}  # Covariance matrices
  π::Array{Float64}           # Mixing coefficients vector : p(z_k = 1)
  γ::Matrix{Float64}          # Responsiblity matrix : [ p(z_nk = 1) ]
  L::Float64                  # Total data Loglikelihood
  K::Int64                    # The Number of clusters we are using to fit
  N::Int64                    # The Number of data points
end

type Points
  X::Array{Float64}
  Y::Array{Float64}
end

# Function definitions
function pdf_MvNormal(μ::Vector{Float64}, Σ::Array{Float64,2}, x::Vector{Float64})
  k = length(μ)
  ((1/sqrt(2*pi))^k)*(1/sqrt(det(Σ)))*exp(-(1/2)*transpose(x-μ)*inv(Σ)*(x-μ))
end

function initialize_EM!(gmm::GMM)
  for k in 1:gmm.K
    gmm.μ[k] = [sum(gmm.X[1,n] for n in 1:length(gmm.X[1,:]))/length(gmm.X[1,:]) , sum(gmm.X[2,n] for n in 1:length(gmm.X[2,:]))/length(gmm.X[1,:])] # SAMPLE MEAN
    gmm.Σ[k] = (5*rand())*eye(2)
    gmm.π[k] = 1/gmm.K
  end
  gmm.L = sum(log(sum(gmm.π[k]*pdf_MvNormal(gmm.μ[k],gmm.Σ[k],gmm.X[:,n])[1] for k in 1:gmm.K)) for n in 1:gmm.N)
end

function E_step!(gmm::GMM)
  for n in 1:gmm.N, k in 1:gmm.K
    gmm.γ[n,k] = (gmm.π[k]*pdf_MvNormal(gmm.μ[k],gmm.Σ[k],gmm.X[:,n])[1]) / sum(gmm.π[j]*pdf_MvNormal(gmm.μ[j],gmm.Σ[j],gmm.X[:,n])[1] for j in 1:gmm.K)
  end
end

function M_step!(gmm::GMM)
  for k in 1:gmm.K
    N_k = sum(gmm.γ[n,k] for n in 1:gmm.N)
    gmm.μ[k] = sum(gmm.γ[n,k]*gmm.X[:,n] for n in 1:gmm.N) / N_k
    gmm.Σ[k] = sum(gmm.γ[n,k]*(gmm.X[:,n]-gmm.μ[k])*transpose(gmm.X[:,n]-gmm.μ[k]) for n in 1:gmm.N) / N_k
    gmm.π[k] = N_k / gmm.N
  end
  gmm.L = sum(log(sum(gmm.π[k]*pdf_MvNormal(gmm.μ[k],gmm.Σ[k],gmm.X[:,n])[1] for k in 1:gmm.K)) for n in 1:gmm.N)
end

function EM!(gmm::GMM, MAX_ITER::Int64, ϵ::Float64)
  initialize_EM!(gmm)
  itr = 0
  improvement = typemax(Float64)
  while itr < MAX_ITER && improvement > ϵ
    itr += 1
    L_old = gmm.L
    E_step!(gmm)
    M_step!(gmm)
    improvement = gmm.L - L_old
    println("Iteration $(itr), Improvement of loglikelihood: $(improvement)")
  end
  println("Finished after $itr iterations.")
end

function Plot_GMM(gmm::GMM)
  figure(figsize = (15,5))
  subplot(1,3,1)
  title("Given Data Points")
  scatter(gmm.X[1,:],gmm.X[2,:]) # without clustering
  subplot(1,3,2)
  title("True Clusters")
  for i in 1:length(mean)
    scatter(gmm.X[1,1+D*(i-1):D*i],gmm.X[2,1+D*(i-1):D*i])
  end
  Z = Points[]
  for k in 1:gmm.K
    push!(Z,Points(Float64[],Float64[]))
  end
  for n in 1:gmm.N
    k = findmax(gmm.γ[n,:])[2]
    push!(Z[k].X, gmm.X[1,n])
    push!(Z[k].Y, gmm.X[2,n])
  end
  subplot(1,3,3)
  title("Clustered Data Points")
  for k in 1:gmm.K
    scatter(Z[k].X,Z[k].Y)
  end
  savefig("GMM clustering result.pdf")
end

# Parameters for generating synthetic data
mean = [[2.0, 6.0] , [7.0, 9.0] , [9.0, 3.0] , [5.0, 5.0]]
covariance = [[1 1.5;1.5 3] , [3 1;1 1] , [2 1;1 1] , [2 0.5;0.5 2]]
D = 100 # the number of data points generated from a parameter set (μ_i,Σ_i)

# Model parameters
K = 4                                                                 # the number of clusters
X = [rand(MvNormal(mean[i],covariance[i]),D) for i in 1:length(mean)] # generating random data points
X = hcat(X...)                                                        # converting the container type (array to matrix). working exactly like cat(2,gmm.X)
N = length(X[1,:])
μ = Array{Vector{Float64}}(K)
Σ = Array{Array{Float64,2}}(K)
π = Array{Float64}(K)
γ = Matrix{Float64}(N,K)                                             
L = 0.0

# Main
gmm = GMM(X,μ,Σ,π,γ,L,K,N)
EM!(gmm, 2000, 0.000000000000001)
Plot_GMM(gmm)
