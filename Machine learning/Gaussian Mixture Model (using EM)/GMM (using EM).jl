# This code is written by Yongkyu Cho.
# Reference: Bishop - Pattern Recognition and Machine Learning

using Distributions, PyPlot

# Synthetic data generation
mean = [[2.0, 6.0] , [7.0, 9.0] , [9.0, 3.0] , [5.0, 5.0]]
covariance = [[1 1.5;1.5 3] , [3 1;1 1] , [2 1;1 1] , [2 0.5;0.5 2]]
D = 100

# Model parameters
K = 4
X = [rand(MvNormal(mean[i],covariance[i]),D) for i in 1:length(mean)]
X = hcat(X...)  # same as cat(2,gmm.X)
N = length(X[1,:])
μ = Array{Vector{Float64}}(K)
Σ = Array{Array{Float64,2}}(K)
π = Array{Float64}(K)
γ = Matrix{Float64}(N,K) #  p(z_nk = 1)
L = 0.0

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
    gmm.μ[k] = [sum(gmm.X[1,n] for n in 1:length(gmm.X[1,:]))/length(gmm.X[1,:]) , sum(gmm.X[2,n] for n in 1:length(gmm.X[2,:]))/length(gmm.X[1,:])] # SAMPgmm.LE MEAgmm.N
    gmm.Σ[k] = (5*rand())*eye(2)
    gmm.π[k] = 1/gmm.K
  end
  for n in 1:gmm.N
    gmm.L += log(sum(gmm.π[k]*pdf_MvNormal(gmm.μ[k],gmm.Σ[k],gmm.X[:,n])[1] for k in 1:gmm.K))
  end
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
  gmm.L = 0.0
  for n in 1:gmm.N
    gmm.L += log(sum(gmm.π[k]*pdf_MvNormal(gmm.μ[k],gmm.Σ[k],gmm.X[:,n])[1] for k in 1:gmm.K))
  end
end

function Plot_GMM(gmm::GMM)
  figure(figsize = (15,5))
  subplot(1,3,1)
  title("Given Data Points")
  scatter(gmm.X[1,:],gmm.X[2,:]) # without clustering
  # savefig("Given Data Points.pdf")
  # figure()
  subplot(1,3,2)
  title("True Clusters")
  for i in 1:length(mean)
    scatter(gmm.X[1,1+D*(i-1):D*i],gmm.X[2,1+D*(i-1):D*i])
  end
  # savefig("True Clusters.pdf")
  Z = Points[]
  for k in 1:gmm.K
    push!(Z,Points(Float64[],Float64[]))
  end
  for n in 1:gmm.N
    k = findmax(gmm.γ[n,:])[2]
    push!(Z[k].X, gmm.X[1,n])
    push!(Z[k].Y, gmm.X[2,n])
  end
  # figure()
  subplot(1,3,3)
  title("Clustered Data Points")
  for k in 1:gmm.K
    scatter(Z[k].X,Z[k].Y)
  end
  # savefig("Clustered Data Points.pdf")
  savefig("GMM clustering result.pdf")
end

function EM!(gmm::GMM)
  initialize_EM!(gmm)
  MAX_ITER = 1000
  itr = 0
  improvement = 10000000
  ϵ = 0.00000000000001
  while itr < MAX_ITER && improvement > ϵ
    itr += 1
    L_old = gmm.L
    E_step!(gmm)
    M_step!(gmm)
    improvement = gmm.L - L_old
    println("Iteration $(itr), Improvement: $(improvement)")
  end
  println("Finished after $itr iterations.")
end

# main
gmm = GMM(X,μ,Σ,π,γ,L,K,N)
EM!(gmm)
Plot_GMM(gmm)
