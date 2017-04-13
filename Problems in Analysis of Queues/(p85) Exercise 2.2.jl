# Written by Yongkyu Cho.
# Reference: Analysis of Queues, Natarajan Gautam.

using PyPlot

# Define parameters
λ = 10.0
μ = 1.25
s = [5,10]
color = ["blue","red"]
what_to_plot = [("P_K", "Blocking probability"), ("W","Average number of customers in the system")]

# Define functions
function p(i::Int64, s::Int64, K::Int64)  # i stands for the number of customers in the system
  temp = 0.0
  for j in 0:s
    temp += (1/factorial(j))*((λ/μ)^j)
  end
  for j in s+1:K
    temp += ((λ/(s*μ))^(j-s))*(1/factorial(s))*((λ/μ)^s)
  end

  if i == 0
    return 1/temp
  elseif 1 <= i <= s
    return (1/factorial(i))*((λ/μ)^i)*(1/temp)
  elseif s+1 <= i <= K
    return ((λ/(s*μ))^(i-s))*(1/factorial(s))*((λ/μ)^s)*(1/temp)
  end
end

function W(s::Int64, K::Int64)
  L = 0.0
  for i in 0:K
    L += i*p(i,s,K)
  end
  return L/(λ*(1-p(K,s,K)))
end

# Plot results
plt = PyPlot
plt.figure(figsize=(25,15))
for i in 1:length(what_to_plot)
  #plt.subplot(length(what_to_plot),1,i)
  plt.subplot(1,length(what_to_plot),i)
  plt.xlabel(L"$K$",fontsize=30)
  plt.ylabel("$(what_to_plot[i][1])",fontsize=30)
  x = linspace(minimum(s),maximum(s)+19,maximum(s)+19-minimum(s)+1)
  y = Float64[]

  if what_to_plot[i][1] == "P_K"
    for j in 1:length(s)
      x = linspace(s[j],s[j]+19,20)
      empty!(y)
      for K in s[j]:s[j]+19
        push!(y, p(K, s[j], K))
      end
      plt.plot(x,y,linewidth=1.0,linestyle="-",color=color[j], label="s = $(s[j])")
    end
  elseif what_to_plot[i][1] == "W"
    for j in 1:length(s)
      x = linspace(s[j],s[j]+19,20)
      empty!(y)
      for K in s[j]:s[j]+19
        push!(y, W(s[j], K))
      end
      plt.plot(x,y,linewidth=1.0,linestyle="-",color=color[j], label="s = $(s[j])")
    end
  end
  plt.legend(fontsize=20)
  plt.tick_params(labelsize=30)

end

plt.savefig("hw2-2_plot.pdf")
