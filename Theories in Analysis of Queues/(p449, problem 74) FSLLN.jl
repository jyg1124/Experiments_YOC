# Written by Yongkyu Cho.

using Distributions
using Plots
pyplot()
plt = PyPlot
const TIME = 10.0

function plot_arrival_process(n::Int64)
  x = [0.0]
  y = [0.0]
  while x[length(x)] < TIME
    push!(x,  x[length(y)] + rand(Pareto(2.095,0.5201))/n  )
    push!(y, y[length(y)] + 1/n)
  end
  plt.xlim(0,10)
  plt.ylim(0,12)
  plt.step(x,y)
end

for i in 1:5
  plot_arrival_process(100)
end
plt.savefig("(p449, problem 74) FSLLN.pdf")
show()
