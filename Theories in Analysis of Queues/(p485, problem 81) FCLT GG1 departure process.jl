# Written by Yongkyu Cho.

using Distributions
using Plots
pyplot()
f = open("sim_record.txt" , "w")
# f2 = open("departure_process.txt.","w")

const λ = 1.0
const μ = 1.25

type customer
  remaining_service_time::Float64
  function customer(_remaining_service_time::Float64)
    remaining_service_time = _remaining_service_time
    new(remaining_service_time)
  end
end

type GG1_queue
  number_of_customers::Int64
  sim_time::Float64
  next_completion_time::Float64
  next_arrival_time::Float64
  customers_in_system::Array{customer}
  function GG1_queue()
    number_of_customers = 0
    sim_time = 0.0
    next_completion_time = 99999999.9
    next_arrival_time = rand(Pareto(2.095,0.5201))
    customers_in_system = customer[]
    new(number_of_customers, sim_time, next_completion_time, next_arrival_time, customers_in_system)
  end
end

function next_event(system::GG1_queue, departure_time::Array{Float64})
  if system.next_arrival_time < system.next_completion_time
    duration = system.next_arrival_time - system.sim_time
    if system.number_of_customers >= 1
      system.customers_in_system[1].remaining_service_time -= duration
    end
    push!(system.customers_in_system, customer(rand(Pareto(2.183,0.432))))
    system.number_of_customers += 1
    system.sim_time = system.next_arrival_time
    system.next_arrival_time = system.sim_time + rand(Pareto(2.095,0.5201))
    system.next_completion_time = system.sim_time + system.customers_in_system[1].remaining_service_time
    println(f,"(Time: $(system.sim_time)) Current event: New arrival (workload: $(system.customers_in_system[length(system.customers_in_system)].remaining_service_time))")
  else
    push!(departure_time, system.next_completion_time) # record departure times
    duration = system.next_completion_time - system.sim_time
    deleteat!(system.customers_in_system, 1)
    system.number_of_customers -= 1
    system.sim_time = system.next_completion_time
    if system.number_of_customers >= 1
      system.next_completion_time = system.sim_time + system.customers_in_system[1].remaining_service_time
    else
      system.next_completion_time = 99999999.9
    end
    println(f,"(Time: $(system.sim_time)) Current event: Service completion (remaining number of customers: $(system.number_of_customers))")
  end
end

function run_simulation(t::Float64) # t: simulation period (should be larger than n*t)
  departure_time = Float64[]
  system = GG1_queue()
  while system.sim_time < t
    next_event(system, departure_time)
  end
  return departure_time
end

function save_diffusion_scaled_process_points!(x::Array{Float64}, y::Array{Float64}, departure_time::Array{Float64}, sampling_period::Float64, scaling_parameter::Int64)
  n = scaling_parameter
  push!(x, 0.0)
  push!(y, 0.0)
  i = 1
  j = 2
  for i in 1:length(departure_time)
    push!(x, departure_time[i]/n)
    if i == 1
      push!(y, y[j-1]-λ*sqrt(n)*((departure_time[i]-0.0)/n) )
    else
      push!(y, y[j-1]-λ*sqrt(n)*((departure_time[i]-departure_time[i-1])/n) )
    end
    push!(x, departure_time[i]/n)
    push!(y, y[j] + (1/sqrt(n)) )
    i += 1
    j += 2
  end
end

function plot_sample_paths(k::Int64, n::Int64, t::Float64) # k: number of paths, n: scaling parameter, t: sampling period
  plot()
  # plot!(xlims = (0, TIME), ylims = (-10,10))
  for i in 1:k
    (x,y) = (Float64[],Float64[])
    departure_time = run_simulation(t*n)
    save_diffusion_scaled_process_points!(x, y, departure_time, t, n)
    plot!(x,y)
  end
  gui()
  savefig("(p485, problem 81) FCLT GG1 departure process.pdf")
end

plot_sample_paths( 3 , 100, 10.0)

 close(f)
# close(f2)
