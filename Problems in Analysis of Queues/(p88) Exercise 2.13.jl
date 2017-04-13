# Written by Yongkyu Cho.
# Reference: Analysis of Queues, Natarajan Gautam.

using Distributions

C=10.0
r=500.0
λ=10.0
μ=20.0
time = 10000.0

function baik(C::Float64, r::Float64, λ::Float64, μ::Float64, time::Float64)
  (  (2*(       μ*((μ+λ)^2)*r  -  μ*λ*(μ^2 + μ*λ + λ^2)*C        )  )   /     (   (  2*(λ^2)+(2*λ*μ)+(μ^2)  )*(λ+μ)  )   )*time
end

function cho(C::Float64, r::Float64, λ::Float64, μ::Float64, time::Float64)
  ( (2*r*(μ^2) + 2*r*λ*μ - C*μ*(2*λ*μ) - 2*C*μ*(λ^2) )/(μ^2 + 2*μ*λ + 2*(λ^2)) )*time
end

function revenue(C::Float64, r::Float64, λ::Float64, μ::Float64, time::Float64)
  ((2*r*(μ^2)+2*r*λ*μ) / (μ^2 + 2*μ*λ + 2*(λ^2)) )*time
end

function baik_cost(C::Float64, r::Float64, λ::Float64, μ::Float64, time::Float64)
  (((((C*(μ^2))/(λ+μ))*2*λ*μ)+2*C*μ*(λ^2))  /  (μ^2 + 2*μ*λ + 2*(λ^2)))*time
end

function cho_cost(C::Float64, r::Float64, λ::Float64, μ::Float64, time::Float64)
 ((C*μ*(2*λ*μ) + 2*C*μ*(λ^2))/(μ^2 + 2*μ*λ + 2*(λ^2)))*time
end

baik(C,r,λ,μ,time)
cho(C,r,λ,μ,time)

const WARMUP_LENGTH = 0.5*time
const REPLICATION_LENGTH = WARMUP_LENGTH + time
const NUM_REPLICATION = 10
f = open("sim_record.txt" , "w")

type Machine_Shop
  C::Float64
  r::Float64
  λ::Float64
  μ::Float64
  revenue::Float64
  cost::Float64
  profit::Float64 # revenue - cost
  num_working_machines::Int64
  life_time::Float64
  repair_time::Float64
  interval_time::Float64

  function Machine_Shop(C::Float64, r::Float64, λ::Float64, μ::Float64)
    revenue = 0.0
    cost = 0.0
    profit = 0.0
    num_working_machines = 2
    life_time = 0.0
    repair_time = 0.0
    interval_time = 0.0
    new(C, r, λ, μ, revenue, cost, profit, num_working_machines, life_time, repair_time, interval_time)
  end
end

type Simulation_Setting
  current_time::Float64
  previous_time::Float64
  function Simulation_Setting()
    current_time = 0.0
    previous_time = 0.0
    new(current_time, previous_time)
  end
end

MS = Machine_Shop(C,r,λ,μ)
SS = Simulation_Setting()

function run_to_end(MS::Machine_Shop, SS::Simulation_Setting)
  while SS.current_time < REPLICATION_LENGTH
    if MS.num_working_machines == 2
      MS.interval_time = rand(Exponential(1/(2*MS.λ)))
      SS.previous_time = SS.current_time
      SS.current_time += MS.interval_time
      MS.num_working_machines = 1
      if SS.current_time > WARMUP_LENGTH
        MS.revenue += 2*MS.r*(SS.current_time-SS.previous_time)
        println(f,"(Time: $(SS.current_time)) Current Event: Fail, #Working machines: $(MS.num_working_machines)")
      end
    elseif MS.num_working_machines == 1
      MS.life_time = rand(Exponential(1/MS.λ))
      MS.repair_time = rand(Exponential(1/MS.μ))
      MS.interval_time = min(MS.life_time, MS.repair_time)

      SS.previous_time = SS.current_time
      SS.current_time += MS.interval_time
      if MS.life_time < MS.repair_time # 고장이 먼저남
        MS.num_working_machines = 0
        if SS.current_time > WARMUP_LENGTH
          MS.revenue += 1*MS.r*(SS.current_time-SS.previous_time)
          println(f,"(Time: $(SS.current_time)) Current Event: Fail, #Working machines: $(MS.num_working_machines)")
        end
      else # 수리가 먼저 됨
        MS.num_working_machines = 2
        if SS.current_time > WARMUP_LENGTH
          MS.revenue += 1*MS.r*(SS.current_time-SS.previous_time)
          MS.cost += MS.C
          println(f,"(Time: $(SS.current_time)) Current Event: Repaired, #Working machines: $(MS.num_working_machines)")
        end
      end
    elseif MS.num_working_machines == 0
      MS.interval_time = rand(Exponential(1/MS.μ))
      SS.previous_time = SS.current_time
      SS.current_time += MS.interval_time
      MS.num_working_machines = 1
      if SS.current_time > WARMUP_LENGTH
        MS.cost += MS.C
        println(f,"(Time: $(SS.current_time)) Current Event: Repaired, #Working machines: $(MS.num_working_machines)")
      end
    end
  end
  MS.profit = MS.revenue - MS.cost
  println(f,"Simulation finished at time $(SS.current_time) with total revenue $(MS.revenue), total cost $(MS.cost), total profit $(MS.profit).")
  println(f,"Warm-up period: $(WARMUP_LENGTH)")
  println(f,"Baik's answer: $(baik(C,r,λ,μ,time))")
  println(f,"Cho's answer: $(cho(C,r,λ,μ,time))")
  println(f, "Computed Revenue: $(revenue(C,r,λ,μ,time))")
  println(f, "Computed Cost (Baik): $(baik_cost(C,r,λ,μ,time))")
  println(f, "Computed Cost (Cho): $(cho_cost(C,r,λ,μ,time))")
  println("Simulation finished at time $(SS.current_time) with total revenue $(MS.revenue), total cost $(MS.cost), total profit $(MS.profit).")
end


profit_sum = 0.0
for i in 1:NUM_REPLICATION
  println("$i th replication")
  MS = Machine_Shop(C,r,λ,μ)
  SS = Simulation_Setting()
  run_to_end(MS,SS)
  profit_sum += MS.profit
end
println("========== Result ==========")
println("Average profit: $(profit_sum/NUM_REPLICATION)")
println("Baik's answer: $(baik(C,r,λ,μ,time))")
println("Cho's answer: $(cho(C,r,λ,μ,time))")
println("Computed Revenue: $(revenue(C,r,λ,μ,time))")
println("Computed Cost (Baik): $(baik_cost(C,r,λ,μ,time))")
println("Computed Cost (Cho): $(cho_cost(C,r,λ,μ,time))")

close(f)
