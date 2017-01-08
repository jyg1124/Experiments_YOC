using PyPlot
using Debug

const ITERATION = 5000
# For plotting
speed_array = Array{Float64}[]
price_array = Array{Float64}[]
obj_array = Float64[]

# type definitions
type Server_Setting
  γ::Float64
  Γ::Float64
  K::Float64
  α::Float64
  n::Int64
  δ::Float64
  ϵ::Float64
  x_0::Float64
  p_0::Float64
  Apps::Tuple
end

type Workload_Setting
  instant_demand::Float64
  inter_arrival_distribution::ASCIIString
  mean_inter_arrival::Float64
  scv_inter_arrival::Float64
  std_inter_arrival::Float64
  workload_distribution::ASCIIString
  mean_workload::Float64
  scv_workload::Float64
  std_workload::Float64
end

type Server
  previous_speed::Float64  # 기존 속도
  current_speed::Float64   # 현재 속도

  previous_price::Float64  # 기존 price
  current_price::Float64   # 현재 price

  previous_remaining_workload::Float64 # 기존 remaining workload 총합 (t 시점)
  current_remaining_workload::Float64  # 현재 remaining workload 총합 (t+1 시점)

  κ::Float64                         # 서버의 버퍼 크기

  function Server(previous_speed::Float64,
                  current_speed::Float64,
                  previous_price::Float64,
                  current_price::Float64)
    previous_remaining_workload = 0.0
    current_remaining_workload = 0.0
    κ = 0.0

    new(previous_speed,
        current_speed,
        previous_price,
        current_price,
        previous_remaining_workload,
        current_remaining_workload,
        κ)
  end
end

type VirtualDataCenter
  ## These variables are set directly by the creator
  SS::Array{Server_Setting}        # Server_Setting을 담는 배열
  S::Array{Server}                 # type Server 담는 배열
  # constructor definition
  function VirtualDataCenter(SS::Array{Server_Setting}, S::Array{Server})
    new(SS,S)
  end
end

# function definitions
function server_power(j::Int64, SS::Array{Server_Setting}, S::Array{Server})
  return SS[j].K + (SS[j].α)*(S[j].current_speed^SS[j].n)
end

function server_power_1st_diff(j::Int64, SS::Array{Server_Setting}, S::Array{Server})
  return (SS[j].α)*(SS[j].n)*(S[j].current_speed^((SS[j].n)-1))
end

function server_power_2nd_diff(j::Int64, SS::Array{Server_Setting}, S::Array{Server})
  return (SS[j].α)*(SS[j].n)*((SS[j].n)-1)*(S[j].current_speed^((SS[j].n)-2))
end

function x_dot(j::Int64, SS::Array{Server_Setting}, S::Array{Server})
  if SS[j].γ < S[j].current_speed < SS[j].Γ
    return S[j].current_price-server_power_1st_diff(j,SS,S)
  elseif S[j].current_speed == SS[j].Γ
    return min( S[j].current_price - server_power_1st_diff(j,SS,S), 0.0 )
  elseif S[j].current_speed == SS[j].γ
    return max(S[j].current_price - server_power_1st_diff(j,SS,S),0.0)
  end
end

function p_dot(j::Int64, S::Array{Server})
  if S[j].current_price >= 0.0
    return S[j].κ + S[j].current_remaining_workload - S[j].current_speed
  else
    return max(S[j].κ + S[j].current_remaining_workload - S[j].current_speed, 0.0)
  end
end

function find_min_price_server(app_type::Int64, SS::Array{Server_Setting}, S::Array{Server})
  server_index = 0
  temp_price = typemax(Float64)

  for j in 1:length(SS)
    if in(app_type, SS[j].Apps) == true
      if temp_price > S[j].current_price
        temp_price = S[j].current_price
        server_index = j
      end
    end
  end

  return server_index
end

function iterate(vdc::VirtualDataCenter, WS::Array{Workload_Setting}, num_iter::Int64)
  iter = 0
  obj = 0.0

  while iter < num_iter
    server_dispatched = Int64[]
    for i in 1:length(WS)
      push!(server_dispatched, find_min_price_server(i, vdc.SS, vdc.S))
      vdc.S[server_dispatched[i]].current_remaining_workload = WS[i].instant_demand
    end
    println("$(server_dispatched[1]), $(server_dispatched[2]), $(server_dispatched[3]), $(server_dispatched[4]), $(server_dispatched[5])")

    for j in 1:length(vdc.S)
      # For plotting
      push!(speed_array[j], vdc.S[j].current_speed)
      push!(price_array[j], vdc.S[j].current_price)
      obj += server_power(j, vdc.SS, vdc.S)

      # updating speeds
      vdc.S[j].previous_speed = vdc.S[j].current_speed
      vdc.S[j].current_speed = vdc.S[j].previous_speed + (1/server_power_2nd_diff(j, vdc.SS, vdc.S))*(x_dot(j, vdc.SS, vdc.S))

      # updating prices
      vdc.S[j].previous_price = vdc.S[j].current_price
      vdc.S[j].current_price = vdc.S[j].previous_price + p_dot(j, vdc.S)
    end
#    println("$obj")
    push!(obj_array, obj)
    obj = 0.0
    iter += 1
  end
end

function workload_setter()
  WS = Workload_Setting[]
  push!(WS, Workload_Setting(20.0, "LogNormal",0.25,2.0,sqrt((0.25^2)*2.0),"LogNormal",5.0,1.5,sqrt((5.0^2)*1.5)))
  push!(WS, Workload_Setting(20.0, "LogNormal", 0.5, 1.5, sqrt((0.5^2)*1.5), "LogNormal", 10.0, 2.0, sqrt((10.0^2)*2.0)))
  push!(WS, Workload_Setting(20.0, "Exponential", 0.25, 1.0, sqrt((0.25^2)*1.0), "LogNormal", 5.0, 1.0, sqrt((5.0^2)*1.0)))
  push!(WS, Workload_Setting(20.0, "LogNormal", 0.1, 0.8, sqrt((0.1^2)*0.8), "LogNormal", 2.0, 0.8, sqrt((2.0^2)*0.8)))
  push!(WS, Workload_Setting(15.0, "LogNormal", 0.2, 2.0, sqrt((0.2^2)*2.0),"LogNormal", 3.0,0.5,sqrt((3.0^2)*0.5)))
  return WS
end

# 서버별 정보를 생성해서 Server_Setting array를 리턴하는 함수
function server_setter()
  SS = Server_Setting[]
  push!(SS, Server_Setting(5.0, 100.0, 150.0, 0.3333, 3, 3.0, 0.001, 100.0, 1000.0, (1,)))
  push!(SS, Server_Setting(7.0, 102.0, 250.0, 0.2, 3, 3.0, 0.001, 102.0, 2000.0, (1,)))
  push!(SS, Server_Setting(6.0, 99.0, 220.0, 1.0, 3, 3.0, 0.001, 99.0, 3000.0, (1,2)))
  push!(SS, Server_Setting(5.0, 105.0, 150.0, 0.6667, 3, 3.0, 0.001, 105.0, 1000.0, (1,2,3)))
  push!(SS, Server_Setting(7.0, 100.0, 300.0, 0.8, 3, 3.0, 0.001, 100.0, 2000.0, (2,3)))
  push!(SS, Server_Setting(8.0, 102.0, 350.0, 0.4, 3, 3.0, 0.001, 102.0, 3000.0, (2,3)))
  push!(SS, Server_Setting(6.0, 100.0, 220.0, 0.4286, 3, 3.0, 0.001, 100.0, 1000.0, (3,)))
  push!(SS, Server_Setting(7.0, 105.0, 350.0, 0.5, 3, 3.0, 0.001, 105.0, 2000.0, (4,5)))
  push!(SS, Server_Setting(8.0, 102.0, 400.0, 0.6, 3, 3.0, 0.001, 102.0, 3000.0, (4,5)))
  push!(SS, Server_Setting(10.0, 105.0, 700.0, 0.4444, 3, 3.0, 0.001, 105.0, 1000.0, (5,)))
  return SS
end

# 서버 객체를 만들어서 Server array를 리턴하는 함수
function server_creater(SS::Array{Server_Setting}, WS::Array{Workload_Setting})
  #aggreated scv를 먼저 계산
  tempv = Float64[]
  for i in 1:length(WS)
    push!(tempv,WS[i].mean_inter_arrival)
  end
  μ_min = minimum(tempv)

  num = 0.0
  denom = 0.0
  for i in 1:length(WS)
    num += (WS[i].std_inter_arrival)^2
    denom += (WS[i].mean_inter_arrival)
  end
  agg_scv = num/((denom)^2)

  #서버 객체 생성 및 초기값 설정
  S = Server[]
  for j in 1:length(SS)
    push!(S, Server(SS[j].x_0, SS[j].x_0, SS[j].p_0, SS[j].p_0))
    # with κ
     S[j].κ = (-log(SS[j].ϵ)*max(1,agg_scv))/(μ_min*SS[j].δ)
    # without κ
    # S[j].κ = 0.0
  end

  return S
end


# main part #

# Initialization
WS = workload_setter()
SS = server_setter()
S = server_creater(SS,WS)
vdc = VirtualDataCenter(SS,S)

# Run

# For plotting
for j = 1:length(vdc.S)
  push!(speed_array, Float64[])
  push!(price_array, Float64[])
end
iterate(vdc, WS, ITERATION)

# Plotting Speeds and Prices and Energy consumption
x = linspace(0,ITERATION,ITERATION)
y = Float64[]
for i in 1:ITERATION
  push!(y, vdc.S[1].κ)
end
plt = PyPlot

plt.figure()
for j in 1:length(S)
  plt.subplot(2,5,j)
  plt.title("Server $j", fontsize=10)
  plt.xlabel("Iteration",fontsize=6)
  if j == 1 || j == 6
    plt.ylabel("Speed (workload/time)",fontsize=10)
  end
  plt.yticks(linspace(0,100,11))
  plt.xticks([0,ITERATION])
  plt.ylim(0,100)
  plt.plot(x,speed_array[j][1:ITERATION],linewidth=1.0,linestyle="-",color="red")
  plt.plot(x,y,linewidth=1.0,linestyle="--", color = "blue")
  plt.tick_params(labelsize=6)

  #plt.plot(x,y,linewidth=2.0,linestyle="--",color="black")
end
plt.savefig("Instantaneous server speeds.pdf")

plt.figure()
for j in 1:length(S)
  plt.subplot(2,5,j)
  plt.title("Server $j", fontsize=10)
  plt.xlabel("Iteration",fontsize=6)
  if j == 1 || j == 6
    plt.ylabel("Price",fontsize=10)
  end
  plt.xticks([0,ITERATION])
  plt.plot(x,price_array[j][1:ITERATION],linewidth=1.0,linestyle="-",color="red")
  plt.tick_params(labelsize=6)
end
plt.savefig("Instantaneous server prices.pdf")

plt.figure()
plt.title("Energy Consumption")
x = linspace(1,5000,5000)
# plt.xticks([1:100])
plt.plot(x, obj_array[1:5000], linewidth=1.0, linestyle="-",color="red")
plt.xlabel("Iterations",fontsize=10)
plt.ylabel("Objective value",fontsize=10)
plt.savefig("Energy consumption.pdf")
