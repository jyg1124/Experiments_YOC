using PyPlot
using Distributions
using Gallium
plt = PyPlot
f = open("sim_record (7.3).txt" , "w")


function graph(f::Function, x_range::LinSpace{Float64})
  plt.plot(x_range,map(f,x_range))
end

function λ(t::Float64)
  if 0.0 <= t%10 < 5.0
    return 3*(t%10) + 4
  else
    return -3*(t%10) + 34
  end
end

# graph(λ,linspace(0.0,100.0,10000))

function NHPP(f::Function, T::Float64)
  x = Float64[] # arrival times
  λ = maximum(map(f,linspace(0.0,T,T*10)))
  t = 0.0
  while t < T
    t -= (1/λ)*log(rand())
    if rand() <= f(t)/λ
      push!(x, t)
    end
  end
  return x
end

type customer
  remaining_service_time::Float64
  function customer()
    remaining_service_time = 0.0 # 사실 거의 쓸일이 없다 (service time이 exponential이므로 볼 때마다 새로워지기 때문)
    new(remaining_service_time)
  end
end

type queue
  arrival_times::Array{Float64}
  number_of_customers::Int64
  number_of_servers::Int64
  sim_time::Float64
  next_completion_time::Float64
  next_arrival_time::Float64
  next_server_return_time::Float64
  customers_in_system::Array{customer}

  function queue(_arrival_times::Array{Float64})
    arrival_times = _arrival_times
    number_of_customers = 0
    number_of_servers = 1
    sim_time = 0.0
    next_completion_time = typemax(Float64)
    next_arrival_time = 0.0
    next_server_return_time = typemax(Float64)
    customers_in_system = customer[]
    new(arrival_times, number_of_customers, number_of_servers, sim_time, next_completion_time, arrival_times[1], next_server_return_time, customers_in_system)
  end
end

type record
  event_times::Array{Float64}
  number_in_system::Array{Int64}
  server_vacation_time::Float64
  function record()
    event_times = [0.0]
    number_in_system = [0]
    server_vacation_time = 0.0
    new(event_times, number_in_system, server_vacation_time)
  end
end

function next_event(system::queue, recording::record)
  if system.next_completion_time == min(system.next_completion_time, system.next_arrival_time, system.next_server_return_time) # 이벤트: 서비스 완료
    inter_event_time = system.next_completion_time - system.sim_time # inter event time 저장
    deleteat!(system.customers_in_system, 1) # 서비스 끝난 고객을 시스템에서 제거
    system.number_of_customers -= 1 # 시스템 내 고객 수 1 감소
    system.sim_time = system.next_completion_time # 시스템의 현재 시각 변경
    if system.number_of_customers >= 1 # 서비스가 끝나서 고객이 나간 직후에, 시스템에 고객이 남아 있으면
      system.next_completion_time = system.sim_time + rand(Exponential(1/25)) # 다음 completion time은 지금시간 + 서비스시간
    else # 서비스가 끝나서 고객이 나간 직후에, 시스템이 텅텅 비어있으면
      system.next_completion_time = typemax(Float64) # 다음 completion time은 infinity로 설정
      system.next_server_return_time = system.sim_time + rand(Uniform(0,0.3)) # 서버가 휴가를 떠나고, Unif(0,0.3)후에 돌아옴 (return time 설정)
      system.number_of_servers = 0 # 시스템 안의 서버의 수는 0명이 됨
      recording.server_vacation_time += system.next_server_return_time - system.sim_time # 서버의 휴가기간 축적
    end

    push!(recording.event_times, system.sim_time)
    push!(recording.number_in_system, system.number_of_customers)
    # println(f, "(Time: $(system.sim_time)) Current event: Service completion (remaining # of customers: $(system.number_of_customers), remaining # of servers: $(system.number_of_servers))")

  elseif system.next_arrival_time == min(system.next_completion_time, system.next_arrival_time, system.next_server_return_time) # 이벤트: 고객 도착
    inter_event_time = system.next_arrival_time - system.sim_time # inter event time 저장
    push!(system.customers_in_system, customer()) # 새로 도착한 고객을 대기열에 삽입
    system.number_of_customers += 1 # 시스템 내 고객 수 1 증가
    system.sim_time = system.next_arrival_time # 시스템의 현재 시각 변경
    deleteat!(system.arrival_times, 1) # arrival times에서 맨 앞에 있는거 하나 제거
    system.next_arrival_time = system.arrival_times[1] # 다음 도착 시각 설정 (arrival times에서 가져옴)
    if system.number_of_servers == 1 # 고객이 도착했을 때, 서버가 있으면,
      system.next_completion_time = system.sim_time + rand(Exponential(1/25)) # 다음 서비스 완료 시각 재할당 (memoryless 하므로 새로운 exponenrial 추가)
    else
      system.next_completion_time = typemax(Float64) # 다음 서비스 완료 시각은 infinity로 할당
    end

    push!(recording.event_times, system.sim_time)
    push!(recording.number_in_system, system.number_of_customers)
    # println(f, "(Time: $(system.sim_time)) Current event: Customer arrival (remaining # of customers: $(system.number_of_customers), remaining # of servers: $(system.number_of_servers))")

  else system.next_server_return_time == min(system.next_completion_time, system.next_arrival_time, system.next_server_return_time)
    inter_event_time = system.next_server_return_time - system.sim_time # inter event time 저장
    system.sim_time = system.next_server_return_time # 시스템의 현재 시각 변경
    if system.number_of_customers == 0 # 만약 서버가 돌아왔을 때, queue가 비어있으면
      system.next_server_return_time = system.sim_time + rand(Uniform(0,0.3)) # 서버가 또 휴가를 떠나고, Unif(0,0.3)후에 돌아옴 (return time 설정)
      recording.server_vacation_time += (system.next_server_return_time - system.sim_time) # 서버의 휴가 기간 축적
      # println(f, "(Time: $(system.sim_time)) Current event: Server return but leave again (remaining # of customers: $(system.number_of_customers), remaining # of servers: $(system.number_of_servers))")
    else # 서버가 돌아왔을 때, 대기 중인 고객들이 있으면
      system.number_of_servers = 1 # 서버는 다시 일을 시작
      system.next_server_return_time = typemax(Float64) # 다음 서버 복귀 시간은 infinity로
      system.next_completion_time = system.sim_time + rand(Exponential(1/25)) # 다음 completion time 설정
      # println(f, "(Time: $(system.sim_time)) Current event: Server return and start working (remaining # of customers: $(system.number_of_customers), remaining # of servers: $(system.number_of_servers))")
    end
  end
end

function run_simulation(T::Float64)
  system = queue(NHPP(λ,T*2))
  X = record()
  while system.sim_time < T
    next_event(system, X)
  end
  return X
end

sample_path = run_simulation(100.0)

# plotting X(t) process
plt.step(sample_path.event_times,sample_path.number_in_system)

println("Amount of time that the servier is on break: $(sample_path.server_vacation_time)")

function do_replication(n::Int64)
  cusum = 0.0
  for i in 1:n
    cusum += run_simulation(100.0).server_vacation_time
  end
  println("Expected amount of time that the servier is on break: $(cusum/n)")
end

do_replication(500)

close(f)
