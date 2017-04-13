# Written by Yongkyu Cho.
# Reference: Simulation, Sheldon M. Ross.

using PyPlot
using Distributions

plt = PyPlot
f = open("sim_record (7.6).txt" , "w")

function PP(λ::Float64, T::Float64)
  x = Float64[] # arrival times
  t = 0.0
  while t < T
    t -= (1/λ)*log(rand())
    push!(x, t)
    end
  return x
end

type customer
  _reneging_time::Float64
  function customer(reneging_time::Float64)
    _reneging_time = reneging_time
    new(_reneging_time)
  end
end

type queue
  _arrival_times::Array{Float64}
  _number_of_customers::Int64
  _number_of_reneged_customers::Int64 # 리니깅한 고객의 수 counting
  _sim_time::Float64
  _next_completion_time::Float64 # event: completion
  _next_arrival_time::Float64 # event: arrival
  _next_reneging_time::Float64 # event: reneging
  _customers_in_system::Array{customer}
  _next_reneging_customer_index::Int64 # 다음에 reneging할 고객의 queue.customers_in_system에서의 인덱스 넘버를 저장
  function queue(arrival_times::Array{Float64})
    _arrival_times = arrival_times
    _number_of_customers = 0
    _number_of_servers = 1
    _sim_time = 0.0
    _next_completion_time = typemax(Float64)
    _next_arrival_time = 0.0
    _next_reneging_time = typemax(Float64)
    _customers_in_system = customer[]
    _next_reneging_customer_index = 0
    new(_arrival_times, _number_of_customers, _number_of_servers, _sim_time, _next_completion_time, _arrival_times[1], _next_reneging_time, _customers_in_system, _next_reneging_customer_index)
  end
end

type record
  _event_times::Array{Float64} # 리니깅이 발생하는 시점을 기록
  _number_in_system::Array{Int64} # X(t)
  _number_of_reneging::Int64 # 누적 리니깅 횟수를 기록
  function record()
    _event_times = [0.0]
    _number_in_system = [0]
    _number_of_reneging = 0
    new(_event_times, _number_in_system, _number_of_reneging)
  end
end

function next_event(system::queue, recording::record)
  if system._next_completion_time == min(system._next_completion_time, system._next_arrival_time, system._next_reneging_time) # 이벤트: 서비스 완료
    inter_event_time = system._next_completion_time - system._sim_time # inter event time 저장
    deleteat!(system._customers_in_system, 1) # 서비스 끝난 고객을 시스템에서 제거
    system._number_of_customers -= 1 # 시스템 내 고객 수 1 감소
    system._sim_time = system._next_completion_time # 시스템의 현재 시각 변경
    if system._number_of_customers >= 2 # 서비스가 끝나서 고객이 나간 직후에, 시스템에 고객이 2명 이상 남아 있으면
      system._next_completion_time = system._sim_time + rand(Exponential(1/4)) # 다음 completion time은 지금시간 + 서비스시간
      # next_reneging_time 을 바꿔야한다
      x = Float64[]
      for i = 2:system._number_of_customers
        push!(x, system._customers_in_system[i]._reneging_time)
      end
      y = findmin(x)
      system._next_reneging_time = y[1] # next reneging time 갱신
      system._next_reneging_customer_index = y[2] + 1 # next reneging customer 갱신
    elseif system._number_of_customers == 1 # 서비스가 끝나서 고객이 나간 직후에, 시스템에 남아있는 고객이 딱 한 명이면,
      system._next_completion_time = system._sim_time + rand(Exponential(1/4)) # 다음 completion time은 지금시간 + 서비스시간
      system._next_reneging_time = typemax(Float64) # 다음 reneging time은 ∞
      system._next_reneging_customer_index = 0
    else # 서비스 끝나서 고객이 나간 직후에, 시스템이 텅 비어버리면,
      system._next_completion_time = typemax(Float64) # 다음 completion time은 ∞
      system._next_reneging_time = typemax(Float64)
      system._next_reneging_customer_index = 0
    end
    println(f, "(Time: $(system._sim_time)) Current event: Service completion (remaining # of customers: $(system._number_of_customers))")
    push!(recording._event_times, system._sim_time)
    push!(recording._number_in_system, system._number_of_customers)

  elseif system._next_arrival_time == min(system._next_completion_time, system._next_arrival_time, system._next_reneging_time) # 이벤트: 고객 도착
    inter_event_time = system._next_arrival_time - system._sim_time # inter event time 저장
    system._sim_time = system._next_arrival_time # 시스템의 현재 시각 변경
    push!(system._customers_in_system, customer(system._sim_time + rand(Uniform(0.0,5.0)))) # 새로 도착한 고객을 대기열에 삽입 (Uniform(0.0,5.0)의 patience time을 갖고 태어남)
    system._number_of_customers += 1 # 시스템 내 고객 수 1 증가
    deleteat!(system._arrival_times, 1) # arrival times에서 맨 앞에 있는거 하나 제거
    system._next_arrival_time = system._arrival_times[1] # 다음 도착 시각 설정 (arrival times에서 가져옴)
    system._next_completion_time = system._sim_time + rand(Exponential(1/4)) # 다음 completion 시각 설정 (exponential이므로 재설정)

    if system._number_of_customers >= 2 # 또한, 시스템내 고객수가 2명 이상이면,
      if system._customers_in_system[system._number_of_customers]._reneging_time < system._next_reneging_time # 만약 새로온 고객의 reneging time이 next reneging time 보다 가까우면,
        system._next_reneging_time = system._customers_in_system[system._number_of_customers]._reneging_time # next reneging time 교체
        system._next_reneging_customer_index = system._number_of_customers # next reneging customer index 교체
      end
    else
      system._next_reneging_time = typemax(Float64)
      system._next_reneging_customer_index = 0
    end
    println(f, "(Time: $(system._sim_time)) Current event: Customer arrival (remaining # of customers: $(system._number_of_customers))")
    push!(recording._event_times, system._sim_time)
    push!(recording._number_in_system, system._number_of_customers)

  else system._next_reneging_time == min(system._next_completion_time, system._next_arrival_time, system._next_reneging_time) # 이벤트: 고객 reneging
    inter_event_time = system._next_reneging_time - system._sim_time # inter event time 저장
    system._sim_time = system._next_reneging_time # 시스템의 현재 시각 변경
    deleteat!(system._customers_in_system, system._next_reneging_customer_index) # 리니깅하는 고객을 system에서 제거
    system._number_of_customers -= 1 # 시스템 내 고객 수 1 감소
    if system._number_of_customers >= 2 # 리니깅한 직후에, 시스템에 고객이 2명 이상 남아 있으면
      # next_reneging_time 을 바꿔야한다
      x = Float64[]
      for i = 2:system._number_of_customers
        push!(x, system._customers_in_system[i]._reneging_time)
      end
      y = findmin(x)
      system._next_reneging_time = y[1] # next reneging time 갱신
      system._next_reneging_customer_index = y[2] + 1 # next reneging customer 갱신
    else # 서비스가 끝나서 고객이 나간 직후에, 시스템에 남아있는 고객이 딱 한 명이면,
      system._next_reneging_time = typemax(Float64) # 다음 reneging time은 ∞
      system._next_reneging_customer_index = 0
    end
    system._next_completion_time = system._sim_time + rand(Exponential(1/4)) # 다음 completion 시각 설정 (exponential이므로 재설정)
    println(f, "(Time: $(system._sim_time)) Current event: Reneging (remaining # of customers: $(system._number_of_customers))")
    recording._number_of_reneging += 1
    push!(recording._event_times, system._sim_time)
    push!(recording._number_in_system, system._number_of_customers)
  end
end

function run_simulation(T::Float64)
  system = queue(PP(5.0,T*2))
  X = record()
  while system._sim_time < T
    next_event(system, X)
  end
  return X
end

sample_path = run_simulation(100.0)
close(f)

# plotting X(t) process

plt.step(sample_path._event_times,sample_path._number_in_system)
plt.xlim(0.0,10.0)
plt.savefig("MM1 queue with impatient customers.pdf")
function do_replication(n::Int64)
  cusum = 0.0
  for i in 1:n
    cusum += run_simulation(100.0)._number_of_reneging
  end
  println("Expected number of lost customers by time 100: $(cusum/n)")
end

do_replication(500)

close(f)
