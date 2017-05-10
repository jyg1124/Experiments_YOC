using Distributions, PyPlot, JuMP, Ipopt
importall Types
importall Functions

const MAX_ARRIVALS = 200000
const WARM_UP_ARRIVALS = 50000

const REPLICATION_TIME = 5000.0
const WARM_UP_TIME = 0.3*REPLICATION_TIME

const REGULAR_UPDATE_INTERVAL = 0.01


# main part #

# Initialization
file_sim_record = open("sim_record.txt" , "w")
file_summarization = open("summarization.txt" , "w")
WS = workload_setter()
SS = server_setter()
S = server_creater(SS,WS)
AI = arrival_generator(WS,REPLICATION_TIME)
vdc = VirtualDataCenter(WS, AI, SS, WARM_UP_ARRIVALS, MAX_ARRIVALS, WARM_UP_TIME, REPLICATION_TIME, REGULAR_UPDATE_INTERVAL, S)
PI = Plot_Information(S,file_sim_record,file_summarization)

# Run
run_to_end(vdc, PI, REPLICATION_TIME, WARM_UP_TIME)      # until a certain number of services are completed

# Plotting Speeds and Prices
x = PI.time_array
plt = PyPlot

plt.figure(figsize = (15,10))
for j in 1:length(S)
  plt.subplot(2,5,j)
  plt.title("Server $j", fontsize=10)
  plt.xlabel("Time",fontsize=6)
  if j == 1 || j == 6
    plt.ylabel("Speed (workload/time)",fontsize=10)
  end
  plt.ylim(0,120)
  plt.plot(x,PI.speed_array[j][1:length(x)],linewidth=1.0,linestyle="-",color="red")
  plt.plot(x,PI.buffer_array[j][1:length(x)],linewidth=1.0,linestyle="--", color = "blue")
  plt.tick_params(labelsize=6)

end
plt.savefig("Long-run server speeds.pdf")

plt.figure(figsize = (15,10))
for j in 1:length(S)
  plt.subplot(2,5,j)
  plt.title("Server $j", fontsize=10)
  plt.xlabel("Time",fontsize=6)
  if j == 1 || j == 6
    plt.ylabel("Price",fontsize=10)
  end
  plt.plot(x,PI.price_array[j][1:length(x)],linewidth=1.0,linestyle="-",color="red")
  plt.tick_params(labelsize=6)
end
plt.savefig("Long-run server prices.pdf")

# Write summarization
for j = 1:10
  println(file_summarization, "P[W_$j>=δ_$j]: $(sum(PI.sojourn_time_array[j])/length(PI.sojourn_time_array[j]))")
end
sum(PI.sojourn_time_array[2])

# Closing IOStreams
close(file_summarization)
close(file_sim_record)
