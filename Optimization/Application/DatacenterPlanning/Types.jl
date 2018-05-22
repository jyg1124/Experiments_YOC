
type Input
    nTimeSlot::Int64 # the number of time slots
    nTypeServer::Int64 # the number of types of servers
    nServer::Int64 # the number of servers in total
    nTypeApp::Int64 # the number of types of applications
    TargetLoad::Float64 # target load
    nForEachTypeServer::Array{Int64} # the number of servers for each type
    nFreqForEachTypeServer::Array{Int64} # the number of frequency for each type of server
    maxServer::Array{Int64} # the number of servers that can host each application
    maxApp::Array{Int64} # the maximum number of apps that can be installed for each type of server
    Capacity::Array{Array{Float64}} # capacity that matches each frequency
    Cost::Array{Array{Float64}} # cost that matches each frequency
    Workload::Array{Array{Float64}} # hourly recorded workload for each application
    function Input(input_file_path::String)
        file = open(input_file_path)
        lines = readlines(file)
        close(file)
        nline = 2
        nTimeSlot = parse(Int64, lines[nline])
        nline += 2
        nTypeServer = parse(Int64, lines[nline])
        nline += 2
        nServer = parse(Int64, lines[nline])
        nline += 2
        nTypeApp = parse(Int64, lines[nline])
        nline += 2
        TargetLoad = parse(Float64, lines[nline])
        nline += 2
        nForEachTypeServer = readdlm( IOBuffer(lines[nline]) , Int64)
        nline += 2
        nFreqForEachTypeServer = readdlm( IOBuffer(lines[nline]) , Int64)
        nline += 2
        maxServer = readdlm( IOBuffer(lines[nline]) , Int64)
        nline += 2
        maxApp = readdlm( IOBuffer(lines[nline]) , Int64)
        nline += 2
        Capacity = Array{Float64}[]
        for l in 1:nTypeServer
            push!(Capacity, readdlm( IOBuffer( lines[nline]) , Float64 ) )
            nline += 1
        end
        nline += 1
        Cost = Array{Float64}[]
        for l in 1:nTypeServer
            push!(Cost, readdlm(IOBuffer(lines[nline]),Float64))
            nline += 1
        end
        nline += 1
        Workload = Array{Float64}[]
        for l in 1:nTypeApp
            push!(Workload, readdlm(IOBuffer(lines[nline]),Float64))
            nline += 1
        end
        new(nTimeSlot, nTypeServer, nServer, nTypeApp, TargetLoad, nForEachTypeServer,
            nFreqForEachTypeServer, maxServer, maxApp, Capacity, Cost, Workload)
    end
end
