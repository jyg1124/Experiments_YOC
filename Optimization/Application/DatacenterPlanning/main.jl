cd(dirname(Base.source_path()))
include("./Types.jl")
# for test
using CPLEX, JuMP

input_file_dir = "./"
input_file_name = "input(20).txt"
input_file_path = string(input_file_dir,input_file_name)

I = Input(input_file_path)

## Declare Model
m = Model(solver = CplexSolver(CPXPARAM_TimeLimit = 30))

## Set Variables
@variables m begin
    A[i = 1:I.nTypeApp , j = 1:I.nServer], Bin
    X[j = 1:I.nServer, f = 0:I.nFreqForEachTypeServer[1]-1, t = 1:I.nTimeSlot], Bin
    R[i = 1:I.nTypeApp , j = 1:I.nServer , t = 1:I.nTimeSlot], (Cont, lowerbound = 0, upperbound = 1)
end

## Set Constraints
U = I.maxApp[1]
Λ = vcat(I.Workload...)
C = vec(I.Capacity[1])
ρ = I.TargetLoad
@constraint(m, constraint2[j = 1:I.nServer], sum(A[:,j]) <= U)
@constraint(m, constraint3[i = 1:I.nTypeApp , t = 1:I.nTimeSlot], sum(R[i,:,t]) == 1)
@constraint(m, constraint4[i = 1:I.nTypeApp , j = 1:I.nServer , t = 1:I.nTimeSlot], R[i,j,t] <= A[i,j])
@constraint(m, constraint5[j = 1:I.nServer , t = 1:I.nTimeSlot], AffExpr(R[:,j,t], Λ[:,t], 0.0) <= ρ*AffExpr(X[j,:,t], C, 0.0))
@constraint(m, constraint6[j = 1:I.nServer , t = 1:I.nTimeSlot], sum(X[j,:,t]) == 1)
@constraint(m, constraint7[j = 1:I.nServer , t = 1:I.nTimeSlot], sum(R[:,j,t]) <= U*AffExpr([X[j,0,t]], [-1], 1))

## Set Objective
Β = vec(I.Cost[1])
@objective(m, Min, sum(Β[f+1]*X[j,f,t] for j = 1:I.nServer, f = 0:I.nFreqForEachTypeServer[1]-1, t = 1:I.nTimeSlot) )
## Solve
solve(m)
println("Obj = $(getobjectivevalue(m))")
