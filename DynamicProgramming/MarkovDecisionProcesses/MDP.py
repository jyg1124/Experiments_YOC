# This python code is written by Yongkyu Cho.

import numpy as np
import random as rd
from scipy import linalg

# Define functions
def g(s,a):
    if s == "Standing":
        if a == "Slow":
            return 1
        elif a == "Fast":
            return 0.8
    elif s == "Moving":
        if a == "Slow":
            return 1
        elif a == "Fast":
            return 1.4
    else:
        if a == "Slow":
            return -0.2

def g_d(S,d):
    temp_list = []
    for s in S:
        temp_list.append(g(s,d[s]))
    return np.array(temp_list)

def p(j,s,a):
    if s == "Standing":
        if a == "Fast":
            if j == "Fallen":
                return 0.4
            elif j == "Moving":
                return 0.6
            else:
                return 0.0
        elif a == "Slow":
            if j == "Moving":
                return 1.0
            else:
                return 0.0
    elif s == "Moving":
        if a == "Fast":
            if j == "Fallen":
                return 0.2
            elif j == "Moving":
                return 0.8
            else:
                return 0.0
        elif a == "Slow":
            if j == "Moving":
                return 1.0
            else:
                return 0.0
    elif s == "Fallen":
        if a == "Slow":
            if j == "Standing":
                return 0.4
            elif j == "Fallen":
                return 0.6
            else:
                return 0.0
        else:
            return 0.0

def P_d(S, d):
    temp_list_1 = []
    for s_1 in S:
        temp_list_2 = []
        for s_2 in S:
            temp_list_2.append(p(s_2,s_1,d[s_1]))
        temp_list_1.append(temp_list_2)
    return np.array(temp_list_1)

def stationary_policies(S,A):
    D = []
    for i in range(len(A["Standing"])):
        for j in range(len(A["Moving"])):
            for k in range(len(A["Fallen"])):
                D.append({"Standing":A["Standing"][i],"Moving":A["Moving"][j],"Fallen":A["Fallen"][k]})
    return D

def sup_norm(V):
    temp = 0.0
    for v in V:
        if abs(v) >= temp:
            temp = abs(v)
    return temp

def policy_evaluation(S,d,δ):
    num_states = len(S)
    return linalg.solve((np.identity(num_states)-δ*P_d(S,d)),g_d(S,d))

def policy_iteration(S,A,δ):
    D = stationary_policies(S,A)
    num_states = len(S)
    n = 1
    old_d = D[0]
    while True:
        # print("Iteration ", n)
        V = policy_evaluation(S,old_d,δ)
        new_d = {}
        temp_values = np.array([-float("inf") for i in range(num_states)])
        for d in D:
            current_values = g_d(S,d)+δ*np.dot(P_d(S,d),V.T)
            if (current_values >= temp_values).all() == True:
                if (current_values == temp_values).all() == False:
                    new_d = d
                    temp_values = current_values
        if old_d == new_d:
            break
        else:
            old_d = new_d
        n += 1
    #print("Number of Iteration: ", n)
    return (policy_evaluation(S,old_d,δ).tolist() , old_d, n)


def value_iteration(S,A,δ,ϵ):
    n = 1
    num_states = len(S)
    V_old = [0.0 for i in range(num_states)]
    d_epsilon = {}
    while True:
        # print("Iteration ", n)
        V_new = []
        for i in range(num_states):
            temp_value = -float("inf")
            for a in A[S[i]]:
                current_value = g(S[i],a) + sum(δ*p(S[j],S[i],a)*V_old[j] for j in range(num_states))
                if current_value >= temp_value:
                    temp_value = current_value
            V_new.append(temp_value)
        if sup_norm(np.array(V_old)-np.array(V_new)) < (ϵ*(1-δ))/(2*δ):
            d = {}
            for s in S:
                temp_value = -float("inf")
                temp_action = ""
                for a in A[s]:
                    current_value = g(s,a) + sum(δ*p(S[j],s,a)*V_new[j] for j in range(num_states))
                    if current_value >= temp_value:
                        temp_value = current_value
                        temp_action = a
                d[s] = temp_action
            d_epsilon = d
            break
        else:
            V_old = V_new
        n += 1
    #print("Number of Iteration: ", n)
    return (V_new, d_epsilon, n)

def sample_next_state(S,state,action):
    rand = rd.random()
    d = {}
    for s in S:
        d[s]=p(s,state,action)
    for s in sorted(d,key = d.get):
        if rand < p(s,state,action):
            return s
        else:
            rand = rand - p(s,state,action)

def dict_to_list(d):
    return [item for item in d.values()]

def dict_to_array(d):
    return np.array([item for item in d.values()])

def max_Q(Q,A,state):
    temp_value = -float("inf")
    opt_action = ""
    for a in A[state]:
        if Q[(state,a)] >= temp_value:
            temp_value = Q[(state,a)]
            opt_action = a
    return (temp_value, opt_action)

def Q_learning_exploitation(S,A,δ,α,N):
    n = 1
    num_states = len(S)
    V_tilde = {s:0.0 for s in S}
    policy = {s:"" for s in S}
    Q = {(s,a):0.0 for s in S for a in A[s]}
    current_state = S[rd.randrange(0,num_states)]
    while n < N:
        #print("Iteration ", n)
        action = max_Q(Q,A,current_state)[1]
        next_state = sample_next_state(S, current_state, action)
        q_tilde = g(current_state,action) + δ*max_Q(Q,A,next_state)[0]
        Q[(current_state,action)] = Q[(current_state,action)] + α*(q_tilde - Q[(current_state,action)])
        V_tilde[current_state], policy[current_state] = max_Q(Q,A,current_state)
        current_state = next_state
        n += 1
    #print("Number of Iteration: ", n)
    return (dict_to_list(V_tilde), policy, n)

def Q_learning_epsilon_greedy(S,A,δ,α,N,ϵ):
    n = 1
    num_states = len(S)
    V_tilde = {s:0.0 for s in S}
    policy = {s:"" for s in S}
    Q = {(s,a):0.0 for s in S for a in A[s]}
    current_state = S[rd.randrange(0,num_states)]
    while n < N:
        #print("Iteration ", n)
        if rd.random() < ϵ: # w.p. ϵ, choose action by exploitation
            action = max_Q(Q,A,current_state)[1]
        else: # w.p. 1-ϵ, choose action by exploration
            action = rd.choice(A[current_state])
        next_state = sample_next_state(S, current_state, action)
        q_tilde = g(current_state,action) + δ*max_Q(Q,A,next_state)[0]
        Q[(current_state,action)] = Q[(current_state,action)] + α*(q_tilde - Q[(current_state,action)])
        V_tilde[current_state], policy[current_state] = max_Q(Q,A,current_state)
        current_state = next_state
        n += 1
    #print("Number of Iteration: ", n)
    return (dict_to_list(V_tilde), policy)

def policy_to_num(d):
    for i in range(len(D)):
        if d == D[i]:
            return i+1

# Set parameters
S = ["Standing","Moving","Fallen"] # state space
A = {"Standing":["Slow","Fast"], "Moving":["Slow","Fast"], "Fallen":["Slow"]} # action space
δ = 0.7 # discount factor
ϵ1 = 0.001 # for value iteration
ϵ2 = 0.5 # for epsilon greedy Q-learning
α = 0.05 # learning rate
N = 10000 # number of iteration

# Run algorithms
policy_iteration(S,A,δ)
value_iteration(S,A,δ,ϵ1)
Q_learning_exploitation(S,A,δ,α,N)
Q_learning_epsilon_greedy(S,A,δ,α,N,ϵ2)
