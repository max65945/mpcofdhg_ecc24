using Revise, LinearAlgebra, JuMP, Ipopt
include("_functions/get_thermalsteadystate.jl")
#=====================================================================
    Steady State Optimization
=====================================================================#

# initialize set point arrays for 2 set points
x_s = zeros(n_v+n_e, 2)
u_s = zeros(n_pr, 2)
d_s = zeros(n_d+1, 2)

# define nominal operating point and 2nd set point as ideal set points
# 1st ideal set point
x_s[:,1] = [T_sup_n_l,T_sup_n_l,T_ret_n_l,T_ret_n_l,T_sup_n_l,T_sup_n_h,T_ret_n_h,T_sup_n_l,T_sup_n_l,T_ret_n_l,T_ret_n_l,T_sup_n_l,T_ret_n_l,T_sup_n_h,T_ret_n_h,T_ret_n_h] # desired steady state
u_s[:,1] = copy(P_pr_n) 
d_s[:,1] = [P_d_n; T_a]
# 2nd ideal set point
x_s[:,2] = [T_sup_l_2,T_sup_l_2,T_ret_l_2,T_ret_l_2,T_sup_l_2,T_sup_h_2,T_ret_h_2,T_sup_l_2,T_sup_l_2,T_ret_l_2,T_ret_l_2,T_sup_l_2,T_ret_l_2,T_sup_h_2,T_ret_h_2,T_ret_h_2] # desired steady state
u_s[:,2] = copy(P_pr_2)
d_s[:,2] = [P_d_2; T_a]

# define box constraints
x_max = [1.03 * m_v_n[1]; 1.03 * m_v_n[6]; T_sup_max_l; T_sup_max_l; T_ret_max_l; T_ret_max_l; T_sup_max_l; T_sup_max_h; T_ret_max_h; T_sup_max_l; T_sup_max_l; T_ret_max_l; T_ret_max_l; T_sup_max_l; T_ret_max_l; T_sup_max_h; T_ret_max_h; T_ret_max_h]
x_min = 0.0 * x_max

u_max = [q_e_max; P_pr_max]
u_min = [q_e_min; 0.0 * P_pr_max]

# initialize steady state arrays for 2 steady states
global x_ss = zeros(n_v+n_e, 2)
global u_ss = zeros(n_pr, 2)

Q_ss = 10.0^3 .* diagm(1 ./ x_s[:,1]) # weight for state deviation
R_ss = 10.0^0 .* diagm(1 ./ u_s[:,1]) # weight for control input deviation

p = n_x, n_tes, n_u, n_c, Δt, Q_ss, R_ss, x_max, x_min, u_max, u_min
x_ss[:,1], u_ss[:,1] = get_thermalsteadystate(x_s[:,1], u_s[:,1], d_s[:,1], A_thermal_n, [B_n E_n], p)
x_ss[:,2], u_ss[:,2] = get_thermalsteadystate(x_s[:,2], u_s[:,2], d_s[:,2], A_thermal_2, [B_2 E_2], p)


u_ss = [q_n q_2; u_ss]
x_ss = [m_v_n[1] m_v_2[1]; m_v_n[6] m_v_2[6]; x_ss]

for i in eachindex(x_ss[1,:])
    if maximum(x_ss[:,i]-x_max)>=0 || maximum(x_ss[:,i]-x_min)<=0
        error((i,"th steady state not reachable!"))
    else
        display((i,"th steady state reachable!")) 
    end
end

for i in eachindex(u_ss[1,:])
    if maximum(u_ss[:,i]-u_max)>=0 || maximum(u_ss[:,i]-u_min)<=0
        error((i,"th steady input not reachable!"))
    else
        display((i,"th steady input reachable!")) 
    end
end
#TODO: Get 2nd steady state closer to bounds
#TODO: Check if 2nd steady state is feasible