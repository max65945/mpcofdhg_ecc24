using Revise, LinearAlgebra, JuMP, Ipopt, ForwardDiff, ControlSystems, Statistics
include("./_functions/mpc_optimization_models.jl")
#=
-----------------------------------------------------------------------------------------------------------------------------------------------------------
    Begin: Closed-loop Simulation of MPC
-----------------------------------------------------------------------------------------------------------------------------------------------------------
=#
delta_t_MPC = 60.   # time step of MPC [s]
k_step = 90         # discrete-time step for set point change
k_end = 180         # amount of discrete-time steps simulated

T_MPC = 60*60
N_MPC = Int(round(T_MPC/delta_t_MPC))         # discrete-time prediction horizon of MPC

T_step = k_step*60
N_step = Int(round(T_step/delta_t_MPC))

T_sim = k_end*60;
N_sim = Int(round(T_sim/delta_t_MPC));         # discrete-time prediction horizon of MPC

x_MPC_1 = zeros(n_x, N_MPC-1, N_sim);
u_MPC_1 = zeros(n_e + n_pr, N_MPC-1, N_sim);
exec_time_MPC_1 = zeros(N_sim)

x_MPC_2 = zeros(n_x, N_MPC-1, N_sim);
u_MPC_2 = zeros(n_e + n_pr, N_MPC-1, N_sim);
exec_time_MPC_2 = zeros(N_sim)

global sol_cl = zeros(n_x,N_sim);

global x0 = copy(x_ss[:,1]);


global x0[1] = x0[1] - 1.0;
global x0[2] = x0[2] - 0.5;
global x0[n_tes+1] = x0[n_tes+1] - 1.0;
global x0[n_tes+6] = x0[n_tes+6] - 1.0;
global x0[n_tes+n_v+1] = x0[n_tes+n_v+1] - 0.0;

# MPC 1
for j = 1:N_sim
    
    p = n_x, n_pr, n_e, N_MPC, x_min, x_max, u_min, u_max, x_ss, u_ss, d_s, Q, R, P, α, x0, sol_cl, N_step
    global model_MPC, x_MPC, u_MPC = getMPCModel_1(j, p)
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Run Optimization
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    start_time = time()
    optimize!(model_MPC); # Run optimization
    end_time = time()

    exec_time_MPC_1[j] = end_time - start_time

    u_MPC_1[:,:,j] = value.(u_MPC[:,1:N_MPC-1]);
    x_MPC_1[:,:,j] = value.(x_MPC[:,1:N_MPC-1]);
    
    global sol_cl[:,j] = x_MPC_1[:,2,j];
    GC.gc()
end

avg_exec_time_MPC_1 = mean(exec_time_MPC_1)
max_exec_time_MPC_1 = maximum(exec_time_MPC_1)

println("############################################################### !!!!!!!!! ###############################################################")
println("############################################################### !!!!!!!!! ###############################################################")
println("############################################################### MPC1 done ###############################################################")
println("############################################################### !!!!!!!!! ###############################################################")
println("############################################################### !!!!!!!!! ###############################################################")

# MPC 2
for j = 1:N_sim
    
    p = n_x, n_pr, n_e, N_MPC, x_min, x_max, u_min, u_max, x_ss, u_ss, d_s, Q, R, P, α, x0, sol_cl, N_step
    global model_MPC, x_MPC, u_MPC = getMPCModel_2(j, p)
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Run Optimization
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    start_time = time()
    optimize!(model_MPC); # Run optimization
    end_time = time()

    exec_time_MPC_2[j] = end_time - start_time

    u_MPC_2[:,:,j] = value.(u_MPC[:,1:N_MPC-1]);
    x_MPC_2[:,:,j] = value.(x_MPC[:,1:N_MPC-1]);
    
    global sol_cl[:,j] = x_MPC_2[:,2,j];
    GC.gc()
end

avg_exec_time_MPC_2 = mean(exec_time_MPC_2)
max_exec_time_MPC_2 = maximum(exec_time_MPC_2)

println("############################################################### !!!!!!!!! ###############################################################")
println("############################################################### !!!!!!!!! ###############################################################")
println("############################################################### MPC2 done ###############################################################")
println("############################################################### !!!!!!!!! ###############################################################")
println("############################################################### !!!!!!!!! ###############################################################")