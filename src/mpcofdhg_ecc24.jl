module mpcofdhg_ecc24

include("./setup.jl")
include("./steadystate_calc.jl")
include("./linearization.jl")
include("./terminalingredients.jl")
include("./simulation_MPC.jl")
include("./postProcessing.jl")

for i in eachindex(x_ss[1,:])
    println("heat losses at ", i, "th steady state: ", round((sum(round.(u_ss[end-n_pr+1:end,i]; digits = 3))+sum(d_s[1:n_d,i])) / sum(round.(u_ss[end-n_pr+1:end,i]; digits = 3)) * 100.0, digits=3), " %")
end

println("Avg. exec. time MPC 1 = ", avg_exec_time_MPC_1, "s")
println("Max. exec. time MPC 1 = ", max_exec_time_MPC_1, "s")
println("Avg. exec. time MPC 2 = ", avg_exec_time_MPC_2, "s")
println("Max. exec. time MPC 2 = ", max_exec_time_MPC_2, "s")

println("Cost reduction by MPC 2 = ",100*(J_cl_MPC_1-J_cl_MPC_2)/J_cl_MPC_1, "%")
println("mass flow input step reduction by MPC 2 = ", 100*(Δu_MPC_1_cl_max[1:n_e]-Δu_MPC_2_cl_max[1:n_e])./Δu_MPC_1_cl_max[1:n_e], "%")
println("heat input step reduction by MPC 2 = ", 100*(Δu_MPC_1_cl_max[end-n_pr+1:end]-Δu_MPC_2_cl_max[end-n_pr+1:end])./Δu_MPC_1_cl_max[end-n_pr+1:end], "%")

end # module mpcofdhg_ecc24