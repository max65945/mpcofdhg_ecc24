using Plots, CSV, DataFrames

#TODO: Add limits and steady states, change step time or lengthen sim
plot(layout = (4,2), size = (800,1000))

plot!(subplot = 1, x_MPC_1[1,1,:], seriestype = :steppost)
plot!(subplot = 1, x_MPC_2[1,1,:], seriestype = :steppost)
plot!(subplot = 1, [0;k_end], [x_max[1]; x_max[1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 1, [0;k_step], [x_ss[1,1]; x_ss[1,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 1, [k_step;k_end], [x_ss[1,2]; x_ss[1,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 3, u_MPC_1[1,1,:], seriestype = :steppost)
plot!(subplot = 3, u_MPC_1[2,1,:], seriestype = :steppost)
plot!(subplot = 3, u_MPC_2[1,1,:], seriestype = :steppost)
plot!(subplot = 3, u_MPC_2[2,1,:], seriestype = :steppost)
plot!(subplot = 3, [0;k_end], [u_max[1]; u_max[1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 3, [0;k_step], [u_ss[1,1]; u_ss[1,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 3, [k_step;k_end], [u_ss[1,2]; u_ss[1,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 5, x_MPC_1[n_tes+1,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+2,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+5,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+1,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+2,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+5,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+9,1,:], seriestype = :steppost)
plot!(subplot = 5, x_MPC_2[n_tes+1,1,:], seriestype = :steppost)
plot!(subplot = 5, [0;k_end], [x_max[n_tes+1]; x_max[n_tes+1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 5, [0;k_step], [x_ss[n_tes+1,1]; x_ss[n_tes+1,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 5, [k_step;k_end], [x_ss[n_tes+1,2]; x_ss[n_tes+1,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 7, u_MPC_1[10,1,:], seriestype = :steppost)
plot!(subplot = 7, u_MPC_2[10,1,:], seriestype = :steppost)
plot!(subplot = 7, [0;k_end], [u_max[10]; u_max[10]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 7, [0;k_step], [u_ss[10,1]; u_ss[10,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 7, [k_step;k_end], [u_ss[10,2]; u_ss[10,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 2, x_MPC_1[2,1,:], seriestype = :steppost)
plot!(subplot = 2, x_MPC_2[2,1,:], seriestype = :steppost)
plot!(subplot = 2, [0;k_end], [x_max[2]; x_max[2]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 2, [0;k_step], [x_ss[2,1]; x_ss[2,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 2, [k_step;k_end], [x_ss[2,2]; x_ss[2,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 4, u_MPC_1[7,1,:], seriestype = :steppost)
plot!(subplot = 4, u_MPC_1[8,1,:], seriestype = :steppost)
plot!(subplot = 4, u_MPC_2[7,1,:], seriestype = :steppost)
plot!(subplot = 4, u_MPC_2[8,1,:], seriestype = :steppost)
plot!(subplot = 4, [0;k_end], [u_max[7]; u_max[7]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 4, [0;k_step], [u_ss[7,1]; u_ss[7,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 4, [k_step;k_end], [u_ss[7,2]; u_ss[7,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 6, x_MPC_1[n_tes+6,1,:], seriestype = :steppost)
plot!(subplot = 6, x_MPC_2[n_tes+6,1,:], seriestype = :steppost)
plot!(subplot = 6, [0;k_end], [x_max[n_tes+6]; x_max[n_tes+6]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 6, [0;k_step], [x_ss[n_tes+6,1]; x_ss[n_tes+6,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 6, [k_step;k_end], [x_ss[n_tes+6,2]; x_ss[n_tes+6,2]], seriestype = :steppost, linestyle = :dash, color = :black)

plot!(subplot = 8, u_MPC_1[11,1,:], seriestype = :steppost)
plot!(subplot = 8, u_MPC_2[11,1,:], seriestype = :steppost)
plot!(subplot = 8, [0;k_end], [u_max[11]; u_max[11]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 8, [0;k_step], [u_ss[11,1]; u_ss[11,1]], seriestype = :steppost, linestyle = :dash, color = :black)
plot!(subplot = 8, [k_step;k_end], [u_ss[11,2]; u_ss[11,2]], seriestype = :steppost, linestyle = :dash, color = :black)

# prepare data vectors for csv export
x_MPC_1_cl = zeros(n_x,N_sim);
x_MPC_2_cl = zeros(n_x,N_sim);
x_ss_cl = zeros(n_x,N_sim);

x_max_cl = zeros(n_x,N_sim);
x_min_cl = zeros(n_x,N_sim);

u_MPC_1_cl = zeros(n_e + n_pr,N_sim);
u_MPC_2_cl = zeros(n_e + n_pr,N_sim);
u_ss_cl = zeros(n_e + n_pr,N_sim);

u_max_cl = zeros(n_e + n_pr,N_sim);
u_min_cl = zeros(n_e + n_pr,N_sim);

for i = 1:N_sim
    x_MPC_1_cl[:,i] = x_MPC_1[:,1,i];
    u_MPC_1_cl[:,i] = u_MPC_1[:,1,i];
    x_MPC_2_cl[:,i] = x_MPC_2[:,1,i];
    u_MPC_2_cl[:,i] = u_MPC_2[:,1,i];
    
    x_max_cl[:,i] = x_max
    x_min_cl[:,i] = x_min

    u_max_cl[:,i] = u_max
    u_min_cl[:,i] = u_min

    if i <= N_step-1
        x_ss_cl[:,i] = x_ss[:, 1];
        u_ss_cl[:,i] = u_ss[:, 1];
    else
        x_ss_cl[:,i] = x_ss[:, 2];
        u_ss_cl[:,i] = u_ss[:, 2];
    end
end

Δu_MPC_1_cl_max = zeros(length(u_MPC_1_cl[:,1]))
for k = 1:length(u_MPC_1_cl[1,:])-1
    for i = 1:length(u_MPC_1_cl[:,1])
        if abs(u_MPC_1_cl[i,k+1]-u_MPC_1_cl[i,k]) > Δu_MPC_1_cl_max[i]
            Δu_MPC_1_cl_max[i] = abs(u_MPC_1_cl[i,k+1]-u_MPC_1_cl[i,k])
        end
    end
end

Δu_MPC_2_cl_max = zeros(length(u_MPC_2_cl[:,1]))
for k = 1:length(u_MPC_2_cl[1,:])-1
    for i = 1:length(u_MPC_2_cl[:,1])
        if abs(u_MPC_2_cl[i,k+1]-u_MPC_2_cl[i,k]) > Δu_MPC_2_cl_max[i]
            Δu_MPC_2_cl_max[i] = abs(u_MPC_2_cl[i,k+1]-u_MPC_2_cl[i,k])
        end
    end
end

global J_cl_MPC_1 = 0
global J_cl_MPC_2 = 0
for k = 1:N_sim
    global J_cl_MPC_1 += J_cl_MPC_1 + (x_MPC_1_cl[:,k] - x_ss_cl[:,k])' * Q * (x_MPC_1_cl[:,k] - x_ss_cl[:,k]) + (u_MPC_1_cl[:,k] - u_ss_cl[:,k])' * R * (u_MPC_1_cl[:,k] - u_ss_cl[:,k])
    global J_cl_MPC_2 += J_cl_MPC_2 + (x_MPC_2_cl[:,k] - x_ss_cl[:,k])' * Q * (x_MPC_2_cl[:,k] - x_ss_cl[:,k]) + (u_MPC_2_cl[:,k] - u_ss_cl[:,k])' * R * (u_MPC_2_cl[:,k] - u_ss_cl[:,k])
end

t = (1:Int(N_sim)) ./ 60

df = DataFrame(time = t, 
            x1_sol_1 = x_MPC_1_cl[1,1:Int(N_sim)],
            x1_sol_2 = x_MPC_1_cl[2,1:Int(N_sim)],
            x1_sol_3 = x_MPC_1_cl[3,1:Int(N_sim)],
            x1_sol_8 = x_MPC_1_cl[8,1:Int(N_sim)],
            x1_sol_15 = x_MPC_1_cl[15,1:Int(N_sim)],

            x2_sol_1 = x_MPC_2_cl[1,1:Int(N_sim)],
            x2_sol_2 = x_MPC_2_cl[2,1:Int(N_sim)],
            x2_sol_3 = x_MPC_2_cl[3,1:Int(N_sim)],
            x2_sol_8 = x_MPC_2_cl[8,1:Int(N_sim)],

            x_ss_1 = x_ss_cl[1,1:Int(N_sim)],
            x_ss_2 = x_ss_cl[2,1:Int(N_sim)],
            x_ss_3 = x_ss_cl[3,1:Int(N_sim)],
            x_ss_8 = x_ss_cl[8,1:Int(N_sim)],

            x_max_1 = x_max_cl[1,1:Int(N_sim)],
            x_max_2 = x_max_cl[2,1:Int(N_sim)],
            x_max_3 = x_max_cl[3,1:Int(N_sim)],
            x_max_8 = x_max_cl[8,1:Int(N_sim)],

            x_min_1 = x_min_cl[1,1:Int(N_sim)],
            x_min_2 = x_min_cl[2,1:Int(N_sim)],
            x_min_3 = x_min_cl[3,1:Int(N_sim)],
            x_min_8 = x_min_cl[8,1:Int(N_sim)],

            u1_sol_1 = u_MPC_1_cl[1,1:Int(N_sim)],
            u1_sol_2 = u_MPC_1_cl[2,1:Int(N_sim)],
            u1_sol_7 = u_MPC_1_cl[7,1:Int(N_sim)],
            u1_sol_8 = u_MPC_1_cl[8,1:Int(N_sim)],
            u1_sol_10 = u_MPC_1_cl[10,1:Int(N_sim)],
            u1_sol_11 = u_MPC_1_cl[11,1:Int(N_sim)],

            u2_sol_1 = u_MPC_2_cl[1,1:Int(N_sim)],
            u2_sol_2 = u_MPC_2_cl[2,1:Int(N_sim)],
            u2_sol_7 = u_MPC_2_cl[7,1:Int(N_sim)],
            u2_sol_8 = u_MPC_2_cl[8,1:Int(N_sim)],
            u2_sol_10 = u_MPC_2_cl[10,1:Int(N_sim)],
            u2_sol_11 = u_MPC_2_cl[11,1:Int(N_sim)],

            u_ss_1 = u_ss_cl[1,1:Int(N_sim)],
            u_ss_7 = u_ss_cl[7,1:Int(N_sim)],
            u_ss_10 = u_ss_cl[10,1:Int(N_sim)],
            u_ss_11 = u_ss_cl[11,1:Int(N_sim)],

            u_max_1 = u_max_cl[1,1:Int(N_sim)],
            u_max_7 = u_max_cl[7,1:Int(N_sim)],
            u_max_10 = u_max_cl[10,1:Int(N_sim)],
            u_max_11 = u_max_cl[11,1:Int(N_sim)],
            u_min_1 = u_min_cl[1,1:Int(N_sim)],
            u_min_7 = u_min_cl[7,1:Int(N_sim)],
            u_min_10 = u_min_cl[10,1:Int(N_sim)],
            u_min_11 = u_min_cl[11,1:Int(N_sim)],

        );

CSV.write("./rawData/results_QIHMPC_250324.csv", df);