using Plots

plot(layout = (6,1), size = (800,1000))

plot!(subplot = 1, x_MPC_1[1,1,:], seriestype = :steppost)
plot!(subplot = 1, x_MPC_2[1,1,:], seriestype = :steppost)

plot!(subplot = 2, u_MPC_1[1,1,:], seriestype = :steppost)
plot!(subplot = 2, u_MPC_1[2,1,:], seriestype = :steppost)
plot!(subplot = 2, u_MPC_2[1,1,:], seriestype = :steppost)
plot!(subplot = 2, u_MPC_2[2,1,:], seriestype = :steppost)

plot!(subplot = 3, x_MPC_1[n_tes+1,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+2,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+5,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+1,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+2,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+5,1,:], seriestype = :steppost)
#plot!(subplot = 3, x_MPC_1[n_tes+n_v+9,1,:], seriestype = :steppost)
plot!(subplot = 3, x_MPC_2[n_tes+1,1,:], seriestype = :steppost)

plot!(subplot = 4, u_MPC_1[10,1,:], seriestype = :steppost)
plot!(subplot = 4, u_MPC_2[10,1,:], seriestype = :steppost)

plot!(subplot = 5, x_MPC_1[n_tes+6,1,:], seriestype = :steppost)
plot!(subplot = 5, x_MPC_2[n_tes+6,1,:], seriestype = :steppost)

plot!(subplot = 6, u_MPC_1[11,1,:], seriestype = :steppost)
plot!(subplot = 6, u_MPC_2[11,1,:], seriestype = :steppost)