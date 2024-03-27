using Revise, LinearAlgebra, JuMP, Ipopt, ForwardDiff, ControlSystems
include("_functions/terminalingredients_functions.jl")

α_init = 10.0^8 # initial α-value 
p = n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E

α_1 = zeros(2)
α = zeros(2)
α[1], α_1[1] = get_alpha(P[:,:,1], κ_MPC[1], α_init, K[:,:,1], AK[:,:,1], u_ss[:,1], u_min, u_max, x_ss[:,1], x_min, x_max, d_s[:,1], p)
α[2], α_1[2] = get_alpha(P[:,:,2], κ_MPC[2], α_init, K[:,:,2], AK[:,:,2], u_ss[:,2], u_min, u_max, x_ss[:,2], x_min, x_max, d_s[:,2], p)

for k in eachindex(α)
    for j = 1:n_x
        global model_max_x = Model(Ipopt.Optimizer);     # define abstract opptimization model
        set_silent(model_max_x);

        @variables(model_max_x, begin
        max_x[1:n_x];
        end);

        for i = 1:n_x
            set_start_value(max_x[i], 0.1);
        end

        @constraint(model_max_x, (max_x - x_ss[:,k])' * P[:,:,k] * (max_x - x_ss[:,k]) <= α[k]);

        for i = 1:n_x
            set_lower_bound(max_x[i], x_min[i]);
            set_upper_bound(max_x[i], x_max[i]);
        end

        J = max_x[j];

        @objective(model_max_x, Max, J);

        optimize!(model_max_x)
        max_x_sol = value.(max_x);

        display(("Set point ", k, "; Δx(", j, ")", maximum(abs.((max_x_sol[j] - x_ss[:,1][j])))))
    end
end
