using Revise, LinearAlgebra, JuMP, Ipopt, ForwardDiff, ControlSystems
include("_functions/linearization_functions.jl")
#= 
-----------------------------------------------------------------------------------------------------------------------------------------------------------
    Calculate Terminal Ingredients
-----------------------------------------------------------------------------------------------------------------------------------------------------------
=#
A_lin = zeros(n_x, n_x, 2)
B_lin = zeros(n_x, n_e+n_pr, 2)

p = n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E
A_lin[:, :, 1] = ForwardDiff.jacobian(x -> dynamic_DHG(x, u_ss[:,1], d_s[:,1], p), x_ss[:,1])
B_lin[:, :, 1] = ForwardDiff.jacobian(u -> dynamic_DHG(x_ss[:,1], u, d_s[:,1], p), u_ss[:,1])
checkStabilizability(A_lin[:, :, 1], B_lin[:, :, 1])

A_lin[:, :, 2] = ForwardDiff.jacobian(x -> dynamic_DHG(x, u_ss[:,2], d_s[:,2], p), x_ss[:,2])
B_lin[:, :, 2] = ForwardDiff.jacobian(u -> dynamic_DHG(x_ss[:,2], u, d_s[:,2], p), u_ss[:,2])
checkStabilizability(A_lin[:, :, 2], B_lin[:, :, 2])

Q = 10.0^0 ./ x_ss[:,1];
Q[1] = Q[1]*10.0^1
Q[2] = Q[2]*10.0^1
Q = diagm(Q);
R = 10.0^0 ./ u_ss[:,1];
R = diagm(R);

K = zeros(n_e+n_pr, n_x, 2)
K[:, :, 1] = -lqr(Continuous, A_lin[:, :, 1], B_lin[:, :, 1], Q, R);
K[:, :, 2] = -lqr(Continuous, A_lin[:, :, 2], B_lin[:, :, 2], Q, R);

AK = zeros(n_x, n_x, 2)
AK[:, :, 1] = A_lin[:, :, 1] + B_lin[:, :, 1] * K[:, :, 1];
AK[:, :, 2] = A_lin[:, :, 2] + B_lin[:, :, 2] * K[:, :, 2];

# Choose feasible κ
κ_MPC = zeros(2)
κ_MPC[1] = -0.95 * maximum(real(eigvals(AK[:, :, 1])));
checkKappa(AK[:, :, 1], κ_MPC[1])
κ_MPC[2] = -0.95 * maximum(real(eigvals(AK[:, :, 2])));
checkKappa(AK[:, :, 2], κ_MPC[2])

P = zeros(n_x, n_x, 2)
P[:, :, 1] = getP_lyap(AK[:, :, 1], K[:, :, 1], R, Q, κ_MPC[1])
P[:, :, 2] = getP_lyap(AK[:, :, 2], K[:, :, 2], R, Q, κ_MPC[2])