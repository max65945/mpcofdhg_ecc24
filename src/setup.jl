using Revise, LinearAlgebra
#=====================================================================
    User Input & Model Setup
=====================================================================#
# physical parameters
c_p = 4.180      # specific heat capacity of water [kJ/kg/K]
λ   = 0.02       # pipe friction coefficient [-]
ρ   = 998.05     # density of water [kg/m^3]
κ   = 2. * 10^-1 # heat loss coefficient [kW/K]
T_a = 8.0        # ambient temperature [°C]

# design parameters
L           = 500.0      # pipe length [m]
Δp_max_spec = 300.0      # maximal specific pressure drop [N/m^3]
m_tes_tot   = [50., 30.] # total mass of TESs [t]
m_tes_h_n   = [25., 15.] # nominal mass of TESs hot layer [t]
v_e_min     = 0.8        # minimal velocity of water [m/s]
v_e_max     = 3.0        # minimal velocity of water [m/s]

# graph parameters
n_d     = 3         # number of consumer
e_d     = [3, 6, 8] # edges associates with consumer
n_pr    = 2         # number of producer
e_pr    = [1, 7]    # edges associated with producer
n_v     = 7         # number of vertices
n_tes   = 2         # number of TES
v_tes_h = [1, 6]    # vertices associated with TESs hot layer
v_tes_c = [4, 7]    # vertices associated with TESs cold layer
n_e     = 9         # number of edges
n_c     = 5         # number of chords
n_x     = n_tes + n_v + n_e # number of states
n_u     = n_c + n_pr        # number of control inputs    

# vertexode-edge incidence matrix
B_0 = [ 1 -1  0  0  0  0  0  0  0; 
        0  1 -1  0 -1  0 -1  0  0;
        0  0  1 -1  0  1  0  0  0;
       -1  0  0  1  0  0  0  0  0;
        0  0  0  0  1 -1  0  0  1;
        0  0  0  0  0  0  1 -1  0;
        0  0  0  0  0  0  0  1 -1];
B_0_plus = 0.5 * (abs.(B_0) + B_0)
B_0_minus = 0.5 * (abs.(B_0) - B_0)

# Fundamental loop matrix for steady state with set of chords \mc{F}_{\mrm{ss}} = {q_3, q_6, q_8} of original graph \mc{G}
F_ss   = [1 1 1 1 0 0 0 0 0; 
          1 1 0 1 1 1 0 0 0;
          1 1 0 1 0 1 1 1 1];

F = [1 0 0 0 0 0 0 0 0; # Fundamental loop matrix
     0 1 1 1 0 0 0 0 0;
     0 1 0 1 1 1 0 0 0;
     0 1 0 1 0 1 1 0 1;
     0 0 0 0 0 0 0 1 0];
     
# rows of vertex-edge incidence matrix associated with hot layer of TES
B_tes = zeros(n_tes,n_e)
[B_tes[i,:] = B_0[v_tes_h[i],:] for i in eachindex(v_tes_h)]

# control input matrix
B = [zeros(n_v, n_pr); zeros(n_e, n_pr)] # initialize empty control input matrix
[B[n_v + e_pr[i], i] = 1.0 for i in eachindex(e_pr)] # assign row of B matrix to corresponding producer
B = 10.0^3 .* c_p^-1 .* B # convert B matrix to unit

# disturbance input matrix
E = [zeros(n_v, n_d+1); zeros(n_e, n_d+1)] # initialize empty disturbance input matrix
[E[n_v + e_d[i], i] = 1.0 for i in eachindex(e_d)] # assign row of E matrix to corresponding consumer
E = 10.0^3 .* c_p^-1 .* E # convert E matrix to unit
E[:,n_d+1] .= c_p^-1 * κ  # add heat loss entry

# nominal operation parameters
T_sup_n_l = 45.0    # nominal supply temperature of low temperature layer [°C]
T_ret_n_l = 20.0    # nominal return temperature of low temperature layer [°C]
ΔT_n_l    = T_sup_n_l - T_ret_n_l # nominal temperature difference of low temperature layer [K]
T_sup_max_l = 1.045 * T_sup_n_l    # maximal supply temperature of low temperature layer [°C]
T_ret_max_l = 1.3 * T_ret_n_l    # minimal return temperature of low temperature layer [°C]

T_sup_n_h = 75.0            # nominal supply temperature of high temperature layer [°C]
T_ret_n_h = copy(T_sup_n_l) # nominal return temperature of high temperature layer [°C]
ΔT_n_h    = T_sup_n_h - T_ret_n_h # nominal temperature difference of high temperature layer [K]
T_sup_max_h = 1.06 * T_sup_n_h    # maximal supply temperature of high temperature layer [°C]
T_ret_max_h = 1.3 * T_ret_n_h    # minimal return temperature of high temperature layer [°C]

ΔT_n   = [ΔT_n_l; ΔT_n_l; ΔT_n_h] # nominal temperature difference at consumers [K]
P_d_n  = [- 1.5; -2.5; -1.];      # nominal heat demand of consumers [MW]

q_c_n  = round.(-10.0^3 * P_d_n ./ (c_p .* ΔT_n); digits = 3) # nominal mass flows over chosen set of chords \mc{F}_{\mrm{ss}} = {q_3, q_6, q_8} of original graph \mc{G} [kg/s]
q_n    = round.(F_ss' * q_c_n; digits = 3)                                 # nominal mass flows over all edges [kg/s]

P_pr_n = round.(10^-3 .* c_p .* q_n[[1,7]] .* [ΔT_n_l;ΔT_n_h]; digits = 3) # nominal produced heat flow [MW]
P_pr_max = [1.1 * P_pr_n[1]; 1.45 * P_pr_n[2]]

Δp_max = Δp_max_spec * L # maximal pressure drop per pipe [N/m^2]
D      = round.(((8.0 * λ * L * q_n.^2) / (ρ * π^2 * Δp_max)).^(1/5); digits = 3) # pipe diameter designed for nominal conditions [m]
A      = π * D.^2 / 4    # pipe's cross sectional area [m^2]

m_p = round.(10^-3 .* A .* L * ρ; digits = 3) # pipe mass [t]

q_e_min = A * v_e_min * ρ # minimal mass flow [kg/s]
q_e_max = A * v_e_max * ρ # maximal mass flow [kg/s]
v       = q_n ./ (ρ .* A) # velocity of water [m/s]

# assign mass of pipe to corresponding edge and sink vertex via
m_v = round.(0.5 * B_0_plus * m_p; digits = 3) # [t]
m_e = round.(0.5 * m_p; digits = 3)            # [t]

# add TES mass to corresponding vertices to obtain vertex mass in nominal operation
m_v_n = copy(m_v)
for i in eachindex(v_tes_h)
    m_v_n[v_tes_h[i]] = m_tes_h_n[i]
    m_v_n[v_tes_c[i]] = m_tes_tot[i] - m_tes_h_n[i]
end

# thermal flow matrix for nominal operation
A_thermal_n = [-diagm(B_0_plus*F_ss'*q_c_n .+ κ) B_0_plus*diagm(F_ss'*q_c_n);
                diagm(F_ss'*q_c_n)*B_0_minus' -diagm(F_ss'*q_c_n .+ κ)]

Θ_n            = 10.0^3 * diagm([m_v_n;m_e]) # thermal inertia matrix [kg]
A_thermal_n = inv(Θ_n) * A_thermal_n     

B_n = inv(Θ_n) * B # control input matrix for nominal operation
E_n = inv(Θ_n) * E # disturbance input matrix for nominal operation

# discretization via implicit midpoint rule
Δt = 60.0 # step time [s]
#A_thermal_n_d = inv(I(n_v+n_e)-Δt/2*A_thermal_n) * (I(n_v+n_e)+Δt/2*A_thermal_n)
#E_thermal_n_d = sqrt(Δt) * inv(I(n_v+n_e)-Δt/2*A_thermal_n) * [B_n E_n]

#=====================================================================
    Derivation of 2nd Set Point
=====================================================================#
T_sup_l_2 = 46.0    # nominal supply temperature of low temperature layer [°C]
T_ret_l_2 = 20.0    # nominal return temperature of low temperature layer [°C]
ΔT_l_2    = T_sup_l_2 - T_ret_l_2 # nominal temperature difference of low temperature layer [K]

T_sup_h_2 = 77.0      # nominal supply temperature of high temperature layer [°C]
T_ret_h_2 = copy(T_sup_l_2) # nominal return temperature of high temperature layer [°C]
ΔT_h_2    = T_sup_h_2 - T_ret_h_2 # nominal temperature difference of high temperature layer [K]

ΔT_2   = [ΔT_l_2; ΔT_l_2; ΔT_h_2] # nominal temperature difference at consumers [K]
#P_d_2  = [-1.5; -2.; -1.];        # nominal heat demand of consumers [MW]
P_d_2  = [-2.35; -2.5; -1.22];        # nominal heat demand of consumers [MW]
q_c_2  = round.(-10.0^3 * P_d_2 ./ (c_p .* ΔT_2); digits = 3)              # nominal mass flows over chosen set of chords \mc{F}_{\mrm{ss}} = {q_3, q_6, q_8} of original graph \mc{G} [kg/s]
q_2    = round.(F_ss' * q_c_2; digits = 3)                                 # nominal mass flows over all edges [kg/s]
P_pr_2 = round.(10^-3 .* c_p .* q_2[[1,7]] .* [ΔT_l_2;ΔT_h_2]; digits = 3) # nominal produced heat flow [MW]

m_v_2 = copy(m_v)
#m_tes_h_2   = [22.5, 16.]; # mass of TESs hot layer for 2nd set point [t]
m_tes_h_2   = [25.5, 15.3]; # mass of TESs hot layer for 2nd set point [t]
for i in eachindex(v_tes_h)
    m_v_2[v_tes_h[i]] = m_tes_h_2[i]
    m_v_2[v_tes_c[i]] = m_tes_tot[i] - m_tes_h_2[i]
end

# thermal flow matrix for 2nd set point
A_thermal_2 = [-diagm(B_0_plus*F_ss'*q_c_2 .+ κ) B_0_plus*diagm(F_ss'*q_c_2);
                diagm(F_ss'*q_c_2)*B_0_minus' -diagm(F_ss'*q_c_2 .+ κ)]
Θ_2            = 10.0^3 * diagm([m_v_2;m_e]) # thermal inertia matrix [kg]
A_thermal_2 = inv(Θ_2) * A_thermal_2  

B_2 = inv(Θ_2) * B
E_2 = inv(Θ_2) * E

# discretization via implicit midpoint rule
#A_thermal_d_2 = inv(I(n_v+n_e)-Δt/2*A_thermal_2) * (I(n_v+n_e)+Δt/2*A_thermal_2)
#E_thermal_d_2 = sqrt(Δt) * inv(I(n_v+n_e)-Δt/2*A_thermal_2) * [B_2 E_2]