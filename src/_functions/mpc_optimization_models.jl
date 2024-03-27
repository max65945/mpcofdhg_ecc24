function getMPCModel_1(j, p)
    
    n_x, n_pr, n_e, N_MPC, x_min, x_max, u_min, u_max, x_ss, u_ss, d_s, Q, R, P, α, x0, sol_cl, N_step = p
    
    global model_MPC = Model(Ipopt.Optimizer);     # define abstract opptimization model_MPC
    set_silent(model_MPC)
    
    @variables(model_MPC, begin
    x_MPC[1:n_x,1:N_MPC];
    u_MPC[1:(n_pr + n_e),1:N_MPC-1];
    end);

    for k = 1:N_MPC-1
        for i = 1:n_x
            set_lower_bound(x_MPC[i,k], x_min[i])
            set_upper_bound(x_MPC[i,k], x_max[i])
        end
        
        for i = 1:(n_pr + n_e)
            set_lower_bound(u_MPC[i,k], u_min[i])
            set_upper_bound(u_MPC[i,k], u_max[i])
        end
    end

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Terminal state constraints
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------   
    if j >= N_step
        @constraint(model_MPC, (x_MPC[:,N_MPC]-x_ss[:,2])' * P[:, :, 2] * (x_MPC[:,N_MPC]-x_ss[:,2]) <= α[2]);
    else
        @constraint(model_MPC, (x_MPC[:,N_MPC]-x_ss[:,1])' * P[:, :, 1] * (x_MPC[:,N_MPC]-x_ss[:,1]) <= α[1]);
    end
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Cost function
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    global J_MPC = 0;
    # integral over stage costs
    for k=1:N_MPC-1
        if j >= N_step
            global J_MPC = J_MPC + (x_MPC[:,k]-x_ss[:,2])' * Q * (x_MPC[:,k]-x_ss[:,2]) + (u_MPC[:,k]-u_ss[:,2])' * R * (u_MPC[:,k]-u_ss[:,2]);
        else
            global J_MPC = J_MPC + (x_MPC[:,k]-x_ss[:,1])' * Q * (x_MPC[:,k]-x_ss[:,1]) + (u_MPC[:,k]-u_ss[:,1])' * R * (u_MPC[:,k]-u_ss[:,1])
        end
    end
    # terminal costs
    if j >= N_step
        global J_MPC = J_MPC + (x_MPC[:,N_MPC]-x_ss[:,2])' * P[:, :, 2] * (x_MPC[:,N_MPC]-x_ss[:,2]);
    else
        global J_MPC = J_MPC + (x_MPC[:,N_MPC]-x_ss[:,1])' * P[:, :, 1] * (x_MPC[:,N_MPC]-x_ss[:,1]);
    end
    
    @objective(model_MPC, Min, J_MPC);
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Start values as first guess for optimization
    #=
    Start values for all time steps equal values of initial steady state
    =#
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    for i = 1:n_x
        for k = 1:N_MPC
            if j >= N_step
                set_start_value(x_MPC[i,k],x_ss[i,2]);
            else
                set_start_value(x_MPC[i,k],x_ss[i,1]);
            end
        end
    end

    for i = 1:(n_pr + n_e)
        for k = 1:N_MPC-1
            if j >= N_step
                set_start_value(u_MPC[i,k],u_ss[i,2]);
            else
                set_start_value(u_MPC[i,k],u_ss[i,1]);
            end
        end
    end
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Fix initial values of state variables
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------   
    if j==1 #, i.e., for initial simulation step:
        for i = 1:n_x
            fix(x_MPC[i],x0[i]; force=true);
        end
    else #. i.e., for all simulation steps except the first one:
        for i = 1:n_x
            fix(x_MPC[i,1],sol_cl[i,j-1]; force=true);
            
        end
    end

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Physical Modeling Constraints
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # ODE Constraints
    for k = 1:N_MPC-1
        if j >= N_step
            p = n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E
            @constraint(model_MPC, x_MPC[:,k+1] == x_MPC[:,k] + delta_t_MPC * dynamic_DHG(x_MPC[:,k], u_MPC[:,k], d_s[:,2], p));
        else
            p = n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E
            @constraint(model_MPC, x_MPC[:,k+1] == x_MPC[:,k] + delta_t_MPC * dynamic_DHG(x_MPC[:,k], u_MPC[:,k], d_s[:,1], p));
        end
        end

    return model_MPC, x_MPC, u_MPC;
end

function getMPCModel_2(j, p)

    n_x, n_pr, n_e, N_MPC, x_min, x_max, u_min, u_max, x_ss, u_ss, d_s, Q, R, P, α, x0, sol_cl, N_step = p

    global model_MPC = Model(Ipopt.Optimizer);     # define abstract opptimization model_MPC
    set_silent(model_MPC)
    
    @variables(model_MPC, begin
    x_MPC[1:n_x,1:N_MPC];
    u_MPC[1:(n_pr + n_e),1:N_MPC-1];
    end);

    for k = 1:N_MPC-1
        for i = 1:n_x
            set_lower_bound(x_MPC[i,k], x_min[i]);
            set_upper_bound(x_MPC[i,k], x_max[i]);
        end
        
        for i = 1:(n_pr + n_e)
            set_lower_bound(u_MPC[i,k], u_min[i]);
            set_upper_bound(u_MPC[i,k], u_max[i]);
        end
    end

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Terminal state constraints
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------   
    if j+N_MPC >= N_step
        @constraint(model_MPC, (x_MPC[:,N_MPC]-x_ss[:,2])' * P[:, :, 2] * (x_MPC[:,N_MPC]-x_ss[:,2]) <= α[2]);
    else
        @constraint(model_MPC, (x_MPC[:,N_MPC]-x_ss[:,1])' * P[:, :, 1] * (x_MPC[:,N_MPC]-x_ss[:,1]) <= α[1]);
    end
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Cost function
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    global J_MPC = 0;
    # integral over stage costs
    for k=1:N_MPC-1
        if j+k >= N_step
            global J_MPC = J_MPC + (x_MPC[:,k]-x_ss[:,2])' * Q * (x_MPC[:,k]-x_ss[:,2]) + (u_MPC[:,k]-u_ss[:,2])' * R * (u_MPC[:,k]-u_ss[:,2]);
        else
            global J_MPC = J_MPC + (x_MPC[:,k]-x_ss[:,1])' * Q * (x_MPC[:,k]-x_ss[:,1]) + (u_MPC[:,k]-u_ss[:,1])' * R * (u_MPC[:,k]-u_ss[:,1])
        end
    end
    # terminal costs
    if j+N_MPC >= N_step
        global J_MPC = J_MPC + (x_MPC[:,N_MPC]-x_ss[:,2])' * P[:, :, 2] * (x_MPC[:,N_MPC]-x_ss[:,2]);
    else
        global J_MPC = J_MPC + (x_MPC[:,N_MPC]-x_ss[:,1])' * P[:, :, 1] * (x_MPC[:,N_MPC]-x_ss[:,1]);
    end
    
    @objective(model_MPC, Min, J_MPC);
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Start values as first guess for optimization
    #=
    Start values for all time steps equal values of initial steady state
    =#
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    for i = 1:n_x
        for k = 1:N_MPC
            if j+k >= N_step
                set_start_value(x_MPC[i,k],x_ss[i,2]);
            else
                set_start_value(x_MPC[i,k],x_ss[i,1]);
            end
        end
    end

    for i = 1:(n_pr + n_e)
        for k = 1:N_MPC-1
            if j+k >= N_step
                set_start_value(u_MPC[i,k],u_ss[i,2]);
            else
                set_start_value(u_MPC[i,k],u_ss[i,1]);
            end
        end
    end
    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Fix initial values of state variables
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------   
    if j==1 #, i.e., for initial simulation step:
        for i = 1:n_x
            fix(x_MPC[i],x0[i]; force=true);
        end
    else #. i.e., for all simulation steps except the first one:
        for i = 1:n_x
            fix(x_MPC[i,1],sol_cl[i,j-1]; force=true);
            
        end
    end

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    ## Physical Modeling Constraints
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # ODE Constraints
    for k = 1:N_MPC-1
        if j+k >= N_step
            p = n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E
            @constraint(model_MPC, x_MPC[:,k+1] == x_MPC[:,k] + delta_t_MPC * dynamic_DHG(x_MPC[:,k], u_MPC[:,k], d_s[:,2], p))
        else
            p = n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E
            @constraint(model_MPC, x_MPC[:,k+1] == x_MPC[:,k] + delta_t_MPC * dynamic_DHG(x_MPC[:,k], u_MPC[:,k], d_s[:,1], p));
        end
        end

    return model_MPC, x_MPC, u_MPC;
end