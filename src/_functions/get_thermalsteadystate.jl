using JuMP, Ipopt

function get_thermalsteadystate(x_s, u_s, d_s, A_thermal, E_thermal, p)

    n_x, n_tes, n_u, n_c, Δt, Q, R, x_max, x_min, u_max, u_min = p

    model = Model(Ipopt.Optimizer)

    @variable(model, x[1:(n_x-n_tes)])
    @variable(model, u[1:(n_u-n_c)])
    
    @objective(model, Min,  (x-x_s)' * Q * (x-x_s) + (u-u_s)' * R * (u-u_s))
    
    #@constraint(model, x .== A_thermal * x + E_thermal * sqrt(Δt) * vcat(u, d_s)) # steady state constraint
    @constraint(model, x .== x + Δt * (A_thermal * x + E_thermal * vcat(u, d_s))) # steady state constraint

    #@constraint(model, x <= 0.99*x_max[n_tes+1:end])
    #@constraint(model, x >= 1.01*x_min[n_tes+1:end])
    #@constraint(model, x[n_tes+1] >= 0.985*x_max[n_tes+1])
    #@constraint(model, x[n_tes+6] >= 0.985*x_max[n_tes+6])
    #@constraint(model, u <= 0.99*u_max[n_e+1:end])
    #@constraint(model, u >= 1.01*u_min[n_e+1:end])

    optimize!(model)

    # store values of optimized steady state
    x_ss = round.(value.(x); digits = 3)
    u_ss = round.(value.(u); digits = 3)
    
    println("x_ss = ", round.(value.(x); digits = 3))
    println("u_ss = ", round.(value.(u); digits = 3))
    println("heat losses at steady state: ", round((sum(round.(value.(u); digits = 3))+sum(d_s[1:n_d,1])) / sum(round.(value.(u); digits = 3)) * 100.0, digits=3), " %")

    return (x_ss, u_ss)
end