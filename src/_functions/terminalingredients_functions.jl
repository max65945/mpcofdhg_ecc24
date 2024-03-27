function get_alpha(P, κ_MPC, α_init, K, AK, u_ss, u_min, u_max, x_ss, x_min, x_max, d_s, p)

    n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E = p

    α = copy(α_init)
    Δx_sol = zeros(n_x)

    check = 0
    i_precision = 0
    α_precicion = 1
    accuracy = 1

    i_α1_max       = 50
    i_α1_rand_max  = 50

    i_rand = 0

    while i_precision <= 2
        for i = 1:i_α1_max

            check += 1

            if i_precision == 0
                α_precicion = 2
            elseif i_precision == 1
                α_precicion = 1.2
            elseif i_precision == 2
                α_precicion = 1.05
            elseif i_precision >= 3
                break
            end

            i_rand = 0

            for k = 1:i_α1_rand_max
                i_rand += 1
                model_alpha_1 = Model(Ipopt.Optimizer)
                set_silent(model_alpha_1)
                @variable(model_alpha_1, Δx[1:n_x])
                J = (2.0 .* rand(Float64, n_x) .- 1.0)'*I(n_x)*Δx 
                @objective(model_alpha_1, Max, J)
                @constraint(model_alpha_1, Δx' * P * Δx == α)
                optimize!(model_alpha_1)
                
                Δx_sol = round.(value.(Δx);digits=accuracy)

                if K * Δx_sol <= u_min-u_ss || K * Δx_sol >= u_max-u_ss || Δx_sol <= x_min-x_ss || Δx_sol >= x_max-x_ss
                    break
                end
            end

            if K * Δx_sol <= u_min-u_ss || K * Δx_sol >= u_max-u_ss || Δx_sol <= x_min-x_ss || Δx_sol >= x_max-x_ss
                α = α/α_precicion
            elseif K * Δx_sol >= u_min-u_ss && K * Δx_sol <= u_max-u_ss && Δx_sol >= x_min-x_ss && Δx_sol <= x_max-x_ss
                α = α_precicion * α
                i_precision += 1
                break
            end
        end
    end
    α = α/α_precicion
    println(("α_init = ", α_init))
    println(("i_rand:",i_rand))
    println(("i_precision:",i_precision))
    println(("check:",check))
    println(("α:",α))

    α_1 = copy(α)

    i_precision = 0
    check = 0

    i_α_max = 20
    i_α_rand_max = 100
    while i_precision <= 2
        for i = 1:i_α_max
            check += 1

            if i_precision == 0
                α_precicion = 4
            elseif i_precision == 1
                α_precicion = 2
            elseif i_precision == 2
                α_precicion = 1.02
            elseif i_precision >= 3
                break
            end

            i_rand = 0
            for k = 1:i_α_rand_max
                i_rand += 1
                model_alpha = Model(Ipopt.Optimizer)
                set_silent(model_alpha)
        
                @variable(model_alpha, Δx[1:n_x])
        
                J = Δx' * P * (dynamic_DHG(Δx + x_ss, (K * Δx) + u_ss, d_s, p) - AK * Δx) - κ_MPC * Δx' * P * Δx
                @objective(model_alpha,Max, J)
        
                for j =1:n_x
                    set_start_value(Δx[j], 10.0*randn(Float64))
                end
        
                @constraint(model_alpha, K*Δx <= u_max-u_ss)
                @constraint(model_alpha, K*Δx >= u_min-u_ss)
                @constraint(model_alpha, Δx <= x_max-x_ss)
                @constraint(model_alpha, Δx >= x_min-x_ss)
                @constraint(model_alpha, Δx' * P * Δx <= α)
        
                optimize!(model_alpha)

                if k == 1
                    Δx_sol = value.(Δx)
                elseif round(value.(Δx)' * P * (dynamic_DHG(value.(Δx) + x_ss, (K * value.(Δx)) + u_ss, d_s, p) - AK * value.(Δx)) - κ_MPC * Δx_sol' * P * value.(Δx); digits=accuracy) >= round(Δx_sol' * P * (dynamic_DHG(Δx_sol + x_ss, (K * Δx_sol) + u_ss, d_s, p) - AK * Δx_sol) - κ_MPC * Δx_sol' * P * Δx_sol; digits=accuracy) 
                    Δx_sol = value.(Δx)
                end
            end

            if round(Δx_sol' * P * (dynamic_DHG(Δx_sol + x_ss, (K * Δx_sol) + u_ss, d_s, p) - AK * Δx_sol) - κ_MPC * Δx_sol' * P * Δx_sol; digits=accuracy) <= 0
                α = α_precicion * α
                i_precision += 1
                break
            else
                α = α/α_precicion
            end

        end
    end
    α = α/α_precicion

    println(("i_rand:",i_rand))
    println(("i_precision:",i_precision))
    println(("check:",check))
    println(("α:",α))
    
    return α, α_1
end