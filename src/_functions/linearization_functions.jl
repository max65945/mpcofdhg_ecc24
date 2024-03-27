function dynamic_DHG(x, u, d, p)

    n_x, n_tes, n_e, n_pr, n_d, m_v, m_tes_tot, m_e, κ, B_tes, B_0_plus, B_0_minus, B, E = p

    m   = x[1:n_tes]
    T_v = x[n_tes+1:n_tes+n_v]
    T_e = x[n_tes+n_v+1:n_x]
    q   = u[1:n_e]
    Q   = u[n_e+1:(n_e+n_pr)]

    Θ_inv = 10.0^-3 .* [1/m[1]; 1/m_v[2]; 1/m_v[3]; 1/(m_tes_tot[1]-m[1]); 1/m_v[5]; 1/m[2]; 1/(m_tes_tot[2]-m[2])]

    return [10.0^-3 * B_tes*q;
            diagm(Θ_inv) * (-diagm(B_0_plus*q.+κ) * T_v + B_0_plus * diagm(q) * T_e .+ E[1:n_v,n_d+1] * d[n_d+1]);
            inv(1000.0 .* diagm(m_e)) * (diagm(q)*B_0_minus' * T_v - diagm(q.+κ) * T_e .+ E[n_v+1:n_v+n_e,n_d+1] * d[n_d+1] + B[n_v+1:n_v+n_e,:] * Q + E[n_v+1:n_v+n_e,1:n_d] * d[1:n_d])]

end

function checkStabilizability(A, B)
    p = eigvals(A);
    p = round.(p; digits = 5);
    for i = eachindex(p)
        if real(p[i]) >= 0 && rank([p[i]*I-A B])==n_x
            display(("Eigenvalue", p[i], "is stabilizable!"))
        elseif real(p[i]) >= 0 && rank([p[i]*I-A B])!=n_x
            error(("Eigenvalue", p[i], "is NOT stabilizable!"))
        end
    end
end

function checkKappa(AK, kappa)
    if kappa + maximum(real(eigvals(AK))) <= 0
        display(("feasible kappa, κ = ", kappa))
    else
        error("infeasible(!) kappa")
    end
end

function getP_lyap(AK, K, R, Q, κ_MPC)
    
    Ak_lyap = AK + κ_MPC * 1. * I(length(AK[:,1]));
    Q_star = Q + transpose(K) *R*K; 
    P = round.(lyap(Ak_lyap',Q_star); digits = 3)
    if isposdef(P) && issymmetric(P)
        display("feasible P, i.e., P>0")
    else
        error("infeasible(!) P, i.e., P<0")
    end

    return P;
end