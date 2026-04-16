


function outflow_speed(mu, VA, a2, S, K, Q)

    C = (1 - mu^2 - K/3 - (Q/4)*(1 - mu^4))*(1 + 2*a2/3 - S/6)
    if VA == 1
        return 1.0
    else
        U2 = (1 - mu^2)*C^2 * VA^2 / (1 - VA^2)
        return sqrt(U2/(1 + U2))
    end

end


function rate(mu, VA, a2, S, K, Q)

    C = (1 - mu^2 - K/3 - (Q/4)*(1 - mu^4))*(1 + 2*a2/3 - S/6)
    return mu*C*outflow_speed(mu, VA, a2, S, K, Q)
end


function upstream_drift(mu, VA)

    C = (1 - mu^2)/(1 + mu^2)
    return mu*C*outflow_speed(mu, VA, 0.0, 0.0, 0.0, 0.0)
end


function static_lambda(r, L, N)

    a2 = N[1]^2 / (4*(r - 1))

    LN12 = L[1]*N[2] - L[2]*N[1]
    LN13 = L[1]*N[3] - L[3]*N[1]
    LN23 = L[2]*N[3] - L[3]*N[2]

    S = -(0.5 - (3/2)*N[1]^2)
    K = LN23 - 0.5*LN12 - 0.5*LN13

    return a2 - 1.25*S - 0.75*K
end


function orbit_basis(r)
    
    ξ = r >= 3.0 ? acosh(sqrt((1 - 1/r)/(1 - (3/2)/r))) : acosh(2*sqrt(2)/(3*sqrt(1 - 1/r)))
    ω = r >= 3.0 ? 0.0 : atan((r/(3*sqrt(3)))*sqrt((3/r-1)^3/(1 - 1/r)))

    T = [cosh(ξ), -sinh(ξ)*sin(ω), sinh(ξ)*cos(ω), 0.0]
    epar = [sinh(ξ), -cosh(ξ)*sin(ω), cosh(ξ)*cos(ω), 0.0]
    eperp = [0.0, cosh(ω), sin(ω), 0.0]
    e3 = [0.0, 0.0, 0.0, 1.0]

    return T, epar, eperp, e3
end


function orbit_basis_cartesian(r, θ)

    T, epar, eperp, e3 = orbit_basis(r)
    xhat = cos(θ)*eperp[2:3] - sin(θ)*epar[2:3]
    yhat = sin(θ)*eperp[2:3] + cos(θ)*epar[2:3]

    return xhat, yhat
end



function orbit_lambda(r, L, N)

    ξ = r >= 3.0 ? acosh(sqrt((1 - 1/r)/(1 - (3/2)/r))) : acosh(2*sqrt(2)/(3*sqrt(1 - 1/r)))
    ω = r >= 3.0 ? 0.0 : atan((r/(3*sqrt(3)))*sqrt((3/r-1)^3/(1 - 1/r)))

    T, epar, eperp, e3 = orbit_basis(r)

    L_4vec = L[1]*epar + L[2]*eperp + L[3]*e3
    N_4vec = N[1]*epar + N[2]*eperp + N[3]*e3

    TN2 = (kron(T', N_4vec) - transpose(kron(T', N_4vec))) .^2
    LN2 = (kron(L_4vec', N_4vec) - transpose(kron(L_4vec', N_4vec))) .^ 2

    S = TN2[1,2] - 0.5*TN2[1,3] - 0.5*TN2[1,4] + 0.5*TN2[2,3] + 0.5*TN2[2,4] - TN2[3,4]
    K = -LN2[1,2] + 0.5*LN2[1,3] + 0.5*LN2[1,4] - 0.5*LN2[2,3] - 0.5*LN2[2,4] + LN2[3,4]

    return -1.25*S - 0.75*K
end


function drop_lambda(L, N)

    LN12 = L[1]*N[2] - L[2]*N[1]
    LN13 = L[1]*N[3] - L[3]*N[1]
    LN23 = L[2]*N[3] - L[3]*N[2]

    S = -(0.5 - (3/2)*N[1]^2)
    K = LN23 - 0.5*LN12 - 0.5*LN13

    return -1.25*S - 0.75*K
end


function lambda_0(α, β)

    λ = 1 - 2*cos(α)^2 - sin(α)^2 * cos(β)^2

    if (0.0 < α < π)
        return λ
    else    
        return -λ - 2
    end
end


function avg_rate(VA, a2, S, K, Q)

    N = 100
    μs = LinRange(0, 1, N)
    Ravg = 0.0
    for i=1:N
        Ravg += rate(μs[i], VA, a2, S, K, Q)
    end

    return Ravg / N
end


function f_actual(VA, α, β)

    Mr = cos(α)
    Nr = sin(α)*cos(β)

    κ = 0.1
    S = 0.5*(3*Nr^2 - 1) * κ
    K = 0.5*(3*Mr^2 - 1) * κ
    Q = S + K

    R = avg_rate(VA, 0.0, S, K, Q)
    Rsr = avg_rate(VA, 0.0, 0.0, 0.0, 0.0)

    return (R/Rsr - 1) / κ
end


function f_error()

    l0(α, β) = 1.2444444 - 2.3166666*cos(α)^2 - 1.41666666*sin(α)^2 * cos(β)^2
    l1(α, β) = 0.75 - 1.5*cos(α)^2 - 0.75*sin(α)^2 * cos(β)^2
    f(α, β) = 1 - 2*cos(α)^2 - sin(α)^2 * cos(β)^2

    x = LinRange(0, π/2, 100)
    y = LinRange(0, π/2, 100)
    err = []

    for i=1:100
        for j=1:100
            err1 = abs((f(x[i], y[j]) - l0(x[i], y[j])))
            err2 = abs((f(x[i], y[j]) - l1(x[i], y[j])))

            if err1 !== NaN
                push!(err, err1)
            end
            if err2 !== NaN
                push!(err, err2)
            end
        end
    end

    return maximum(err)
end

