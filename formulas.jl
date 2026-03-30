


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



function lambda(r, chi_a, chi_v, N, M)

    L = Vector([N[2]*M[3]-N[3]*M[2], N[3]*M[1] - N[1]*M[3], N[1]*M[2] - N[2]*M[1]])

    v2 = chi_v/(2*r)
    g2v2 = v2/(1 - v2)
    a2 = chi_a / r

    LN12 = L[1]*N[2] - L[2]*N[1]
    LN13 = L[1]*N[3] - L[3]*N[1]
    LN23 = L[2]*N[3] - L[3]*N[2]

    A = 1 + (3/2)*g2v2
    B = 0.5*(1 + 3*g2v2)

    S = A*N[1] - B*N[2] - 0.5*N[3]
    K = -0.5*LN12 -B*LN13 + A*LN23

    return a2 - 1.25*K - 0.75*S
end