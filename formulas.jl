


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