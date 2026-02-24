function tau(y, a, S)

    return sqrt((1 + y*a)^2 - S*y^2)
end

function zeta(y, h, Q)

    return sqrt((1 - y*h)^2 + Q*y^2)
end


function V(y, R, a, h, S, Q)

    tu = tau(1.0, a, S)
    zu = zeta(1.0, h, Q)

    E_H = (R / tu) / sqrt(1 - (1-1/(tu*zu)^2)*R^2)
    return tau(y, a, S) * E_H / zeta(y, h, Q)
end


function psi_contours(E, R, a, h, S, Q)

    N = 100
    ys = LinRange(0,1,N+1)
    dy = ys[2]
    psi_vec_1 = zeros(N+1)
    psi_vec_2 = zeros(N+1)
    psi = zeros(N, 2*N+1)

    for i=1:N
        y = ys[i+1]
        V1 = V(y, R, a, h, S, Q)
        V2 = V(-y, R, a, h, S, Q)
        psi_vec_1[i+1] = psi_vec_1[i] + dy*(E/V1)
        psi_vec_2[i+1] = psi_vec_2[i] + dy*(E/V2)
    end
    psi[1,1:N] = reverse(psi_vec_2[2:end])
    psi[1,N+2:end] = psi_vec_1[2:end]

    for i=2:N
        psi[i,:] = psi[i-1,:] .+ E*dy
    end

    return psi
end

