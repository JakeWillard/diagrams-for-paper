

function ybar(y, a, S)

    return y - a*(1 + a*y/2)*y^2 + S*y^3/6
end


function rate(dl, a, K, S, VA)

    Bm_B = (1 - dl^2 - K/3)*(1 - (a/2)*(1 - a/3) - S/6)

    if VA == 1.0
        Vo = 1.0
    else
        Uo2 = (1 - dl^2)*Bm_B^2 * VA^2/(1 - VA^2)
        Vo = sqrt(Uo2 / (1 + Uo2))
    end

    return dl*Bm_B * Vo

end

function psi_for_conformal(y, t)

    return y - t
end


function plot_psi_conformal(a, S)

    N = 100

    psi = zeros(N, N)
    ybars = LinRange(0, 2, N)
    ylines = [ybar(y, a, S) for y in LinRange(0, 1, 6)]
    ts = LinRange(0, 2, N)

    for i=1:N
        for j=1:N
            psi[j, i] = ybars[i] +  ts[j]
        end
    end

    psi_lines = LinRange(0.0, ybar(1.0, a, S), 5)

    p = contour(ybars, ts, psi, levels=psi_lines, cbar=false, color=:black, xtickfontcolor=:white, ytickfontcolor=:white)
    for y in ylines
        plot!(p, y*ones(N), ts, color=:black, legend=false, linestyle=:dash)
    end
    xlims!(p, 0.0, 2.0)
    ylims!(p, 0.0, 2.0)
    xlabel!(p, L"\bar{y}")
    ylabel!(p, L"t")

    return p
end


function plot_rate(VA)

    dls = LinRange(0, 1, 100)

    p1 = plot(dls, [rate(dl, 0.0, 0.0, 0.0, VA) for dl in dls], title=L"\mathcal{\bar{K}}=\bar{S}=0", label=L"\bar{a}=0.0")
    plot!(p1, dls, [rate(dl, 0.1, 0.0, 0.0, VA) for dl in dls], label=L"\bar{a}=0.1")
    plot!(p1, dls, [rate(dl, 0.2, 0.0, 0.0, VA) for dl in dls], label=L"\bar{a}=0.2")
    plot!(p1, dls, [rate(dl, 0.3, 0.0, 0.0, VA) for dl in dls], label=L"\bar{a}=0.3")

    p2 = plot(dls, [rate(dl, 0.0, -0.5, 0.0, VA) for dl in dls], title=L"\bar{a}=\bar{S}=0", label=L"\mathcal{\bar{K}}=-0.5")
    plot!(p2, dls, [rate(dl, 0.0, -0.3, 0.0, VA) for dl in dls], label=L"\mathcal{\bar{K}}=-0.3")
    plot!(p2, dls, [rate(dl, 0.0, 0.0, 0.0, VA) for dl in dls], label=L"\mathcal{\bar{K}}=0.0")
    plot!(p2, dls, [rate(dl, 0.0, 0.3, 0.0, VA) for dl in dls], label=L"\mathcal{\bar{K}}=0.3")
    plot!(p2, dls, [rate(dl, 0.0, 0.5, 0.0, VA) for dl in dls], label=L"\mathcal{\bar{K}}=0.5")

    p3 = plot(dls, [rate(dl, 0.0, 0.0, -0.5, VA) for dl in dls], title=L"\bar{a}=\bar{\mathcal{K}}=0", label=L"\bar{S}=-0.5")
    plot!(p3, dls, [rate(dl, 0.0, 0.0, -0.3, VA) for dl in dls], label=L"\bar{S}=-0.3")
    plot!(p3, dls, [rate(dl, 0.0, 0.0, 0.0, VA) for dl in dls], label=L"\bar{S}=0.0")
    plot!(p3, dls, [rate(dl, 0.0, 0.0, 0.3, VA) for dl in dls], label=L"\bar{S}=0.3")
    plot!(p3, dls, [rate(dl, 0.0, 0.0, 0.5, VA) for dl in dls], label=L"\bar{S}=0.5")

    return plot(p1, p2, p3, layout=(3, 1))


end
