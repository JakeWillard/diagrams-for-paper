



function plot_rate(VA)

    μs = LinRange(0, 1, 100)

    C_pos = [0.0, 0.1, 0.2, 0.3, 0.4]
    C_all = [-0.4, -0.2, 0.0, 0.2, 0.4]

    pa = plot(μs, [rate(μ, VA, 0.0, 0.0, 0.0, 0.0) for μ in μs], label=L"\bar{a}^2=%$(C_pos[1])")
    for i=2:5
        plot!(pa, μs, [rate(μ, VA, C_pos[i], 0.0, 0.0, 0.0) for μ in μs], label=L"\bar{a}^2=%$(C_pos[i])")
    end
    
    pS = plot(μs, [rate(μ, VA, 0.0, 0.0, 0.0, 0.0) for μ in μs], label=L"\bar{S}=%$(C_all[1])")
    pK = plot(μs, [rate(μ, VA, 0.0, 0.0, 0.0, 0.0) for μ in μs], label=L"\mathcal{\bar{K}}=%$(C_all[1])")
    pQ = plot(μs, [rate(μ, VA, 0.0, 0.0, 0.0, 0.0) for μ in μs], label=L"\bar{Q}=%$(C_all[1])")
    for i=2:5
        plot!(pS, μs, [rate(μ, VA, 0.0, C_all[i], 0.0, 0.0) for μ in μs], label=L"\bar{S}=%$(C_all[i])")
        plot!(pK, μs, [rate(μ, VA, 0.0, 0.0, C_all[i], 0.0) for μ in μs], label=L"\mathcal{\bar{K}}=%$(C_all[i])")
        plot!(pQ, μs, [rate(μ, VA, 0.0, 0.0, 0.0, C_all[i]) for μ in μs], label=L"\bar{Q}=%$(C_all[i])")
    end

    for p in [pa, pS, pK, pQ]

        xlabel!(p, L"\mu")
        ylabel!(p, L"\mathcal{R}")
        ylims!(p, 0.0, 0.25)
    end

    title!(pa, L"\bar{S}=\mathcal{\bar{K}}=\bar{Q}=0")
    title!(pS, L"\bar{a}^2=\mathcal{\bar{K}}=\bar{Q}=0")
    title!(pK, L"\bar{a}^2=\bar{S}=\bar{Q}=0")
    title!(pQ, L"\bar{a}^2=\bar{S}=\mathcal{\bar{K}}=0")


    l = @layout [a b; c d]
    return plot(pa, pS, pK, pQ, layout=l)

end


function plot_difference()

    μs = LinRange(0, 1, 100)
    VAs = [0.2, 0.4, 0.6, 0.8, 1.0]
    colors = [:blue, :red, :green, :black, :purple]

    p = plot(μs, [rate(μ, VAs[1], 0.0, 0.0, 0.0, 0.0) for μ in μs], label=L"V_A=0.2", color=colors[1])
    for i=2:5
        plot!(p, μs, [rate(μ, VAs[i], 0.0, 0.0, 0.0, 0.0) for μ in μs], label=L"V_A=%$(VAs[i])", color=colors[i])
    end
    for i=1:5
        #plot!(p, μs, [rate(μ, VAs[i], 0.0, 0.0, 0.0, 0.0)/(1 + μ^2) for μ in μs], linestyle=:dash, label="", color=colors[i])
        plot!(p, μs, [upstream_drift(μ, VAs[i]) for μ in μs], linestyle=:dash, label="", color=colors[i])
    end

    xlabel!(p, L"\mu")
    title!(p, L"$\mathcal{R}$ (solid) and $V_d|_{\mathfrak{u}}$ (dashed)")



    return p
end


function plot_lambda_heatmap(r, N, M)

    chi_as = LinRange(0, 1, 100)
    chi_vs = LinRange(0, 1, 100)
    ls = zeros(100, 100)
    for i=1:100
        for j=1:100
            ls[j, i] = lambda(r, chi_as[i], chi_vs[j], N, M)
        end
    end

    p = heatmap(chi_as, chi_vs, ls)
    xlabel!(p, L"\chi_a")
    ylabel!(p, L"\chi_v")

    return p
end


function light_rate_diagram(Δ̄s1, Δ̄s2, xlimits)


    p = vline([Δ̄s1], label=L"s=w_y^+", color=:darkred, aspect_ratio=:equal)
    vline!(p, [-Δ̄s2], label=L"s=w_y^-", color=:darkgreen, aspect_ratio=:equal)
    x1 = [0.0, Δ̄s1, 0.0]
    y1 = [0.0, Δ̄s1, 2*Δ̄s1]
    x2 = [0.0, -Δ̄s2, 0.0]
    y2 = [0.0, Δ̄s2, 2*Δ̄s2]

    plot!(p, x1, y1, label="light signal", color=:black, linestyle=:dash, aspect_ratio=:equal)
    plot!(p, x2, y2, label="", color=:black, linestyle=:dash, aspect_ratio=:equal)

    xlabel!(p, L"\bar{s}")
    ylabel!(p, L"t")
    xlims!(p, xlimits...)

    return p
end


function plot_lambda()

    ϵ = 0.3
    αs = LinRange(0-ϵ, π+ϵ, 100)
    βs = LinRange(0-ϵ, π+ϵ, 100)
    λ1s = zeros(100, 100)
    λ2s = zeros(100, 100)
    fs = zeros(100, 100)

    for i=1:100
        for j=1:100
            λ1s[j,i] = 1.01111 - 1.61666*cos(αs[i])^2 - 1.41666*sin(αs[i])^2 * cos(βs[j])^2
            λ2s[j,i] = 0.58333 - cos(αs[i])^2 - 0.75*sin(αs[i])^2 * cos(βs[j])^2
            fs[j,i] = 1 - 2*cos(αs[i])^2 - sin(αs[i])^2 * cos(βs[j])^2
        end
    end


    #levels = LinRange(-1.0, 1.0, 10)
    ls1 = [round(x, sigdigits=2) for x in LinRange(minimum(λ1s), maximum(λ2s), 10)]
    ls2 = [round(x, sigdigits=2) for x in LinRange(minimum(λ2s), maximum(λ2s), 10)]

    println(ls1)
    println(ls2)
    p1 = heatmap(αs, βs, λ1s, clim=(-0.75, 1.0))
    p2 = heatmap(αs, βs, λ2s, clim=(-0.75, 1.0))
    p3 = contour(αs, βs, λ1s, levels=[0.0])
    p4 = contour(αs, βs, λ2s, levels=[0.0])
    # p1 = contour(αs, βs, λ1s, aspect_ratio = :equal, clabels=true, color=:darkrainbow, levels=ls1)
    # p2 = contour(αs, βs, λ2s, aspect_ratio = :equal, clabels=true, color=:darkrainbow, levels=ls2)
    #p3 = contour(αs, βs, fs, aspect_ratio = :equal, clabels=true, color=:darkrainbow, levels=5)
    # p1 = contour(αs, βs, λ1s, clabels=true, cbar=false, color=:darkrainbow, levels=[-1.0, -0.5, 0.0, 0.5])
    # p2 = contour(αs, βs, λ1s, clabels=true, cbar=false, color=:darkrainbow, levels=[-1.0, -0.5, 0.0, 0.5])
    #scatter!(p, [π/2], [π/2], legend=false)
    #annotate!(p, π/2, π/2, ("  1.0", 8, :black, :left))
    title!(p1, L"\lambda / \bar{\kappa}" * " for " * L"V_A \ll 1")
    xlabel!(p1, L"\alpha")
    ylabel!(p1, L"\beta")

    title!(p2, L"\lambda / \bar{\kappa}" * " as " * L"V_A \rightarrow 1")
    xlabel!(p2, L"\alpha")
    ylabel!(p2, L"\beta")

    # title!(p3, L"1 - 2\cos^2(\alpha) - \sin^2(\alpha) \cos^2(\beta)")
    # xlabel!(p3, L"\alpha")
    # ylabel!(p3, L"\beta")

    #l = @layout [a; b; c d]
    return plot(p1, p2)
end


# function plot_static_lambda(r1, r2)

#     N = 5
#     rs = LinRange(r1, r2, 100)

#     L1s = [[cos(θ), -sin(θ), 0.0] for θ in LinRange(0, π/2, N)]
#     N1s = [[sin(θ), cos(θ), 0.0] for θ in LinRange(0, π/2, N)]
#     L2s = [[cos(θ), 0.0, -sin(θ)] for θ in LinRange(0, π/2, N)]
#     N2s = [[sin(θ), 0.0, cos(θ)] for θ in LinRange(0, π/2, N)]
#     L3s = [[0.0, cos(θ), -sin(θ)] for θ in LinRange(0, π/2, N)]
#     N3s = [[0.0, sin(θ), cos(θ)] for θ in LinRange(0, π/2, N)]

#     N11 = [N1s[i][1] for i=1:N]
#     N21 = [N2s[i][1] for i=1:N]
#     N32 = [N3s[i][2] for i=1:N]

#     λ1 = zeros(100, N)
#     λ2 = zeros(100, N)
#     λ3 = zeros(100, N)

#     for i=1:N
#         λ1[:, i] = [static_lambda(r, L1s[i], N1s[i]) for r in rs]
#         λ2[:, i] = [static_lambda(r, L2s[i], N2s[i]) for r in rs]
#         λ3[:, i] = [static_lambda(r, L3s[i], N3s[i]) for r in rs]
#     end

#     p1 = plot(rs, λ1, label=reshape([L"N_1 = %$(round(N11[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_3")
#     p2 = plot(rs, λ2, label=reshape([L"N_1 = %$(round(N21[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_2")
#     #p3 = plot(rs, λ3, label=reshape([L"N_2 = %$(round(N32[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_1")

#     return plot(p1, p2)
# end


# function plot_orbit_lambda(r1, r2)

#     N = 5
#     rs = LinRange(r1, r2, 100)

#     L1s = [[cos(θ), -sin(θ), 0.0] for θ in LinRange(0, π/2, N)]
#     N1s = [[sin(θ), cos(θ), 0.0] for θ in LinRange(0, π/2, N)]
#     L2s = [[cos(θ), 0.0, -sin(θ)] for θ in LinRange(0, π/2, N)]
#     N2s = [[sin(θ), 0.0, cos(θ)] for θ in LinRange(0, π/2, N)]
#     L3s = [[0.0, cos(θ), -sin(θ)] for θ in LinRange(0, π/2, N)]
#     N3s = [[0.0, sin(θ), cos(θ)] for θ in LinRange(0, π/2, N)]

#     N1par = [N1s[i][1] for i=1:N]
#     N2par = [N2s[i][1] for i=1:N]
#     N3perp = [N3s[i][2] for i=1:N]

#     λ1 = zeros(100, N)
#     λ2 = zeros(100, N)
#     λ3 = zeros(100, N)

#     for i=1:N
#         λ1[:, i] = [orbit_lambda(r, L1s[i], N1s[i]) for r in rs]
#         λ2[:, i] = [orbit_lambda(r, L2s[i], N2s[i]) for r in rs]
#         λ3[:, i] = [orbit_lambda(r, L3s[i], N3s[i]) for r in rs]
#     end

#     p1 = plot(rs, λ1, label=reshape([L"N_{\parallel} = %$(round(N1par[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_3")
#     p2 = plot(rs, λ2, label=reshape([L"N_{\parallel} = %$(round(N2par[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_{\perp}")
#     p3 = plot(rs, λ3, label=reshape([L"N_{\perp} = %$(round(N3perp[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_{\parallel}")

#     return plot(p1, p2, p3)
# end


# function plot_drop_lambda(r1, r2)

#     N = 5
#     rs = LinRange(r1, r2, 100)

#     L1s = [[cos(θ), -sin(θ), 0.0] for θ in LinRange(0, π/2, N)]
#     N1s = [[sin(θ), cos(θ), 0.0] for θ in LinRange(0, π/2, N)]
#     L2s = [[cos(θ), 0.0, -sin(θ)] for θ in LinRange(0, π/2, N)]
#     N2s = [[sin(θ), 0.0, cos(θ)] for θ in LinRange(0, π/2, N)]
#     L3s = [[0.0, cos(θ), -sin(θ)] for θ in LinRange(0, π/2, N)]
#     N3s = [[0.0, sin(θ), cos(θ)] for θ in LinRange(0, π/2, N)]

#     N11 = [N1s[i][1] for i=1:N]
#     N21 = [N2s[i][1] for i=1:N]
#     N32 = [N3s[i][2] for i=1:N]

#     λ1 = zeros(100, N)
#     λ2 = zeros(100, N)
#     λ3 = zeros(100, N)

#     for i=1:N
#         λ1[:, i] = [drop_lambda(r, L1s[i], N1s[i]) for r in rs]
#         λ2[:, i] = [drop_lambda(r, L2s[i], N2s[i]) for r in rs]
#         λ3[:, i] = [drop_lambda(r, L3s[i], N3s[i]) for r in rs]
#     end

#     p1 = plot(rs, λ1, label=reshape([L"N_1 = %$(round(N11[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_3")
#     p2 = plot(rs, λ2, label=reshape([L"N_1 = %$(round(N21[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_2")
#     p3 = plot(rs, λ3, label=reshape([L"N_2 = %$(round(N32[i], digits=2))" for i=1:N], (1,N)), title=L"\vec{M}=\vec{e}_1")

#     return plot(p1, p2, p3)

# end


# function plot_lambda_0()

#     ϵ = 0.3
#     αs = LinRange(0-ϵ, π+ϵ, 100)
#     βs = LinRange(0-ϵ, π+ϵ, 100)
#     λs = zeros(100, 100)

#     for i=1:100
#         for j=1:100
#             λs[j,i] = lambda_0(αs[i], βs[j])
#         end
#     end

#     p = contour(αs, βs, lambda_0, clabels=true, cbar=false, color=:darkrainbow, levels=[-1.0, -0.5, 0.0, 0.5])
#     scatter!(p, [π/2], [π/2], legend=false)
#     annotate!(p, π/2, π/2, ("  1.0", 8, :black, :left))
#     title!(p, L"f(\alpha, \beta)")
#     xlabel!(p, L"\alpha")
#     ylabel!(p, L"\beta")


#     return p
# end



# function plot_orbit_vectors()

#     N = 30
#     x = LinRange(-5.0, 5.0, N)
#     y = LinRange(-5.0, 5.0, N)
#     u = []
#     v = []
#     X = zeros(N, N)
#     Y = zeros(N, N)

#     for i=1:N
#         for j=1:N
#             r = sqrt(x[i]^2 + y[j]^2)
#             θ = atan(y[j]/x[i])
#             X[j,i] = x[i]
#             Y[j,i] = y[j]
#             if r <= 1.0
#                 push!(u, 0.0)
#                 push!(v, 0.0)
#             else
#                 T, epar, eperp, e3 = orbit_basis(r)
#                 epar_x = cos(θ)*epar[2] - sin(θ)*epar[3]
#                 epar_y = sin(θ)*epar[2] + cos(θ)*epar[3]
#                 push!(u, epar_y)
#                 push!(u, epar_x)
#             end
#         end
#     end

#     return quiver(X, Y, quiver=(u, v))
# end



# function plot_f_actual_heatmaps(VA)

#     αs = LinRange(0, π, 100)
#     βs = LinRange(0, π, 100)

    
#     return heatmap(αs, βs, (x,y) -> f_actual(VA, x, y))
# end


# function plot_f_actual_comparisons()

#     VAs = [0.01, 0.5, 1.0]
#     βs = [0.0, π/6, π/3, π/2]
#     αs = LinRange(0, π, 100)

#     p1 = plot(αs, x -> f_actual(VAs[1], x, βs[1]), label=L"V_A=%$(VAs[1])")
#     p2 = plot(αs, x -> f_actual(VAs[1], x, βs[2]), label=L"V_A=%$(VAs[1])")
#     p3 = plot(αs, x -> f_actual(VAs[1], x, βs[3]), label=L"V_A=%$(VAs[1])")
#     p4 = plot(αs, x -> f_actual(VAs[1], x, βs[4]), label=L"V_A=%$(VAs[1])")

#     for VA in VAs[2:end]

#         plot!(p1, αs, x -> f_actual(VA, x, βs[1]), label=L"V_A=%$(VA)")
#         plot!(p2, αs, x -> f_actual(VA, x, βs[2]), label=L"V_A=%$(VA)")
#         plot!(p3, αs, x -> f_actual(VA, x, βs[3]), label=L"V_A=%$(VA)")
#         plot!(p4, αs, x -> f_actual(VA, x, βs[4]), label=L"V_A=%$(VA)")
#     end

#     title!(p1, L"\beta=0")
#     title!(p2, L"\beta=\pi/6")
#     title!(p3, L"\beta=\pi/3")
#     title!(p4, L"\beta=\pi/2")

#     xlabel!(p1, L"\alpha")
#     xlabel!(p2, L"\alpha")
#     xlabel!(p3, L"\alpha")
#     xlabel!(p4, L"\alpha")

#     ylabel!(p1, L"\lambda/\bar{\kappa}")
#     ylabel!(p2, L"\lambda/\bar{\kappa}")
#     ylabel!(p3, L"\lambda/\bar{\kappa}")
#     ylabel!(p4, L"\lambda/\bar{\kappa}")

#     l = @layout [a b; c d]
#     return plot(p1, p2, p3, p4, layout=l)
# end

# function plot_f_actual_shaded()

#     βs = [0.0, π/6, π/3, π/2]
#     αs = LinRange(0, π, 100)

#     l1 = [f_actual(0.0001, α, βs[1]) for α in αs]
#     l2 = [f_actual(0.0001, α, βs[2]) for α in αs]
#     l3 = [f_actual(0.0001, α, βs[3]) for α in αs]
#     l4 = [f_actual(0.0001, α, βs[4]) for α in αs]

#     u1 = [f_actual(1.0, α, βs[1]) for α in αs]
#     u2 = [f_actual(1.0, α, βs[2]) for α in αs]
#     u3 = [f_actual(1.0, α, βs[3]) for α in αs]
#     u4 = [f_actual(1.0, α, βs[4]) for α in αs]

#     p1 = plot(αs, l1, fillrange=u1, fillalpha = 0.35, label=L"\lambda/\bar{\kappa}")
#     p2 = plot(αs, l2, fillrange=u2, fillalpha = 0.35, label=L"\lambda/\bar{\kappa}")
#     p3 = plot(αs, l3, fillrange=u3, fillalpha = 0.35, label=L"\lambda/\bar{\kappa}")
#     p4 = plot(αs, l4, fillrange=u4, fillalpha = 0.35, label=L"\lambda/\bar{\kappa}")

#     plot!(p1, αs, x -> lambda_0(x, βs[1]), label=L"f(\alpha,\beta)")
#     plot!(p2, αs, x -> lambda_0(x, βs[2]), label=L"f(\alpha,\beta)")
#     plot!(p3, αs, x -> lambda_0(x, βs[3]), label=L"f(\alpha,\beta)")
#     plot!(p4, αs, x -> lambda_0(x, βs[4]), label=L"f(\alpha,\beta)")

#     title!(p1, L"\beta=0")
#     title!(p2, L"\beta=\pi/6")
#     title!(p3, L"\beta=\pi/3")
#     title!(p4, L"\beta=\pi/2")

#     xlabel!(p1, L"\alpha")
#     xlabel!(p2, L"\alpha")
#     xlabel!(p3, L"\alpha")
#     xlabel!(p4, L"\alpha")

#     l = @layout [a b; c d]
#     return plot(p1, p2, p3, p4, layout=l)
# end



# function conformal_alt(a, S)

#     sbar(s) = s - 0.5*a*s^2 + (1/6)*(S - a^2)*s^3

#     p = vline([sbar(s) for s in LinRange(-1, 1, 5)], color=:black, linestyle=:dash, legend=false)

#     xlims!(p, -1.0, 1.0)
#     ylims!(p, -1.0, 1.0)

#     # n = 4
#     # pos_s = [y + (1 - Cp)*y^2 for y in LinRange(0, 1, n+1)[2:end]]
#     # neg_s = [-(y + (1 - Cm)*y^2) for y in LinRange(0, 1, n+1)[2:end]]

#     # p = vline(pos_s, color=:black, linestyle=:dash, legend=false)
#     # #vline!(neg_s, color=:black, linestyle=:dash, legend=false)

#     return p
# end




# struct ConformalDiagramBase
#     Cp :: Float64
#     Cm :: Float64
# end

# function conformal_diagram(c::ConformalDiagramBase)

#     x = [-1.0, 0.0, 1.0]

#     p = plot(x, x, aspect_ratio = :equal, linestyle=:dash, color=:black, legend=false, framestyle=:origin)
#     plot!(p, x, -x, aspect_ratio = :equal, linestyle=:dash, color=:black)
#     for x0 in LinRange(0.0, 2.0, 5)[2:end]
#         plot!(p, x, x .+ x0, aspect_ratio = :equal, linestyle=:dash, color=:black)
#         plot!(p, x, -x .+ x0, aspect_ratio = :equal, linestyle=:dash, color=:black)
#     end

#     xlabel!(p, L"\bar{s}")
#     ylabel!(p, L"t")

#     xlims!(-1.0, 1.0)
#     ylims!(0.0, 2.0)
#     return p
# end
