

function field_sheet_iso(V)

    x = LinRange(-1, 1, 40)
    X, Y, Z = mgrid(x, x, x)

    psi = V .* Z + X.^2 - Y.^2

    trace = isosurface(
        x=X[:], y=Y[:], z=Z[:], value=psi[:], 
        isomin=0, isomax=1,
        surface=(count=6, fill=0.7),
        caps=attr(x_show=false, y_show=false, z_show=false),
        opacity=0.6,
        showscale=false
        )

    return plot(trace)
end