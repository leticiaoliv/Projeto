function gradiente(f::Function, grad::Array, z::Array)
    k = 0
    t = 1
    γ = 0.5
    η = 0.1
    while abs(grad(z)) > 0
        d = -grad(z)
        while f(z + t*d) > f(z) + η*t*grad(z)*d
            t = γ*t
        end
        z = z + t*d
        k = k + 1
    end
    return z
end