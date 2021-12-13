function newton(A::Matrix, G:: Matrix, y::Array, g::Array, z::Array)
    k = 0
    t = 1
    γ = 0.5
    η = 0.1
    f(z) = sum((A*z - y).^2)
    grad(z) = G*z - g
    hes = G
    while norm(grad(z)) > 1e-4 && k < 10000
        d = -inv(hes)*grad(z)
        while f(z + t*d) > f(z) + η*t*grad(z)'*d
            t = γ*t
        end
        z = z + t*d
        k = k + 1
    end
    return z
end
