function erro(x::Array, y::Array, α::Real)
    n = length(x)
    erro = 0
    if α = 1
        grad = zeros(2)
        for i = 1:n
            erro = erro + (a₀ + a₁*x[i] - y[i])^2
            grad[1] = grad[1] + 2*(a₀ + a₁*x[i] - y[i]) 
            grad[2] = grad[2] + 2*(a₀ + a₁*x[i] - y[i])*x[i]
        end
    elseif α = 2
        grad = zeros(3)
        for i = 1:n
            erro = erro + (a₀ + a₁*x[i] + a₂*x[i]^2 - y[i])^2
            grad[1] = grad[1] + 2*(a₀ + a₁*x[i] + a₂*x[i]^2 - y[i]) 
            grad[2] = grad[2] + 2*(a₀ + a₁*x[i] + a₂*x[i]^2 - y[i])*x[i]
            grad[3] = grad[3] + 2*(a₀ + a₁*x[i] + a₂*x[i]^2 - y[i])*x[i]^2
        end
    elseif α = 3
        grad = zeros(2)
        for i = 1:n
            erro = erro + (a₀*x[i]^a₁ - y[i])^2
            grad[1] = grad[1] + 2*x[i]^a₁*(a₀*x[i]^a₁ - y[i])
            grad[2] = grad[2] + 2*a₀*x[i]^a₁*log(x[1])*(a₀*x[i]^a₁ - y[i])
        end
    elseif α = 4
        grad = zeros(2)
        for i = 1:n
            erro = erro + (a₀*exp(a₁*x[i]) - y[i])^2
            grad[1] = grad[1] + 2*exp(a₁*x[i])*(a₀*exp(a₁*x[i]) - y[i])
            grad[2] = grad[2] + 2*a₀*x[i]*exp(a₁*x[i])*(a₀*exp(a₁*x[i]) - y[i])
        end
    elseif α = 5
        grad = zeros(3)
        for i = 1:n
            erro = erro + (exp(a₀ + a₁*x[i] + a₃*x[i]^2) - y[i])^2
            grad[1] = grad[1] + 2*exp(a₀ + a₁*x[i] + a₃*x[i]^2)*(exp(a₀ + a₁*x[i] + a₃*x[i]^2) - y[i])
            grad[2] = grad[2] + 2*x[i]*exp(a₀ + a₁*x[i] + a₃*x[i]^2)*(exp(a₀ + a₁*x[i] + a₃*x[i]^2) - y[i])
            grad[3] = grad[3] + 2*x[i]^2*exp(a₀ + a₁*x[i] + a₃*x[i]^2)*(exp(a₀ + a₁*x[i] + a₃*x[i]^2) - y[i])
        end
    end
    return erro    
end
