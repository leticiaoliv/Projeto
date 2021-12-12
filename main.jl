using Plots
using LinearAlgebra
using DelimitedFiles

function main()
    x = readdlm("dados_x.csv");
    y = readdlm("dados_y.csv");
    A = [ones(52) x x.^2];
    G = [52 sum(x) sum(x.^2); sum(x) sum(x.^2) sum(x.^3); sum(x.^2) sum(x.^3) sum(x.^4)];
    g = [sum(y); sum(x.*y); sum(x.^2 .*y)];
    z = ones(3);

    z = gradiente(A, G, y, g, z)
end
