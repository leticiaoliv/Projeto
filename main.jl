using Plots
using LinearAlgebra
using DelimitedFiles

function main()
    x = readdlm("dados_x.csv");
    y = readdlm("dados_y.csv");

    A = [ones(52) x];
    G = [52 sum(x); sum(x) sum(x.^2)];
    g = [sum(y); sum(x.*y)];
    z = ones(2);

    gr = gradiente(A, G, y, g, z)
    nt = newton(A, G, y, g, z)
    errograd = sqrt(sum((A*gr - y).^2))/norm(y);
    erronewt = sqrt(sum((A*nt - y).^2))/norm(y);
    print(" Aproximação Linear")
    print(" Coeficientes Gradiente = ", gr, " Coeficientes Newton = ", nt)
    print(" Erro Absoluto Gradiente = ", sqrt(sum((A*gr - y).^2)),
          " Erro Absoluto Newton =", sqrt(sum((A*nt - y).^2)))
    print(" Erro Relativo Gradiente = ", errograd, " Erro Relativo Newton = ", erronewt)


    scatter(x, y, title = "Contaminados COVID-19",label = "Casos por Semana", lw = 3)
    plot!(x -> gr[1] + gr[2]*x, label = "Gradiente", c=:red, lw=3)
    plot!(x -> nt[1] + nt[2]*x, label = "Newton", c=:green, lw=3)
    xlabel!("Semanas")
    ylabel!("Número de Casos")
    png("AproxLin")

    A = [ones(52) x x.^2];
    G = [52 sum(x) sum(x.^2); sum(x) sum(x.^2) sum(x.^3); sum(x.^2) sum(x.^3) sum(x.^4)];
    g = [sum(y); sum(x.*y); sum(x.^2 .*y)];
    z = ones(3);

    gr = gradiente(A, G, y, g, z)
    nt = newton(A, G, y, g, z)
    errograd = sqrt(sum((A*gr - y).^2))/norm(y);
    erronewt = sqrt(sum((A*nt - y).^2))/norm(y);
    print(" Aproximação Quadrática")
    print(" Coeficientes Gradiente = ", gr, " Coeficientes Newton = ", nt)
    print(" Erro Absoluto Gradiente = ", sqrt(sum((A*gr - y).^2)),
          " Erro Absoluto Newton =", sqrt(sum((A*nt - y).^2)))
    print(" Erro Relativo Gradiente = ", errograd, " Erro Relativo Newton = ", erronewt)

    scatter(x, y, title = "Contaminados COVID-19",label = "Casos por Semana", lw = 3)
    plot!(x -> gr[1] + gr[2]*x + gr[3]*x.^2, label = "Gradiente", c=:red, lw=3)
    plot!(x -> nt[1] + nt[2]*x + nt[3]*x.^2, label = "Newton", c=:green, lw=3)
    xlabel!("Semanas")
    ylabel!("Número de Casos")
    png("AproxQuad")

    A = [ones(52) x x.^2 x.^3];
    G = [52 sum(x) sum(x.^2) sum(x.^3);
        sum(x) sum(x.^2) sum(x.^3) sum(x.^4);
        sum(x.^2) sum(x.^3) sum(x.^4) sum(x.^5);
        sum(x.^3) sum(x.^4) sum(x.^5) sum(x.^6)];
    g = [sum(y); sum(x.*y); sum(x.^2 .*y); sum(x.^3 .*y)];
    z = ones(4);

    gr = gradiente(A, G, y, g, z)
    nt = newton(A, G, y, g, z)
    errograd = sqrt(sum((A*gr - y).^2))/norm(y);
    erronewt = sqrt(sum((A*nt - y).^2))/norm(y);
    print(" Aproximação Cúbica")
    print(" Coeficientes Gradiente = ", gr, " Coeficientes Newton = ", nt)
    print(" Erro Absoluto Gradiente = ", sqrt(sum((A*gr - y).^2)),
          " Erro Absoluto Newton =", sqrt(sum((A*nt - y).^2)))
    print(" Erro Relativo Gradiente = ", errograd, " Erro Relativo Newton = ", erronewt)


    scatter(x, y, title = "Contaminados COVID-19",label = "Casos por Semana", lw = 3)
    plot!(x -> gr[1] + gr[2]*x + gr[3]*x.^2 + gr[4]*x.^3, label = "Gradiente", c=:red, lw=3)
    plot!(x -> nt[1] + nt[2]*x + nt[3]*x.^2 + nt[4]*x.^3, label = "Newton", c=:green, lw=3)
    xlabel!("Semanas")
    ylabel!("Número de Casos")
    png("AproxCub")

end
