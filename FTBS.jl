using Plots

L=1
T=2
c=1

errors = Float64[]
deltas = Float64[]
solutions = []

for dx in [0.01,0.005,0.0025,0.002,0.001]

 N=Int(L/dx)+1
 x=LinRange(0,1,N)

dt=0.5*dx/c
 M = Int(T / dt)

 u = [sin(2π * xi) for xi in x]

 for j in 1:M
    uold=copy(u)
    for i in 2:N
        u[i]=uold[i]-c*(dt/dx)*(uold[i]-uold[i-1])
    end
    u[1]=u[N]
 end 
 u_exact = [sin(2π * (xi - c*T)) for xi in x]
 error = sqrt(sum((u .- u_exact).^2) * dx)
 push!(errors, error)
 push!(deltas, dx)
 push!(solutions, (x, copy(u)))
end

println("Order of convergence:")
for i in 1:length(errors)-1
    p = log2(errors[i] / errors[i+1])/log2(deltas[i] / deltas[i+1])
    println(p)
end

plot(title="Solutions at final time T=$T", xlabel="x", ylabel="u(x,T)", legend=:bottomleft)

for (x, u) in solutions
    plot!(x, u, label="dx=$(round(x[2]-x[1], digits=5))")
end
 
x_fine = LinRange(0, 1, 1000)
u_exact_fine = [sin(2π * (xi - c * T)) for xi in x_fine]
plot!(x_fine, u_exact_fine, lw=2, linestyle=:dash, color=:black, label="Exact")
display(current())