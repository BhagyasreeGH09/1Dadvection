using Plots

L=1
T=1
c=0.2

solutions = []

for dx in [0.002,0.001,0.0005,0.00025]

 N=Int(L/dx)+1
 x=LinRange(0,1,N)

dt=0.5*dx/c
 M = Int(T / dt)

 
v=c*(dt/dx)
u = [0.0 for xi in x]
for i in 1:N
    if 0.4 <= x[i] <= 0.6
        u[i] = 1.0
    end
end

 for j in 1:M
    uold=copy(u)
    for i in 2:N-1
        u[i]=(uold[i+1]+uold[i-1])/2 -v*(uold[i+1] - uold[i-1])/2
    end
    u[1]=0
    u[N]=0
 end 
 
 push!(solutions, (x, copy(u)))
end

plot(title="Solutions at final time T=$T", xlabel="x", ylabel="u(x,T)", legend=:bottomleft)

for (x, u) in solutions
    plot!(x, u, label="dx=$(round(x[2]-x[1], digits=5))")
end
 
x_fine = LinRange(0, 1, 1000)
u_exact_fine = u_exact_fine = [if  0.4<=xi - c*T <= 0.6 1.0 else 0.0 end for xi in x_fine]
plot!(x_fine, u_exact_fine, lw=2, linestyle=:dash, color=:black, label="Exact")
display(current())
