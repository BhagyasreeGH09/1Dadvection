using Plots

L=1
T=1
c=1

L1_errors=Float64[]
L2_errors=Float64[]
Linf_errors=Float64[]
deltas = Float64[]
solutions = []
dt=0.001
M = Int(T / dt)
for dx in [0.01,0.005,0.0025,0.002]

 N=Int(L/dx)+1
 x=LinRange(0,1,N)

v=c*(dt/dx)
 u = [sin(2π * xi) for xi in x]

 for j in 1:M
   uold=copy(u)
   for i in 2:N-1
      u[i]=uold[i]-v*(uold[i+1]-uold[i-1])/2+ (v^2)*(uold[i+1] - 2*uold[i] + uold[i-1])/2
end
    u[N]= uold[N] - v*(uold[1] - uold[N-1])/2 + (v^2)*(uold[1] - 2*uold[N] + uold[N-1])/2
    #u[1] = uold[1] -v*(uold[2] - uold[N])/2 + (v^2)*(uold[2] - 2*uold[1] + uold[N])/2
    u[1]=u[N]
end


  u_exact = [sin(2π * (xi - c*T)) for xi in x]

  error = [abs(u_exact[i] - u[i]) for i in 1:N]
  push!(L1_errors, dx * sum(error))
  push!(L2_errors, sqrt(dx * sum(e^2 for e in error)))
  push!(Linf_errors, maximum(error))
  push!(deltas, dx)
  push!(solutions, (x, copy(u)))
end
L1_cgs=Float64[]
L2_cgs=Float64[]
Linf_cgs=Float64[]

for i in 1: length(deltas)-1
     push!(L1_cgs, log2(L1_errors[i] / L1_errors[i+1]) / log2(deltas[i] / deltas[i+1]))
    push!(L2_cgs, log2(L2_errors[i] / L2_errors[i+1]) / log2(deltas[i] / deltas[i+1]))
    push!(Linf_cgs, log2(Linf_errors[i] / Linf_errors[i+1]) / log2(deltas[i] / deltas[i+1]))
end
 println(L1_errors)
 println(L2_errors)
 println(Linf_errors)
 println(L1_cgs)
 println(L2_cgs)
 println(Linf_cgs)

plot(title="Solutions at final time T=$T", xlabel="x", ylabel="u(x,T)", legend=:bottomleft)

for (x, u) in solutions
    plot!(x, u, label="dx=$(round(x[2]-x[1], digits=5))")
end
 
x_fine = LinRange(0, 1, 1000)
u_exact_fine = [sin(2π * (xi - c * T)) for xi in x_fine]
plot!(x_fine, u_exact_fine, lw=2, linestyle=:dash, color=:black, label="Exact")
display(current())
