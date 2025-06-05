using Trixi
using LinearAlgebra
using Plots
using OrdinaryDiffEqLowStorageRK

a=0
b=2
T=1
c=1
polydeg=4

initial_condition_sine_wave(x)=sin(2π*x)
function advection_solve(dx;T,polydeg,c)
    
 n=Int((b-a)/dx)

 basis= LobattoLegendreBasis(polydeg)
 nodes= basis.nodes
 weights = basis.weights

 x= Matrix{Float64}(undef, length(nodes), n)
 for i in 1:n
     x_l = a +(i-1)*dx +dx/2
     for j in eachindex(nodes)
         e= nodes[j]
         x[j,i]=x_l+ dx/2 * e
     end
 end

 u0= initial_condition_sine_wave.(x)

 M=diagm(weights)
 B=diagm([-1; zeros(polydeg - 1); 1])
 D = basis.derivative_matrix

 surface_flux = flux_lax_friedrichs
 function rhs!(du,u,x,t)
     du.=zero(eltype(du))
     numerical_flux=copy(du)
     
     equation=LinearScalarAdvectionEquation1D(c)

     for i in 2:(n-1)
         numerical_flux[1,i] = surface_flux(u[end,i-1],u[1,i],1,equation)
         numerical_flux[end,i-1]= numerical_flux[1,i]

         numerical_flux[end,i]=surface_flux(u[end,i],u[1,i+1],1,equation)
         numerical_flux[1,i+1]=numerical_flux[end,i]
     end

     numerical_flux[1,1]=surface_flux(u[end,end],u[1,1],1,equation)
     numerical_flux[end,end]=numerical_flux[1,1]

     for i in 1:n
         du[:,i]= (2/dx) *(du[:,i]+ (M\ transpose(D))*M* u[:,i] - (M\B)*numerical_flux[:,i])
     end
     return nothing
 end
 dt=(0.5*(x[2]-x[1]))/((2*polydeg+1)*c)
 t = (0.0, T)
 ode = ODEProblem(rhs!, u0, t, x)

 soln = solve(ode, RDPK3SpFSAL49(); dt=dt,abstol = 1.0e-9, reltol = 1.0e-9,
             ode_default_options()...)

 return x, soln.u[end]
end

deltas=Float64[]

L2_errors = Float64[]

for dx in [0.25,0.2,0.125,0.1,0.05]   
 x,u = advection_solve(dx;T,polydeg,c)
  u_exact = sin.(2π .* (x .- c * T))
  error = abs.(u_exact .- u)
  push!(L2_errors, sqrt((x[2]-x[1]) * sum(error.^2)))
  push!(deltas, x[2]-x[1])
 if dx==0.05
  plot(title="Solutions at final time T=$T", xlabel="x", ylabel="u(x,T)", legend=:bottomleft)
  plot!(vec(x),vec(u),label=false)
  display(current())
 end
end

L2_cgs=Float64[]

for i in 1: length(deltas)-1
    push!(L2_cgs, log2(L2_errors[i] / L2_errors[i+1]) / log2(deltas[i] / deltas[i+1]))
end
 println(L2_errors)
 println(L2_cgs)

