using Trixi
using LinearAlgebra
using Plots
using OrdinaryDiffEqLowStorageRK

a=0
b=2
T=1
c=1
polydeg=5
dx=0.125

initial_condition_step(x) = (0.2 ≤ x ≤ 0.5) ? 1.0 : 0.0
    
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

 u0= initial_condition_step.(x)

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

     numerical_flux[1,1]=surface_flux(0,u[1,1],1,equation)
     numerical_flux[end,end]=surface_flux(u[end,end],0,1,equation)

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

 plot(vec(x), vec(soln.u[end]), label = "Numerical solution at t=$(t[2])", legend = :topleft,
     lw = 3)

x_fine = LinRange(0, 2, 1000)
u_exact_fine = [if  0.2<=xi - c*T <= 0.5 1.0 else 0.0 end for xi in x_fine]
plot!(x_fine, u_exact_fine, lw=2, linestyle=:dash, color=:black, label="Exact solution at t=$(t[2])")
display(current())



