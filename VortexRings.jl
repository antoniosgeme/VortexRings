using PotentialFlow
using Plots
gr() #Plotting Backend
default(grid = false) #Plotting default options
tkfont = Plots.font("Times New Roman",15) #Plotting default options

#Function to compute the velocities that the system induces on itself
function compute_ẋ!(ẋ, x, t)
    reset_velocity!(ẋ, x) # Allocates an array for velocity computation
    self_induce_velocity!(ẋ, x, t) # Computes induced velocity of each vortex element
    
end

#Function to create a vortex patch from vortices
function vortex_patch!(vort,zc,Γ,radius,nring::Int;δ=0)
    Δr = radius/(nring-1/2) # Distance between rings
    dΓ = Γ/(1+8*nring*(nring-1)/2) # Differential circulation per vortex
    push!(vort,Vortex.Blob(zc,dΓ,δ)) # Create centroid vortex
    for ir in 1:nring-1
        nθ = 8*ir
        for j = 0:nθ-1
            push!(vort,Vortex.Blob(zc + ir*Δr*exp(im*2π*j/nθ),dΓ,δ)) # Create other vorticies
        end 
    end
    return vort
end

#Convenience constructor
vortex_patch(zc,Γ,radius,nring::Int;δ=0) = vortex_patch!(Vortex.Blob[],zc,Γ,radius,nring,δ=δ)

# Function to plot flow field 
function plotme(final_sys)
    p = plot(ratio=1,legend=:none,#=xlims = (-8,3),=# ylims = (-2,2), markerstrokewidth=0, markersize=2, tickfont=tkfont,dpi=300)
    for i in 1:length(final_sys)
        marker_z = i ÷ 2 == 0 ? 1 : -1
        plot!(final_sys[i],markerstrokewidth=0, markersize=2, marker_z = marker_z, color= cgrad(:RdBu))
    end 
    display(p)
end 


# function to time march for specified time using RK4. 
function time_march(time)
    T = 0:Δt:time-Δt
    for (iter,tloc) in enumerate(T)
        println("Current time: $t of $time")
        global sys₊,sys,t, plot_array# Declaring that variables to be modified are outside the function namescope 
        TimeMarching.rk4!(sys₊, sys, t, Δt, compute_ẋ!, advect!, ẋs) # March forward in time
        sys₊, sys = sys, sys₊ # save new system
        t += Δt 
        x,y = extract_coordinates(sys)
        #final_x,final_y,final_z = assign_z(x,y,30)
        final_x,final_y,final_z = reorder(x,y,30)
        if mod(iter,1) == 0
            push!(plot_array, scatter(final_x,final_y,final_z,camera=(60,60),markersize=2,size=(800,800),markercolor=:blue,legend=false,zlim=(-1.2,1.2),ylim=(-1.2,1.2)))
        end 
    end
end 

function reorder(x,y,n)
    final_x = Array{Float64,2}(undef,n,length(x))
    final_y = Array{Float64,2}(undef,n,length(y))
    final_z = Array{Float64,2}(undef,n,length(y))

    θ = 0: 2π/n :2π - 2π/n

    for i in 1:length(x)
        for k in 1:n
            final_x[k,i] = x[i] 
            final_y[k,i] = y[i] * cos(θ[k])
            final_z[k,i] = y[i] * sin(θ[k])
        end 
    end 
    return final_x,final_y,final_z
end

function extract_coordinates(sys)
    x = Float64[]
    y = Float64[]
    for patch in sys
        append!(x,real(Elements.position(patch)))
        append!(y,imag(Elements.position(patch)))
    end
    
    return x,y
end


function makegif(plot_array)
    anim = animate(plot_array)
    #gif(anim,"VortexRings.gif",fps=1)
end


r0 = 0.3 # initial radius of the vortex patch
d0 = 1    # initial distance between patch centroids
h0 = 2 # Initial horizontal distance between patches
Γ0 = 1.0 # strength of patch.
nring = 5   # number of rings in each patch.
Δt = 0.01*π^2*d0^2/abs(Γ0) # set the time step
δ = 0.05 # Vortex regularization radius 
t = 0.0  # Inital time 

# All coordinates are in the complex plane
sys = (vortex_patch(h0 + 0.5im*d0,-1.25Γ0,r0,nring,δ=δ), # Top right patch
    vortex_patch(h0 + -0.5im*d0,1.25Γ0,r0,nring,δ=δ), # Bottom right
    vortex_patch(0.5im*d0,-0.75Γ0,r0,nring,δ=δ), #Top left
    vortex_patch(-0.5im*d0,0.75Γ0,r0,nring,δ=δ)) # Bottom left

npatch = length(sys) # Total number of patches
sys₊ = deepcopy(sys) # Preallocating for next timestep
ẋs = [allocate_velocity(sys) for k = 1:4] # Preallocating for RK4
plot_array = []
time_march(50) 
makegif(plot_array)














