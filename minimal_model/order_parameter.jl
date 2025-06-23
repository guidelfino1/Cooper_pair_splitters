include("/Users/guilhermedelfino/Documents/Y_junction/minimal_model/core.jl")




dim=6

res = 50; 

#reduced phase diagram
# phi_ar = [i for i=3/2-1/4:0.125/4:3/2+1/4]
# res_p=1

# delta_ar = [i for i=0.001:0.3:2]

#fine grid
# delta_ar = [i for i=0.001+0.01:0.02:1.2+0.01]
# phi_ar = [i for i=-1:0.125/8:1]

delta_ar = [0.001]
phi_ar = [i for i=-1:0.125/16:1]

#extended phase diagram
# phi_ar = [i for i=-3:0.125/2:-1.0]
# delta_ar = [i for i=0.001:0.02:2]

res_p= size(phi_ar)[1]
res_d = size(delta_ar)[1]


phi_grid =  [i for i=-1:1/12:1]
grid_phi = size(phi_grid)[1]

cou = 0

for j=1:res_d,i=1:res_p

    cou = cou+1

    en_phi = zeros(grid_phi, grid_phi)
    phi = phi_ar[i]*pi
    delta = delta_ar[j]

    for m=1:grid_phi, n=1:grid_phi
        phim = phi_grid[m]*pi
        phin = phi_grid[n]*pi
        en_phi[m,n] = energy(res, phi, delta, phim, phin)
    end

    Emin =minimum(en_phi)
    coord = findfirst(item -> item == Emin, en_phi)
    sav = [phi_grid[coord[1]] phi_grid[coord[2]]] #no factors of pi
    path = "/Users/guilhermedelfino/Documents/Y_junction/minimal_model/phase_boundary/"
    ti = path*"Emin_phi"*string(phi_ar[i])*"delta"*string(delta)*".csv"

    println("(Δ, ϕ) = (", delta_ar[j]," ,",phi_ar[i],")")
    println(cou, " / ",res_p*res_d)

    open(ti, "w") do io
         writedlm(io, sav)
    end
end


