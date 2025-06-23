include("/Users/guilhermedelfino/Documents/Y_junction/minimal_model/core.jl")

dim=6;
res=200


#reduced phase diagram
# phi_ar = [i for i=3/2-1/4:0.125/4:3/2+1/4]
# res_p=1

# phi_ar = [i for i=1/2-0.125:0.125/4:1/2+0.125]

#fine grid shifted
# delta_ar = [i for i=0.001+0.01:0.02:1.2+0.01]
# phi_ar = [i for i=-1:0.125/8:1]

delta_ar = [i for i=0.001:0.02:1.2]
phi_ar = [i for i=-1:0.125/8:1]

# delta_ar = [i for i=0.001:0.02:2]

res_p= size(phi_ar)[1]
res_d = size(delta_ar)[1]


path_write = "/Users/guilhermedelfino/Documents/Y_junction/minimal_model/chern_number_fine_grid/"
path_read = "/Users/guilhermedelfino/Documents/Y_junction/minimal_model/order_parameter_fine_grid/"

cou = 0

for i=1:res_p, j=1:res_d
    cou = cou+1

    phi = phi_ar[i]*pi
    delta = delta_ar[j]

    title_read = path_read*"Emin_phi"*string(phi_ar[i])*"delta"*string(delta)*".csv"

    phi12 = readdlm(title_read, '\t', Float64, '\n')

    phi1 = phi12[1]*pi #does not contain pi factors
    phi2 = phi12[2]*pi #does contains pi factors
    # println(phi1, phi2)

    # println(phi1)
    # println(phi2)
    
    cn = chernFHS(res, phi, delta, phi1, phi2)
    println("(Δ, ϕ) = (", delta_ar[j]," ,",phi_ar[i],")")
    println(cou, " / ",res_p*res_d)

    #  println(cn)

    title_write = path_write*"Cn"*string(phi_ar[i])*"delta"*string(delta)*".csv"
    open(title_write, "w") do io
             writedlm(io, cn)
     end
end
