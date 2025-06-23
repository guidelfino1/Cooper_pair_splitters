include("/Users/guilhermedelfino/Documents/Y_junction/molecule_model/core_molecule.jl")





res = 50

N = 2 #only even values for N
M=Int(N/2)

t = 1.0
gamma = 0.7*t


phi_arr = [i for i=-1:0.25:1]
delta_arr = [i for i=0.001:0.3:2]

res_p= size(phi_arr)[1]
res_d = size(delta_arr)[1]

count = 0
for j=1:res_d,i=1:res_p
    phi = phi_arr[i]*pi
    delta = delta_arr[j]
    positions = order_parameter(res, N, M, phi,gamma,t, delta)
    path = "/Users/guilhermedelfino/Documents/Y_junction/molecule_model/order_parameter_N"*string(N)*"_negative/"
    ti = path*"Emin_phi"*string(phi_arr[i])*"delta"*string(delta)*".csv"

    open(ti, "w") do io
         writedlm(io, positions)
    end
    count = count+1
    println(count, " / 63:   " , "(Δ, Φ) =  ("*string(delta)*", "*string(phi_arr[i])*")")
end

