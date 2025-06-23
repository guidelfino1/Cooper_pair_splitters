include("/Users/guilhermedelfino/Documents/Y_junction/molecule_model/core_molecule.jl")









res=200

N = 2 #only even values for N
M=Int(N/2)

t = 1.0
gamma = 0.7*t


#New standard choice
# phi_arr = [i for i=-1:0.125/2:1]
# delta_arr = [i for i=0.001:0.02:2]

#Reduced phase diagram
phi_arr = [i for i=-1:0.25:1]
delta_arr = [i for i=0.001:0.3:2]

res_p= size(phi_arr)[1]
res_d = size(delta_arr)[1]


path_write = "/Users/guilhermedelfino/Documents/Y_junction/molecule_model/chern_number_N"*string(N)*"_negative/"
path_read = "/Users/guilhermedelfino/Documents/Y_junction/molecule_model/order_parameter_N"*string(N)*"_negative/"



for i=1:res_p, j=1:res_d

    phi = phi_arr[i]*pi
    delta = delta_arr[j]

    title_read = path_read*"Emin_phi"*string(phi_arr[i])*"delta"*string(delta)*".csv"

    theta12 = readdlm(title_read, '\t', Float64, '\n')

    theta1 = theta12[1]*pi #does not contain pi factors
    theta2 = theta12[2]*pi #does contains pi factors

    
    cn = chernFHS(res, N, M, phi,gamma,t, delta, theta1, theta2)
    
    #  println(cn)

    title_write = path_write*"Cn"*string(phi_arr[i])*"delta"*string(delta)*".csv"
    open(title_write, "w") do io
             writedlm(io, cn)
     end

     println("(Δ, Φ) =  ("*string(delta)*", "*string(phi_arr[i])*")"*" ----> c3 = "*string(cn))
end
