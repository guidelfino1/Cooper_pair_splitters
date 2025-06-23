include("/Users/guilhermedelfino/Documents/Y_junction/molecule_model/core_molecule.jl")


res1 = 50
res2 = 200

M=2

N=Int(2*M+1)


delta=1.0
#large t (test)
alp = 6
t = N*delta*alp


phi_arr = [i for i=-1:0.25:1]


#new gamma array
step = 7
gamma_min = 1.0
gamma_max = ((2*M+1)/(3*M+1))*alp
delta_gamma = (gamma_max - gamma_min)/step
gamma_arr = [i for i=gamma_min:delta_gamma:gamma_max]

res_p= size(phi_arr)[1]
res_g= size(gamma_arr)[1]

total_str = string(res_p*res_g)

path_write = "/Users/guilhermedelfino/Documents/Y_junction/molecule_model/even_M_v4/SC_chern_M"*string(M)*"_alp"*string(alp)*"/"
path_read = "/Users/guilhermedelfino/Documents/Y_junction/molecule_model/even_M_v4/SC_phase_M"*string(M)*"_alp"*string(alp)*"/"

count = 0
for j=1:res_g,i=1:res_p
    phi = phi_arr[i]*pi
    gamma = gamma_arr[j]
    gamma0 = (3*M+1)*gamma
    positions = order_parameter(res1, N, M, phi,gamma0, t, delta)
    ti = path_read*"phase"*string(phi_arr[i])*"gamma"*string(round(gamma, digits=4))*".csv"

    open(ti, "w") do io
         writedlm(io, positions)
    end


    count = count+1
    println(count, " / "*total_str*":   " , "(Γ, Φ) =  ("*string(round(gamma, digits=4))*", "*string(phi_arr[i])*")")

    theta1 = positions[1]*pi
    theta2 = positions[2]*pi

    cn = chernFHS(res2, N, M, phi,gamma0,t, delta, theta1, theta2)

    title_write = path_write*"Cn"*string(phi_arr[i])*"gamma"*string(round(gamma, digits=4))*".csv"

    open(title_write, "w") do io
             writedlm(io, cn)
     end


    println("----> c3 = "*string(cn))
end

