using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")

using LinearAlgebra
using DelimitedFiles
using Plots



function ham_molecule(N,M, phi,gamma,t, delta, theta1, theta2,kx, ky)
    #returns a dim = 6(4M+N) Hamiltonian. 
    #Here we are using phi_+ - phi_- := phi

    dim_ham = 4*M+N #dimension of Hamiltonian in a section of wire
    
    a = exp(1im*phi/3)
    d12e = -ky-(sqrt(3)*kx/2+ky/2)
    d13e = -ky-(-sqrt(3)*kx/2+ky/2)
    d23e = sqrt(3)*kx
    d12 = exp(1im*d12e)
    d13 = exp(1im*d13e)
    d21 = exp(-1im*d12e)
    d31 = exp(-1im*d13e)
    d23 = exp(1im*d23e)
    d32 = exp(-1im*d23e)

    D = zeros(ComplexF64,(dim_ham,dim_ham))
    alpha_plus = zeros(ComplexF64,(dim_ham,dim_ham))
    alpha_minus = zeros(ComplexF64,(dim_ham,dim_ham))

    for i=1:2*M, j=1:2*M
        if abs(i-j)==2
            D[i,j] = -t
        end
    end

    if M!=0
        D[2*M-1,2*M+1]=-t
        D[2*M+1,2*M-1]=-t

        D[2*M,2*M+1]=-t
        D[2*M+1,2*M]=-t


        D[2*M+N,2*M+N+1] = -t
        D[2*M+N+1,2*M+N] = -t

        D[2*M+N,2*M+N+2] = -t
        D[2*M+N+2,2*M+N] = -t

    end

    for i=2*M+1:2*M+N, j=2*M+1:2*M+N
        if abs(i-j)==1
            D[i,j] = -t
        end
    end

    

    for i=2*M+N+1:4*M+N, j=2*M+N+1:4*M+N
        if abs(i-j)==2
            D[i,j]=-t
        end
    end
    
    alpha_plus[1,2] = -gamma*a
    alpha_minus[dim_ham-1,dim_ham] = -gamma*a

    delta_mat = delta*[Matrix(I, dim_ham, dim_ham) zeros((dim_ham,dim_ham)) zeros((dim_ham,dim_ham));
    zeros((dim_ham,dim_ham)) exp(1im*theta1)*Matrix(I, dim_ham, dim_ham) zeros((dim_ham,dim_ham));
    zeros((dim_ham,dim_ham)) zeros((dim_ham,dim_ham)) exp(1im*theta2)*Matrix(I, dim_ham, dim_ham)]

    hk = [D alpha_plus+alpha_minus*d12 transpose(conj(alpha_plus))+transpose(conj(alpha_minus))d13;
    transpose(conj(alpha_plus))+transpose(conj(alpha_minus))d21 D alpha_plus+alpha_minus*d23;
    alpha_plus+alpha_minus*d31 transpose(conj(alpha_plus))+transpose(conj(alpha_minus))d32 D]

    h_minus_k = [D alpha_plus+alpha_minus*d21 transpose(conj(alpha_plus))+transpose(conj(alpha_minus))d31;
    transpose(conj(alpha_plus))+transpose(conj(alpha_minus))d12 D alpha_plus+alpha_minus*d32;
    alpha_plus+alpha_minus*d13 transpose(conj(alpha_plus))+transpose(conj(alpha_minus))d23 D]

    mat = [hk delta_mat;
    conj(delta_mat) -conj(h_minus_k)]

    return mat
end




function eigsys(N,M, phi,gamma,t, delta, theta1, theta2,kx, ky)
    mat = ham_molecule(N,M, phi,gamma,t, delta, theta1, theta2,kx, ky)

    vals, vec = eigen(mat)

    return vals,vec

end

function brillouin(kx,ky)
    if ((kx/sqrt(3) + ky/3 < 4*pi/9) &&
        (-kx/sqrt(3) + ky/3 < 4*pi/9) && 
        (2*ky/3 < 4*pi/9) &&
        (-kx/sqrt(3) - ky/3 < 4*pi/9) &&
        (kx/sqrt(3)- ky/3 < 4*pi/9) &&
        (-2*ky/3 < 4*pi/9))
        ans = true
    else
        ans=false
    end

    return ans
end


function energy(res, N, M, phi,gamma,t, delta, theta1, theta2)
    dim_ham = 4*M+N

    k_arr = LinRange(-pi,pi,res)
    delta_k = k_arr[2] - k_arr[1]

    ans = 0

    for kx in k_arr, ky in k_arr
        if brillouin(kx, ky)
            en = eigsys(N,M, phi,gamma,t, delta, theta1, theta2,kx, ky)[1]
            #cumulative = 0
            for j=1:3*dim_ham
                #cumulative = cumulative + en[j]
                ans = ans+en[j]
            end
            #ans = ans+cumulative
        end

    end
    return ans*delta_k^2

end

function order_parameter(res, N,M, phi,gamma,t, delta)
    # grid_theta =  [i for i=-1:1/12:1]
    # size_theta = size(grid_theta)[1]

    grid_theta = [-0.6666666666666666 0 0.6666666666666666] #reduced grid, to save time
    size_theta = 3
    

    en_theta = zeros(size_theta, size_theta)

    for m=1:size_theta, n=1:size_theta
        thetam = grid_theta[m]*pi
        thetan = grid_theta[n]*pi
        en_theta[m,n] = energy(res, N, M, phi, gamma,t, delta, thetam, thetan)
    end

    Emin =minimum(en_theta)
    println(Emin)
    coord = findfirst(item -> item == Emin, en_theta)
    sav = [grid_theta[coord[1]] grid_theta[coord[2]]] #no factors of pi

    return sav
end






function chern_number_half_filling(res, N, M, phi,gamma,t, delta, theta1, theta2, kx0, ky0)
    #input BZ point kx0 and ky0; N, resolution, and Hamiltonian paramters

    deltak = 2*pi/res

    dim_ham = 4*M+N

    uk0 = eigsys(N,M, phi,gamma,t, delta, theta1, theta2,kx0, ky0)[2]
    uk1 = eigsys(N,M, phi,gamma,t, delta, theta1, theta2,kx0+deltak,ky0)[2]
    uk2 = eigsys(N,M, phi,gamma,t, delta, theta1, theta2,kx0,ky0+deltak)[2]
    uk3 = eigsys(N,M, phi,gamma,t, delta, theta1, theta2,kx0+deltak,ky0+deltak)[2]
    ans=1

    for a=1:3*dim_ham

        U1 = transpose(conj(uk0[:,a]))*uk1[:,a]
        norm1 = abs(U1)
        U1 = U1/norm1

        U2 = transpose(conj(uk0[:,a]))*uk2[:,a]
        norm2 = abs(U2)
        U2 = U2/norm2

        U3 = transpose(conj(uk1[:,a]))*uk3[:,a]
        norm3 = abs(U3)
        U3 = U3/norm3

        U4 = transpose(conj(uk2[:,a]))*uk3[:,a]
        norm4 = abs(U4)
        U4 = U4/norm4


        ans = ans*U1*U3*conj(U4)*conj(U2)

    end

    logans = log(ans)

    return logans

end



function chernFHS(res, N,M, phi,gamma,t, delta, theta1, theta2)
    k_arr= LinRange(-pi,pi,res)

        numb3 = 0
        for kx0 in k_arr, ky0 in k_arr
            if brillouin(kx0, ky0)
                numb3 +=chern_number_half_filling(res, N,M, phi,gamma,t, delta, theta1, theta2, kx0, ky0)
            end
        end

    return real((-1*im/(2*pi))*numb3)

end





function overprint(str)  
    print("\u1b[1F")
    #Moves cursor to beginning of the line n (default 1) lines up   
    print(str)   #prints the new line
   print("\u1b[0K") 
   # clears  part of the line.
   #If n is 0 (or missing), clear from cursor to the end of the line. 
   #If n is 1, clear from cursor to beginning of the line. 
   #If n is 2, clear entire line. 
   #Cursor position does not change. 

    println() #prints a new line, i really don't like this arcane codes
end