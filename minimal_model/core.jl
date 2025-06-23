using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")

using LinearAlgebra
using DelimitedFiles


function ham(kx,ky, phi, delta, theta1, theta2)
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

    mat = [0 (a+a*d12) (conj(a)+conj(a)*d13) (delta) 0 0; 
    (conj(a)+conj(a)*d21) 0 (a+a*d23) 0 (exp(1im*theta1)*delta) 0; 
    (a+a*d31) (conj(a)+conj(a)*d32) 0 0 0 (exp(1im*theta2)*delta); 
    (conj(delta)) 0 0 0 (-conj(a)-conj(a)*d12) (-a-a*d13); 
    0 conj(exp(1im*theta1)*delta) 0 -(a+a*d21) 0 -(conj(a)+conj(a)*d23); 
    0 0 conj(exp(1im*theta2)*delta) (-conj(a)-conj(a)*d31) -(a+a*d32) 0]

    return mat
end


function eigsys(kx,ky,phi, delta,theta1, theta2)
   
    mat = ham(kx,ky, phi, delta, theta1, theta2)

    vals, vec = eigen(mat)

    return vals,vec

end




function claudiobz(kx,ky)
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


function energy(res, phi, delta, theta1, theta2)
    kx_ar = LinRange(-pi,pi,res)
    ky_ar = LinRange(-pi,pi,res)
    delta_k = kx_ar[2] - kx_ar[1]

    ans = 0

    for kx in kx_ar, ky in ky_ar
        if claudiobz(kx, ky)
            en = eigsys(kx,ky, phi,delta, theta1, theta2)[1]
            # ans = ans+ en[1]+en[2]+en[3]
            for a=1:3
                if en[a]<0
                    ans = ans+ en[a]
                end
            end
        end

    end
    return ans*delta_k^2

end




function chern_number_half_filling(kx0,ky0, res, phi, delta, theta1, theta2)
    #input BZ point kx0 and ky0; resolution, band number 1<a<dim and SC pairing delta
    #outputs the Berry curvature F for the a band at point (kx0, ky0)

    deltak = 2*pi/res

    uk0 = eigsys(kx0,ky0,phi,delta,theta1, theta2)[2]
    uk1 = eigsys(kx0+deltak,ky0,phi, delta,theta1, theta2)[2]
    uk2 = eigsys(kx0,ky0+deltak,phi, delta,theta1, theta2)[2]
    uk3 = eigsys(kx0+deltak,ky0+deltak,phi, delta,theta1, theta2)[2]
    ans=1

    for a=1:3

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

        #berr = zeros(res,res,dim,dim)

        ans = ans*U1*U3*conj(U4)*conj(U2)

    end

    logans = log(ans)


    return logans

end

function FHS2(kx0,ky0, res, phi, delta, theta1, theta2)
    #input BZ point kx0 and ky0; resolution, band number 1<a<dim and SC pairing delta
    #outputs the Berry curvature F for the a band at point (kx0, ky0)

    deltak = 2*pi/res

    uk0 = eigsys(kx0,ky0,phi,delta,theta1, theta2)[2]
    uk1 = eigsys(kx0+deltak,ky0,phi, delta,theta1, theta2)[2]
    uk2 = eigsys(kx0,ky0+deltak,phi, delta,theta1, theta2)[2]
    uk3 = eigsys(kx0+deltak,ky0+deltak,phi, delta,theta1, theta2)[2]
    ans=1

    for a=1:2

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

        #berr = zeros(res,res,dim,dim)

        ans = ans*U1*U3*conj(U4)*conj(U2)

    end

    logans = log(ans)


    return logans

end



function chernFHS(res, phi, delta, theta1, theta2)

    k_arr= LinRange(-pi,pi,res)

        numb3 = 0
        for kx in k_arr, ky in k_arr
            if claudiobz(kx, ky)
                # numb2 +=FHS2(kx,ky, res, phi, delta, theta1, theta2)
                numb3 +=chern_number_half_filling(kx,ky, res, phi, delta, theta1, theta2)
                
            end
        end

    return real((-1*im/(2*pi))*numb3)

end

function chernFHS2(res, phi, delta, theta1, theta2)

    k_arr= LinRange(-pi,pi,res)

        numb2 = 0
        for kx in k_arr, ky in k_arr
            if claudiobz(kx, ky)
                numb2 +=FHS2(kx,ky, res, phi, delta, theta1, theta2)
                # numb3 +=chern_number_half_filling(kx,ky, res, phi, delta, theta1, theta2)
                
            end
        end

    return real((-1*im/(2*pi))*numb2)

end





