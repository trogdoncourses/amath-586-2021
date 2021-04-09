module amath586

using LinearAlgebra, Plots, LaTeXStrings

ρ = (α,z) -> (z.^(length(α)-1:-1:0))'*α
σ = (β,z) -> (z.^(length(β)-1:-1:0))'*β
R = (α,β,z) -> ρ(α,z)/σ(β,z)

export convergence_stability, run_newton, ChebyshevT, ChebyshevU, DChebyshevT


function ChebyshevT(n::Int64,x)
    if n < 0
        return 0.
    elseif n == 0
        return 1.
    elseif n == 1
        return x
    else
        u1 = x
        u0 = 1
        u2 = 0.
        for j = 2:n
            u2 = 2*x*u1-u0
            u0 = u1
            u1 = u2
        end
        return u2
    end 
end

function ChebyshevU(n::Int64,x)
    if n < 0
        return 0.
    elseif n == 0
        return 1.
    elseif n == 1
        return 2*x
    else
        u1 = 2*x
        u0 = 1
        u2 = 0.
        for j = 2:n
            u2 = 2*x*u1-u0
            u0 = u1
            u1 = u2
        end
        return u2
    end 
end

DChebyshevT(n::Int64,x) = n*ChebyshevU(n-1,x)

function find_roots(c) # supposing that the leading order coefficient is 1
    # c contains the rememing coefficients
    r = length(c)
    A = zeros(Complex{Float64},r,r)
    A[1,:] = -c
    A[2:end,1:end-1] = A[2:end,1:end-1] + I # add identity matrix to lower-left block
    return eigvals(A)
end

function check_condition(λ)
    if maximum(abs.(λ)) > 1
        return 0
    else
        for i = 1:length(λ)
            if abs(λ[i]) ≈ 1. && sum(map(t -> λ[i] ≈ t,λ)) > 1
                return 0
            end
        end
    end
    return 1
end

function root_condition(α,β,z)
    r = length(α)-1
    c = α-z*β
    if α[1]-z*β[1] ≈ 0.
        λ = find_roots(c[3:end]/c[2])
    else
        λ = find_roots(c[2:end]/c[1]) # let's suppose that first and second coefficients don't vanish simultaneously
    end
    return check_condition(λ)
end

function check_convergence(α,β) #supposing that α[1] = 1
    r = length(α)-1
    if sum(α) ≈ 0. && sum(β) ≈ sum(α.*(r:-1:0))
        println("Method is consistent")
    else
        println("Method is inconsistent")
    end
    if root_condition(α,β,0.) == 1
        println("Method is zero-stable")
    else
        println("Method is unstable")
    end
end

function run_newton(g,Dg,u₀,tol,max_iter)
    Unew = u₀
    Uold = u₀
    for j = 1:max_iter
        Unew = Uold - (Dg(Uold)\g(Uold))
        if maximum(abs.(Unew-Uold)) < tol 
            return Unew
        end
        Uold = Unew
    end
    println("Newton didn't terminate")
    return Unew
end

function convergence_stability(α,β)
    check_convergence(α,β)
    θ = 0:0.01(1+rand()/10):2*π # random perturbation to avoid singularities
    z = map(t -> R(α,β,exp(1im*t)),θ);
    
    if abs(minimum(real(z)) - maximum(real(z))) < .1
        xrange = [-4.,4.]
    else 
        xrange = [minimum(real(z))-1,maximum(real(z))+1]
    end
    
    if abs(minimum(imag(z)) - maximum(imag(z))) < .1
        yrange = [-4.,4.]
    else 
        yrange = [minimum(imag(z))-1,maximum(imag(z))+1]
    end
    
    if xrange[1] < -10
        xrange[1] = -4
    end
    if yrange[1] < -10
        yrange[1] = -4
    end
    if xrange[2] > 10
        xrange[2] = 4
    end
    if yrange[2] > 10
        yrange[2] = 4
    end
    
    contourf(xrange[1]:0.01:xrange[2],yrange[1]:0.01(1+rand()/10):yrange[2],(x,y)-> root_condition(α,β,x+1im*y),colorbar=false)
    plot!(real(z),imag(z),xlim=xrange,ylim=yrange,aspectratio=1,legend=false,lw=4,linecolor=:orange)
end

end