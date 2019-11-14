__precompile__()

module PAPA

export papa_reconstruction, papa_reconstruction_debug

using LinearAlgebra, LeastSquaresOptim

greet() = print("Welcome to PAPA world!")
include("Process_Solver.jl")
include("SupportFunctions.jl")

"""
    papa_reconstruction(N, Npairs, pair_order, initial_chilist, sigmaVec[, CP_penalty, f_tol, iter])

Function for PAPA bootstrapping from a list of χ matrices to a larger χ matrix.

# Arguments
- `N::Integer`: the number of qubits
- `Npairs::Integer`: the number of qubit pairs

# Examples
```jldoctest
julia>

```
"""
function papa_reconstruction(N::Int,Npairs::Int,pair_order::Array{Int,2},initial_chilist::Vector{Float64},sigmaVec::Vector{Float64};CP_penalty = 1.0,x_tol=1.0e-7,f_tol=1.0e-7,iter=1000000,alg=LevenbergMarquardt(),flag=true)

    setup_tuple = setup_problem(N,Npairs,pair_order)

    Nel_tot = convert(Int,2*Npairs*16*17/2 + 16*2*Npairs + Npairs)

    function min_fun!(F::Vector{Float64},x::Vector{Float64})
        process_solver!(F,x,sigmaVec,setup_tuple,CP_penalty)
    end

    # function min_fun(x)
    #     return process_solver(x,sigmaVec,setup_tuple,CP_penalty)
    # end

    # result  = nlsolve(min_fun!,initial_chilist;xtol=x_tol,ftol=f_tol,iterations=iter)
    # chi_list_final = result.zero
    if flag
        print("Optimization Beginning\n")
    end
    lsa_papa = LeastSquaresProblem(x = initial_chilist,f! = min_fun!,output_length = Nel_tot)
    result = optimize!(lsa_papa,alg;x_tol=x_tol,f_tol=f_tol,iterations=iter)
    chilist = result.minimizer
    # optimize(min_fun,initial_chilist)

    dim2q = 16
    chi_PAPA = zeros(ComplexF64,Npairs,dim2q,dim2q)
    chitemp = zeros(ComplexF64,16,16)
    Nel = 136
    for nn = 1:1:Npairs
        fill!(chitemp,0.)
        start = 1 + (nn-1)*Nel
        stop = nn*Nel;
        chitemp[triu(trues(dim2q,dim2q))] = chilist[start:1:stop] + 1im*chilist[(start+Npairs*Nel):1:(stop+Npairs*Nel)]
        chitemp[:,:] = chitemp+chitemp'
        chitemp[:,:] = chitemp - diagm(0 => diag(chitemp))/2
        chi_PAPA[nn,:,:] = chitemp[:,:]
    end

    return (chi_PAPA, result.ssr)
end

function papa_reconstruction_debug(N::Int,Npairs::Int,pair_order::Array{Int,2},initial_chilist::Vector{Float64},sigmaVec::Vector{Float64};CP_penalty = 1.0,x_tol=1.0e-7,f_tol=1.0e-7,iter=1000000,alg=LevenbergMarquardt(),flag=true)

    setup_tuple = setup_problem(N,Npairs,pair_order)

    Nel_tot = convert(Int,2*Npairs*16*17/2 + 16*2*Npairs + Npairs)

    function min_fun!(F::Vector{Float64},x::Vector{Float64})
        process_solver!(F,x,sigmaVec,setup_tuple,CP_penalty)
    end

    # function min_fun(x)
    #     return process_solver(x,sigmaVec,setup_tuple,CP_penalty)
    # end

    # result  = nlsolve(min_fun!,initial_chilist;xtol=x_tol,ftol=f_tol,iterations=iter)
    # chi_list_final = result.zero
    if flag
        print("Optimization Beginning\n")
    end
    lsa_papa = LeastSquaresProblem(x = initial_chilist,f! = min_fun!,output_length = Nel_tot)
    result = optimize!(lsa_papa,alg;x_tol=x_tol,f_tol=f_tol,iterations=iter)
    chilist = result.minimizer
    # optimize(min_fun,initial_chilist)

    dim2q = 16
    chi_PAPA = zeros(ComplexF64,Npairs,dim2q,dim2q)
    chitemp = zeros(ComplexF64,16,16)
    Nel = 136
    for nn = 1:1:Npairs
        fill!(chitemp,0.)
        start = 1 + (nn-1)*Nel
        stop = nn*Nel;
        chitemp[triu(trues(dim2q,dim2q))] = chilist[start:1:stop] + 1im*chilist[(start+Npairs*Nel):1:(stop+Npairs*Nel)]
        chitemp[:,:] = chitemp+chitemp'
        chitemp[:,:] = chitemp - diagm(0 => diag(chitemp))/2
        chi_PAPA[nn,:,:] = chitemp[:,:]
    end

    return (chi_PAPA, result)
end

pair_order_1 = [1 2]

# Experimental data
sigma12_I = zeros(16,16)
max_ent_state!(4;mrho=sigma12_I)
sigma12_I = convert(Array{ComplexF64,2},sigma12_I)

# Conver to vector of data

sigmaVec_I = [real(sigma12_I[triu(trues(16,16))]);imag(sigma12_I[triu(trues(16,16))])]

# Initial state
chi12_idI = zeros(ComplexF64,16,16)

for i1 = 1:1:4
    for i2 = 1:1:4
        for i3 = 1:1:4
            for i4 = 1:1:4
                chi12_idI[4*(i1-1)+i2,4*(i3-1)+i4] = 4*sigma12_I[4*(i2-1)+i1,4*(i4-1)+i3]
            end
        end
    end
end

initial_chilist_I = [real(chi12_idI[triu(trues(16,16))]);imag(chi12_idI[triu(trues(16,16))])]

function __init__()
    print("Precompiling for speed\n")
    papa_reconstruction(2,1,pair_order_1,initial_chilist_I,sigmaVec_I,CP_penalty = 1.0,x_tol=1.0e-7,f_tol=1.0e-7,iter=1000000,alg=LevenbergMarquardt(),flag=false);
    papa_reconstruction(2,1,pair_order_1,initial_chilist_I,sigmaVec_I,CP_penalty = 1.0,x_tol=1.0e-7,f_tol=1.0e-7,iter=1000000,alg=Dogleg(),flag=false);
    nothing
end


end # module
