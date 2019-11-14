export partial_trace, trace_distance, post_CPTP!, max_ent_state!, iden, bases_create!, full_process

"""
    post_CPTP!(ρ)

Applies a CPTP contraction to the input Hermitian operator
Input arguments:
    ρ: Hermitian operator to be transformed to a positive semi-definite version
"""
function post_CPTP!(rho::Array{T,2}) where {T <: Number}

    A = eigen(rho)
    A.values[real(A.values) .<0.] .= 0.
    rho[:,:] = A.vectors*Diagonal(A.values)*A.vectors'
    nothing
end

"""
    trace_distance(A, B)

Calculates the trace distance between two matrices
Input arguments:
    A, B: matrices of the same dimensions
Output arguments:
    : trace distance
"""
function trace_distance(A::Array{T1,2},B::Array{T2,2}) where {T1 <: Number, T2 <: Number}

    tr(sqrt((A-B)'*(A-B)))/2
end

"""
    partial_trace(rho, qs, out, dim, num)

This function calculates the partial trace of a density matrix
Input arguments:
    rho: full density matrix
    qs: a vector of Ints of length num, with a 1 if the qudit is to be traced out and a 0 otherwise
    dim: vector of Ints of the dimensions of the each qudit (same order as qs)
    num: total number of subsystems
Output:
       rhoPT: the partial trace density matrix
"""
function partial_trace(rho::Array{T,2},qs::Array{Int64,1},dim::Array{Int64,1},num::Int64) where {T <: Number}

    rhoPT = rho
    for ii = 1:1:num
        tempV = zeros(dim[ii])

        if qs[ii] == 1
            dim1 = dim[1:ii]
            qs1 = qs[1:ii]
            tdim = div(prod(dim),prod(dim1[qs1.==1]))
            temp = zeros(tdim,tdim)

            for jj = 1:1:dim[ii]
                dimt = dim[1:(ii-1)]
                qst = qs[1:(ii-1)]
                tdim2 = prod(dimt[qst.==0])
                tdim3 = prod(dim[(ii+1):end])

                unitV(tempV,jj)
                tempOp = kron(Matrix{ComplexF64}(I,tdim2,tdim2),kron(tempV,Matrix{ComplexF64}(I,tdim3,tdim3)))
                temp = temp + tempOp'*rhoPT*tempOp
            end
            rhoPT = temp
        end
    end
    return rhoPT
end

"""In place version of the `partial_trace` function. Note the explicit type required for out."""
function partial_trace(out::Array{ComplexF64,2},rho::Array{T,2},qs::Array{Int64,1},dim::Array{Int64,1},num::Int64) where {T <: Number}

    rhoPT = rho
    for ii = 1:1:num
        tempV = zeros(dim[ii])

        if qs[ii] == 1
            dim1 = dim[1:ii]
            qs1 = qs[1:ii]
            tdim = div(prod(dim),prod(dim1[qs1.==1]))
            temp = zeros(tdim,tdim)

            for jj = 1:1:dim[ii]
                dimt = dim[1:(ii-1)]
                qst = qs[1:(ii-1)]
                tdim2 = prod(dimt[qst.==0])
                tdim3 = prod(dim[(ii+1):end])

                unitV(tempV,jj)
                tempOp = kron(Matrix{ComplexF64}(I,tdim2,tdim2),kron(tempV,Matrix{ComplexF64}(I,tdim3,tdim3)))
                temp = temp + tempOp'*rhoPT*tempOp
            end
            rhoPT = temp
        end
    end
    out[:,:] = rhoPT
    nothing
end

"""
    unitV(tempV, el)

This function creates a unit vector
Input arguments:
    tempV: vector of zeros that will be overwritten by the unit vector
    el: location of the element to set to 1.
"""
function unitV(tempV::Array{Float64,1},el::Int64)
    fill!(tempV,0.)
    tempV[el] = 1.
    nothing
end

"""
    full_process(out, ρ, χ, I)

This function calculates the output state for a given input state and set of two-qubit quantum processes
Input:
    out: output state
    rho: input state
    Choi: [Nprocesses, dim^2, dim^2]cell array of Choi states describing the processes
    iden: identity the same size as rho
"""
function full_process(out::Array{ComplexF64,2},rho::Array{ComplexF64,2},Choi::Array{ComplexF64,3},iden::Array{ComplexF64,2})

    out[:,:] = rho
    for ii = 1:1:size(Choi,1)
        out[:,:] = partial_trace(kron(transpose(out),iden)*Choi[ii,:,:],[1,0],[size(rho,1),size(rho,1)],2).*size(rho,1)
    end
end

"""
    iden(dim)

This function creates an identity matrix of complex floats
Input:
    dim: dimension of the identity matrix
Output:
    identity matrix of complex floats of dimension dim
"""
function iden(dim::Int)
    return Matrix{ComplexF64}(I,dim,dim)
end

"""
    bases_create!(Bases, N, Npairs, pair_order)

This function creates a two-qubit operator basis dilated to N qubits, for each pair of qubits in the full set
Input:
    Bases: [Npairs,16,2^Nqubits,2^Nqubits] size array that stores the bases
    N: number of qubits
    Npairs: number of qubit pairs
    pair_order: ordering of the qubit pairs, qubit number counting starts at 1
Output:
    nothing (modifies the elements of Bases)
"""
function bases_create!(Bases::Array{T,4},N::Int,Npairs::Int,pair_order::Array{Int,2}) where {T <: Number}

    Bases[:,:,:,:] = convert(Array{ComplexF64,4},Bases)

    tempV1 = zeros(2)
    tempV2 = zeros(2)

    for nn = 1:1:Npairs
        for aa = 1:1:2
            for bb = 1:1:2
                for cc = 1:1:2
                    for dd = 1:1:2

                        unitV(tempV1,aa)
                        unitV(tempV2,bb)
                        tempab = kron(iden(2^(pair_order[nn,1]-1)),kron(tempV1,kron(iden(2^(pair_order[nn,2]-pair_order[nn,1]-1)),kron(tempV2,iden(2^(N-pair_order[nn,2]))))))

                        unitV(tempV1,cc)
                        unitV(tempV2,dd)
                        tempcd = kron(iden(2^(pair_order[nn,1]-1)),kron(tempV1,kron(iden(2^(pair_order[nn,2]-pair_order[nn,1]-1)),kron(tempV2,iden(2^(N-pair_order[nn,2]))))));

                        Bases[nn,4*((2*(aa-1)+bb)-1)+2*(cc-1)+dd,:,:] = tempab*tempcd'
                    end
                end
            end
        end
    end
    nothing
end

"""
    max_ent_state!(dim, mvec, mrho)

This function creates a maximally entangled state. It can return either a state vector or density matrix (or both). Note the syntax used!
Input:
    dim: Hilbert space dimension of the system
Keyword optional input:
    mvec: empty vector to store the maximally entangled state vector
    mrho: empty matrix to store the maximally entangled density matrix
Output:
    nothing
"""
function max_ent_state!(dim::Int;mvec::Vector{Float64}=zeros(dim*dim),mrho::Array{Float64,2}=zeros(dim*dim,dim*dim))

    tempV = zeros(dim)
    for mm = 1:1:dim
        unitV(tempV,mm)
        mvec[:] = mvec[:] + kron(tempV,tempV)
    end
    mvec[:] = mvec[:]./sqrt(dim)

    mrho[:,:] = mvec*mvec'
    nothing
end
