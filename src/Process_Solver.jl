"""
    setup_problem(N, Npairs, pair_order)

Create the list of objects necessary for the PAPA reconstruction.
"""
function setup_problem(N::Int,Npairs::Int,pair_order::Array{Int,2})

    # Dimensions and parameter numbers
    dim = 2^N
    Choi_dim = dim^2
    dim2q = 16
    Nel = convert(Int,dim2q*(dim2q+1)/2)

    # Two-qubit bases
    Bases = zeros(Npairs,dim2q,dim,dim)
    B2 = zeros(1,16,4,4)

    bases_create!(Bases,N,Npairs,pair_order)
    Bases = convert(Array{ComplexF64,4},Bases)

    p_temp = zeros(Int,1,2)
    p_temp[1,:] = [1 2]
    # p_temp = convert(Array{Int,2},p_temp)
    bases_create!(B2,2,1,p_temp)
    B2 = convert(Array{ComplexF64,4},B2)

    # Maximally entangled states
    mrho = zeros(dim^2,dim^2)
    mrho2 = zeros(16,16)
    max_ent_state!(dim;mrho=mrho)
    max_ent_state!(4,mrho=mrho2)

    # Identities
    IN = iden(dim)
    I4 = iden(4)

    # Choi and chi matrices
    Choi = zeros(ComplexF64,Npairs,Choi_dim,Choi_dim)
    Choi2q = zeros(ComplexF64,dim2q,dim2q)
    chitemp = zeros(ComplexF64,dim2q,dim2q)
    rhotemp = zeros(ComplexF64,dim,dim)
    tempE = zeros(ComplexF64,16)
    tempF = zeros(ComplexF64,Nel)

    # Counters and qubit labels
    start = 0
    stop = 0
    qstemp = ones(Int,N)

    return (Bases,B2,N,Npairs,pair_order,mrho,mrho2,IN,I4,Choi_dim,dim2q,Nel,Choi,Choi2q,chitemp,rhotemp,tempE,tempF,start,stop,qstemp)
end

"""
Functionalized version of the full process solver
"""
function process_solver(chilist::Vector{Float64},sigmaVec::Vector{Float64},setup_tuple::Tuple,CP_penalty::Float64)
    return process_solver_full(chilist,sigmaVec,CP_penalty,setup_tuple...)
end

function process_solver!(F::Vector{Float64},chilist::Vector{Float64},sigmaVec::Vector{Float64},setup_tuple::Tuple,CP_penalty::Float64)
    F[:] = process_solver_full(chilist,sigmaVec,CP_penalty,setup_tuple...)
    nothing
end


"""
Full process solver that implements the PAPA_BEAR bootstrapping
"""
function process_solver_full(chilist::Vector{Float64},sigmaVec::Vector{Float64},CP_penalty::Float64,Bases::Array{ComplexF64,4},B2::Array{ComplexF64,4},N::Int,Npairs::Int,pair_order::Array{Int,2},mrho::Array{Float64,2},mrho2::Array{Float64,2},IN::Array{ComplexF64,2},I4::Array{ComplexF64,2},Choi_dim::Int,dim2q::Int,Nel::Int,Choi::Array{ComplexF64,3},Choi2q::Array{ComplexF64,2},chitemp::Array{ComplexF64,2},rhotemp::Array{ComplexF64,2},tempE::Array{ComplexF64,1},tempF::Array{ComplexF64,1},start::Int,stop::Int,qstemp::Vector{Int})::Vector{Float64}

    F = []

    fill!(Choi,0.)

    for nn = 1:1:Npairs
        fill!(chitemp,0.)
        start = 1 + (nn-1)*Nel
        stop = nn*Nel;
        chitemp[triu(trues(dim2q,dim2q))] = chilist[start:1:stop] + 1im*chilist[(start+Npairs*Nel):1:(stop+Npairs*Nel)]
        chitemp[:,:] = chitemp+chitemp'
        chitemp[:,:] = chitemp - diagm(0 => diag(chitemp))/2

        fill!(Choi2q,0.)
        for aa = 1:1:dim2q
            for bb = 1:1:dim2q
                # Choi state for the 2-qubit process on all N qubits
                Choi[nn,:,:] = Choi[nn,:,:] + chitemp[aa,bb]*kron(IN,Bases[nn,aa,:,:])*mrho*kron(IN,Bases[nn,bb,:,:]')

                # Choi state for the 2-qubit process on the 2-qubits only (for constraints use)
                Choi2q[:,:] = Choi2q[:,:] + chitemp[aa,bb]*kron(I4,B2[1,aa,:,:])*mrho2*kron(I4,B2[1,bb,:,:]')
            end
        end

        tempE = eigvals(Choi2q)

        F = vcat(F,[real(4*partial_trace(Choi2q,[0,1],[4,4],2)[:]-I4[:]);imag(4*partial_trace(Choi2q,[0,1],[4,4],2)[:]-I4[:]);CP_penalty*abs(sum(tempE[tempE.<0]))])

    end

    for pp = 1:1:Npairs
        fill!(chitemp,0.)  #ignore the name, used as a temp array

        fill!(qstemp,1)
        qstemp[pair_order[pp,1]] = 0
        qstemp[pair_order[pp,2]] = 0

        # Calculate the reduced process
        for dd = 1:1:dim2q
            fill!(rhotemp,0.) #ignore the name, used as a temp array
            full_process(rhotemp,Bases[pp,dd,:,:],Choi,IN)
            chitemp[:,:] = chitemp[:,:] + kron(B2[1,dd,:,:],partial_trace(rhotemp,qstemp,2*ones(Int,N),N));
        end
        chitemp[:,:] = chitemp[:,:]/2^N

        # Compare to input data
        start = 1 + (pp-1)*Nel
        stop = pp*Nel
        tempF[:] = chitemp[triu(trues(dim2q,dim2q))]

        F = vcat(F,[real(tempF)-sigmaVec[start:1:stop];imag(tempF)-sigmaVec[(start+Npairs*Nel):1:(stop+Npairs*Nel)]])
    end

    return abs.(F)

end
