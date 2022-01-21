#!/usr/bin/env julia


function orcaHess_reader(filename::String)
    #= Reads .hess file from ORCA calculation.
        Returns matrix C with normal modes as rows (including translation and rotation).
    =#
    open(filename) do file
        while true
            line = readline(file)
            if startswith(line , "\$vibrational_frequencies")
                break
            end
        end
        dof = parse(Int, readline(file))
        freqs = Array{Float64}(undef, dof)
        for i in 1:dof
            line =  split(readline(file))
            freq = parse(Float64, line[2])
            freqs[i] = freq
        end
        while true
            line = readline(file)
            if startswith(line, "\$normal_modes")
                break
            end
        end
        readline(file)
        #dof = parse(Int, split(readline(file))[1])
        C = Matrix{Float64}(undef, dof, dof)
        i_block = 1
        while i_block <= fld(dof, 5) + 1
            readline(file)
            i_dof = 1
            while i_dof <= dof
                c_ij = parse.(Float64, split(readline(file))[2:end])
                for (j, cij) in enumerate(c_ij)
                    C[ (i_block-1)*5 + j, i_dof ] = cij
                end
                i_dof += 1
            end
            i_block += 1
        end
        return (freqs, C)
    end
end

