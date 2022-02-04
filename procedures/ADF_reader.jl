#!/usr/bin/env julia

#f = "/home/adam/work/shared_ruthenium_stratocaster/frq_benzene.1094055.out"

function ADF_reader(filepath::String)
	#= Reads ADF output file with frequency calculation.
		Returns	- xyz coordinates as N×3 Matrix{Float64}
				- atom types as Array{String}
				- frequencies as Array{String}
				- normal modes as 3*N×3*N Matrix{Float64} =#
	open(filepath) do file
		# Skip to the definition of atomic coordinates
        while true
			line = readline(file)
			if startswith(line, " atom")
				break
			end
		end
		readline(file)
		xyzs = []
		atoms = []
		# Read xyz coordinates
		while true
			line = split(readline(file))
			if length(line) == 0
                break
            end
			atom = line[1]
			append!(atoms, atom)
			xyz = parse.(Float64, line[6:8])
			append!(xyzs, [xyz])	
		end
		xyzs = vcat(map(x->x', xyzs)...)
		atoms = string.(atoms)
		natoms = length(atoms)
		dof = 3*natoms
		# Skip to frequency calculation section in output
		while true
			line = readline(file)
			if startswith(line, " Vibrations and Normal Modes")
				break
			end
		end
		for i in 1:4
			readline(file)
		end
		# Read normal modes
		C = Matrix{Float64}(undef, dof, dof)
		freqs = Array{Float64}(undef, dof)
		i_block = 1
		while i_block < fld(dof-6, 3) + 1
			readline(file)
			readline(file)
			line = readline(file)
			frqs_block = parse.(Float64, split(line))
			for (j, f) in enumerate(frqs_block)
				freqs[ 6 + j + (i_block-1)*3 ] = f
			end
			readline(file)
			for i in 1:natoms
				line = split(readline(file))
				cij = parse.(Float64, line[2:end])
				cij = reshape(cij, (3,3))'
				for (a_count, atom_disp) in enumerate(eachrow(cij))
					for (j, c) in enumerate(atom_disp)
						C[ 6 + (i_block-1)*3 + a_count, (i-1)*3 + j ] = c
					end
				end
			end
 			i_block += 1
		end
		return (xyzs, atoms, freqs, C)
	end
end

#xyzs, atoms, freqs, C = ADF_reader(f)
