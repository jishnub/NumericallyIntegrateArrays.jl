###########################################################################################

# Regular grid for 1D array
function simps(y::AbstractVector{<:Number},dx::Real=1;even="avg")

	if iseven(size(y,1))
		if !(even ∈ ("avg","first","last"))
			error("even has to be one of avg, first or last."*
				"Specified even is $(even)")
		end

		T = promote_type(eltype(y), Float64)
		result = zero(T)
		val = zero(T)
		first_ind = first(axes(y,1))
		last_ind = last(axes(y,1))

		if even ∈ ["avg","first"]
			result += simpsRegular1D(view(y,first_ind:last_ind-1),dx)
			val += trapzRegular1D(view(y,last_ind-1:last_ind),dx)
		end
		if even ∈ ["avg","last"]
			result += simpsRegular1D(view(y,first_ind+1:last_ind),dx)
			val += trapzRegular1D(view(y,first_ind:first_ind+1),dx)
		end
		if even == "avg"
			val /= 2
			result /=2
		end
		result += val
	else
		result = simpsRegular1D(y,dx)
	end
	return result 
end

# Regular grid for nD array
function simps(y::AbstractArray{<:Number},dx::Real=1;even="avg",axis::Integer=1)
	
	T = promote_type(eltype(y), Float64)

	if axis != 1
		perm = [axis]
		append!(perm,[i for i in 1:axis-1])
		append!(perm,[i for i in axis+1:ndims(y)])
		y = permutedims(y,perm)
	end

	trailing_axes,inds = trailing_indices(y)
	result = zeros(T,inds)

	if iseven(size(y,1))
		if !(even ∈ ("avg","first","last"))
			error("even has to be one of avg, first or last."*
				"Specified even is $(even)")
		end

		val = zeros(T,inds)
		first_ind = first(axes(y,1))
		last_ind = last(axes(y,1))

		if even ∈ ["avg","first"]
			simpsRegular1DAdd!(view(y,first_ind:last_ind-1,trailing_axes),
						result,dx)
			trapzRegular1DAdd!(view(y,last_ind-1:last_ind,trailing_axes),
						val,dx)
		end
		if even ∈ ["avg","last"]
			simpsRegular1DAdd!(view(y,first_ind+1:last_ind,trailing_axes),
						result,dx)
			trapzRegular1DAdd!(view(y,first_ind:first_ind+1,trailing_axes),
						val,dx)
		end
		if even == "avg"
			val ./= 2
			result ./=2
		end
		result .+= val
	else
		simpsRegular1D!(y,result,dx)
	end
	return result 
end

# Irregular grid for 1D array
function simps(y::AbstractVector{<:Number},x::AbstractVector{<:Real};even="avg")

	if iseven(size(y,1))
		if !(even ∈ ("avg","first","last"))
			error("even has to be one of avg, first or last."*
				"Specified even is $(even)")
		end

		T = promote_type(eltype(y), Float64)
		result = zero(T)
		val = zero(T)
		first_ind = first(axes(y,1))
		last_ind = last(axes(y,1))
		first_ind_x = first(axes(x,1))
		last_ind_x = last(axes(x,1))

		if even ∈ ["avg","first"]
			result += simpsIrregular1D(view(y,first_ind:last_ind-1),
				view(x,first_ind_x:last_ind_x-1))
			val += trapzRegular1D(view(y,last_ind-1:last_ind),
				x[last_ind_x]-x[last_ind_x-1])
		end
		if even ∈ ["avg","last"]
			result += simpsIrregular1D(view(y,first_ind+1:last_ind),
				view(x,first_ind_x+1:last_ind_x))
			val += trapzRegular1D(view(y,first_ind:first_ind+1),
				x[first_ind_x+1]-x[first_ind_x])
		end
		if even == "avg"
			val /= 2.
			result /=2.
		end
		result += val
	else
		result = simpsIrregular1D(y,x)
	end
	return result
end

# Irregular grid for nD array
function simps(y::AbstractArray{<:Number},x::AbstractVector{<:Real};even="avg",axis::Integer=1)

	T = promote_type(eltype(y), Float64)

	if axis != 1
		perm = [axis]
		append!(perm,[i for i in 1:axis-1])
		append!(perm,[i for i in axis+1:ndims(y)])
		y = permutedims(y,perm)
	end

	trailing_axes,inds = trailing_indices(y)
	result = zeros(T,inds)

	if iseven(size(y,1))
		if !(even ∈ ("avg","first","last"))
			error("even has to be one of avg, first or last."*
				"Specified even is $(even)")
		end

		val = zeros(T,inds)

		first_ind = first(axes(y,1))
		last_ind = last(axes(y,1))

		first_ind_x = first(axes(x,1))
		last_ind_x = last(axes(x,1))

		if even ∈ ["avg","first"]
			simpsIrregular1DAdd!(view(y,first_ind:last_ind-1,trailing_axes),
				result,view(x,first_ind_x:last_ind_x-1))

			trapzRegular1DAdd!(view(y,last_ind-1:last_ind,trailing_axes),
				val,x[last_ind_x]-x[last_ind_x-1])
		end
		if even ∈ ["avg","last"]
			simpsIrregular1DAdd!(view(y,first_ind+1:last_ind,trailing_axes),
				result,view(x,first_ind_x+1:last_ind_x))

			trapzRegular1DAdd!(view(y,first_ind:first_ind+1,trailing_axes),
				val,x[first_ind_x+1]-x[first_ind_x])
		end
		if even == "avg"
			val ./= 2.
			result ./=2.
		end
		result .+= val
	else
		simpsIrregular1D!(y,result,x)
	end
	return result
end

#############################################################################

function simpsweight(Nelems,index)
	@assert(1<=index<=Nelems,"index not in range")
	if index==1 || index==Nelems 
		return 1
	elseif iseven(index)
		return 4
	end
	return 2
end

function quadratic_fit_integrate!(coeffs,y,x)
	# Find the lagrange polynomial fit to y(x) (quadratic)
	# integrate this between first(x) and last(x)
	@assert(length(y)==length(x)==length(coeffs)==3)

	coeffs[1] = -(x[1]-x[3])*(2x[1]-3x[2]+x[3])/(6*(x[1]-x[2]))

	coeffs[2] = -(x[1]-x[3])^3/(6*(x[1]-x[2])*(x[2]-x[3]))

	coeffs[3] = (x[1]-x[3])*(x[1]-3x[2]+2x[3])/(6*(x[2]-x[3]))

	dot(coeffs,y)
end

#############################################################################
## Simpson's 1/3 rule
#############################################################################

# 1D array with uniform grid spacing for an odd number of points
function simpsRegular1D(y::AbstractVector{<:Number},dx::Real=1)
	
	@assert(isodd(length(y)),"Number of elements must be odd to apply Simpson's rule")
	T = promote_type(eltype(y),Float64)
	int_y = zero(T)

	N = size(y,1)

	for (ind,y_i) in enumerate(y)
		w = simpsweight(N,ind)
		int_y += w * y_i
	end
	
	return int_y*dx/3
end

function simpsRegular1D!(y::AbstractVector{<:Number},
	int_y::AbstractArray{T},dx::Real=1) where {T<:Number}
	
	@assert(isodd(length(y)),"Number of elements must be odd to apply Simpson's rule")
	
	N = size(y,1)

	fill!(int_y,zero(T))

	for (ind,y_i) in enumerate(y)
		w = simpsweight(N,ind)
		@. int_y += w * y_i
	end

	@. int_y *= dx/3

	return first(int_y)
end

# UnitRanges can be integrated exactly
function simpsRegular1D(y::AbstractUnitRange{<:Real},dx::Real=1)
	(last(y)^2 - first(y)^2) * dx/2
end

# N dimensional array, integrate along the first axis
function simpsRegular1D(y::AbstractArray{<:Number},dx::Real=1)

	T = promote_type(eltype(y),Float64)
	trailing_axes = trailing_indices(y)
	int_y = zeros(T,inds)

	for ind_t in trailing_axes
		int_y[ind_t] = simpsRegular1D(view(y,:,ind_t),dx)
	end
	return int_y
end

function simpsRegular1D!(y::AbstractArray{<:Number},
	int_y::AbstractArray{<:Number},dx::Real=1)

	trailing_axes,_ = trailing_indices(y)

	for ind_t in trailing_axes
		int_y[ind_t] = simpsRegular1D(view(y,:,ind_t),dx)
	end

	return int_y
end

function simpsRegular1DAdd!(y::AbstractArray{<:Number},
	int_y::AbstractArray{<:Number},dx::Real=1)

	trailing_axes,_ = trailing_indices(y)

	for ind_t in trailing_axes
		int_y[ind_t] += simpsRegular1D(view(y,:,ind_t),dx)
	end

	return int_y
end

#############################################################################
## Simpson's 1/3 rule for irregular grids
#############################################################################

# 1D array with non-uniform grid spacing, integral along one axis
function simpsIrregular1D(y::AbstractVector{<:Number},x::AbstractVector{<:Real})

	@assert(isodd(length(y)),"The array should have an odd number of elements")
	@assert(length(x)==length(y),"y and x need to have the same number of elements")

	T = promote_type(eltype(y), Float64)
	int_y = zero(T)

	y_no_offset = OffsetArrays.no_offset_view(y)
	x_no_offset = OffsetArrays.no_offset_view(x)

	# The strategy is to fit a quadratic to three points, and compute its integral
	coeffs = zeros(T,3)

	for ind in 1:2:length(x)-2

		int_y += quadratic_fit_integrate!(coeffs,view(y_no_offset,ind:ind+2),
						view(x_no_offset,ind:ind+2))
	end

	return int_y
end

function simpsIrregular1D!(y::AbstractVector{<:Number},
	int_y::AbstractArray{T},x::AbstractVector{<:Real}) where {T<:Number}

	@assert(isodd(length(y)),"The array should have an odd number of elements")
	@assert(length(x)==length(y),"y and x need to have the same number of elements")

	fill!(int_y,zero(T))

	y_no_offset = OffsetArrays.no_offset_view(y)
	x_no_offset = OffsetArrays.no_offset_view(x)

	# The strategy is to fit a quadratic to three points, and compute its integral
	coeffs = zeros(T,3)

	for ind in 1:2:length(x)-2

		int_y .+= quadratic_fit_integrate!(coeffs,view(y_no_offset,ind:ind+2),
						view(x_no_offset,ind:ind+2))
	end

	return first(int_y)
end

# N dimensional array, integrate along the first axis
function simpsIrregular1D(y::AbstractArray{<:Number},x::AbstractVector{<:Real})

	T = promote_type(eltype(y),Float64)
	trailing_axes,inds = trailing_indices(y)
	int_y = zeros(T,inds)

	for ind_t in trailing_axes
		int_y[ind_t] = simpsIrregular1D(view(y,:,ind_t),x)
	end

	return int_y
end

function simpsIrregular1D!(y::AbstractArray{<:Number},
	int_y::AbstractArray{T},x::AbstractVector{<:Real}) where {T<:Number}

	trailing_axes,_ = trailing_indices(y)

	for ind_t in trailing_axes
		int_y[ind_t] = simpsIrregular1D(view(y,:,ind_t),x)
	end

	return int_y
end

function simpsIrregular1DAdd!(y::AbstractArray{<:Number},
	int_y::AbstractArray{T},x::AbstractVector{<:Real}) where {T<:Number}

	trailing_axes,_ = trailing_indices(y)

	for ind_t in trailing_axes
		int_y[ind_t] += simpsIrregular1D(view(y,:,ind_t),x)
	end

	return int_y
end

function simpsCircle(y::AbstractVector{<:Number},
	ϕ::AbstractArray{<:Real,1}=LinRange(0,2π,size(y,1));
	even::String="avg",ϕ2π::Bool=true)::Array{<:Number,0}

	if ϕ2π
		z = simps(y,x=ϕ,even=even)
	else

		info("Assuming that 2π is not included,"*
			" will append it to the array before integrating.\n"*
			"If 2π is included set ϕ2π=true in the function call")

		check_last_value(ϕ,2π) # check if the last value of ϕ is 2π

		push!(y,y[1])

		z = simps(y,x=push!(ϕ,2π),even=even)

	end

	return z 
end

function simpsSpherical2D(y::AbstractArray{<:Number},
	θ::AbstractArray{<:Real,1}=LinRange(0,π,size(y,1)),
	ϕ::AbstractArray{<:Real,1}=LinRange(0,2π,size(y,2));
	even="avg",θint="clenshaw_quadrature",GLweights=nothing,
	ϕ2π=true)::Array{<:Number,0}

	if θint=="gauss_quadrature"
		
		if GLweights != nothing
			z = gauss_quadrature(y,GLweights)
		end

	elseif θint=="clenshaw_quadrature"
		
		z = clenshaw_curtis_quadrature(y)
		
	else
		z = simps(sin.(θ).*y,x=θ,even=even)
	end

	simpsCircle(z,ϕ,even=even,ϕ2π=ϕ2π) 
end

function simpsSpherical3D(y::AbstractArray{<:Number},
	r::Union{Real,AbstractArray{<:Real,1}}=1,
	θ::AbstractArray{<:Real,1}=LinRange(0,π,size(y,2)),
	ϕ::AbstractArray{<:Real,1}=LinRange(0,π,size(y,3));
	even="avg",θint="clenshaw_quadrature",GLweights = nothing,
	ϕ2π=true)
	
	if typeof(r) <: AbstractArray
		z = simps(r.^2 .* y,x=r,even=even)
	else
		dr = r
		z = dr^2*simps(@. ((1:size(y,1))^2 * y),dx=dr,even=even)
	end

	simpsSpherical2D(z,θ,ϕ,θint=θint,
		GLweights=GLweights,even=even,ϕ2π=ϕ2π)
end