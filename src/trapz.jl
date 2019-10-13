####################################################################################
####################################################################################
# Trapezoidal rule
####################################################################################
####################################################################################

# Regular grid for 1D array
trapz(y::AbstractArray,dx::Real=1;kwargs...)  = trapzRegular1D(
												OffsetArrays.no_offset_view(y),dx;
												kwargs...)

# Irregular grid for 1D array
trapz(y::AbstractArray,x::AbstractVector;kwargs...) = trapzIrregular1D(
												OffsetArrays.no_offset_view(y),
												OffsetArrays.no_offset_view(x);kwargs...)


####################################################################################

function trapzweight(Nelems,index)
	@assert(1<=index<=Nelems,"index not in range")
	if index==1 || index==Nelems
		return 1
	end
	return 2
end

######################################################################################

# Regular grid for 1D array
function trapzRegular1D(y::AbstractVector,dx::Real=1)
	N = size(y,1)
	T = promote_type(eltype(y),Float64)
	int_y = zero(T)
	for (ind,y_i) in enumerate(y)
		w = trapzweight(N,ind)
		int_y += w*y_i
	end
	return dx/2*int_y
end

# UnitRanges can be integrated exactly
function trapzRegular1D(y::AbstractUnitRange,dx::Real=1)
	(last(y)^2 - first(y)^2)/2 * dx
end

# N dimensional array, integrate along any axis
function trapzRegular1D(y::AbstractArray,dx::Real=1;axis=1)

	leading_axes = CartesianIndices(axes(y)[1:axis-1])
	trailing_axes = CartesianIndices(axes(y)[axis+1:end])
	T = promote_type(eltype(y),Float64)
	inds_type = NTuple{ndims(y)-1,<:AbstractUnitRange}
	inds = (leading_axes.indices...,trailing_axes.indices...) :: inds_type
	int_y = zeros(T,inds)
	for ind_t in trailing_axes,ind_l in leading_axes
		int_y[ind_l,ind_t] = trapzRegular1D(view(y,ind_l,:,ind_t),dx)
	end
	return int_y

end

####################################################################################

# Irregular grid for 1D array
function trapzIrregular1D(y::AbstractVector,x::AbstractVector)
	first_ind = first(axes(y,1))
	last_ind = last(axes(y,1))

	first_ind_x = first(axes(x,1))
	last_ind_x = last(axes(x,1))

	T = promote_type(eltype(y),Float64)
	int_y = zero(T)

	for (ind_x,ind_y) in zip(first_ind_x:last_ind_x-1,first_ind:last_ind-1)
		dx = x[ind_x+1] - x[ind_x]
		int_y += (y[ind_y]+y[ind_y+1])/2 * dx
	end
	return int_y
end

# N dimensional array, integrate along any axis
function trapzIrregular1D(y::AbstractArray,x::AbstractVector;axis=1)
	
	leading_axes = CartesianIndices(axes(y)[1:axis-1])
	trailing_axes = CartesianIndices(axes(y)[axis+1:end])
	T = promote_type(eltype(y),Float64)
	inds_type = NTuple{ndims(y)-1,<:AbstractUnitRange}
	inds = (leading_axes.indices...,trailing_axes.indices...) :: inds_type
	int_y = zeros(T,inds)
	for ind_t in trailing_axes,ind_l in leading_axes
		int_y[ind_l,ind_t] = trapzIrregular1D(view(y,ind_l,:,ind_t),x) 
	end
	return int_y
end

####################################################################################
# Error in trapz
#####################################################################################

function trapzerror(y::AbstractVector,dx::Number)::Number
	y_2nddiff = diff(diff(y))
	x_range = (length(y)-1)*dx
	return max(eps(),abs(x_range/12*mean(y_2nddiff)))
end

#####################################################################################
# trapz on a sphere
#####################################################################################

function trapzCircle(y::AbstractVector,
	ϕ::AbstractArray{<:Real,1}=linspace(0,2π,size(y,1));
	ϕ2π::Bool=true):: Array{<:Number,0}

	if ϕ2π
		z = trapz(y,x=ϕ)
	else

		info("Assuming that 2π is not included,"*
		" will append it to the array before integrating.\n"*
		"If 2π is included set ϕ2π=true in the function call")
		
		check_last_value(ϕ,2π)

		push!(y,y[1])

		z = trapz(y,x=push!(ϕ,2π))
		
	end

	return z
end

function trapzSpherical2D(y::AbstractArray{<:Number,2},
	θ::AbstractArray{<:Real,1}=linspace(0,π,size(y,1)),
	ϕ::AbstractArray{<:Real,1}=linspace(0,2π,size(y,2));
	θint="clenshaw_quadrature",GLweights = nothing,
	ϕ2π=true)

	if θint=="clenshaw_quadrature"

		z = clenshaw_curtis_quadrature(y) :: Vector{<:Number}
		
	elseif θint=="gauss_quadrature"
		
		if GLweights != nothing
			z = gauss_quadrature(y,GLweights) :: Vector{<:Number}
		end

	else
		
		z = trapz(sin.(θ).*y,x=θ) :: Vector{<:Number}

	end

	trapzCircle(z,ϕ,ϕ2π=ϕ2π) :: Array{<:Number,0}
end

function trapzSpherical3D(y::AbstractArray{<:Number,3},
	r::Union{Real,AbstractArray{<:Real,1}}=1,
	θ::AbstractArray{<:Real,1}=linspace(0,π,size(y,2)),
	ϕ::AbstractArray{<:Real,1}=linspace(0,2π,size(y,3));
	θint="clenshaw_quadrature",GLweights = nothing,
	ϕ2π=true):: Array{<:Number,0}
	
	if typeof(r) <: AbstractArray
		z = trapz(r.^2 .* y,x=r)
	else
		dr =r
		z = dr^2 .* trapz((1:size(y,1)).^2 .* y,dx=dr)
	end

	trapzSpherical2D(z,θ,ϕ,GLweights=GLweights,ϕ2π=ϕ2π,θint=θint)
end