module NumericallyIntegrateArrays
using OffsetArrays

# package code goes here
export simps,trapz
# export sphericalIntegral,cylindricalIntegral

############################################################################################
## helper functions
############################################################################################

function check_last_value(a::AbstractArray,val=2π;var="ϕ",switch="ϕ2π")
	if a[end]≈val
		warn("$(var) appears to have $(val) as the last element, "*
			"but the keyword argument $(switch) is set to false.\n"*
			"To fix this, add $(switch)=true to your function call")
	end
end

scalar_to_0dim_array(z::Number)::Array{<:Number,0} = setindex!(Array{eltype(z)}(),z)
scalar_to_0dim_array(T::Type{<:Number},z::Number)::Array{<:Number,0} = setindex!(Array{T}(),z)

##############################################################################################

function leading_trailing_axes(y::AbstractArray{<:Number},dim::Integer)
	ax = axes(y)
	leading_axes = CartesianIndices(ax[1:dim-1])
	trailing_axes = CartesianIndices(ax[dim+1:end])
	inds_type = NTuple{ndims(y)-1,<:AbstractUnitRange}
	inds = (leading_axes.indices...,trailing_axes.indices...) :: inds_type
	return leading_axes,trailing_axes,inds
end

include("trapz.jl")
include("simps.jl")

#####################################################################################

end # module
