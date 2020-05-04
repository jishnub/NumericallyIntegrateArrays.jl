using Test
using NumericallyIntegrateArrays

@testset "simps" begin
    x = range(0,π,length=10); dx = step(x)

    @testset "dx" begin
	    @testset "1D" begin
		    f = sin.(x)
		    @test isapprox(simps(f,step(x)),2,rtol=1e-3)
	    end
	    @testset "2D" begin
	        f = sin.(x) * cos.(x)'
	       	@test isapprox(simps(f,step(x)),2cos.(x),rtol=1e-3) 
	    end
    end
    @testset "x array" begin
	    @testset "1D" begin
		    f = sin.(x)
		    @test isapprox(simps(f,x),2,rtol=1e-3)
	    end
	    @testset "2D" begin
	        f = sin.(x) * cos.(x)'
	       	@test isapprox(simps(f,x),2cos.(x),rtol=1e-3) 
	    end
    end
end

@testset "trapz" begin
    x = range(0,π,length=10); dx = step(x)

    @testset "dx" begin
	    @testset "1D" begin
		    f = sin.(x)
		    @test isapprox(trapz(f,step(x)),2,rtol=1e-1)
	    end
	    @testset "2D" begin
	        f = sin.(x) * cos.(x)'
	       	@test isapprox(trapz(f,step(x)),2cos.(x),rtol=1e-1) 
	    end
    end
    @testset "x array" begin
	    @testset "1D" begin
		    f = sin.(x)
		    @test isapprox(trapz(f,x),2,rtol=1e-1)
	    end
	    @testset "2D" begin
	        f = sin.(x) * cos.(x)'
	       	@test isapprox(trapz(f,x),2cos.(x),rtol=1e-1) 
	    end
    end
end