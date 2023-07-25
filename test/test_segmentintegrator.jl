function run_segmentintegrator_tests()
    
    @testset "SegmentIntegrator" begin
        println("\n")
        println("=========================")
        println("Testing SegmentIntegrator")
        println("=========================")
        
        ## initial grid
        xgrid = grid_unitsquare(Triangle2D)
        
        ## Taylor--Hood FESpace
        FES = FESpace{H1P2{2,2}}(xgrid)
        
        ## Hagen-Poiseuille flow
        function u(result, qpinfo)
            x = qpinfo.x
            result[1] = x[2]*(1.0-x[2])
            result[2] = 0.0
        end
        
        ## interpolate
        uh = FEVector(FES)
        interpolate!(uh[1], u)
        
        ## init segment integrator
        SI = SegmentIntegrator(Edge1D, [(1, Identity)])
        initialize!(SI, uh)

        ## integrate along line [1/4,1/4] to [3/4,1/4] in first triangle
        ## exact integral should be [3//32,0]
        result = zeros(Float64, 2)
        world = Array{Array{Float64,1},1}([[1//4,1//4], [3//4,1//4]])
        bary = Array{Array{Float64,1},1}([[1//4,1//2], [3//4, 1//2]])
        integrate_segment!(result, SI, world, bary, 1)
        error1 = sqrt((result[1] - 3//32)^2 + result[2]^2)
        
        ## integrate along line [1/2, 0] to [1/2, 1/2]
        ## exact integral should be [1//12, 0]
        world = Array{Array{Float64,1},1}([[1//2,0], [1//2, 1//2]])
        bary = Array{Array{Float64,1},1}([[1//2,0], [1//2, 1//1]])
        integrate_segment!(result, SI, world, bary, 1)
        error2 = sqrt((result[1] - 1//12)^2 + result[2]^2)
        
        @test max(error1,error2) ≈ 0
    end
end