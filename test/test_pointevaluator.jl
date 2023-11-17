

function run_pointevaluator_tests()
    @testset "PointEvaluator" begin
        println("\n")
        println("======================")
        println("Testing PointEvaluator")
        println("======================")
        test_pointevaluation()
    end
end

function test_pointevaluation()
    ## check if a linear function is interpolated correctly
    ## and can be evaluated correctly within the cell
    xgrid = grid_triangle([-1.0 0.0; 1.0 0.0; 0.0 1.0]')
    FES = FESpace{H1P1{1}}(xgrid)
    Iu = FEVector(FES)
    interpolate!(Iu[1], (result,qpinfo) -> (result[1] = qpinfo.x[1] + qpinfo.x[2];))
    PE = PointEvaluator([(1, Identity)])
    initialize!(PE, Iu)
    CF = CellFinder(xgrid)
    eval = zeros(Float64,1)
    x = [0.15234,0.2234] # point inside the cell
    xref = [0.0,0.0]
    cell = gFindLocal!(xref, CF, x; icellstart = 1)
    evaluate_bary!(eval, PE, xref, cell)
    @test abs(eval[1] - sum(x)) < 1e-15

    evaluate!(eval, PE, x)
    @test abs(eval[1] - sum(x)) < 1e-15
end