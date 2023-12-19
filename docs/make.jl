using Documenter
using ExtendableFEMBase
using ExampleJuggler
using CairoMakie


function make_all(; with_examples::Bool = true)

	module_examples = []
    pluto_examples = []

    if with_examples
        
		DocMeta.setdocmeta!(ExampleJuggler, :DocTestSetup, :(using ExampleJuggler); recursive = true)

		example_dir = joinpath(@__DIR__, "..", "examples")
		pluto_example_dir = joinpath(@__DIR__, "..", "pluto-examples")

		modules = [
			"Example200_LowLevelPoisson.jl",
			"Example210_LowLevelNavierStokes.jl",
			"Example280_BasisPlotter.jl",
			"Example281_DiscontinuousPlot.jl",
			"Example290_InterpolationBetweenMeshes.jl",
		]
        
        notebooks = [
            "Low level Poisson" => "LowLevelPoisson.jl"
            "Low level Navier-Stokes" => "LowLevelNavierStokes.jl"
        ]

		cleanexamples()

		module_examples = @docmodules(example_dir, modules, Plotter = CairoMakie)
        pluto_examples = @docplutonotebooks(pluto_example_dir, notebooks, iframe=true)
        pushfirst!(module_examples, "Introduction" => "examples_intro.md")

    end

    makedocs(
        modules=[ExtendableFEMBase],
        sitename="ExtendableFEMBase.jl",
        authors="Christian Merdon",
        repo = "github.com/chmerdon/ExtendableFEMBase.jl",
        clean = false,
        format = Documenter.HTML(size_threshold = 250000),
        checkdocs = :all,
        warnonly = true,
        doctest = true,
        pages = [
            "Home" => "index.md",
            "Index" => "package_index.md",
            "List of Finite Elements" => "fems.md",
            "Base Structures" => Any[
                    "fespace.md",
                    "fevector.md",
                    "fematrix.md",
                    "functionoperators.md",
                    "feevaluator.md",
                    "interpolations.md",
                    "quadrature.md"
                ],
            "Advanced Stuff" => Any[
                "pointevaluators.md",
                "segmentintegrators.md",
                "plots.md"
            ],
            "Tutorial Notebooks" => pluto_examples,
            "Examples" => module_examples,
        ]
    )

	cleanexamples()
    
end

make_all(; with_examples = true)

deploydocs(
    repo = "github.com/chmerdon/ExtendableFEMBase.jl",
)