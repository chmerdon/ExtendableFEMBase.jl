using Documenter
using Literate
using ExtendableSparse
using ExtendableGrids
using PlutoSliderServer
using ExtendableFEMBase


# Turn block comments starting in the first column into "normal" hash comments
# as block comments currently are not handled by Literate.jl.
function hashify_block_comments(input)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out=IOBuffer()
    line_number=0
    in_block_comment_region=false
    for line in lines_in
        line_number+=1
        if occursin(r"^#=", line)
            if in_block_comment_region
                error("line $(line_number): already in block comment region\n$(line)")
            end
            println(lines_out,replace(line,r"^#=" => "#"))
            in_block_comment_region=true
        elseif occursin(r"^=#", line)
            if !in_block_comment_region
                error("line $(line_number): not in block comment region\n$(line)")
            end
            println(lines_out,replace(line,r"^=#" => "#"))
            in_block_comment_region=false
        else
            if in_block_comment_region
                println(lines_out,"# "*line)
            else
                println(lines_out,line)
            end
        end
    end
    return String(take!(lines_out))
end

#
# Replace SOURCE_URL marker with github url of source
#
function replace_source_url(input,source_url)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out=IOBuffer()
    for line in lines_in
        println(lines_out,replace(line,"SOURCE_URL" => source_url))
    end
    return String(take!(lines_out))
end

function make_all(; with_examples::Bool = true, run_examples = true, run_notebooks::Bool = true)

    generated_examples = []
    notebooks = []

    if run_notebooks
        #
        # Run notebooks
        #
        notebooks = [
            "Low level Poisson" => "LowLevelPoisson.jl"
            "Low level Navier-Stokes" => "LowLevelNavierStokes.jl"
        ]

        notebookjl = last.(notebooks)
        notebookmd = []

        # function rendernotebook(name)
        #     base=split(name,".")[1]
        #     input=joinpath(@__DIR__,"..","pluto-examples",base*".jl")
        #     output=joinpath(@__DIR__,"src","nbhtml",base*".html")
        #     session = Pluto.ServerSession();
        #     html_contents=PlutoStaticHTML.notebook2html(input;session)
        #     write(output, html_contents)
        # end


        # for notebook in notebookjl
        #     @info "Converting $(notebook)"
        #     rendernotebook(notebook)
        # end

        # Use sliderserver to generate html
        notebook_html_dir = joinpath(@__DIR__, "src", "nbhtml")
        export_directory(
            joinpath(@__DIR__, "..", "pluto-examples"),
            notebook_paths = notebookjl,
            Export_output_dir = joinpath(notebook_html_dir),
            Export_offer_binder = false,
        )

        # generate frame markdown for each notebook
        for notebook in notebookjl
            base = split(notebook, ".")[1]
            mdstring = """
                       ##### [$(base).jl](@id $(base))
                       [Download](https://github.com/chmerdon/ExtendableFEMBase.jl/blob/master/examples/pluto/$(notebook))
                       this [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook.
                       ```@raw html
                       <iframe style="height:20000px" width="100%" src="../$(base).html"> </iframe>
                       ```
                       """
            mdname = base * ".md"
            push!(notebookmd, joinpath("nbhtml", mdname))
            io = open(joinpath(notebook_html_dir, mdname), "w")
            write(io, mdstring)
            close(io)
        end

        @show notebookmd
        notebooks = first.(notebooks) .=> notebookmd
        @show notebooks
        pushfirst!(notebooks, "About the notebooks" => "notebooks_intro.md")
        @show notebooks
    end


    #
    # Generate Markdown pages from examples
    #
    if with_examples
        example_jl_dir = joinpath(@__DIR__,"..","examples")
        example_md_dir  = joinpath(@__DIR__,"src","examples")
        excluded_examples = ["XXX","A05","231","260","401","402"] # excludes just the run of these examples
        image_dir = joinpath(@__DIR__,"src","images")

        for example_source in readdir(example_jl_dir)
            base,ext=splitext(example_source)
            if example_source == "pluto"
                break;
            end
            if example_source[1:7] == "Example" && ext==".jl"
                number = example_source[8:10]
                source_url="https://github.com/chmerdon/ExtendableFEMBase.jl/raw/master/examples/"*example_source
                preprocess(buffer)=replace_source_url(buffer,source_url)|>hashify_block_comments
                    Literate.markdown(joinpath(@__DIR__,"..","examples",example_source),
                                example_md_dir,
                                documenter=false,
                                execute=false,
                                info=false,
                                preprocess=preprocess)

                filename = example_md_dir * "/" * base * ".md"
                if (run_examples) && !(number in excluded_examples) # exclude these examples for now (because they take long or require extra packages)
                    # generate default main run output file 
                    include(example_jl_dir * "/" * example_source)
                    @time open(filename, "a") do io
                        redirect_stdout(io) do
                            println("**Default output:**")
                            println("```")
                            println("julia> $base.main()")
                            eval(Meta.parse("$base.main()"))
                            println("```")
                        end
                    end
                end
                for k = 1 : 4
                    imgfile = "../images/" * base * "_$k.png"
                    if isfile(image_dir * "/" * base * "_$k.png")
                        open(filename, "a") do io
                            redirect_stdout(io) do
                                println("![]($imgfile)")
                            end
                        end
                    end
                end
            end
        end

        generated_examples=joinpath.("examples",readdir(example_md_dir))

        pushfirst!(generated_examples, "examples_intro.md")
    end

    makedocs(
        modules=[ExtendableFEMBase],
        sitename="ExtendableFEMBase.jl",
        authors="Christian Merdon",
        repo = "github.com/chmerdon/ExtendableFEMBase.jl",
        clean = false,
        checkdocs = :all,
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
                    "quadrature.md",
                    "pointevaluators.md",
                    "segmentintegrators.md",
                ],
            "Tutorial Notebooks" => notebooks,
            "Examples" => generated_examples,
        ]
    )

    with_examples && rm(example_md_dir, recursive = true)
    run_notebooks && rm(notebook_html_dir, recursive = true)
    
end

#make_all(; with_examples = true, run_examples = true, run_notebooks = true)
make_all(; with_examples = true, run_examples = false, run_notebooks = false)

deploydocs(
    repo = "github.com/chmerdon/ExtendableFEMBase.jl",
)