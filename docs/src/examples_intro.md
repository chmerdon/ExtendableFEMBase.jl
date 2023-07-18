# About the examples

The examples have been designed with the following issues in mind:
- they run from the Julia REPL
- each example is a Julia module named similar to the basename of the example file.
- an example can be used as the starting point for a project 
- some examples define test cases for the test suite
- ExampleXYZ with X = A can be considered advanced and uses low-level structures
  and/or demonstrates customisation features or experimental features
- the default output of the main function is printed on the website and can be
  used to check if the code runs as expected (unfortunately REPL messages are not recorded)
- printed assembly and solving times (especially in a first iteration) can be much larger due to first-run compilation times


## Running the examples

In order to run `ExampleXXX`, peform the following steps:

- Download the example file (e.g. via the source code link at the top)
- Make sure all used packages are installed in your Julia environment
- In the REPL: 
```
julia> include("ExampleXXX.jl")`

julia> ExampleXXX.main()
```
- Some examples offer visual output via the optional argument Plotter = PyPlot or Plotter = GLMakie
(provided the package PyPlot/GLMakie is installed and loaded)