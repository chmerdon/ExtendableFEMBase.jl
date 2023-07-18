
# Meshing

Meshes are stored as an ExtendableGrid, see [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl) for details and constructors.
Grid generators for simplex grids can be found e.g. in the external module [SimplexGridFactory.jl](https://github.com/j-fu/SimplexGridFactory.jl)

Cells, faces and edges of the mesh are associated to AbstractElementGeometries (defined by [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl)) that are used to dispatch functionality (local/global transformation, enumeration rules, set of basis functions, volume calculation, refinements etc.). See further below for a list of recognized element geometries.


## Recognized Geometries and Reference Domains

The following list contains all subtypes of ExtendableGrids.AbstractElementGeometries and their reference domains for which the package offers finite elements on them.

----
##### Edge1D <: AbstractElementGeometry1D

    [1]-----[2]               [1] = [0]
                              [2] = [1]
                              
----
##### Triangle2D
    
    [3]                 
     | \   
     |   \                    [1] = [0,0]
     |     \                  [2] = [1,0]
     |       \                [3] = [0,1]
     |         \ 
    [1]--------[2]
            
----
##### Parallelogram2D <: Quadrilateral2D

    [4]--------[3]               
     |          |             [1] = [0,0]
     |          |             [2] = [1,0]
     |          |             [3] = [1,1]
     |          |             [4] = [0,1]
    [1]--------[2]

    Note: most finite elements only work as intended on Parallelogram2D
          since the local<>global map stays affine in this case


----
##### Tetrahedron3D

    [4]                 
     |\\   
     | \ \                    [1] = [0,0,0]
     |  \  \                  [2] = [1,0,0]
     |   \   \                [3] = [0,1,0]
     | _-[3]-_ \              [4] = [0,0,1]
    [1]--------[2]


----
##### Parallelepiped3D <: Hexahedron3D
                         
        [8]--------[7]        [1] = [0,0,0]
       / |        / |         [2] = [1,0,0]
    [5]--------[6]  |         [3] = [1,1,0]
     |   |      |   |         [4] = [0,1,0]
     |   |      |   |         [5] = [0,0,1]
     |  [4]-----|--[3]        [6] = [1,0,1]
     | /        | /           [7] = [1,1,1]
    [1]--------[2]            [8] = [0,1,1]

    Note: most finite elements only work as intended on Parallelepiped3D
          since the local<>global map stays affine in this case
