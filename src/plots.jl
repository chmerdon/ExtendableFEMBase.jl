
import UnicodePlots: lines!, points!, pixel!

function unicode_gridplot(
        xgrid;
        title = "gridplot",
        resolution = (40,20),
        color = (200,200,200),
        bface_color = (255,0,0),
        kwargs...
        )
    coords = xgrid[Coordinates]
    ex = extrema(view(coords,1,:))
    ey = extrema(view(coords,2,:))
    canvas = BrailleCanvas(resolution[2], resolution[1],            # number of rows and columns (characters)
                       origin_y=ey[1], origin_x=ex[1],              # position in virtual space
                       height=ey[2], width=ex[2]; blend = false, kwargs...)    # size of the virtual space
    
    bfacenodes = xgrid[BFaceNodes]
    bfaceregions = xgrid[BFaceRegions]
    ebfr = extrema(bfaceregions)
    facenodes = xgrid[FaceNodes]

    nfaces = size(facenodes,2)
    for j in 1 : nfaces
       lines!(canvas, coords[1,facenodes[1,j]], coords[2,facenodes[1,j]], coords[1,facenodes[2,j]], coords[2,facenodes[2,j]]; color = color)
    end

    nbfaces = size(bfacenodes,2)
    for j in 1 : nbfaces
        cscale = (bfaceregions[j] - ebfr[1])/(ebfr[2] - ebfr[1])
        c = Int.(round.(bface_color .* cscale))
        lines!(canvas, coords[1,bfacenodes[1,j]], coords[2,bfacenodes[1,j]], coords[1,bfacenodes[2,j]], coords[2,bfacenodes[2,j]]; color = c)                 # pixel space
    end
    plot = Plot(canvas; title = title, kwargs...)
    return plot
end


function unicode_scalarplot(
        u::FEVectorBlock; 
        component = 1,
        resolution = (40,40),
        colormap = :inferno,
        title = u.name,
        kwargs...)
        
    xgrid = u.FES.xgrid
    coords = xgrid[Coordinates]
    ex = extrema(view(coords,1,:))
    ey = extrema(view(coords,2,:))

	X = LinRange(0,1,resolution[1])
	Y = LinRange(0,1,resolution[2])
	xgrid_plot = simplexgrid(X,Y)

    I = FEVector(FESpace{H1P1{1}}(xgrid_plot))
    lazy_interpolate!(I[1], [u], [(1, IdentityComponent{component})])
    
    return heatmap(reshape(I.entries, (resolution[1], resolution[2]))', xfact=(ex[2]-ex[1])/(resolution[1]-1), yfact=(ey[2]-ey[1])/(resolution[2]-1), border = border, canvas = canvas, xoffset=ex[1], yoffset=ey[1], title = title, colormap=colormap, kwargs...)
end