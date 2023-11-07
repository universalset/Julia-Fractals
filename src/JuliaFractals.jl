module JuliaFractals

using Images: colorview, RGB

export complexcircle, complexsegment
export juliaset, mandelbrotset
export findinterestingjulia, findmandelbrotedge
export juliamovie, zoomjuliamovie
export showimg

const quadraticfunc(z::ComplexF64, c::ComplexF64)::ComplexF64 = z^2 + c

# keep these three in sync (xdim:ydim = xdomain:ydomain)
const defaultimagedims::Tuple{Int,Int} = (700,500)
const defaultxdomain::Tuple{Float64,Float64} = (-1.75, 1.75)
const defaultydomain::Tuple{Float64,Float64} = (-1.25, 1.25)

const defaultmandelbrotimagedims::Tuple{Int,Int} = (800,800)
const defaultmandelbrotxdomain::Tuple{Float64,Float64} = (-2.0, 2.0)
const defaultmandelbrotydomain::Tuple{Float64,Float64} = (-2.0, 2.0)

# good for higher max iteration counts, like the default 800
const defaultcolors::Vector{Vector{Float64}} = [
    [0.0,0.0,0.0], # black
    [1.0, 0.4, 0.4], # Red
    [0.4, 1.0, 0.4], # Blue
    [0.4, 0.7, 0.7], # Teal
    [0.4, 0.55, 0.85], #Teal/Green
    [0.4,0.4,1.0] # Green
]

const defaultmaxiters::Int = 800 # can be slow, but gives best results

const defaultmandelbrotmaxiters::Int = 400

const defaultescapelimit::Float64 = 10.0 # keep higher than max absolute value of f(z) in domain


function pixeltocomplex(
    a::Int,
    b::Int,
    minx::Float64,
    maxx::Float64,
    miny::Float64,
    maxy::Float64,
    scalex::Int,
    scaley::Int,
)::ComplexF64
    """Find the complex number corresponding to a given pixel in an image.
    
    a,b: x and y pixel coordinates in the image, with (0,0) at top left.
    {min|max}{x|y}: the min/max values of the real/imaginary part for the targen region on the complex plane.
    scalex, scaley: the number of pixels in each row / columm of the image.

    returns: a Complex{Float64} corresponding to the given pixel (a,b).
    """
    minx + (maxx-minx)*a/scalex + (maxy + (miny-maxy)*b/scaley)*im
end


function interpolatecolors(val::Float64, colors::Vector{Vector{Float64}})::Vector{Float64}
    """Interpolates a color between adjacent colors in a given array.

    val: a floating point between 0 and 1, representing a distance along a color range.
    colors: a vector of triples (each a vector of Floats in the range 0-1), each the RGB value of a color.
            The first value is the color for 0.0, the last is the color for 1.0, and the returns
            are evenly spaced in between.
    
    returns: a Vector of three Floats, giving the color at the given point in RGB. The color in assigned
            by interpolation between the neighboring colors in the given color array.
    """
    scaledval = val*(length(colors)-1)
    major = Int64(floor(scaledval))
    frac = scaledval-major
    if major+1 == length(colors)
        colors[major+1]
    else
        colors[major+1] + (colors[major+2] - colors[major+1])*frac
    end
end


function complexcircle(
    center::ComplexF64,
    radius::Float64,
    points::Int,
)::Vector{ComplexF64}
    """Defines a sequence of equidistant points around a complex circle.
    Useful for making looping animations.

    center: the center of the circle
    radius: the radius of the circle
    points: the number of points to generate

    returns: a vector of `points` complex numbers equally spaced around the circumference of the
        given circle in the complex plane.
    """

    [center + radius * exp(2*pi*i*im/points) for i in 1:points]
end


function complexsegment(
    vstart::ComplexF64,
    vend::ComplexF64,
    points::Int,
)::Vector{ComplexF64}
    """Defines a sequence of equidistant points along a complex segment.
    Useful for making animations.

    vstart, vend: the endpoints of the segment.
    points: the number of points to generate

    returns: a vector of `points` complex numbers equally spaced along the
        given segment in the complex plane.
    """
    [vstart + (vend-vstart)*i/(points-1) for i in 0:(points-1)]
end


function juliaset(
    c::ComplexF64;
    func=quadraticfunc, 
    dims::Tuple{Int,Int}=defaultimagedims, 
    xdomain::Tuple{Float64,Float64}=defaultxdomain, 
    ydomain::Tuple{Float64,Float64}=defaultydomain, 
    colors::Vector{Vector{Float64}}=defaultcolors,
    maxiters::Int=defaultmaxiters, 
    escapelimit::Float64=defaultescapelimit, 
    smooth=true,
)
    """Generates an array of colors corresponding to a colorization of the complement of the filled
    Julia set of the given function within the given domain.

    Parameters:
    c: Complex parameter for defining the function.

    Named parameters:
    func: A function of two variables f(z,c).  The function g(z) = f(z,c) is the function we'll be
        finding the Julia set for.  Functions of the form f(z,c) = h(z) + c where h is a polynomial
        work best, and play nicely with other helper functions.
    dims: the (x,y) pixel dimensions of the image we want.
    xdomain, ydomain: the min/max real and imaginary parts (resp) for the rectangular view in the image.
    colors: a vector of RGB triples colors to use (with interpolation) for colorizing the complement
        of the Julia set.  RGB components are in the range [0.0, 1.0].
    maxiters: the maximum number of iterations to apply before assuming that the orbit does not escape
        Larger values yield higher quality images when many orbits barely escape, but run more slowly.
    escapelimit: the value which, if an iterate's norm exceeds it, that iteration is considered to
        have escaped.  Keep large enough to not have non-escaping orbits be mistakenly labeled, but
        otherwise shouldn't have much effect.  Usually the default is fine.
    smooth: should we smooth colors between the values assigned to consecutive numbers of iterations
        based on the norm of the first value in the orbit that exceeds escapelimit?
        Generally should leave this on, as turning it off results in unsightly color banding in many
        images.
    
    Notes:
    - Keep dims and xdomain / ydomain in the same ratios, or your image will be stretched in one dimension
    - changing the value of maxiters may require changing the colors array, because colors are assigned
        based on the ratio (number of iterations to escape) / maxiters and so changing maxiters
        will shift colors.  A good idea is to have colors change more slowly at higher iterations.
        At some point I may implement adaptive coloration logic that takes this into account.
    - values whose orbits don't escape in maxiters iterations are assumed to be in the filled
        Julia set and will be colored black.  For some functions/parameters c, this means that
        changing maxiters will change the points which are colored black.  In such cases larger
        maxiters are usually preferred as the results are generally more correct.

        Usage: call showimg() on the result to generate an image (and show it in a notebook).
    """
    llescape = log(log(escapelimit))
    logorder = log(log(abs(func(func(Complex(escapelimit), c), c))) / log(abs(func(Complex(escapelimit), c))))
    A = zeros(Float64,3,dims[2],dims[1])
    for a in 1:dims[1]
        for b in 1:dims[2]
            z = pixeltocomplex(a,b,xdomain[1],xdomain[2],ydomain[1],ydomain[2],dims[1],dims[2])
            numiters = 0
            while numiters < maxiters && abs(z)<=escapelimit
                z = func(z,c)
                numiters += 1
            end
            if smooth && (0 < numiters) && (numiters < maxiters)
                numiters -= (log(log(abs(z))) - llescape)/logorder
                # in unusual cases, the smoothing computation could result in numiters < 0 or an out of domain error.
                # it's safe to just set numiters = 0 in such cases
                if isnan(numiters) || numiters < 0
                    numiters = 0
                end
            end
            A[:,b,a] = if (numiters >= maxiters)
                [0.0,0.0,0.0]
            else
                interpolatecolors(numiters/maxiters, colors)
            end
        end
    end
    A
end


function mandelbrotset(;
    func=quadraticfunc, 
    dims::Tuple{Int,Int}=defaultmandelbrotimagedims, 
    xdomain::Tuple{Float64,Float64}=defaultmandelbrotxdomain, 
    ydomain::Tuple{Float64,Float64}=defaultmandelbrotydomain, 
    colors::Vector{Vector{Float64}}=defaultcolors,
    maxiters::Int=defaultmandelbrotmaxiters, 
    escapelimit::Float64=defaultescapelimit,
)
    """Creates a colorization of (the complement of) a Mandelbrot set for the given function family.

    Named parameters:
    func: A function of two variables f(z,c).  For each complex c in the domain, the function
        g(z) = f(z,c) will be iterated starting at 0; i.e the forward orbit of 0 is computed.
        For best results, f(z,c) should have the form f(z,c) = h(z) + c where h(z) is a rational
        function with degree >= 2; polynomials tend to be nicest (as for Julia sets). 
    dims: the (x,y) pixel dimensions of the image we want.
    xdomain, ydomain: the min/max real and imaginary parts (resp) for the rectangular view in the image.
    colors: a vector of RGB triples colors to use (with interpolation) for colorizing the complement
        of the Julia set.  RGB components are in the range [0.0, 1.0].
    maxiters: the maximum number of iterations to apply before assuming that the orbit does not escape
        Larger values yield higher quality images when many orbits barely escape, but run more slowly.
    escapelimit: the value which, if an iterate's norm exceeds it, that iteration is considered to
        have escaped.  Keep large enough to not have non-escaping orbits be mistakenly labeled, but
        otherwise shouldn't have much effect.  Usually the default is fine.

    Notes:
    - Keep dims and xdomain / ydomain in the same ratios, or your image will be stretched in one dimension
    - changing the value of maxiters may require changing the colors array, because colors are assigned
        based on the ratio (number of iterations to escape) / maxiters and so changing maxiters
        will shift colors.  A good idea is to have colors change more slowly at higher iterations.
        At some point I may implement adaptive coloration logic that takes this into account.
    - values whose orbits don't escape in maxiters iterations are assumed to be in the Mandelbrot
        set and will be colored black.  For some functions and windows, this means that
        changing maxiters will change the points which are colored black.  In such cases larger
        maxiters are usually preferred as the results are generally more correct.
    - Does not currently support color smoothing as the logic doesn't exactly tronslate,
        but it's on the todo list.

        Usage: call showimg() on the result to generate an image (and show it in a notebook).
    """
    M = zeros(Float64,3,dims[2],dims[1])
    for a in 1:dims[1]
        for b in 1:dims[2]
            c = pixeltocomplex(a,b,xdomain[1],xdomain[2],ydomain[1],ydomain[2],dims[1],dims[2])
            numiters = 0
            z = Complex(0.0)
            while numiters < maxiters && abs(z)<=escapelimit
                z = func(z,c)
                numiters += 1
            end
            M[:,b,a] = if (numiters==maxiters) 
                [0,0,0] 
            else
                interpolatecolors(numiters/maxiters, colors)
            end
        end
    end
    M
end


function findinterestingjulia(;
    func=quadraticfunc,
    xdomain::Tuple{Float64,Float64}=defaultmandelbrotxdomain,
    ydomain::Tuple{Float64,Float64}=defaultmandelbrotydomain,
    iterrange=(40,400),
    escapelimit=defaultescapelimit,
)::ComplexF64
    """Find a random value c likely to produce an 'interesting' looking Julia set for the given function.
    
    Named Parameters:
    func: A function of two variables f(z,c), to be used in a subsequent Julia set generation.
        For best results, f(z,c) should have the form f(z,c) = h(z) + c where h(z) is a rational
        function with degree >= 2.
    xdomain, ydomain: the min/max real and imaginary parts (resp) for the constant c
    iterrange: the value c is considered interesting if the orbit of 0 under g(z) = f(z,c) escapes
        in a number of iterations within this range.
    escapelimit: the value which, if an iterate's norm exceeds it, that iteration is considered to
    have escaped.  Keep large enough to not have non-escaping orbits be mistakenly labeled, but
    otherwise shouldn't have much effect.  Usually the default is fine.

    returns: a complex value c such that juliaset(c;func=func) is likely to produce an interesting image.

    Notes:
    - This function works by choosing random values for c and checking if they meet the criteria.
        Thus changes that decrease the fraction of c in the domain which meet the criteria will
        increase expected runtime.
        - in particular, increasing the size of the domains will usually slow this down
        - narrowing the iterrange will slow it down, but if the domain contains many
            values for which the orbit does not escape, then increasing the upper limit may also
            slow it down because individual choices sometimes take longer to evaluate.
    - If few or no values in the domain meet the criteria, the function may run for an unbounded
        amount of time.
            - TODO: include a failsafe to break after a number of failed iterations.
    - This function works on the principle that the Julia set for h(z) + c is connected
        iff the parameter c is in the Mandelbrot set for h(z).  It tries to pick values for which
        the divergence is slow, so the set is 'barely' disconnected, giving an 'interesting' looking
        colorized complement.
    - Thus you may wish to look at the Mandelbrot set in your chosen domain to see if there are many
        eligible parameters c.

    Usage:
    > c = finditerestingjulia(func=f, params...)
    > img = showimg(juliaset(c, func=f, otherparams...))
    """
    numiters = 0
    c = Complex(0.0)
    while(numiters < iterrange[1] || numiters > iterrange[2])
        c = rand()*(xdomain[2]-xdomain[1])+xdomain[1] + (rand()*(ydomain[2]-ydomain[1])+ydomain[1])*im
        numiters = 0
        z = Complex(0.0)
        while numiters < iterrange[2]+1 && abs(z)<=escapelimit
            z = func(z,c)
            numiters += 1
        end
    end
    c
end


function findmandelbrotedge(;
    func=quadraticfunc,
    inset::ComplexF64=Complex(0.0),
    outset::ComplexF64=Complex(0.0),
    maxiters::Int=defaultmaxiters,
    prec::Float64=0.0000000001,
    escapelimit=defaultescapelimit,
    useinside=true,
)::ComplexF64
    """Find a value of c on a given segment (or random radial ray) close to the edge of the
    Mandelbrot set for the give function.  This should be likely to produce an 'interesting'
    looking Julia set for the given function.  Alternative (often better) to findinterestingjulia.

    Named Parameters:
    func: A function of two variables f(z,c), to be used in a subsequent Julia set generation.
        For best results, f(z,c) should have the form f(z,c) = h(z) + c where h(z) is a rational
        function with degree >= 2.
    inset, outset: two complex numbers, required to have the property that inset is in the
        Mandelbrot set for the given function (orbit does not escape) and outset is not.
        Both default to 0.  If outset is zero, it's chosen randomly in a circle, attempting to be
        outside the Mandelbrot set.
    maxiters: the maximum number of iterations to apply before assuming that the orbit does not escape
    prec: precision with which to calculate c before stopping.
    escapelimit: the value which, if an iterate's norm exceeds it, that iteration is considered to
    have escaped.  Keep large enough to not have non-escaping orbits be mistakenly labeled, but
    otherwise shouldn't have much effect.  Usually the default is fine.
    useinside: whether to try to pick a point just *inside* (if true) or just *outside* (if false)
        the Mandelbrot set.

    returns: a complex value c near the boundary of the Mandelbrot set for func,
        such that juliaset(c;func=func) is likely to produce an interesting image
    
    Notes:
    - This function uses binary search to try to find a boundary point of the Mandelbrot set
        along the given segment.  If the segment intersects the Mandelbrot set in multiple
        places, any intersection may be found.
    - If inset is not inside the Mandelbrot set, or outset is not outside, the point found will
        not be any good.  If you are using an unusual function, consider drawing the Mandelbrot
        set first to choose the segment.
        - TODO: better default/random behavior for nonstandard Mandelbrot sets.
    - Values of prec closer to 0 tend to produce points closer to the boundary.
        Lower maxiters can sometimes pick points further outside the Mandelbrot set even with
        useinside=true.  Try different parameters here.
    - This function tends to produce more saturated / closer to connected sets than
        findinterestingjulia and using maxiters in the juliaset function at least equal to the
        maxiters used here is encouraged.
    
    Usage:
    > c = findmandelbrotedge(func=f, params...)
    > img = showimg(juliaset(c, func=f, otherparams...))
    """
    if outset==Complex(0.0)
        θ = 2π*rand()
        outset = 3*(cos(θ) + sin(θ)*im)
    end
    while abs(outset-inset) > prec
        mid = (outset+inset)/2
        numiters = 0
        z = Complex(0.0)
        while numiters < maxiters && abs(z) <= escapelimit
            z = func(z, mid)
            numiters += 1
        end
        if abs(z) > escapelimit
            outset = mid
        else
            inset = mid
        end
    end
    if useinside inset else outset end
end

function juliamovie(
    path::Vector{ComplexF64};
    func=quadraticfunc, 
    dims::Tuple{Int,Int}=defaultimagedims, 
    xdomain::Tuple{Float64,Float64}=defaultxdomain, 
    ydomain::Tuple{Float64,Float64}=defaultydomain, 
    colors::Vector{Vector{Float64}}=defaultcolors,
    maxiters::Int=defaultmaxiters, 
    escapelimit::Float64=defaultescapelimit, 
    smooth=true,
)
    """Creates an animated image (a sequence of images) of Julia sets based on a sequence of Complex
    parameters.

    Parameters:
    path: an array of complex numbers, corresponding to a sequence of c parameter values for
        drawing the Julia set of frame of the animations.
    
    Named Parameters:
    As in juliaset, passed through.

    returns: an animated image, as a sequence of frames.  Save as a gif to view.

    Notes:
    - Consider using complexcircle to generate a path for a standard looping animation.
        A concatenation of complexsegment results can also serve this purpose, though mind
        duplication of endpoints.
    - This can take a long time when many frames must be generated!
    - Doesn't seem to properly display in a Jupyter notebook.
    - TODO: more automation for path generation; more options for frame rate when saving as GIF.
    
    usage: call save (from FileIO) to save as an animated GIF, e.g:
    > path = complexcircle(center, radius, points)
    > mov = juliamovie(path)
    > save("animation.gif", mov)
"""
    MOVIE = zeros(3, dims[2], dims[1], length(path))
    for i in eachindex(path)
        MOVIE[:,:,:,i] = juliaset(
            path[i],
            func=func,
            dims=dims,
            xdomain=xdomain,
            ydomain=ydomain,
            colors=colors,
            maxiters=maxiters,
            escapelimit=escapelimit,
            smooth=smooth,
        )
    end
    showimg(MOVIE)
end

function zoomjuliamovie(
    c;
    zoompoint::ComplexF64=Complex(0.0),
    zoomsteps::Int=240,
    finalzoomratio::Float64=1000.0,
    func=quadraticfunc, 
    dims::Tuple{Int,Int}=defaultmandelbrotimagedims, 
    xdomain_init::Tuple{Float64,Float64}=defaultmandelbrotxdomain, 
    ydomain_init::Tuple{Float64,Float64}=defaultmandelbrotydomain, 
    colors::Vector{Vector{Float64}}=defaultcolors,
    maxiters::Int=defaultmaxiters, 
    escapelimit::Float64=defaultescapelimit, 
    smooth=true,
)
    """Creates an animated image (a sequence of images) of Julia sets by gradually changing the
    domain to zoom in on a given point at a constant rate.

    Parameters: 
    c: Complex parameter for defining the function, as in juliaset.

    Named Parameters:
    zoompoint: the point to zoom in on
    zoomsteps: the number of frames to generate
    finalzoomratio: the amount of zoom in the last frame, relative to the first.
    xdomain_init, ydomain_init: real/imaginary parts of the inital bounding box.
    All other named parameters (ekcept xdomain, ydomain) as in juliaset, passed through.

    returns: an animated image, as a sequence of frames.  Save as a gif to view.

    Notes:
    - This can take a long time when many frames must be generated!
    - Doesn't seem to properly display in a Jupyter notebook.
    - Zoom animation is often a bit choppy / unpleasant.
        - TODO: look into fixing this in some way.
    - TODO: support arbitrary zoom paths and not just zooming in on a single point.
        - This could also be used for non-zooming, or flexibly-zooming, "tours"
        - Need to supply automated path generation in that case.
    - TODO: create a version of this for mandelbrotset as well.

    Usage:
    usage: call save (from FileIO) to save as an animated GIF, e.g:
    > mov = zoomjuliamovie(c, params...)
    > save("animation.gif", mov)
    """
    MOVIE = zeros(3, dims[2], dims[1], zoomsteps)
    for i in 1:zoomsteps
        zoomlevel = finalzoomratio^((i-1)/zoomsteps)
        mincorner = ((xdomain_init[1] + ydomain_init[1]*im) - zoompoint) / zoomlevel + zoompoint
        maxcorner = ((xdomain_init[2] + ydomain_init[2]*im) - zoompoint) / zoomlevel + zoompoint
        MOVIE[:,:,:,i] = juliaset(c,
            func=func,
            dims=dims,
            xdomain=(real(mincorner), real(maxcorner)),
            ydomain=(imag(mincorner), imag(maxcorner)),
            colors=colors,
            maxiters=maxiters,
            escapelimit=escapelimit,
            smooth=smooth,
        )
    end
    showimg(MOVIE)
end

function showimg(img)
    """Turn an array of colors into an image object, enabling it to be shown in a notebook and
    saved in a variety of formats.
    """
    colorview(RGB, img)
end

end # module JuliaFractals
