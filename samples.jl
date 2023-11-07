using JuliaFractals
using FileIO: save

println("Generating and saving sample images.")

# basic invocation with defoult settings
samplejuliabasic = showimg(
    juliaset(0.0+.642im)
)

# save the image
save("sample_julia_basic.gif", samplejuliabasic)


# same julia set, different color scheme
purplecolors = [[1.0,1.0,1.0], [0.9,0.5,1.0], [0.85,0.25,1.0], [0.825,0.125,1.0], [0.8,0.0,1.0]]

samplejuliapurple = showimg(
    juliaset(
        0.0+.642im,
        colors=purplecolors,
    )
)

save("sample_julia_purple.gif", samplejuliapurple)

# We can change the default function to get different patterns

cubicfunc(z,c) = z^3 + c

samplejuliacubic = showimg(
    juliaset(
        -0.582819 + 0.314035im,
        func = cubicfunc,
    )
)

save("sample_julia_cubic.gif", samplejuliacubic)

# while we're at it, let's change the size and viewport

samplejuliacubicbig = showimg(
    juliaset(
        -0.582819 + 0.314035im,
        func = cubicfunc,
        dims=(800,800),
        xdomain=(-1.25,1.25),
        ydomain=(-1.25,1.25),
    )
)

save("sample_julia_cubic_big.gif", samplejuliacubicbig)

# note how colors change and black areas increase when reducing maxiters

samplejuliacubiclowiter = showimg(
    juliaset(
        -0.582819 + 0.314035im,
        func = cubicfunc,
        dims=(800,800),
        xdomain=(-1.25,1.25),
        ydomain=(-1.25,1.25),
        maxiters=200,
    )
)

save("sample_julia_cubic_lowiter.gif", samplejuliacubiclowiter)


# we can also generate the Mandelbrot set

mandelbrotbasic = showimg(mandelbrotset())

save("basic_mandelbrot.gif", mandelbrotbasic)

# and Mandelbrot sets for other functions

mandelbrotcubic = showimg(
    mandelbrotset(func=cubicfunc)
)

save("cubic_mandelbrot.gif", mandelbrotcubic)


println("Generating and saving random Julia sets images.")

# We can use either findinterestingjulia or findmandelbrotedge
# to generate c parameters to try for julia sets.

# note: image generated here are random and change each time!

for i in 1:5
    cval = findinterestingjulia()
    randomjulia = showimg(
        juliaset(cval)
    )
    save("random_interesting_julia_" + string(i) + ".gif", randomjulia)
end

for i in 1:5
    cval = findmandelbrotedge()
    randomjulia = showimg(
        juliaset(cval)
    )
    save("random_mandelbrotedge_julia_" + string(i) + ".gif", randomjulia)
end


println("Generating sample animations, this may take a while.")

# We can also generate animations.  This one uses complexcircle.

path = complexcircle(
    -1.022156 + 0.251445im,
    0.01,
    120,
)

juliaanim = juliamovie(path)

# save as .gif and get an animated gif!
save("sample_animated_julia.gif", juliaanim)


# optional: A zoom animation. uncomment before running if you want it
# this one will take quite a while.

#= cval=0.03894257714098054 - 1.1229861597876438im

juliazoomanim = zoomjuliamovie(
    cval,
    func=cubicfunc,
)

save("sample_julia_zoom_animation.gif", juliazoomanim)
=#