module PairDistribution

using StatsBase
using Interpolations

export rdf, rdf1d, rdf2p, rdfinter2d

"""
    distpbc(v1, v2, box)

Calculate a distance vector wrt accoring to periodic boundary conditions, defined in `box`.
"""
function distpbc(v1, v2, box)
    rij = v2 - v1
    for (k, l) in enumerate(box)
        hl = l/2
        if rij[k] > hl
            rij[k] = rij[k] - l
        elseif rij[k] < -hl
            rij[k] = rij[k] + l
        end
    end
    return rij
end

"""
    rdfhist(func, data::Array{T,3}, bins::Int...) -> Array{N,Int}

Calculate `N`D histogrm with `bins` bins using a function `func` for evaluating
values of an empirical distribution.
"""
function rdfhist(func, data::Array{T,3}, bins::Int...; debug=0,
                 minmax=fill((zero(T), maximum(data)/2),length(bins)) ) where T
    pnum, _, fnum = size(data)
    edges = tuple([range(mn, mx, length=b+1) for ((mn,mx),b) in zip(minmax,bins)]...)
    hist = Histogram(edges)
    Δs = map(step, hist.edges)
    for t in 1:fnum # t = 1
        debug>0 && println("Frame $t")
        for i in 1:pnum # i = 1
            ri = @view data[i,2:end,t]
            for j in i+1:pnum # j = 2
                rj = @view data[j,2:end,t]
                for v in func((Δs, ri, rj))
                    push!(hist, v)
                end
            end
        end
    end
    hist
end

"""
    rdf(func, data::Array{T,3}, box::Vector{T}, bins::Int...; kwargs...)

Calculate a `N`D radial distribution function (RDF) using a function `func` for constructing
an empirical distribution for RDF.
"""
function rdf(func, data::Array{T,3}, box::Vector{T}, bins::Int...;
            volfunc = (up,lw)->(up^3 - lw^3), kwargs...) where T
    hist = rdfhist(func, data, bins...; kwargs...)
    Δs = map(step, hist.edges)
    pnum, _, fnum = size(data)
    N = fnum * pnum
    ρ = pnum / prod(box) #  Density of the particles
    c = 4 / 3 * π * ρ

    idxs = CartesianIndices(tuple((1:i for i in bins)...))
    gr = zeros(T, bins)
    for ci in idxs, i in 1:length(ci)
        if i == 1
            lw = (ci.I[i] - 1)*Δs[i]
            up = lw + Δs[i]
            dVρ = c * volfunc(up, lw)
            gr[ci] = hist.weights[ci] / N / dVρ
        end
    end
    rs = [collect(r)[1:end-1].+dr for (r,dr) in zip(hist.edges, Δs)]
    (gr,rs)
end

"""
    rdf1d(data::Array{T,3}, box::Vector{T}, bins::Int...; kwargs...)

Calculate a 1D radial distribution function (RDF) using distance between particles.
"""
function rdf1d(data::Array{T,3}, box::Vector{T}, bins::Int...; kwargs...) where T
    rdf(data, box, bins...; kwargs...) do (_, ri, rj)
        rij = PairDistribution.distpbc(rj, ri, box)
        r2 = sum(abs2, rij)
        [sqrt(r2), sqrt(r2)]
    end
end

"""
    rdf2p(data::Array{T,3}, box::Vector{T}, bins::Int...; kwargs...)

Calculate a 2D radial distribution function (RDF) using distance between particles
and azimuthal angle.
"""
function rdf2p(data::Array{T,3}, box::Vector{T}, bins::Int...; kwargs...) where T
    rdf(data, box, bins...; kwargs...) do (_, ri, rj)
        rij = PairDistribution.distpbc(rj, ri, box)
        r2 = sum(abs2, rij)
        r = sqrt(r2)
        ϕ = atan(rij[2],rij[3]) # azimuthal angle
        ϕs = atan(-rij[2],rij[3]) # azimuthal angle
        [(r,ϕ),(r,ϕs)]
    end
end

"""
    rdfinter2d(gr, rs; bins=30, maxval=0.1, dist=(0.0,Inf))

Construct interpolated 2D RDF with dimensions `bins`x`bins`.

While constructing interpolation filter all value larger then `maxval`,
and which are not in the distance range `dist`.
"""
function rdfinter2d(gr, rs; bins=30, maxval=0.1, dist=(0.0,Inf))
    z = hcat(copy(gr), gr[:,1])
    # filter distances
    minidx = findfirst(d->d>=dist[1],rs[1])
    minidx = minidx === nothing ? 1 : minidx
    maxidx = findfirst(d->d>=dist[2],rs[1])
    maxidx = maxidx === nothing ? length(rs[1]) : maxidx
    dd = rs[1][minidx:maxidx]
    # convert angles to [0,2π] range
    rr = [0; rs[2].+π]
    itp = interpolate((dd, rr), z[minidx:maxidx,:], Gridded(Linear()))
    xx = range(extrema(rr)..., length=bins) |> collect
    yy = range(extrema(dd)..., length=bins) |> collect
    zz = [
        let v = itp(y,x)
            v > maxval ? 0.0 : v
        end
        for y in yy, x in xx]
    xx, yy, zz # dim(yy) x dim(xx)
end

end
