using PairDistribution
using Test
using Random
using DelimitedFiles

function readlfdem(fname)
    pnum, vf, lx, ly, lz = open(fname, "r") do io
        readline(io)
        parse(Int, split(readline(io), " ")[end]),
        parse(Float32, split(readline(io), " ")[end]),
        parse(Float32, split(readline(io), " ")[end]),
        parse(Float32, split(readline(io), " ")[end]),
        parse(Float32, split(readline(io), " ")[end])
    end
    box = Float32[lx,ly,lz]
    X = readdlm(fname, ' ', Float32, skipstart=23, comments=true)
    X, box, pnum, vf
end

function frames(fname)
    X, box, pnum, _  = readlfdem(fname)
    # make frames
    tot = size(X,1)
    fnum = Int(tot/pnum) # number of frames
    F = Array{Float32,3}(undef, pnum, 4, fnum)
    for i in 1:fnum
        F[:,:,i] .= X[(i-1)*pnum+1:i*pnum,2:5]
    end
    F, box
end

@testset "In-Box Distnace" begin
    vbox = [2;2]
    @test PairDistribution.distpbc([0;0],[0;1],vbox) == [0;1]
    @test PairDistribution.distpbc([0;-1],[0;1],vbox) == [0;0]
end
