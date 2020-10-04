include("Graph.jl")

using Random
using SparseArrays

function powerLaw(n; alp = 2.5, xmin = 1)
    Random.seed!(round(Int, time() * 10000))
    x = rand(n)
    for i = 1 : n
        x[i] = xmin * ((1 - x[i])^(-1.0/(alp - 1.0)))
    end
    xm = x[argmax(x)]
    x ./= xm
    return x
end

function Uniform(n)
    Random.seed!(round(Int, time() * 10000))
    x = rand(n)
    return x
end

function Exponential(n; lmd = 1, xmin = 1)
    Random.seed!(round(Int, time() * 10000))
    x = rand(n)
    for i = 1 : n
        x[i] = xmin - (1.0/lmd)*log(1-x[i])
    end
    xm = x[argmax(x)]
    x ./= xm
    return x
end

function getW(G)
	L = zeros(G.n, G.n)
	for (ID, u, v, w) in G.E
		L[u, u] += w
		L[v, v] += w
		L[u, v] -= w
		L[v, u] -= w
	end
	for i = 1 : G.n
		L[i, i] += 1.0
	end
	return inv(L)
end

function getL(G)
	L = zeros(G.n, G.n)
	for (ID, u, v, w) in G.E
		L[u, u] += w
		L[v, v] += w
		L[u, v] -= w
		L[v, u] -= w
	end
	return L
end

function getSparseIpL(G)
    d = ones(G.n)
    for (ID, u, v, w) in G.E
        d[u] += w
        d[v] += w
    end
    Is = zeros(Int32, G.m*2+G.n)
    Js = zeros(Int32, G.m*2+G.n)
    Vs = zeros(G.m*2+G.n)
    for (ID, u, v, w) in G.E
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = -w
        Is[ID + G.m] = v
        Js[ID + G.m] = u
        Vs[ID + G.m] = -w
    end
    for i = 1 : G.n
        Is[G.m + G.m + i] = i
        Js[G.m + G.m + i] = i
        Vs[G.m + G.m + i] = d[i]
    end
    return sparse(Is, Js, Vs, G.n, G.n)
end

function getSparseL(G)
    d = zeros(G.n)
    for (ID, u, v, w) in G.E
        d[u] += w
        d[v] += w
    end
    Is = zeros(Int32, G.m*2+G.n)
    Js = zeros(Int32, G.m*2+G.n)
    Vs = zeros(G.m*2+G.n)
    for (ID, u, v, w) in G.E
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = -w
        Is[ID + G.m] = v
        Js[ID + G.m] = u
        Vs[ID + G.m] = -w
    end
    for i = 1 : G.n
        Is[G.m + G.m + i] = i
        Js[G.m + G.m + i] = i
        Vs[G.m + G.m + i] = d[i]
    end
    return sparse(Is, Js, Vs, G.n, G.n)
end
