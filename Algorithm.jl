include("Graph.jl")
include("Tools.jl")

using SparseArrays
using Laplacians

function Exact(G, s) # Returns ci, d, p, idc
	T = time()
	L = getSparseL(G)
	W = getW(G)
	avgs = sum(s)/G.n
	s2 = zeros(G.n)
	for i = 1 : G.n
		s2[i] = s[i] - avgs
	end
	z = W * s
	z2 = W * s2
	# calculate C_I(G)
	lz = L * z
	ci = lz' * lz
	# calculate D(G)
	d = z2' * L * z2
	# calculate P(G)
	p = z2' * z2
	# calculate C(G)
	c = z' * z
	# calculate I_dc(G)
	idc = d + c
	# END CALCULATION
	T = time() - T
	return T, ci, d, p, idc
end

function Approx(G, s; eps = 1e-6) # Returns aci, ad, ap, aidc
	T = time()
	IpL = getSparseIpL(G)
	sL = getSparseL(G)
	f = approxchol_sddm(IpL, tol=0.1*eps)
	avgs = sum(s)/G.n
	s2 = zeros(G.n)
	for i = 1 : G.n
		s2[i] = s[i] - avgs
	end
	z = f(s)
	z2 = f(s2)

	# calculate aC_I(G)
	alz = sL * z
	aci = alz' * alz
	# calculate aD(G)
	ad = z2' * sL * z2
	# calculate aP(G)
	ap = z2' * z2
	# calculate aC(G)
	ac = z' * z
	# calculate aI_dc(G)
	aidc = ad + ac
	# END CALCULATION
	T = time() - T
	return T, aci, ad, ap, aidc
end

function doExp(G, s, lg)
	T, ci, d, p, idc = Exact(G, s)
	T2, aci, ad, ap, aidc = Approx(G, s)
	println(lg, "Exact  Time : ", T)
	println(lg, "Approx Time : ", T2)
	println(lg, "ERROR of ci : ", abs(aci-ci)/ci)
	println(lg, "ERROR of d : ", abs(ad-d)/d)
	println(lg, "ERROR of p  : ", abs(ap-p)/p)
	println(lg, "ERROR of idc : ", abs(aidc-idc)/idc)
	println(lg)
end

function doLarge(G, s, lg)
	T2, aci, ad, ap, aidc = Approx(G, s)
	println(lg, "Approx Time : ", T2)
	println(lg)
end
