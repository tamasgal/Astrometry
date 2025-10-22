####    Constants    ####

const KMAX = 10
const PXMIN = 1e-7
const VMAX  = 0.5
const IMAX  = 100

const lr_1998_n  = SMatrix{60, 4}(vcat([l.n' for l in lr_1998]...)...)
const lr_1998_a  = SMatrix{60, 2}(vcat([l.a' for l in lr_1998]...)...)
const b_1998_n   = SMatrix{60, 4}(vcat([b.n' for b in b_1998]...)...)
const b_1998_a   = SMatrix{60, 1}(vcat([b.a' for b in b_1998]...)...)

#  2000A constants
const ϕ0_2000_equinox = SMatrix{33, 8}(vcat([t.n' for t in iau_2000_equinox_0_series]...)...)
const a0_2000_equinox = SMatrix{33, 8}(vcat([t.a' for t in iau_2000_equinox_0_series]...)...)
const ϕ1_2000_equinox = SMatrix{ 1, 8}(vcat([t.n' for t in iau_2000_equinox_1_series]...)...)
const a1_2000_equinox = SMatrix{ 1, 8}(vcat([t.a' for t in iau_2000_equinox_1_series]...)...)

const ln_2000A_nutation = SMatrix{678, 5}(vcat([t.n' for t in iau_2000A_nutation_lunisolar_series]...)...)
const la_2000A_nutation = SMatrix{678, 6}(vcat([t.a' for t in iau_2000A_nutation_lunisolar_series]...)...)
const pn_2000A_nutation = SMatrix{687,13}(vcat([t.n' for t in iau_2000A_nutation_planetary_series]...)...)
const pa_2000A_nutation = SMatrix{687, 4}(vcat([t.a' for t in iau_2000A_nutation_planetary_series]...)...)

const ϕ0_2000As = SMatrix{33, 8}(vcat([t.n' for t in s0_2000A]...)...)
const a0_2000As = SMatrix{33, 2}(vcat([t.a' for t in s0_2000A]...)...)
const ϕ1_2000As = SMatrix{ 3, 8}(vcat([t.n' for t in s1_2000A]...)...)
const a1_2000As = SMatrix{ 3, 2}(vcat([t.a' for t in s1_2000A]...)...)
const ϕ2_2000As = SMatrix{25, 8}(vcat([t.n' for t in s2_2000A]...)...)
const a2_2000As = SMatrix{25, 2}(vcat([t.a' for t in s2_2000A]...)...)
const ϕ3_2000As = SMatrix{ 4, 8}(vcat([t.n' for t in s3_2000A]...)...)
const a3_2000As = SMatrix{ 4, 2}(vcat([t.a' for t in s3_2000A]...)...)
const ϕ4_2000As = SMatrix{ 1, 8}(vcat([t.n' for t in s4_2000A]...)...)
const a4_2000As = SMatrix{ 1, 2}(vcat([t.a' for t in s4_2000A]...)...)

#  2000B constants
const ln_2000B = SMatrix{77, 5}(vcat([t.n' for t in iau_2000B_nutation_lunisolar_series]...)...)
const la_2000B = SMatrix{77, 6}(vcat([t.a' for t in iau_2000B_nutation_lunisolar_series]...)...)

#  2006 constants
const jaxy_2006 = SVector(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1) .+ 1
const jasc_2006 = SVector(0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0) .+ 1
const japt_2006 = SVector(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)

const ϕ0_2006s = SMatrix{33, 8}(vcat([t.n' for t in iau_2006_equinox_0_series]...)...)
const a0_2006s = SMatrix{33, 2}(vcat([t.a' for t in iau_2006_equinox_0_series]...)...)
const ϕ1_2006s = SMatrix{ 3, 8}(vcat([t.n' for t in iau_2006_equinox_1_series]...)...)
const a1_2006s = SMatrix{ 3, 2}(vcat([t.a' for t in iau_2006_equinox_1_series]...)...)
const ϕ2_2006s = SMatrix{25, 8}(vcat([t.n' for t in iau_2006_equinox_2_series]...)...)
const a2_2006s = SMatrix{25, 2}(vcat([t.a' for t in iau_2006_equinox_2_series]...)...)
const ϕ3_2006s = SMatrix{ 4, 8}(vcat([t.n' for t in iau_2006_equinox_3_series]...)...)
const a3_2006s = SMatrix{ 4, 2}(vcat([t.a' for t in iau_2006_equinox_3_series]...)...)
const ϕ4_2006s = SMatrix{ 1, 8}(vcat([t.n' for t in iau_2006_equinox_4_series]...)...)
const a4_2006s = SMatrix{ 1, 2}(vcat([t.a' for t in iau_2006_equinox_4_series]...)...)
