include("pyp.jl")

###################################################################################
# 以下，具体的な HPYP の定義
###################################################################################

struct RULE_type0 <: PYPModel
    t0_rule1::PL_DicY
    t0_rule ::PL_MapY
    t0_mono ::PL_NonY
    t0_uni  ::Uniform
    function RULE_type0(Σ::Int64)
        new(
        PL_DicY( 0x0000ffff, 0xfffff0000 ), # a|wA
        PL_MapY( 0x0000ffff, 0x000ff0000 ), # a|A
        PL_NonY( 0x0000ffff              ), # a|*
        Uniform( Σ ) 
        )
    end
end

@inline function cascade_pq(x::UInt64, y::UInt64, c::RULE_type0)
    xy = y << 2LEN_CHA + x
    p,q = pq(c.t0_rule1, xy)
    p,q = pq(c.t0_rule, xy, p, q)
    p,q = pq(c.t0_mono, x, p, q)
    p,q = pq(c.t0_uni, p, q)
    p + q
end

function cascade_add(x::UInt64, y::UInt64, c::RULE_type0)
    #@show "type0 add" x, y
    xy = y << 2LEN_CHA + x
    p_all = cascade_pq(x,y,c)
    water = p_all * rand()
    water < ((p, q) = add(c.t0_rule1, xy, water      ))[1] && return
    water < ((p, q) = add(c.t0_rule , xy, water, p, q))[1] && return
    water < ((p, q) = add(c.t0_mono , xy, water, p, q))[1]
end

function cascade_del(x::UInt64, y::UInt64, c::RULE_type0)
    xy = y << 2LEN_CHA + x
    del(c.t0_rule1, xy) && return
    del(c.t0_rule,  xy) && return
    del(c.t0_mono, xy)
end


###################################################################################
# cutoffが必要な場合
###################################################################################
# スライスサンプリングの場合の枝刈りつき p,q:
# HPYP全体の確率の上限: p+q*1 = p+q, 下限 p+ q*0 = p
# p < p_slice < p+q でなければ　p または p+q を返す
# (p or p+q のみで判断できる)

function cascade_pq_with_cutoff(x::UInt64, y::UInt64, p_slice::Float64, p_start::Float64, c::RULE_type0)
    xy = y << 2LEN_CHA + x
    p,q = 0.0, p_start
    ((p,q) = pq_with_cutoff(c.t0_rule1, xy, p_slice, p, q))[2] == 0.0 && return (p, 0.0)
    ((p,q) = pq_with_cutoff(c.t0_rule , xy, p_slice, p, q))[2] == 0.0 && return (p, 0.0)
    ((p,q) = pq_with_cutoff(c.t0_mono , xy, p_slice, p, q))[2] == 0.0 && return (p, 0.0)
    (p,q) = pq(c.t0_uni, p, q)
    p, q
end

