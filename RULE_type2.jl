include("pyp.jl")

###################################################################################
# 以下，具体的な HPYP の定義
###################################################################################

struct RULE_type2 <: PYPModel
    rule1::PL_DicY
    rule ::PL_MapY
    ruleL::PL_MapY
    monoL::PL_NonY
    ruleR::PL_MapY
    monoR::PL_NonY
    uni  ::Uniform
    function RULE_type2(NMAX::Int64)
        new(
        PL_DicY( 0x0000ffff, 0xfffff0000 ), # BC|wA
        PL_MapY( 0x0000ffff, 0x000ff0000 ), # BC|A
        PL_MapY( 0x0000ff00, 0x000ff0000 ), # B*|A
        PL_NonY( 0x0000ff00              ), # B*|*
        PL_MapY( 0x000000ff, 0x000ff0000 ), # *C|A
        PL_NonY( 0x000000ff              ), # *C|*
        Uniform( NMAX ) 
        )
    end
end

@inline function cascade_pq(x::UInt64, y::UInt64, c::RULE_type2)
    xy = y << 2LEN_CHA + x
    p,q = pq(c.rule1, xy)
    #@show p,q
    p,q = pq(c.rule, xy, p, q)
    #@show p,q
    p1,q1 = pq(c.ruleL, xy)
    p1,q1 = pq(c.monoL, x, p1, q1)
    p1,q1 = pq(c.uni, p1, q1)
    p2,q2 = pq(c.ruleR, xy)
    p2,q2 = pq(c.monoR, x, p2, q2)
    p2,q2 = pq(c.uni, p2, q2)
    #@show p1,q1, p2,q2
    p + q * (p1+q1)*(p2+q2), (p1+q1), (p2+q2)
end

function cascade_add(x::UInt64, y::UInt64, c::RULE_type2)
    #@show "type2"
    xy = y << 2LEN_CHA + x
    p_all, p1, p2 = cascade_pq(x,y,c)
    water = p_all * rand()
    water <  ((p, q) = add(c.rule1, xy, water,     ))[1] && return
    water <  ((p, q) = add(c.rule , xy, water, p, q))[1] && return
    water = p1 * rand()
    water < ((p, q) = add(c.ruleL, xy, water,     ))[1] ||
    water < ((p, q) = add(c.monoL, xy, water, p, q))[1]
    water = p2 * rand()
    water < ((p, q) = add(c.ruleR, xy, water,     ))[1] ||
    water < ((p, q) = add(c.monoR, xy, water, p, q))[1]
end

function cascade_del(x::UInt64, y::UInt64, c::RULE_type2)
    xy = y << 2LEN_CHA + x
    del(c.rule1, xy) && return
    del(c.rule,  xy) && return
    del(c.ruleL, xy) ||
    del(c.monoL, xy)
    del(c.ruleR, xy) ||
    del(c.monoR, xy)
end


###################################################################################
# cutoffが必要な場合
###################################################################################
# スライスサンプリングの場合の枝刈りつき p,q:
# HPYP全体の確率の上限: p+q*1 = p+q, 下限 p+ q*0 = p
# p < p_slice < p+q でなければ　p を返す (p のみで判断できる)

@inline function cascade_pq_with_cutoff(x::UInt64, y::UInt64, p_slice::Float64, p_start::Float64, c::RULE_type2)
    xy = y << 2LEN_CHA + x
    p,q = 0.0, p_start
    #@show 1, xy, p_slice, p, q
    ((p,q) = pq_with_cutoff(c.rule1, xy, p_slice, p, q))[2] == 0.0 && return (p,q)
    #@show 2, xy, p_slice, p, q
    ((p,q) = pq_with_cutoff(c.rule , xy, p_slice, p, q))[2] == 0.0 && return (p,q)
    #@show 3, xy, p_slice, p, q
    p1,q1 = pq(c.ruleL, xy)
    p1,q1 = pq(c.monoL, x, p1, q1)
    p1,q1 = pq(c.uni, p1, q1)
    p2,q2 = pq(c.ruleR, xy)
    p2,q2 = pq(c.monoR, x, p2, q2)
    p2,q2 = pq(c.uni, p2, q2)
    #@show 4, xy, p_slice, p, q * (p1+q1)*(p2+q2)
    p, q * (p1+q1)*(p2+q2)
end
