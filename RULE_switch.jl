include("pyp.jl")

###################################################################################
# 以下，具体的な HPYP の定義
###################################################################################

struct RULE_switch <: PYPModel
    #sw_rule1::PL_DicY
    sw_rule ::PL_MapY
    sw_mono ::PL_NonY
    sw_uni  ::Uniform
    function RULE_switch()
        new(
        #PL_DicY( 0xffffffff, 0xffff0000 ), # {0,1}|wA
        PL_MapY( 0x0000ffff, 0x00ff0000 ), # {0,1}|A
        PL_NonY( 0x0000ffff             ), # {0,1}|*
        Uniform( 2 ) 
        )
    end
end

@inline function cascade_pq(x::UInt64, y::UInt64, c::RULE_switch)
    xy = y << 2LEN_CHA + x
    p,q = 0.0, 1.0
    #p,q = pq(c.sw_rule1, xy)
    p,q = pq(c.sw_rule, xy, p, q)
    p,q = pq(c.sw_mono, x, p, q)
    p,q = pq(c.sw_uni, p, q)
    p + q
end

function cascade_add(x::UInt64, y::UInt64, c::RULE_switch)
    xy = y << 2LEN_CHA + x
    p_all = cascade_pq(x,y,c)
    water = p_all * rand()
    #@show "switch", x, y, p_all
    p,q = 0.0, 1.0
    #water < ((p, q) = add(c.sw_rule1, xy, water,     ))[1] && return
    water < ((p, q) = add(c.sw_rule , xy, water, p, q))[1] && return
    water < ((p, q) = add(c.sw_mono , xy, water, p, q))[1]
end

function cascade_del(x::UInt64, y::UInt64, c::RULE_switch)
    xy = y << 2LEN_CHA + x
    #del(c.sw_rule1, xy) && return
    del(c.sw_rule,  xy) && return
    del(c.sw_mono, xy)
end


'''
###################################################################################
# cutoffが必要な場合
###################################################################################
# スライスサンプリングの場合の枝刈りつき p,q:
# HPYP全体の確率の上限: p+q*1 = p+q, 下限 p+ q*0 = p
# p < p_slice < p+q でなければ　p または p+q を返す
# (p or p+q のみで判断できる)

function cascade_pq_with_cutoff(x::UInt64, y::UInt64, p_slice::Float64, p_start::Float64, c::RULE_switch)
    y <<= 16
    xy = x+y
    p,q = 0.0, p_start
    #((p,q) = pq_with_cutoff(c.sw_rule1, xy, p_slice      ))[2] == 0.0 && return (p, 0.0)
    ((p,q) = pq_with_cutoff(c.sw_rule , xy, p_slice, p, q))[2] == 0.0 && return (p, 0.0)
    ((p,q) = pq_with_cutoff(c.sw_mono , xy, p_slice, p, q))[2] == 0.0 && return (p, 0.0)
    (p,q) = pq(c.uni, p, q)
    p, q
end
'''


