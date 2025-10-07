###################################################################################
# PCFGの非終端記号を，type0とtype2に分ける場合．
###################################################################################

include("PCFG.jl")
include("RULE_switch.jl")
include("RULE_type0.jl")
include("RULE_type2.jl")
include("hyper_parameters.jl")

const EX0    = cha(0x80) # 0x80 が type0 非終端記号のための any シンボル
const EX2    = cha(0x00) # 0x00 が type0 非終端記号のための any シンボル

###################################################################################
# 以下，具体的な HPYP の定義
###################################################################################

mutable struct PCFG_sep  <: PCFG
    t0::RULE_type0
    t2::RULE_type2
    V0::Vector{cha}
    V2::Vector{cha}
    NEXT_V0::cha
    NEXT_V2::cha
    function PCFG_sep(Σ::Int64)
        obj = new(
            RULE_type0(Σ),
            RULE_type2(NMAX),
            [cha(EX0+1)],
            [cha(EX2+1)],
            cha(EX0+2),
            cha(EX2+2)
        )
        init_mcmc_hyperparam(obj)
        obj
    end
end

@inline _N(  c::PCFG_sep) = (_N0(c), _N2(c))
@inline _N0(  c::PCFG_sep) = length(c.V0)
@inline _N2(  c::PCFG_sep) = length(c.V2)
@inline _Σ(   c::PCFG_sep) = c.t0.t0_uni.N
@inline _NV0(  c::PCFG_sep) = ( _N0(c), c.V0 )
@inline _NV2(  c::PCFG_sep) = ( _N2(c), c.V2 )
@inline _V(c::PCFG_sep) = [c.V0;c.V2]

#ワイルドカードに対応するための，確率の和のための係数の計算．
@inline function p_mulitply_wildcard(pcfg::PCFG_sep, B::cha, C::cha, p::Float64)
    n0_others = NMAX/2-_N0(pcfg)+1 #+1は wildcardが_N0に含まれるから．
    n2_others = NMAX/2-_N2(pcfg)+1 #+1は wildcardが_N2に含まれるから．
    #@show n0_others n2_others
    B==EX0 && (p*=n0_others)
    B==EX2 && (p*=n2_others)
    C==EX0 && (p*=n0_others)
    C==EX2 && (p*=n2_others)
    p
end

###################################################################################
# ノンターミナルの集合の更新: 
#   カウントが 0 のものをリスト V から削除，
#   新しい id を見つける．
###################################################################################

# V{0,2} は，EX{0,2} から始まり，ソートされているように更新する．
function update_V(c::PCFG_sep)
    V0 = ( get_nonzero_eys(c.t0.t0_rule) .>> 2LEN_CHA )
    V2 = ( get_nonzero_eys(c.t2.rule   ) .>> 2LEN_CHA )
    c.V0 = [EX0;V0]
    c.V2 = [EX2;V2]
    # 未使用の非終端記号番号を見つけておく．
    c.NEXT_V0 = sort_find_least_new(c.V0)
    c.NEXT_V2 = sort_find_least_new(c.V2)
end

###################################################################################
# サンプルされたルール中の指定の列([2],[2,3,4])に，EX があれば NEXT_V に置き換える．
###################################################################################
function replace_EX(c, x)
    if typeof(x)==cha
        if x == EX0
            return c.NEXT_V0
        elseif x == EX2 
            return c.NEXT_V2
        end
    end
    return x
end
    
function assign_EX(c::PCFG_sep, rules_t0::Vector{enc},rules_t2::Vector{enc})
    lambda = 
    for i in eachindex(rules_t0)
#        tup = map( x->(typeof(x)==cha && x==EX0 ? c.NEXT_V0 : x ), decode_t0(rules_t0[i]) )
        tup = map( x -> replace_EX(c, x), decode_t0(rules_t0[i]) )

        rules_t0[i] = join_xy(encode(tup...))
    end
    for i in eachindex(rules_t2)
#        tup = map( x->(typeof(x)==cha && x==EX2 ? c.NEXT_V2 : x ), decode_t2(rules_t2[i]) )
        tup = map( x -> replace_EX(c, x), decode_t2(rules_t2[i]) )
        rules_t2[i] = join_xy(encode(tup...))
    end
end

###################################################################################
# 通常の操作(取得，追加，削除)
###################################################################################

@inline function get_p(pcfg::PCFG_sep, w1::ta, A::cha, B::cha, C::cha)
    p = p_mulitply_wildcard(pcfg, B, C, 1.0)
    x, y = encode(w1,A,B,C)
    p*= cascade_pq( x, y, pcfg.t2)[1]
end
@inline function get_p(pcfg::PCFG_sep, w1::ta, A::cha, u::ta)
    x, y = encode(w1,A,u)
    p=1.0
    p*= cascade_pq( x, y, pcfg.t0)
end

function increment(pcfg::PCFG_sep, w1::ta, A::cha, B::cha, C::cha)
    #@assert A!= EX0
    #@assert B!= EX0
    #@assert C!= EX0
    #@assert A!= EX2
    #@assert B!= EX2
    #@assert C!= EX2
    x, y = encode(w1,A,B,C)
    cascade_add( x, y, pcfg.t2)
end
function increment(pcfg::PCFG_sep, w1::ta, A::cha, u::ta)
    #@assert A!= EX0
    #@assert A!= EX2
    x, y = encode(w1,A,u)
    cascade_add( x, y, pcfg.t0)
end

function decrement(pcfg::PCFG_sep, w1::ta, A::cha, B::cha, C::cha)
    x, y = encode(w1,A,B,C)
    cascade_del( x, y, pcfg.t2)
end
function decrement(pcfg::PCFG_sep, w1::ta, A::cha, u::ta)
    x, y = encode(w1,A,u)
    cascade_del( x, y, pcfg.t0)
end

###################################################################################
# p_slice をしきい値とし，1か0のみ返せば良い場合(スライスサンプリング)
###################################################################################
# スライスサンプリングの場合の枝刈りつき p,q:
# HPYP全体の確率の上限: p+q*1 = p+q, 下限 p+ q*0 = p
# p < p_slice < p+q でなければ　(p, 0.0) を返す (p のみで判断できる)

function slice(pcfg::PCFG_sep, w1::ta, A::cha, B::cha, C::cha, p_slice::Float64 )
    p = p_mulitply_wildcard(pcfg, B, C, 1.0)
    x, y = encode(w1,A,B,C)
    p,q = cascade_pq_with_cutoff(x, y, p_slice,  p, pcfg.t2)
    Int(p_slice < p+q)
end
function slice(pcfg::PCFG_sep, w1::ta, A::cha, u::ta, p_slice::Float64 )
    x, y = encode(w1,A,u)
    p,q = cascade_pq_with_cutoff(x, y, p_slice,  p, pcfg.t0)
    Int(p_slice < p+q)
end


#=
make_get_all_layers(T) = "all_lay(c::$T) = [" *
     join( [ "c."*string(lay) for lay in fieldnames(T)], "," ) * "]"
eval(Meta.parse(make_get_all_layers(PCFG_type0)))
eval(Meta.parse(make_get_all_layers(PCFG_type2)))

function get_all_layers(c::PCFG_sep)
    [all_lay(c.t0);all_lay(c.t2)]
end
=#

function check_probsum(c::PCFG_sep)
    Σ = _Σ(c)
    for w1 in ta(0):ta(Σ)
        for A in c.V2
            p2 = 0.0
            for B in cha(1):cha(255) # c.V #
                B == cha(0x80) && continue
                for C in cha(1):cha(255) # c.V #
                    C == cha(0x80) && continue
                    p2 += get_p(c, w1, A, B, C)
                end
            end
            p2 = round(p2,digits=8)
            @show "type2", w1, A,  p2 # p_pos + p_neg #p0, p2, psum, p3,
        end
        for A in c.V0
            p0 = 0.0
            for u in ta(1):ta(Σ)
                p0 += get_p(c, w1, A, u)
            end
            p0 = round(p0,digits=8)
            @show "type0", w1, A,  p0 # p_pos + p_neg #p0, p2, psum, p3,
        end
    end
end

function check_probsum_V(c::PCFG_sep)
    Σ = _Σ(c)
    for w1 in ta(0):ta(Σ)
        for A in c.V2
            p = 0.0
            for B in [c.V0;c.V2] # c.V #
                for C in [c.V0;c.V2] # c.V #
                    p += get_p(c, w1, A, B, C)
                    x, y = encode( w1, A, B, C)
                end
            end
            p = round(p, digits=8)
            #p_more = round(p_more, digits=8)            
            @show "type2", w1,A,  p # p_pos + p_neg #p0, p2, psum, p3,
            #check_consistency(c.sw)
        end
    end
end



###################################################################################
# 小さい値の加算により生じる丸め誤差について
# 小数点以下 12桁で誤差
###################################################################################

function demostrate_round_error()
    c = PCFG_sep(1)
    x = 0.0
    delta = 1.0*(1/254)*(1/254) #1.5500031000061998e-5 
    for i in 1:(127*2)^2
        x += delta
    end
    @show x, "生じる丸め誤差は", 1.0-x
    check_probsum(c)
    check_consistency(c)
end