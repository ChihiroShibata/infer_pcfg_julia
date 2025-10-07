###################################################################################
# PCFGの非終端記号を，type0とtype2で共通とする場合．
###################################################################################
include("PCFG.jl")
include("RULE_switch.jl")
include("RULE_type0.jl")
include("RULE_type2.jl")
include("hyper_parameters.jl")



const EX    = cha(0x00) # 0x00 が 非終端記号のための any シンボル


###################################################################################
# 以下，具体的な HPYP の定義
###################################################################################

mutable struct PCFG_all <: PCFG
    sw::RULE_switch
    t0::RULE_type0
    t2::RULE_type2
    V::Vector{cha}
    NEXT_V::cha
    function PCFG_all(Σ::Int64, initial_V=[cha(EX+1)])
        obj = new(
            RULE_switch(),
            RULE_type0(Σ),
            RULE_type2(NMAX),
            initial_V,
            cha(maximum(initial_V)+1)
        )
        init_mcmc_hyperparam(obj)
        obj
    end
end

@inline _N(   c::PCFG_all) = length(c.V)
@inline _Σ(   c::PCFG_all) = c.t0.t0_uni.N

#ワイルドカードに対応するための，確率の和のための係数の計算．
@inline function p_mulitply_wildcard(pcfg::PCFG_all, B::cha, C::cha, p::Float64)
    n_others = NMAX-_N(pcfg)+1 #+1は wildcardが_Nに含まれるから．
    B==EX && (p*=n_others)
    C==EX && (p*=n_others)
    p
end

###################################################################################
# ノンターミナルの集合の更新: 
#   カウントが 0 のものをリスト V から削除，
#   新しい id を見つける．
###################################################################################

# V は，EX から始まり，ソートされているように更新する．
function update_V(pcfg::PCFG_all)
    V = ( get_nonzero_eys(pcfg.sw.sw_rule) .>> 2LEN_CHA )
    pcfg.V = [EX;V]
        # 未使用の非終端記号番号を見つけておく．
    pcfg.NEXT_V = sort_find_least_new(pcfg.V)
end

###################################################################################
# サンプルされたルール中の指定の列([2],[2,3,4])に，EX があれば NEXT_V に置き換える．
###################################################################################
function assign_EX(c::PCFG_all, rules_t0::Vector{enc},rules_t2::Vector{enc})
    for i in eachindex(rules_t0)
        tup = map( x->(typeof(x)==cha && x==EX ? c.NEXT_V : x ), decode_t0(rules_t0[i]) )
        rules_t0[i] = join_xy(encode(tup...))
    end
    for i in eachindex(rules_t2)
        tup = map( x->(typeof(x)==cha && x==EX ? c.NEXT_V : x ), decode_t2(rules_t2[i]) )
        rules_t2[i] = join_xy(encode(tup...))
    end
end

###################################################################################
# 通常の操作(取得，追加，削除)
###################################################################################

@inline function get_p(pcfg::PCFG_all, w1::ta, A::cha, B::cha, C::cha)
    p = p_mulitply_wildcard(pcfg, B, C, 1.0)
    x, y = encode(w1,A,B,C)
    p*= cascade_pq(enc(1), y, pcfg.sw)
    p*= cascade_pq(        x, y, pcfg.t2)[1]
end
@inline function get_p(pcfg::PCFG_all, w1::ta, A::cha, u::ta)
    x, y = encode(w1,A,u)
    p = cascade_pq(enc(0), y, pcfg.sw)
    p*= cascade_pq(        x, y, pcfg.t0)
end

function increment(pcfg::PCFG_all, w1::ta, A::cha, B::cha, C::cha)
    @assert A!= EX
    @assert B!= EX
    @assert C!= EX
    x, y = encode(w1,A,B,C)
    cascade_add(enc(1), y, pcfg.sw)
    cascade_add(        x, y, pcfg.t2)
end
function increment(pcfg::PCFG_all, w1::ta, A::cha, u::ta)
    @assert A!= EX
    x, y = encode(w1,A,u)
    cascade_add(enc(0), y, pcfg.sw)
    cascade_add(        x, y, pcfg.t0)
end

function decrement(pcfg::PCFG_all, w1::ta, A::cha, B::cha, C::cha)
    x, y = encode(w1,A,B,C)
    cascade_del(enc(1), y, pcfg.sw)
    cascade_del(        x, y, pcfg.t2)
end
function decrement(pcfg::PCFG_all, w1::ta, A::cha, u::ta)
    x, y = encode(w1,A,u)
    cascade_del(enc(0), y, pcfg.sw)
    cascade_del(        x, y, pcfg.t0)
end

###################################################################################
# p_slice をしきい値とし，1か0のみ返せば良い場合(スライスサンプリング)
###################################################################################
# スライスサンプリングの場合の枝刈りつき p,q:
# HPYP全体の確率の上限: p+q*1 = p+q, 下限 p+ q*0 = p
# p < p_slice < p+q でなければ　(p, 0.0) を返す (p のみで判断できる)

function slice(pcfg::PCFG_all, w1::ta, A::cha, B::cha, C::cha, p_slice::Float64 )
    p = p_mulitply_wildcard(pcfg, B, C, 1.0)
    x, y = encode(w1,A,B,C)
    p *= cascade_pq(enc(1), y, pcfg.sw)
    p,q = cascade_pq_with_cutoff(x, y, p_slice,  p, pcfg.t2)
    q == 0.0 ? Int(p_slice < p) : Int(p_slice < p+q)
end
function slice(pcfg::PCFG_all, w1::ta, A::cha, u::ta, p_slice::Float64 )
    x, y = encode(w1,A,u)
    p = cascade_pq(enc(0), y, pcfg.sw) 
    p,q = cascade_pq_with_cutoff(x, y, p_slice,  p, pcfg.t0)
    q == 0.0 ? Int(p_slice < p) : Int(p_slice < p+q)
end

function check_consistency(pcfg::PCFG_all)
    c = pcfg.sw
    check_consistency(c)
    check_consistency(pcfg.t0)
    check_consistency(pcfg.t2)
end


###################################################################################
# すべてのレイヤーのリストを取得
###################################################################################

make_get_all_layers(T) = "all_lay(c::$T) = [" *
     join( [ "c."*string(lay) for lay in fieldnames(T)], "," ) * "]"
eval(Meta.parse(make_get_all_layers(RULE_switch)))
eval(Meta.parse(make_get_all_layers(RULE_type0)))
eval(Meta.parse(make_get_all_layers(RULE_type2)))

function get_all_layers(c::PCFG_all)
    [all_lay(c.sw);all_lay(c.t0);all_lay(c.t2)]
end

###################################################################################
# 確率の総和が1になるか，それぞれの y についてチェックする．
###################################################################################

function check_probsum(c::PCFG_all)
    Σ = _Σ(c)
    for w1 in ta(0):ta(Σ)
        for A in c.V
            p0,p2,p3 = 0.0,0.0,0.0
            for B in cha(1):cha(255) # c.V #
                for C in cha(1):cha(255) # c.V #
                    p2 += get_p(c, w1, A, B, C)
                    x, y = encode( w1, A, B, C)
                    p3 += cascade_pq(        x, y, c.t2)[1]
                end
            end
            for u in ta(1):ta(Σ)
                p0 += get_p(c, w1, A, u)
            end
            #_, y = encode( w1, A, cha(0), cha(0))
            _, y = encode( ta(0), A, cha(0), cha(0))

            p_pos = cascade_pq(enc(1), y, c.sw)
            p_neg = cascade_pq(enc(0), y, c.sw)
            p_more = cascade_pq(enc(2), y, c.sw)

            
            psum = round(p0+p2,digits=8)
            p0 = round(p0, digits=8)
            p2 = round(p2, digits=8)
            p_pos = round(p_pos, digits=8)            
            p_neg = round(p_neg, digits=8)            
            #p_more = round(p_more, digits=8)            
            @show w1,A,  psum, p0, p2 # p_pos + p_neg #p0, p2, psum, p3,
            #check_consistency(c.sw)
        end
    end
end

function check_probsum(c::PCFG_all)
    Σ = _Σ(c)
    for w1 in ta(0):ta(Σ)
        for A in c.V
            p0,p2 = 0.0,0.0
            for B in cha(1):cha(NMAX) # c.V #
                for C in cha(1):cha(NMAX) # c.V #
                    p2 += get_p(c, w1, A, B, C)
                    x, y = encode( w1, A, B, C)
                end
            end
            for u in ta(1):ta(Σ)
                p0 += get_p(c, w1, A, u)
            end
            #_, y = encode( w1, A, cha(0), cha(0))
            _, y = encode( ta(0), A, cha(0), cha(0))
            
            psum = round(p0+p2,digits=8)
            p0 = round(p0, digits=8)
            p2 = round(p2, digits=8)      
            #p_more = round(p_more, digits=8)            
            @show w1,A,  psum, p0, p2 # p_pos + p_neg #p0, p2, psum, p3,
            #check_consistency(c.sw)
        end
    end
end

function check_probsum_V(c::PCFG_all)
    Σ = _Σ(c)
    for w1 in ta(0):ta(Σ)
        for A in c.V
            p0,p2 = 0.0,0.0
            for B in c.V # c.V #
                for C in c.V # c.V #
                    p2 += get_p(c, w1, A, B, C)
                    x, y = encode( w1, A, B, C)
                end
            end
            for u in ta(1):ta(Σ)
                p0 += get_p(c, w1, A, u)
            end
            #_, y = encode( w1, A, cha(0), cha(0))
            _, y = encode( ta(0), A, cha(0), cha(0))
            
            psum = round(p0+p2,digits=8)
            p0 = round(p0, digits=8)
            p2 = round(p2, digits=8)      
            #p_more = round(p_more, digits=8)            
            @show w1,A,  psum, p0, p2 # p_pos + p_neg #p0, p2, psum, p3,
            #check_consistency(c.sw)
        end
    end
end

