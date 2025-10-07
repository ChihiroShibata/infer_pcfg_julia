include("PCFG.jl")

###################################################################################
# Sampling の抽象タイプ・共通の型・定義・関数
###################################################################################
const FMIN = sqrt(floatmin())
const FMAX = sqrt(floatmax())

abstract type Bag end

struct BagSave <: Bag
    sample::Vector{ta}
    sampled_type2::Vector{enc}
    sampled_type0::Vector{enc}
    function BagSave(bag::T) where T <: Bag
        new(bag.sample, copy(bag.sampled_type2), copy(bag.sampled_type0))
    end
end


###################################################################################
# サンプルされたルール(Matrix)中の指定の列([2],[2,3,4])に，Aがあれば X に置き換える．
###################################################################################


###################################################################################
# サンプルされたルールの一括の追加と削除
###################################################################################

#add_sampled_new(bag, c) = add_sampled(bag, c, false)

function add_sampled( bag::T, c::G, need_update=true) where {T <: Bag, G <: PCFG}
    # 0 に割り当てられた場合は，新しいVに割り当て直す
    assign_EX(c, bag.sampled_type0, bag.sampled_type2)
    for e in bag.sampled_type0
        increment(c, decode_t0(e)...)
    end
    for e in bag.sampled_type2
        increment(c, decode_t2(e)...)
    end
    need_update && update_V(c)
    0
end

function del_sampled(bag::T, c::G, need_update=true) where {T <: Bag, G<: PCFG}
    for e in bag.sampled_type0
        #@show "del", decode_t0(e)
        decrement(c, decode_t0(e)...)
    end
    for e in bag.sampled_type2
        decrement(c, decode_t2(e)...)
    end
    need_update && update_V(c)
    0
end

###################################################################################
# ブロックサンプリングに必要な一式を収めた構造体
###################################################################################
mutable struct Bag_block <: Bag
    sample::Vector{ta}
    tab::Array{Float64, 3}
    scaler::Float64
    p_top::Float64
    sampled_type2::Vector{enc}
    sampled_type0::Vector{enc}
    function Bag_block(sample::Vector{ta})
        L = length(sample)
        new(sample, 
            zeros(Float64, L, L, 2), 
            1.0,
            0.0,
            zeros(enc, L-1),
            zeros(enc, L)
        )
    end
    function Bag_block(bag::BagSave)
        b = Bag_block(bag.sample)
        b.sampled_type2 = bag.sampled_type2
        b.sampled_type0 = bag.sampled_type0
        b
    end
end