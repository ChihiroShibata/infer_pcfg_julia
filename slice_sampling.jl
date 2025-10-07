# Making Block sampling
const FMIN = sqrt(floatmin())
const FMAX = sqrt(floatmax())
include("sampling.jl")


###################################################################################
# サンプリングに必要な一式を収めた構造体
###################################################################################
mutable struct Bag_slice <: Bag
    sample::Vector{ta}
    tab::Array{Int64, 3}
    n_path::Int64
    eval_count::Int64
    p_slice_type2::Vector{Float64}
    p_slice_type0::Vector{Float64}
    sampled_type2::Vector{enc}
    sampled_type0::Vector{enc}
    function Bag_slice(sample::Vector{ta}, c::PCFG_all)
        L = length(sample)
        new(sample, 
            zeros(Int64, L, L, _N(c)), 
            0,0,
            zeros(Float64, L-1),
            zeros(Float64, L  ),
            zeros(enc, L-1),
            zeros(enc, L)
        )
    end
    function Bag_slice(bag::BagSave, c::PCFG_all)
        b= Bag_slice(bag.sample, c)
        b.sampled_type2 = bag.sampled_type2
        b.sampled_type0 = bag.sampled_type0
        b
    end
end

###################################################################################
# ルールをモデル確率に基づき，サンプルし直す．
###################################################################################
function resample(bag::Bag_slice, c::PCFG_all)
    L, N = length(bag.sample), _N(c)
    # テーブルサイズを必要に応じて更新
    size(bag.tab)[3]==N || ( bag.tab = zeros(Int64,L,L,N) )
    #スライス値を更新
    for (i,e) in enumerate(bag.sampled_type0)
        bag.p_slice_type0[i] = rand() * get_p(c, decode_t0(e)...)
    end
    for (i,e) in enumerate(bag.sampled_type2)
        bag.p_slice_type2[i] = rand() * get_p(c, decode_t2(e)...)
    end
    #CYKブロック作成
    bag.n_path, bag.eval_count = build_block_main(bag, c)
    # ルールのサンプリング
    sub_sample_rules(1, length(bag.sample), min(ST_ID,_N(c)), 1, bag, c)
end

###################################################################################
# CYK table の作成 (ボトムアップ)
###################################################################################
function build_block_main(bag::Bag_slice, cfg::PCFG_all)
    eval_count = 0
    N= _N(cfg)
    # maxrix を先に作っておく．
    #mat = BitArray(undef, N, N, N)
    V = cfg.V
    w = bag.sample
    L= length(w)
    tab = bag.tab
    
    #CYKブロックの最下辺
    for i in 1:L
        w1 = (i>1 ? w[i-1] : PAD)
        u = w[i]
        p_slice = bag.p_slice_type0[i]
        for a in 1:N
            val = slice(cfg, w1, V[a], u, p_slice)
            tab[i,i,a] = val
        end
    end
    
    # slice[k]の値で B,C に対し，どれかのAで slice count が発生したとき1をたてる．
    need_check_A = trues(L,N,N)

    #CYKブロックを上位へ向かって構築
    for i in L-1:-1:1
        fill!(need_check_A, true)
        w1 = (i>1 ? w[i-1] : PAD)
        for j in i+1:L
            for a in 1:N; tab[i,j,a] = 0; end
            for k in i:j-1
                p_slice = bag.p_slice_type2[k]
                for b in 1:N
                    tab[i,k,b] == 0 && continue
                    for c in 1:N
                        tab[k+1,j,c] == 0 && continue
                        need_check_A[k,b,c] || continue
                        slice_count_A = 0
                        for a in 1:N
                            i==1 && j==L && a!=min(ST_ID,N) && continue
                            eval_count += 1
                            bit = slice(cfg, w1, V[a], V[b], V[c], p_slice)
                            slice_count_A += bit
                            tab[i,j,a] += bit * tab[i,k,b] * tab[k+1,j,c]
                        end
                        need_check_A[k,b,c] = (slice_count_A>0)
                    end
                end
                #@show i,j,n_sum
            end
        end
    end
    
    #display(tab)

    tab[1,L,min(ST_ID,N)], eval_count
end


###################################################################################
# ルールの再帰的なサンプリング (トップダウン)
###################################################################################
function sub_sample_rules(i::Int64, j::Int64, a::Int64, idx::Int64, bag::Bag_slice, cfg::PCFG_all)
    #@show i, j, a
    N = _N(cfg)
    V = cfg.V
    w = bag.sample
    w1 = (i>1 ? w[i-1] : PAD)
    i == j &&( bag.sampled_type0[i] = join_xy(encode(w1, V[a], w[i])); return idx)
    tab = bag.tab
    tab[i,j,a] == 0 && (@show (i,j,a);display(tab))
    water = rand(1:tab[i,j,a])
    for k in i:j-1
        p_slice = bag.p_slice_type2[k]
        for b in 1:N
            tab[i,k,b] == 0 && continue
            for c in 1:N
                tab[k+1,j,c] == 0 && continue
                water -= slice(cfg, w1, V[a], V[b], V[c], p_slice) * 
                         tab[i,k,b] * tab[k+1,j,c] 
                if water <= 0
                    #if (i,j,a) == (1,4,2)
                    #    @show (k,b,c),(i,j,a)
                    #end
                    #@show i,j,k,b,c, tab[i,k,b], tab[k+1,j,c] 
                    bag.sampled_type2[k,:] .= join_xy(encode(w1, V[a], V[b], V[c]))
                    idx = sub_sample_rules(i,  k, b,  idx+1, bag, cfg)
                    idx = sub_sample_rules(k+1,j, c,  idx, bag, cfg)
                    return idx
                end
            end
        end
    end
    @show i,j,a,water
end

