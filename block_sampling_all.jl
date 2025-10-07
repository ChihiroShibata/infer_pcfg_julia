include("sampling.jl")
include("PCFG_all.jl")

###################################################################################
# ルールをモデル確率に基づき，サンプルし直す．
###################################################################################
function resample(bag::Bag_block, c::PCFG_all)
    L, N = length(bag.sample), _N(c)
    # テーブルサイズを必要に応じて更新
    size(bag.tab)[3]==N || ( bag.tab = zeros(Float64,L,L,N) )
    #CYKブロック作成
    bag.p_top = build_block_main(bag, c, Val(N))
    #確率消失した場合，scaler の値を更新して繰り返す． 
    for i in 1:10
        FMIN < bag.p_top < FMAX && break
        bag.p_top < FMIN && (bag.p_top = FMIN)
        bag.scaler /= (bag.p_top^(1.0/(L-1)) )
        build_block_main(bag, c, Val(N))
        @show "re-scaling happen:" bag.scaler bag.p_top
    end
    # ルールのサンプリング
    sub_sample_rules(1, length(bag.sample), min(ST_ID,_N(c)), 1, bag, c)
end

###################################################################################
# 現在の文生成確率を計算する．
###################################################################################
function get_p_top(bag::Bag_block, cfg::PCFG_all)
    L, N = length(bag.sample), _N(cfg)
    # テーブルサイズを必要に応じて更新
    size(bag.tab)[3]==N || ( bag.tab = zeros(Float64,L,L,N) )
    # CYK ブロック作成
    bag.p_top = build_block_main(bag, cfg, Val(N))
    if bag.p_top > 1.0 || bag.p_top < 0.0 
        display(bag.tab)
    end
    return bag.p_top * bag.scaler^(1-length(bag.sample))
end


###################################################################################
# CYK table の作成 (ボトムアップ)
###################################################################################
function build_block_main(bag::Bag_block, cfg::PCFG_all, ::Val{N}) where {N}
    #N= _N(cfg)
    # maxrix を先に作っておく．
    #mat = zero(MArray{Tuple{N,N,N},Float64})
    mat = zeros(N,N,N)
    V = cfg.V
    w = bag.sample
    L= length(w)
    tab = bag.tab
    
    #CYKブロックの最下辺
    for i in 1:L
        w1 = (i>1 ? w[i-1] : PAD)
        u = w[i]
        for a in 1:N
            p = get_p(cfg, w1, V[a], u)
            tab[i,i,a] = p
        end
    end
    #CYKブロックを上位へ向かって構築
    for i in L-1:-1:1
        w1 = (i>1 ? w[i-1] : PAD)
        for a in 1:N; for b in 1:N; for c in 1:N
            mat[a, b, c] = get_p(cfg, w1, V[a], V[b], V[c])
        end;end;end
        for j in i+1:L
#            j = i+d
            for a in 1:N
                i==1 && j==L && a!=min(ST_ID,N) && continue
                psum = 0.0
                for k in i:j-1
                    for b in 1:N
                        for c in 1:N
                            #p = get_p(cfg, w1, V[a], V[b], V[c])
                            p = mat[a,b,c] 
                            psum += tab[i,k,b] * tab[k+1,j,c] * p
                        end
                    end
                end
                tab[i,j,a] = bag.scaler * psum
            end
        end
    end
    tab[1,L,min(ST_ID,N)] 
end


###################################################################################
# ルールの再帰的なサンプリング (トップダウン)
###################################################################################
function sub_sample_rules(i::Int64, j::Int64, a::Int64, idx::Int64, bag::Bag_block, cfg::PCFG_all)
    N = _N(cfg)
    V = cfg.V
    w = bag.sample
    w1 = (i>1 ? w[i-1] : PAD)
    i == j &&( bag.sampled_type0[i] = join_xy(encode(w1, V[a], w[i])); return idx)
    tab = bag.tab
    water = rand() * tab[i,j,a] / bag.scaler
    for k in i:j-1
        for b in 1:N
            for c in 1:N
                p = get_p(cfg, w1, V[a], V[b], V[c])
                water -= tab[i,k,b] * tab[k+1,j,c] * p
                if water < 0.0
                    bag.sampled_type2[k] = join_xy(encode(w1, V[a], V[b], V[c]))
                    idx = sub_sample_rules(i,  k, b,  idx+1, bag, cfg)
                    idx = sub_sample_rules(k+1,j, c,  idx, bag, cfg)
                    return idx
                end
            end
        end
    end
    @show i,j,a,water
end
