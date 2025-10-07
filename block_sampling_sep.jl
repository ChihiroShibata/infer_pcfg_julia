include("sampling.jl")
include("PCFG_sep.jl")

###################################################################################
# ルールをモデル確率に基づき，サンプルし直す．
###################################################################################
function resample(bag::Bag_block, c::PCFG_sep)
    L, N0, N2 = length(bag.sample), _N0(c), _N2(c)
    # テーブルサイズを必要に応じて更新
    size(bag.tab)[3]==max(N0,N2) || ( bag.tab = zeros(Float64,L,L,max(N0,N2)) )
    #CYKブロック作成
    bag.p_top = build_block_main(bag, c)
    #確率消失した場合，scaler の値を更新して繰り返す． 
    for i in 1:10
        FMIN < bag.p_top < FMAX && break
        bag.p_top < FMIN && (bag.p_top = FMIN)
        bag.scaler /= (bag.p_top^(1.0/(L-1)) )
        build_block_main(bag, c)
        @show "re-scaling happen:" bag.scaler bag.p_top
        @assert false
    end
    # ルールのサンプリング
    sub_sample_rules(1, length(bag.sample), min(ST_ID,N2), 1, bag, c)
end

###################################################################################
# 現在の文生成確率を計算する．
###################################################################################
function get_p_top(bag::Bag_block, cfg::PCFG_sep)
    L, N0, N2 = length(bag.sample), _N0(cfg), _N2(cfg)
    # テーブルサイズを必要に応じて更新
    size(bag.tab)[3]==max(N0,N2) || ( bag.tab = zeros(Float64,L,L,max(N0,N2)) )
    # CYK ブロック作成
    bag.p_top = build_block_main(bag, cfg)
    if bag.p_top > 1.0 || bag.p_top < 0.0 
        display(bag.tab)
    end
    return bag.p_top * bag.scaler^(1-length(bag.sample))
end


###################################################################################
# CYK table の作成 (ボトムアップ)
###################################################################################
function build_block_main(bag::Bag_block, cfg::PCFG_sep)
    #N= _N(cfg)
    # maxrix を先に作っておく．
    #mat = zero(MArray{Tuple{N,N,N},Float64})
    w = bag.sample
    L= length(w)
    tab = bag.tab
    
    N0,V0 = _NV0(cfg)
    N2,V2 = _NV2(cfg)
    mat02 = fill(-1.0,N2,N0,N2) # 確率の再利用用
    mat20 = fill(-1.0,N2,N2,N0) # 確率の再利用用
    mat22 = fill(-1.0,N2,N2,N2) # 確率の再利用用

    #CYKブロックの最下辺を計算．
    for i in 1:L
        w1 = (i>1 ? w[i-1] : PAD)
        u = w[i]
        for a in 1:N0
            tab[i,i,a] = get_p(cfg, w1, V0[a], u)
        end
    end
    #CYKブロックの最下辺のひとつ上を計算．
    for i in 1:L-1
        w1 = (i>1 ? w[i-1] : PAD)
        for a in 1:N2
            psum = 0.0
            for b in 1:N0
                for c in 1:N0
                    p = get_p(cfg, w1, V2[a], V0[b], V0[c])
                    psum += tab[i,i,b] * tab[i+1,i+1,c] * p
                end
            end            
            tab[i,i+1,a] = bag.scaler * psum
        end
    end
    
    #CYKブロックを上位へ向かって構築．
    for i in L-1:-1:1
        w1 = (i>1 ? w[i-1] : PAD)
        #複数回出てくる確率を事前に計算しておく．
        for a in 1:N2
            for b in 1:N2
                for c in 1:N0
                    mat20[a, b, c] = get_p(cfg, w1, V2[a], V2[b], V0[c])
                    mat02[a, c, b] = get_p(cfg, w1, V2[a], V0[c], V2[b])
                end
                # 以下は i と j が 4離れていないと起きない．
                #@show i,L
                if i+3 <= L
                    for c in 1:N2
                        mat22[a, b, c] = get_p(cfg, w1, V2[a], V2[b], V2[c])
                    end
                end
            end
        end
        # 2段上からなので，j=i+2よりスタート．
        for j in i+2:L
            for a in 1:N2
                i==1 && j==L && a!=min(ST_ID,N2) && continue
                psum = 0.0
                for k in i:j-1
                    # mat, Nb, Nc のセレクト
                    mat,Nb,Nc = ( i == k ? (mat02,N0,N2) : ( k == j-1 ? (mat20,N2,N0) : (mat22,N2,N2) ) )
                    for b in 1:Nb
                        for c in 1:Nc
                            p = mat[a,b,c] 
                            #@assert p > 0.0 @show i,k,j,a,b,c, p
                            psum += tab[i,k,b] * tab[k+1,j,c] * p
                        end
                    end
                end
                tab[i,j,a] = bag.scaler * psum
            end
        end
    end
    tab[1,L,min(ST_ID,N2)] 
end


###################################################################################
# ルールの再帰的なサンプリング (トップダウン)
###################################################################################
function sub_sample_rules(i::Int64, j::Int64, a::Int64, idx::Int64, bag::Bag_block, cfg::PCFG_sep)
    N0,V0 = _NV0(cfg)
    N2,V2 = _NV2(cfg)
    w = bag.sample
    w1 = (i>1 ? w[i-1] : PAD)
    #最下辺の場合，type0のサンプリング
    i == j &&( bag.sampled_type0[i] = join_xy(encode(w1, V0[a], w[i])); return idx)
    tab = bag.tab
    #中間ノードの場合，type2のサンプリング
    ratio = rand()
    water = ratio * tab[i,j,a] / bag.scaler
    for k in i:j-1
        #Vb, Vc , Nb, Nc のセレクト
        Vb,Nb = ( i == k ? (V0,N0) : (V2,N2) )
        Vc,Nc = ( k == j-1 ? (V0,N0) : (V2,N2) )
        for b in 1:Nb
            for c in 1:Nc
                p = get_p(cfg, w1, V2[a], Vb[b], Vc[c])
                water -= tab[i,k,b] * tab[k+1,j,c] * p
                if water < 0.0
                    bag.sampled_type2[k] = join_xy(encode(w1, V2[a], Vb[b], Vc[c]))
                    idx = sub_sample_rules(i,  k, b,  idx+1, bag, cfg)
                    idx = sub_sample_rules(k+1,j, c,  idx, bag, cfg)
                    return idx
                end
            end
        end
    end
    @show i,j,a,water, ratio
end
