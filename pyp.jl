

NUD = Dict{UInt64,Tuple{Int32,Int32}}
NUV = Vector{Tuple{Int32,Int32}}
BinD = Dict{UInt64, Vector{Int32}}
_NUV() = [(Int32(0),Int32(0)) for i in 1:256]

abstract type PYPLayer end
abstract type PYPModel end

mutable struct PL_DicY <: PYPLayer
    mask::UInt64
    mask_y::UInt64
    θ::Float64
    d::Float64
    bins_dic::BinD
    dic::NUD
    dic_y::NUD
    PL_DicY(mask_x, mask_y) = new(mask_x+mask_y, mask_y, 0.5, 0.5, BinD(), NUD(),  NUD())
end
mutable struct PL_MapY <: PYPLayer
    mask::UInt64
    mask_y::UInt64
    θ::Float64
    d::Float64
    bins_dic::BinD
    dic::NUD
    vec_y::NUV
    PL_MapY(mask_x, mask_y) = new(mask_x+mask_y, mask_y, 0.5, 0.5, BinD(), NUD(), _NUV())

end
mutable struct PL_NonY <: PYPLayer
    mask::UInt64
    mask_y::UInt64
    θ::Float64
    d::Float64
    bins_dic::BinD
    dic::NUD
    ny::Int32
    uy::Int32
    PL_NonY(mask_x) = new(mask_x,    0x0, 0.5, 0.5, BinD(), NUD(), 0, 0  )

end
struct Uniform
    N::UInt64
end
#PL_DicY(mask, mask_y) = PL_DicY(mask, mask_y, 0.5, 0.5, BinD(), NUD(),  NUD())
#PL_MapY(mask, mask_y) = PL_MapY(mask, mask_y, 0.5, 0.5, BinD(), NUD(), _NUV())
#PL_NonY(mask)         = PL_NonY(mask,    0x0, 0.5, 0.5, BinD(), NUD(), 0, 0  )

get_nuy(lay::PL_DicY, ey) = get(lay.dic_y, ey, (Int32(0),Int32(0)) )
get_nuy(lay::PL_MapY, ey) = lay.vec_y[ey >> 2LEN_CHA + 1]
get_nuy(lay::PL_NonY, ey) = (lay.ny, lay.uy)
set_nuy(lay::PL_DicY, ny, uy, ey) = (lay.dic_y[ey] = (ny,uy) )
set_nuy(lay::PL_MapY, ny, uy, ey) = (lay.vec_y[ey >> 2LEN_CHA + 1] = (ny,uy) )
set_nuy(lay::PL_NonY, ny, uy, ey) = (lay.ny = ny; lay.uy = uy)
del_nuy(lay::PL_DicY, ey) = delete!(lay.dic_y, ey)
del_nuy(lay::PL_MapY, ey) = true
del_nuy(lay::PL_NonY, ey) = true

get_nu(lay, e)       = get(lay.dic, e, (Int32(0),Int32(0)) )
set_nu(lay, n, u, e) = (lay.dic[e] = (n,u) )
del_nu(lay, e)       = delete!(lay.dic, e) 
get_bins(lay, e) = get( lay.bins_dic,   e, Vector{Int32}() )
set_bins(lay, bins, e) = (lay.bins_dic[e] = bins )
del_bins(lay, e)       = delete!(lay.bins_dic, e) 

get_nonzero_eys(lay::PL_DicY) = sort!(collect(keys(lay.dic_y)))
get_nonzero_eys(lay::PL_MapY) = UInt64.(findall(y->y[1]>0, lay.vec_y).- 1).<< 2LEN_CHA
get_nonzero_eys(lay::PL_NonY) = isempty(lay.dic) ? [] : [UInt64(0)]

function rename(lay, e, e2, ey, ey2)
    ny, uy= get_nuy( lay, ey); set_nuy( lay, ny, uy, ey2); del_nuy( lay, ey)
    n, u  = get_nu(  lay, e ); set_nu(  lay, n,  u,  e2 ); del_nu(  lay, e )
    bins  = get_bins(lay, e ); set_bins(lay, bins,   e2 ); del_bins(lay, e )
end

# PYPレイヤにおける p, q を計算する．
# p : このレイヤ，またはより上位のレイヤからサンプル x|y を引くときの確率
# q : 下位レイヤ(ベース)から x|y を引くときの確率の係数
@inline function pq(lay::T, xy::UInt64, p=0.0, q=1.0 ) where T<:PYPLayer
    n,u   = get_nu(lay, xy&lay.mask    )
    ny,uy = get_nuy( lay, xy&lay.mask_y  )
    deno = 1.0/( ny + lay.θ )
    cp = (     n - lay.d * u  ) * deno
    cq = ( lay.θ + lay.d * uy ) * deno
    p += q * cp
    q *= cq
    #@show typeof(lay),p,q
    (p, q)
end
@inline function pq( lay::Uniform, p=0.0, q=1.0  )
     p, q/lay.N
end

################################################################################
# 破壊的操作 : sample の追加と削除
################################################################################
function roulette(bins::Vector{Int32}, n::Int32, d::Float64)
    water = rand()*(n-length(bins)*d)    
    for (idx,v) in enumerate(bins)
        water -= (v-d)
        water <= 0.0 && return idx
    end
    #@assert true @show n, sum(bins), bins
    return 0
end

#各レイヤから，x|y をサンプルしてしまう．
# water > 0.0 なら，下位レイヤーからもサンプルが必要．
function add(lay::T, xy::UInt64, water, p=0.0, q=1.0) where T<:PYPLayer
    e, ey = xy .& (lay.mask, lay.mask_y)
    n, u  = get_nu(lay, e )
    ny,uy = get_nuy( lay, ey)
    bins  = get_bins(lay, e )
    #@show ey, xy
    #@assert length(bins) == u @show bins, u
    #@show e, ey, typeof(lay), bins, get_nuy(lay, ey) ,water, p, q
    p, q = pq(lay, xy, p, q)
    #s@show typeof(lay), n,u, p, length(bins)
    if water < p
        idx = roulette(bins, n, lay.d)
        bins[idx] += 1
        set_nu( lay, n+1,  u,  e )
        set_nuy(lay, ny+1, uy, ey)
        set_bins(lay, bins, e)
    else
        append!(bins,1)
        set_nu( lay, n+1,  u+1,  e )
        set_nuy(lay, ny+1, uy+1, ey)
        set_bins(lay, bins, e)
    end
    #@show typeof(lay), bins, get_nuy(lay, ey) ,water, p, q
    return p, q
end

#もし下位レイヤーを削除する必要がない場合は true を返す
function del(lay::T, xy::UInt64) where T <:PYPLayer
    e, ey = xy .& (lay.mask, lay.mask_y)
    n, u  = get_nu(lay, e )
    ny,uy = get_nuy( lay, ey)
    bins  = get_bins(lay, e ) 
    idx = roulette(bins, n, 0.0)
    bins[idx] -= 1
    if bins[idx] > 0
        set_nu( lay, n-1,  u,  e )
        set_nuy(lay, ny-1, uy, ey)
        set_bins(lay, bins, e)
        true
    else
        deleteat!(bins, idx)
        set_nu( lay, n-1,  u-1,  e )
        set_nuy(lay, ny-1, uy-1, ey)
        u -1 > 0 || del_nu(lay, e)
        uy-1 > 0 || del_nuy(lay, ey)
        u-1  > 0 || del_bins(lay, e)
        false
    end
end

################################################################################
# cutoffが必要な場合
################################################################################    
# スライスサンプリングの場合の枝刈りつき p,q:
# HPYP全体の確率の上限: p+q*1 = p+q, 下限 p+ q*0 = p
# p < p_slice < p+q でなければ　(p or p+q, 0.0) を返す 
# (p or p+q のみで判断できる)

@inline function pq_with_cutoff(lay::T, xy::UInt64, p_slice::Float64, 
                                p=0.0, q=1.0 ) where T<:PYPLayer
    #@show p,q
    p, q = pq(lay, xy, p, q)
    #@show p,q
    if p_slice < p
        (p, 0.0)
    elseif p+q < p_slice
        (p+q, 0.0)
    else
        (p, q)
    end
end


################################################################################
# エラーチェック用
################################################################################

to_str(xy::UInt64)=join(["$((xy>>(8*i))&0xff)" for i in 0:Int(floor(log(0x100,xy)))]," ")

get_hyperparam(lay::T) where T<:PYPLayer = (lay.θ, lay.d)

@eval function check_consistency(c::T) where T<:PYPModel
    for name in fieldnames(typeof(c))
        lay = getproperty(c, Meta.parse("$name") )
        check_consistency(lay, name)
    end
end

function check_consistency(lay::T, name="") where T<:PYPLayer
       """
        1. dic_y, dic_xy に登録されている  y の集合が等しいかチェックする．
        2. サンプル Ny の総和， アサイン Uy の総和を出力する．
        3. 全ての y について 次をチェックする：
           a) dic_xy.{N, U} を xについて総和を取ると Ny, Uy と等しいか チェックする．
           b) dic_xy.U が 1 以上かチェックする．
           c) dic_xy.N >= urn_xy.U となっているかチェックする．
        4. 登録されている x の数を出力する．
        """
    print("[$name] ")
    #lay.mask == 0xffffffff && @show lay.bins_dic
    #1. dic_y, dic_xy に登録されている  y の集合が等しいかチェックする．
    eys = get_nonzero_eys(lay)
    nus = [get_nuy(lay,ey) for ey in eys]
    xys = collect(keys(lay.dic))
    eys_agg_xy   = collect(xy & lay.mask_y for xy in xys)
    @assert Set(eys)==Set(eys_agg_xy) @show Set(eys),Set(eys_agg_xy)
    
    #2. サンプル Ny の総和， アサイン Uy の総和を出力する．
    nsum = sum((x->x[1]).(nus))
    usum = sum((x->x[2]).(nus))
    print("N,U = $nsum, $usum, " )
    print("|{(x,y)}| = $(length(xys)), ")
    print("|{y}| = $(length(eys)), ")
    
    x_set = Set{UInt64}() # for 4
    #3. 全ての ey について 次をチェックする：
    for ey in eys
        # a) dic_xy.{N, U} を xについて総和を取ると Ny, Uy と等しいか チェックする．
        select = xys[findall(xy->(xy& lay.mask_y)==ey , xys)]
        nus = [lay.dic[xy] for xy in select]
        nsum = sum((x->x[1]).(nus))
        usum = sum((x->x[2]).(nus))
        @assert nsum == get_nuy(lay,ey)[1] @show nsum, get_nuy(lay,ey)[1]
        @assert usum == get_nuy(lay,ey)[2] @show nsum, get_nuy(lay,ey)[2]
        
        #  b) すべての x でdic_xy.U が 1 以上かチェックする．
        #  c) dic_xy.N >= urn_xy.U となっているかチェックする．
        for xy in select
            #@show select
            push!(x_set, xy-ey)
            (n,u) = lay.dic[xy]
            @assert u>0 @show  u,to_str(xy)
            @assert n>=u @show n,u,to_str(xy)
        end
    end
    print("|{x}| = $(length(x_set))")
    println("")
end

check_consistency(lay::Uniform, name="") = println("[$name] $(lay.N)")