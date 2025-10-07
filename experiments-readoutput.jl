include("slice_sampling.jl")
include("block_sampling.jl")
include("hyper_parameters.jl")

#using BenchmarkTools
using Random
using JLD
using DelimitedFiles
using Dates


###################################################################################
# 初期化，テスト，訓練，ハイパーパラメータの表示，保存
###################################################################################

function init(sentences)
    c = PCFG_all(301,255)
    bags_init = [ Bag_block(smp,c) for smp in sentences ]
    for bag in bags_init
        resample(bag, c)
        add_sampled_new(bag, c)
    end
    @show _N(c)
    c, bags_init
end

function test(c::PCFG_all, bags::Vector{Bag_block}, logging)
    ps_top = [get_p_top(bag, c) for bag in bags]
    score = mean(log.(ps_top))
    push!(logging, [_N(c), score])
    println( Dates.format(now(), "dd-HH:MM:SS "), _N(c), " ", score  )
end

function train(c::PCFG_all, bags::Vector{T}) where T<:Bag
    for (j, bag) in enumerate(bags)
        del_sampled(bag, c)
        resample(bag, c)
        add_sampled(bag, c)
    end
end

function show_hyper(c::PCFG_all)
    for lay in get_all_layers(c)
        if typeof(lay) <: PYPLayer
            @show lay.θ, lay.d
        end
    end
end

function save_model(i, c::PCFG_all, bags_tr::Vector{T}, logging::Vector{T2}) where {T<:Bag, T2}
    save("./output/model-$(typeof(T)-$i.jld", 
        "pcfg", c, 
        "rules", map(BagSave, bags_tr), 
        "logging", logging )
end

function load_model(i)
    obj = load("./output/model.jld")
    c = obj["pcfg"]
    bags = (x->Bag_slice(x,c)).(obj["rules"])
    c, bags, obj["logging"]
end

###################################################################################
# アウトプットを読み込み，テストを行う．
###################################################################################

sentences = [split(l) for l in readlines("brown16.txt")][3:1000]
sentences = [ [ta(parse(ta, a)+1) for a in l] for l in sentences ]
te = collect(1:length(sentences)).%10 .== 0
data_te = sentences[te]
data_tr = sentences[te .== false]

@show length(data_te) length(data_tr)

#sentences = [[1,2,3,4,1,2], [1,2,3,4,1,2,1,2], [3,4,3,4,3,4], [1,2,3,4,1,2]]
#sentences = [ [ta(a) for a in l] for l in sentences ]

Random.seed!(0)

T= Bag_block

logging = []

#初期化
c, bags_init = init(data_tr)
#訓練用バッグ
bags_tr = [ T(BagSave(bag),c) for bag in bags_init ]
#テスト用バッグ
bags_te = [ Bag_block(smp,c) for smp in data_te ]

update_V(c)
test(c, bags_te, logging)
for i in 1:100
    mcmc_hyperparam(c)
    train(c, bags_tr)
    save_model(i, c, bags_tr, logging)
    if i%10==0
        test(c, bags_te, logging)
    end
end



#check_consistency(c)
#check_probsum(c)
