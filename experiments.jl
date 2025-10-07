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

function test(c::PCFG_all, bags_te::Vector{Bag_block}, logging)
    ps_top = [get_p_top(bag, c) for bag in bags_te]
    score = mean(log.(ps_top))
    time_str = Dates.format(now(), "mm/dd-HH:MM:SS")
    push!(logging, [time_str, _N(c), score])
    println(time_str, " " , _N(c), " ", score  )
end

function train(c::PCFG_all, bags::Vector{T}) where T<:Bag
    for (j, bag) in enumerate(bags)
        #@show bag.sampled_type2 bag.sampled_type0 c.V _N(c)
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

function save_model(output_dirname, i, 
        c::PCFG_all, bags_tr::Vector{T}, logging::Vector{T2}) where {T<:Bag, T2}
    save(output_dirname*"model-$i.jld", 
        "pcfg", c, 
        "rules", map(BagSave, bags_tr), 
        "logging", logging )
end

function load_model(output_file_name, T)
    obj = load(output_file_name)
    c = obj["pcfg"]
    #@show obj["rules"]
    bags_tr = (x->T(x,c)).(obj["rules"])
    c, bags_tr, obj["logging"]
end

function read_brown()
    sentences = [split(l) for l in readlines("brown16.txt")][3:1000]
    sentences = [ [ta(parse(ta, a)+1) for a in l] for l in sentences ]
    te = collect(1:length(sentences)).%10 .== 0
    data_te = sentences[te]
    data_tr = sentences[te .== false]
    @show length(data_te) length(data_tr)
    data_tr, data_te
end    

###################################################################################
# データを読み込み，学習を行う．
###################################################################################
function main_train(T, output_dirname)
    mkpath(output_dirname)
    data_tr, data_te = read_brown()
    
    #sentences = [[1,2,3,4,1,2], [1,2,3,4,1,2,1,2], [3,4,3,4,3,4], [1,2,3,4,1,2]]
    #sentences = [ [ta(a) for a in l] for l in sentences ]

    Random.seed!(0)

    logging = []

    #初期化
    c, bags_init = init(data_tr)
    #訓練用バッグ
    bags_tr = [ T(BagSave(bag),c) for bag in bags_init ]
    #テスト用バッグ
    bags_te = [ Bag_block(smp,c) for smp in data_te ]

    update_V(c)
    test(c, bags_te, logging)
    save_model(output_dirname, 0, c, bags_tr, logging)
    for i in 1:100
        mcmc_hyperparam(c)
        train(c, bags_tr)
        if i%10==0
            test(c, bags_te, logging)
            save_model(output_dirname, i, c, bags_tr, logging)
        end
    end
end

###################################################################################
# メインパート の切り替え：学習とログの確認
###################################################################################
if length(ARGS)>=2 && ARGS[1] == "train" && ARGS[2] == "block"

    output_dirname =  "./output/block"*Dates.format(now(), "-mm.dd-HH:MM/")
    main_train(Bag_block, output_dirname)

elseif length(ARGS)>=2 && ARGS[1] == "train" && ARGS[2] == "slice"

    output_dirname =  "./output/slice"*Dates.format(now(), "-mm.dd-HH:MM/")
    main_train(Bag_slice, output_dirname)

elseif length(ARGS)>=2 && ARGS[1] == "readlog" 

    output_file_name = ARGS[2]
    c, bags_tr, logging = load_model(output_file_name)
    println("last logging:",  logging[end])
    data_tr, data_te = read_brown()
    bags_te = [ Bag_block(smp,c) for smp in data_te ]
    test(c, bags_te, logging)

elseif length(ARGS)>=2 && ARGS[1] == "resume" 

    output_file_name = ARGS[2]
    c, bags_tr, logging = load_model(output_file_name, Bag_block)
    println("last logging:",  logging[end])
    output_dirname =  "./output/slice"*Dates.format(now(), "-mm.dd-HH:MM/")
    mkpath(output_dirname)
    data_tr, data_te = read_brown()
    Random.seed!(0)
    logging = []
    bags_te = [ Bag_block(smp,c) for smp in data_te ]
    #test(c, bags_te, logging)
    save_model(output_dirname, 0, c, bags_te, logging)
    for i in 1:100
        mcmc_hyperparam(c)
        train(c, bags_tr)
        if i%10==0
            test(c, bags_te, logging)
            save_model(output_dirname, i, c, bags_te, logging)
        end
    end
end


#check_consistency(c)
#check_probsum(c)
