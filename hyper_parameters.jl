include("PCFG.jl")
const FMIN = sqrt(floatmin())
const FMAX = sqrt(floatmax())
###################################################################################
# ハイパーパラメータ (θとd) のサンプリング
###################################################################################
using Random
import Distributions

beta( a, b) = rand(Distributions.Beta(a, b))
gamma(a, b) = rand(Distributions.Gamma(a, b))
bernoulli(r_0, r_1) = Int( rand() * (r_0 + r_1) < r_1 )

function sample_cμ(u::Int32, θ::Float64, d::Float64)
    cμ = 0
    #for n in 1: u-1
    for n in 0: u-1
        cμ += bernoulli( θ, n*d )
    end
    cμ
end

function sample_cνi(bins::Vector{Int32}, d::Float64)
    cνi = 0
    for ci in bins
        for n in 1:ci-1
            cνi += bernoulli( 1-d, n-1 )
        end
    end
    cνi
end

const γa = 1.0
const γb = 1.0
const βa = 0.0
const βb = 0.0

function mcmc_θ_d(lay::T) where T<:PYPLayer
    sum_cμ = 0
    sum_cμ_bar = 0
    sum_cν = 0
    sum_nlogz = 0.0
    for ey in get_nonzero_eys(lay)
        m,u = get_nuy(lay,ey)
        cμ  = sample_cμ(u, lay.θ, lay.d)
        sum_cμ     += cμ
        sum_cμ_bar += u - cμ
        #z = ( m == 1 ? 1.0 : beta(lay.θ+1, m-1) )
        z =  beta(lay.θ, m)
        z =  (z<FMIN ? FMIN : z)
        sum_nlogz  += -log(z)
    end
    for bins in sort!(collect(values(lay.bins_dic)))
        sum_cν += sample_cνi(bins, lay.d)
    end
    lay.θ = gamma(1.0 + γa + sum_cμ,  1.0/(γb + sum_nlogz) )
    lay.d = beta( sum_cμ_bar + 1.0 + βa, sum_cν + 1.0 + βb )
    lay.θ = max(lay.θ, 0.01 )
    lay.θ = min(lay.θ, 1.0)
    lay.d = min(lay.d, 0.5  )
    lay.θ, lay.d
end

mcmc_θ_d(lay::Uniform) = 0


###################################################################################
# 以下，メタプログラミング．
# PCFG 複合型 の コンストラクタ 内で呼び出すと，
# mcmc_hyperparam() の関数が使えるようになる． 
###################################################################################
function init_mcmc_hyperparam(c::T) where T<:PCFG
    #println(str_show_hyperparam(c))
    eval(Meta.parse(str_mcmc_hyperparam(c)))
    eval(Meta.parse(str_show_hyperparam(c)))
end

function str_mcmc_hyperparam(c::T) where T<:PCFG  
    lines=[]
    push!(lines, "function mcmc_hyperparam(pcfg::$T)")
    for modelname in fieldnames(T)
        model = getproperty(c, Meta.parse("$modelname") )
        for layname in fieldnames(typeof(model))
            push!(lines, "    mcmc_θ_d(pcfg.$modelname.$layname)")
        end
    end
    push!(lines,"end")
    join(lines, "\n")
end

function str_show_hyperparam(c::T) where T<:PCFG  
    lines=[]
    push!(lines, "function show_hyperparam(pcfg::$T)")
    for modelname in fieldnames(T)
        model = getproperty(c, Meta.parse("$modelname") )
        for layname in fieldnames(typeof(model))
            lay = getproperty(model, Meta.parse("$layname") )
            if typeof(lay) <: PYPLayer
                push!(lines, "    println(\"[$layname]\", get_hyperparam(pcfg.$modelname.$layname))")
            end
        end
    end
    push!(lines,"end")
    join(lines, "\n")
end
#= 下記のような文字列を生成する．
function mcmc_hyperparam(pcfg::PCFG_all)
    mcmc_θ_d(pcfg.sw.sw_rule)
    mcmc_θ_d(pcfg.sw.sw_mono)
    mcmc_θ_d(pcfg.sw.sw_uni)
    mcmc_θ_d(pcfg.t0.t0_rule1)
    mcmc_θ_d(pcfg.t0.t0_rule)
    mcmc_θ_d(pcfg.t0.t0_mono)
    mcmc_θ_d(pcfg.t0.t0_uni)
    mcmc_θ_d(pcfg.t2.rule1)
    mcmc_θ_d(pcfg.t2.rule)
    mcmc_θ_d(pcfg.t2.ruleL)
    mcmc_θ_d(pcfg.t2.monoL)
    mcmc_θ_d(pcfg.t2.ruleR)
    mcmc_θ_d(pcfg.t2.monoR)
    mcmc_θ_d(pcfg.t2.uni)
end
function show_hyperparam(pcfg::PCFG_sep)
    println("[t0_rule1]", get_hyperparam(pcfg.t0.t0_rule1))
    println("[t0_rule]", get_hyperparam(pcfg.t0.t0_rule))
    println("[t0_mono]", get_hyperparam(pcfg.t0.t0_mono))
    println("[rule1]", get_hyperparam(pcfg.t2.rule1))
    println("[rule]", get_hyperparam(pcfg.t2.rule))
    println("[ruleL]", get_hyperparam(pcfg.t2.ruleL))
    println("[monoL]", get_hyperparam(pcfg.t2.monoL))
    println("[ruleR]", get_hyperparam(pcfg.t2.ruleR))
    println("[monoR]", get_hyperparam(pcfg.t2.monoR))
end
=#
