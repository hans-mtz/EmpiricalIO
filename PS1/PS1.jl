## PS1 Empirical IO
# Hans Martinez

## Packages
# Pkg.add("Distributions")
# Pkg.add("Parameters")
# Pkg.add(["TexTables","StatsModels","RDatasets"])
# # Pkg.add("Plots")
# Pkg.add("JLD2")
# Pkg.add("Latexify")

# cd("/Volumes/SSD Hans/Github/Emp IO")
using Random, Distributions
using Parameters
using DataFrames, Statistics
using TexTables, StatsModels, GLM, RDatasets
using LinearAlgebra
using Plots
using JLD2
using Latexify



Random.seed!(3331)
N = Normal()
G = Gumbel()



## Question 1 ##
# Generate prices, market shares and observed characteristics 100 products
@with_kw struct P
    j=100 #100 products
    n=1 # consumers
    aux=rand(N,(j,15))
    x = aux[:,1] #+ aux[:,2] + aux[:,3]  # observed product characteristics
    ξ = aux[:,4] #+ aux[:,5] + aux[:,6] # mean utility of product j
    w = aux[:,7] #+ aux[:,8] + aux[:,9]# observed instrument cost shifter
    ω = aux[:,10] # + aux[:,11] + aux[:,12] # unobserved cost shifter
    v = rand(N,n) #
    β = 2.5
    α = 4.0
    η = 1.5
    σ = 2.0
    # ε = rand(G,(j,n)) # idiosyncratic consumer error
    p = η.*x .+ ξ .+ w .+ ω # prices
    δ = β.*x .- α.*p .+ ξ # mean utilities
    μ = σ.*(x*v')
    u = δ' .+ μ' #.+ ε
    # z = aux[:,7] + aux[:,8] + aux[:,10] + aux[:,11] + aux[:,13] + aux[:,14]
    # z1 = aux[:,4] + aux[:,7] + aux[:,10] + aux[:,11]
    # z2 = aux[:,7] + aux[:,10] + aux[:,15]
    df= DataFrame(x=x,ξ=ξ,w=w,ω=ω,p=p,δ=δ)
end

p = P()

# @with_kw struct Par
#     j=100 #100 products
#     n=10000 # consumers
#     aux=rand(N,(j,15))
#     x = aux[:,1] # + aux[:,2] + aux[:,3]  # observed product characteristics
#     ξ = aux[:,4] #+ aux[:,5] + aux[:,6] # mean utility of product j
#     w = aux[:,7] #+ aux[:,8] + aux[:,9]# observed instrument cost shifter
#     ω = aux[:,10] #+ aux[:,11] + aux[:,12] # unobserved cost shifter
#     v = rand(N,n) #
#     β = 1.5
#     α = 0.08
#     η = 0.4
#     σ = 2.0
#     ε = rand(G,(j,n)) # idiosyncratic consumer error
#     p = η.*x .+ ξ .+ w .+ ω # prices
#     δ = β.*x .- α.*p .+ ξ # mean utilities
#     μ = σ.*x*v'
#     uij = δ' .+ μ' # .+ ε
#     z = aux[:,7] + aux[:,8] + aux[:,10] + aux[:,11] + aux[:,13] + aux[:,14]
#     z1 = aux[:,4] + aux[:,7] + aux[:,10] + aux[:,11]
#     z2 = aux[:,7] + aux[:,10] + aux[:,15]
#     df= DataFrame(x=x,ξ=ξ,w=w,ω=ω,p=p,δ=δ,z=z)
# end





## Market shares

function mkt_shr(p::P=p,δ0=0)
    @unpack δ = p
    den = exp(δ0)+sum(exp.(δ))
    share = exp.(δ)/den
    return share
end

𝓈 = mkt_shr()

sum(𝓈)

den = exp(0)+sum(exp.(p.δ))
𝓈0= 1/den

## Normalize price of outside good to 20%

δ0= log(0.2) - log(0.8) + log(sum(exp.(p.δ)))

s = mkt_shr(p,δ0)

sum(s)

den = exp(δ0)+sum(exp.(p.δ))
𝓈0= exp(δ0)/den

## Description of the simulated Data

# @unpack x, p, w, ω, ξ, δ = p,
df= DataFrame(x=p.x,ξ=p.ξ,w=p.w,ω=p.ω,p=p.p,δ=p.δ)
# df=p.df

# sumtab=describe(p.df)

sumtab=summarize(p.df)

to_tex(sumtab) |> print



## Estimate α and β ignoring endogeneity of Prices

y = log.(s).-log(0.2)
y1 = log.(𝓈).-log(𝓈0)
# y0 = log.(s).-log(0.2)

data= DataFrame( y = log.(s).-log(0.2),y1 = log.(𝓈).-log(𝓈0),x=p.x,w=p.w,p=p.p)
ols=lm(@formula(y ~ x + p ), data)
z=predict(lm(@formula(p ~ w + x), data))
data= DataFrame( y = y,y1 = y1,x=p.x,w=p.w,p=p.p,z=z)
iv2s=lm(@formula(y ~  x + z), data)
iv=lm(@formula(y ~x + p + w), data)
iv2s2=lm(@formula(y ~  0 + x + z ), data)

ols2=lm(@formula(y ~  0 + x + p ), data)
z=predict(lm(@formula(p ~ 0 + w + x), data))
data= DataFrame( y = y,y1 = y1,x=p.x,w=p.w,p=p.p,z=z)
iv2s3=lm(@formula(y ~  0 + x + z ), data)
iv2=lm(@formula(y ~  0 + x + p + w), data)
iv2s4=lm(@formula(y ~  x + z ), data)

table = regtable( "OLS" => (ols, ols2),
                  "IV Reduced Form" => (iv, iv2),
                  "IV2S" => (iv2s,iv2s2,iv2s3,iv2s4))

copy_to_clipboard(true)
to_tex(table) |> print
## Instrument correlation

# cor(p.p,p.z)
# cor(p.x,p.z)
# cor(p.ξ,p.z)

println(describe(p.df),"\n",
    "cor(p,w)= ",cor(p.p,p.w),"\n",
    "cor(ξ,w)= ",cor(p.ξ,p.w),"\n",
    table, "\n",
    iv2, "\n",
    iv)

# mdtable(IE, side=fmet, head=["n=4","n=6","n=11","n=21"],fmt = FancyNumberFormatter(4))|> print

ins = [ cor(p.p, p.x) cor(p.p, p.w) cor(p.p, p.p);
        cor(p.ξ, p.x) cor(p.ξ, p.w) cor(p.ξ, p.p);
        cor(p.ω, p.x) cor(p.ω, p.w) cor(p.ω, p.p)]

cortab=latextabular(ins, side=["p","ξ","ω"], head=["x","w","p"], fmt=FancyNumberFormatter(4)) |> print

## Question 2

m = P(n=10000)

# m.u
# m.μ
# m.v

## Description of the simulated Data

Atab=summarize(m.df)

to_tex(Atab) |> print

write_tex("t4.tex",Atab)

v = [mean(m.v), std(m.v), minimum(m.v), maximum(m.v)]

ins_m = [ cor(m.p, m.x) cor(m.p, m.w) cor(m.p , m.p);
          cor(m.ξ, m.x) cor(m.ξ, m.w) cor(m.ξ , m.p);
          cor(m.ω, m.x) cor(m.ω, m.w) cor(m.ω , m.p)]

# to_tex(ins_m) |> print

latextabular(ins_m, side=["p","ξ","ω"], head=["x","w","p"], fmt=FancyNumberFormatter(4)) |> print
# write_tex("t5.tex",ins_m)
## mkt shares
# m.δ
# function mkt_shr_BLP(δ=m.δ, σ=m.σ, n=m.n, m::Par=m)
#     mkt = Par(δ=δ,σ=σ,n=n)
#     @unpack uij,n = mkt
#     den = 1 .+ sum(exp.(uij), dims=2)
#     f =exp.(uij)./den
#     m_s=sum(f,dims=1)/n
#     # sum(m_s)
#     m_s0=sum(1 ./den)/n
#     return m_s', m_s0
# end
#
# m_s, m_s0 = mkt_shr_BLP()
# sum(m_s)


function mkt_shr_BLP(δ0=0, σ=m.σ, δ=m.δ, m::P=m)
    # mkt = Par(δ=δ,σ=σ,n=n)
    @unpack n, x, v = m

    μ = σ.*(x*v')
    # μ = σ.*(x*v')
    uij =  δ' .+ μ'
    # maximum(isnan.(uij))>0 && error("   Got Nan at uij")
    den = exp(δ0) .+ sum(exp.(uij), dims=2)
    # maximum(isnan.(den))>0 && error("   Got Nan at den")
    f = exp.(uij)./den
    # maximum(isnan.(f))>0 && error(" Got Nan at f")
    m_s=sum(f,dims=1)/n
    # maximum(isnan.(m_s))>0 && error("   Got Nan at m_s")
    # sum(m_s)
    m_s0=sum(exp(δ0) ./den)/n
    return m_s', m_s0
end

m_s, m_s0 = mkt_shr_BLP()
sum(m_s)


## Assume consumer have homogeneous taste parameters
# m_s=m_s'
m_s
y=log.(m_s[:,1]).-log(m_s0)
@unpack x, p, w = m
df = DataFrame(y=vec(y),x=x,p=p,w=w,)
olsq2=lm(@formula(y ~ x + p + w), df)
olsq21=lm(@formula(y ~ x + p ), df)
olsq2_2=lm(@formula(y ~ 0 + x + p), df)
olsq2_21=lm(@formula(y ~ 0 + x + p + w), df)

table_2 = regtable(olsq21,olsq2,olsq2_2,olsq2_21)

write_tex("ols_blp.tex",table_2)

## CM (fixed point) algorithm to invert δ

#Initial guess
# δ_old=vec(log.(m_s) .- log(m_s0))

# BLP CM

function BLP_CM(σ0;m::P=m)
    # Initializing parameters
    @unpack j, σ, n, δ = m

    dist=1

    # if max(abs(σ0-σ)) < 0.01
    #     tol=1E-9
    #     flag=0
    # else
    #     tol=1E-6
    #     flag=1
    # end

    norm = 1
    avgnorm = 1
    max_iter=3000
    prev_old=zeros(j)
    tol=1E-6

    #For any guess σ0 ≠ σ
    # σ0=σ
    # n=n
    # s=s
    # δ_old=vec(log.(m_s) .- log(m_s0))

    # guess δ, first time δ=ln(s_j)-ln(s_0) observed shares
    # setting outside good greater than 20%
    x=0.0025
    exp_δ_20 = (x/(1-x)).*mean(sum(exp.(m.u), dims=2))
    δ0_20=log(exp_δ_20)

    # m1, m0 = mkt_shr_BLP(δ0_20)
    # sum(m1)

    s_j, s_0 = mkt_shr_BLP(δ0_20)
    # δ_old=log.(s_j) .- log(s_0)
    δ_old=s_j./s_0
    # δ_old=log.(δ_old)

    println("   Starting CM Loop \n")
    for iter=1:max_iter


        # compute predicted shares with this value of δ
        # 𝓈, 𝓈0 = mkt_shr_BLP(vec(δ_old),σ0,n)

        𝓈, 𝓈0 = mkt_shr_BLP(0,σ0,log.(δ_old))

        if isnan(sum(𝓈))
            δ_old = big.(δ_old)
            𝓈, 𝓈0 = mkt_shr_BLP(0,σ0,log.(δ_old))
            if isnan(sum(𝓈))
                println("Got NaN after big. iter = $iter \n shares = $𝓈")
                return log.(prev_old)
                error("     BLP CM - Got Nan's from market shares \n")
            end
        end

        # println("Sum of shares=",sum(𝓈), "; iter = $iter \n")
        # CM : compute new value of δ
        # δ_new = δ_old .+ log.(s_j) .- log.(𝓈)
        δ_new = δ_old .* s_j ./ 𝓈
        # println("New delta is = $δ_new, iter = $iter \n")
        # check convergence, repeat
        # Convergence criteria taken from Rasmus, taken from Nevo
        # http://www.rasmusen.org/zg604/lectures/blp/frontpage.htm
        # consulted february 26, 2021
        t=abs.(δ_new.-δ_old)
        norm=maximum(t)
        avgnorm=mean(t)
        dist=maximum(abs.(δ_new./δ_old.-1))
        # println("New norm = $norm, avgnorm = $avgnorm, iter = $iter \n")

        if dist< tol #*10^(flag*floor(iter/50))) & (avgnorm < 1e-3*tol*10^(flag*floor(iter/50))) #dist <= tol
            println("δ CM converged")
            println("Iterations = $iter and Distance = ",100*dist,"%")
            # println("Iterations = $iter and Norm = ",norm," avgnorm = ", avgnorm)
            println("------------------------")
            println(" ")
            return log.(δ_new)
        end

        #
        # if (norm < tol*10^(flag*floor(iter/50))) & (avgnorm < 1e-3*tol*10^(flag*floor(iter/50))) #dist <= tol
        #     println("δ CM converged")
        #     # println("Iterations = $iter and Distance = ",100*dist,"%")
        #     println("Iterations = $iter and Norm = ",norm," avgnorm = ", avgnorm)
        #     println("------------------------")
        #     println(" ")
        #     return log(δ_new)
        # end

        prev_old = δ_old

        # update δ value to repeat
        δ_old = δ_new
        iter += 1
        mod(iter,100)==0 && println("   BLP CM Loop: iter=$iter, dist=",100*dist,"%")
        # mod(iter,100)==0 && println("   BLP CM Loop: iter=$iter, norm = ",norm,", avgnorm = ", avgnorm)
    end
    # If loop ends there was no convergence -> Error!
    error("     BLP CM - Solution not found")
end

## CM
# σ0=1.5,n=10000,δ_old=δ_old
δ_hat = BLP_CM(1.98)

# 𝓈, 𝓈0 = mkt_shr_BLP(P(δ=vec(log.(δ_hat)),σ=1.98))
# sum(𝓈)
histogram(δ_hat, nbins=25, label="δ", color=:purple)
title!("Estimated mean utilities with σ=1.98")
savefig("figs/meanu.png")

## Estimate α_hat and β_hat via IV using δ_hat
@unpack x, p, w = m
data= DataFrame( y = vec(δ_hat),x=x,w=w,p=p)
ols_20=lm(@formula(y ~ x + p ), data)
z=predict(lm(@formula(p ~ w + x), data))
data= DataFrame( y = vec(δ_hat),x=x,w=w,p=p,z=z)
iv2s_20=lm(@formula(y ~  x + z), data)
iv_20=lm(@formula(y ~x + p + w), data)
iv2s2=lm(@formula(y ~  0 + x + z ), data)

ols=lm(@formula(y ~  0 + x + p ), data)
z=predict(lm(@formula(p ~ 0 + w + x), data))
data= DataFrame( y = vec(δ_hat),x=x,w=w,p=p,z=z)
iv2s=lm(@formula(y ~  0 + x + z ), data)
iv=lm(@formula(y ~  0 + x + p + w), data)
iv2s3=lm(@formula(y ~  x + z ), data)

table_2 = regtable( "OLS" => (ols_20,ols),
                  "IV Reduced Form" => (iv_20,iv),
                  "IV2S" => (iv2s_20,iv2s2,iv2s3,iv2s))

println(table_2)

write_tex("t7.tex",table_2)

coef(iv2s_20)


## Obtain ξ

ξ_hat = residuals(iv2s_20)

histogram(ξ_hat, nbins=25, label="ξ", color=:purple)
title!("Implied unobserved product attributes with σ=1.98")
savefig("figs/xi.png")

## Instruments

# E[ξ|x,w]=0, then we could unobserved

##

# instrument sum of the characteristics of the other xi

z_rk = [ sum(m.x) - m.x[i] for i in 1:length(m.x)]
w_rk = [ sum(m.w) - m.w[i] for i in 1:length(m.w)]



##

Z = [ones(m.j) z_rk w_rk]
W = diagm([1, 1, 1])
G = (Z'ξ_hat)'W*(Z'ξ_hat)

## repeat chaging σ0 at several point around σ

σ_grid = collect(m.σ-0.5:0.01:m.σ+0.5)
length(σ_grid)
# Z = [ones(m.j) z_rk w_rk]
Z = [ones(m.j) z_rk.^2 w_rk.^2]
W = diagm([1, 1, 1])
coefs=zeros(length(σ_grid),3)
se=zeros(length(σ_grid),3)
G_i=zeros(length(σ_grid))

for i=1:length(σ_grid)
    @unpack x, p, w = m
    y = BLP_CM(σ_grid[i])
    z=predict(lm(@formula(p ~ w + x), DataFrame(x=x,w=w,p=p)))
    df = DataFrame(y=vec(y), α=z, β=x)
    IV = lm(@formula(y ~ β + α), df)
    coefs[i,:] = coef(IV)
    se[i,:] = stderror(IV)
    xi = residuals(IV)
    G_i[i] = (Z'xi)'W*(Z'xi)
end

## Getting  the values for the min
val, index=findmin(G_i)

println(coefs[index,:],"\n",
            se[index,:],"\n",
            σ_grid[index],"\n",
            val)
σ_grid[index]

blp_r = DataFrame( coefs=[coefs[index,1], coefs[index,2 ], coefs[index,3], σ_grid[index]],
                    se=[se[index,1], se[index,2], se[index,3], NaN ])


blp_tab2=latextabular(blp_r, side=["Const","β","α","σ"], head=["Coef","SE"],fmt = FancyNumberFormatter(4)) |> print
## graphs

gr()
plot(σ_grid,G_i, label="Criterion function value")
savefig("figs/cfun2.png")

# plot!(xax,fn,linewidth=3,label = "Function: title",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
# plot!(xax,interp[:,1],linewidth=3,label=fmet[2], linestyle=:dash)
# plot!(xax,interp[:,2],linewidth=3,label=fmet[3],linestyle=:dot)
# plot!(xax,interp[:,3],linewidth=3,label=fmet[4], linestyle=:dashdot)
# plot!(xi,yi,linetype=:scatter,marker=(:circle,6),markercolor=:blue,label = "Data")
# savefig("graphs/$plotname $n.png")
## Saving data to load in Markdown filter
@save "ps1.jld2" sumtab table cortab table_2 G_i blp_tab

# histogram(G_i, nbins=50)
