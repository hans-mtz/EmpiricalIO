Pkg.add(["Ipopt","ForwardDiff"])

using JuMP, Ipopt, ForwardDiff

##

model = Model(Ipopt.Optimizer)

function mkt_shr_BLP(δ0=0,m::P=m)
    # mkt = Par(δ=δ,σ=σ,n=n)
    @unpack n, δ, μ = m
    uij =  δ' .+ μ'
    # maximum(isnan.(uij))>0 && error("   Got Nan at uij")
    sum = 
    den = 1 .+ sum(exp.(uij), dims=2)
    # maximum(isnan.(den))>0 && error("   Got Nan at den")
    f = exp.(uij)./den
    # maximum(isnan.(f))>0 && error(" Got Nan at f")
    m_s=sum(f,dims=1)/n
    # maximum(isnan.(m_s))>0 && error("   Got Nan at m_s")
    # sum(m_s)
    m_s0=sum(1 ./den)/n
    return m_s', m_s0
end

m_s, m_s0 = mkt_shr_BLP()
sum(m_s)
