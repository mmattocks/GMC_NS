@testset "GMC_NS_Ensemble tests..." begin
    sample_dist=Normal(100.,5.)
    samples=rand(sample_dist, 10000)

    n_chains=10

    max_μ=1000
    min_μ=eps()
    max_λ=100
    min_λ=1e-4

    ln_max_μ=log(max_μ^2/sqrt(max_μ^2 + inv(max_λ)))
    ln_min_μ=log(min_μ^2/sqrt(min_μ^2 + inv(min_λ)))
    ln_max_λ=1/log(1+(inv(max_λ)/max_μ^2))
    ln_min_λ=1/log(1+(inv(min_λ)/min_μ^2))

    n_priors=[Uniform(min_μ,max_μ),Uniform(min_λ,max_λ)]
    n_box=[min_μ max_μ;min_λ max_λ]
    ln_priors=[Uniform(ln_min_μ,ln_max_μ),Uniform(ln_min_λ,ln_max_λ)]
    ln_box=[ln_min_μ ln_max_μ;ln_min_λ ln_max_λ]

    gmc=GMC_DEFAULTS
    gmc[1]=3
    gmc[2]=1e-5
    gmc[end]=1e4

    ne=Normal_Ensemble("NE_test", n_chains, samples, n_priors, n_box, gmc...)
    lne=LogNormal_Ensemble("LNE_test", n_chains, samples, ln_priors, ln_box, gmc...)

    uds=Vector{Vector{Function}}()

    lds=Vector{Vector{Function}}()

    ne_ev = converge_ensemble!(ne,backup=(true,100),upper_displays=uds, lower_displays=lds)
    lne_ev = converge_ensemble!(lne,backup=(true,100),upper_displays=uds, lower_displays=lds)

    @test ne_ev > lne_ev

    for e in [ne, lne]
        @info "Testing detailed balance at $(e.path)..."

        @test issorted(e.log_Li)
        @test all(e.log_Li.<0.)
        @test issorted(e.log_Xi,rev=true)
        @test all(e.log_Xi.<=0.)
        @test all(e.log_wi.<0.)
        @test all(e.log_Liwi.<0.)
        @test all(e.log_Zi.<0.)

        Nvec=Vector{Int64}()
        for i in 2:length(e.log_Li)
            Xi=e.log_Xi[i]
            N=round(-(i-1)/Xi)
            push!(Nvec, N)
        end

        @test all(e.GMC_Nmin.<=Nvec.<=n_chains)

        for chain in 1:e.t_counter-1
            ch_pfx=e.path*'/'*string(chain)*'.'
            total_spls=0
            while isfile(ch_pfx*string(total_spls+1))
                total_spls+=1
            end
    
            lastm=deserialize(ch_pfx*"1")
    
            for spl in 2:total_spls
                splm=deserialize(ch_pfx*string(spl))

                @test !(splm.θ!=lastm.θ && splm.v!=lastm.v)
                @test !(splm.θ==lastm.θ && splm.v==lastm.v)
                if splm.θ==lastm.θ
                    @test splm.v!=lastm.v
                    @test splm.log_Li==lastm.log_Li
                elseif splm.v==lastm.v
                    @test splm.θ!=lastm.θ
                    @test splm.log_Li >= lastm.log_Li
                else
                    @error "Galilean trajectory generation is broken!"
                end
                lastm=splm
            end
        end
    end

    rm("NE_test", recursive=true)
    rm("LNE_test", recursive=true)
end    
    