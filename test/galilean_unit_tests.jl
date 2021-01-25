@testset "Galilean trajectory unit tests..." begin
    #test box movement
    pos=[.5,.5]
    box=[0 1.
         0 1.]
    allow_d=[.1,.1]
    box_d=[1.,1.]
    
    allowed_pos, allowed_d=box_move(pos,allow_d,box)
    @test allowed_pos==pos+allow_d
    @test allowed_d == allow_d

    boxed_pos, boxed_d=box_move(pos, box_d, box)
    @test boxed_pos==[1.,1.]
    @test boxed_d==[.5,.5]

    box_reflect_v=box_reflect([.5,.99999999],box,[.1,.1])
    @test box_reflect_v==[.1,-.1]

    #test boundary norm finding
    norm=1000
    v=[.5,.5]
    bn=boundary_norm(v, norm)
    @test bn == normalize(norm./v)

    #test reflection given norm
    rv=reflect(v,bn,0.)
    @test isapprox(rv, -v)
    
    #test reflection perturbaton
    pv=r_perturb(v,.01)
    @test v!=pv
end