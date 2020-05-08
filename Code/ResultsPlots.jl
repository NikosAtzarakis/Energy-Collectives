function ResultsPlots(p,str)
    stringb=string(str)
    str1="by Prosumer "
    pq = Array{Plot}(undef, number_of_prosumers)
    for i in 1:number_of_prosumers

        ib=string(i)
        pq[i]=Gadfly.plot((12cm,8cm),x=1:time_lapses,y=p[1:time_lapses,i],Guide.xticks(ticks=[1:1:24;]),Guide.yticks(ticks=[-0.5:0.05:0.5;]),Geom.point,Geom.line,Guide.xlabel("Time(h)"),
        Guide.ylabel(  stringb * str1 * ib * " (kW)"))
        draw(PNG(joinpath(pwd(),datapath,"results\\" * str * "\\ Prosumer " * ib * ".png"),20cm,20cm),pq[i])


    end
end
