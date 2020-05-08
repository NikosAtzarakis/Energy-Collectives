function SetPointsPlots(p,str)
    pu = Array{Plot}(undef, number_of_prosumers)
    if str=="Set-Points for generation"
        for j in 1:number_of_prosumers
                jb=string(j)
                pu[j]=Gadfly.plot((12cm,8cm),x=1:time_lapses,y=p[2:2:time_lapses*2,j],Geom.point,Geom.line,Guide.xlabel("Time(h)"),Guide.ylabel("u(kW)"))
                draw(PDF(joinpath(pwd(),"Energy_C.optimization","Energy Collectives 7_11_2019","t=24,pr=6","results\\" * str * "\\Prosumer " * jb * ".pdf"),20cm,20cm),pu[j])

        end
   else
           for j in 1:number_of_prosumers
                   jb=string(j)
                   pu[j]=Gadfly.plot((12cm,8cm),x=1:time_lapses,y=p[1:2:time_lapses*2,j],Geom.point,Geom.line,Guide.xlabel("Time(h)"),Guide.ylabel("u(kW)"))
                   draw(PDF(joinpath(pwd(),"Energy_C.optimization","Energy Collectives 7_11_2019","t=24,pr=6","results\\" * str * "\\Prosumer " * jb * ".pdf"),20cm,20cm),pu[j])

           end
   end
end
