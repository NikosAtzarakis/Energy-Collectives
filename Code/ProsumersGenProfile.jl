function ProsumersGenProfile(Lower_bounds,Upper_bounds)
    # using GraphPlot
    pGen = Array{Plot}(undef, number_of_prosumers)
    Gen_upper_diag= Array{Float64}(undef, time_lapses)
    Gen_lower_diag= Array{Float64}(undef, time_lapses)
    for j in 1:number_of_prosumers
        for i in 2:2:(time_lapses*2)
            Gen_upper_diag[div(i,2)]=Upper_bounds[i,j]
        end
        for i in 2:2:(time_lapses*2)
            Gen_lower_diag[div(i,2)]=Lower_bounds[i,j]
        end
        jb=string(j)

        pGen=Gadfly.plot(set_default_plot_size(8inch,5inch),layer(x=1:time_lapses,y=Gen_upper_diag,Geom.point,Geom.line,Theme(default_color=color("green"))), layer(x=1:time_lapses,y=Gen_lower_diag,Geom.point,Geom.line,Theme(major_label_font_size=35pt,minor_label_font_size=20pt,default_color=color("red"))),Guide.XLabel("Ώρες (h)"), Guide.xticks(ticks=[1:1:24;]),Guide.yticks(ticks=[0:0.5:10;]) ,Guide.YLabel("Πρόβλεψη παραγωγής Φ/Β του μέλους " * jb * " (kW)"))
        draw(PNG(joinpath(pwd(),datapath,"Prosumers Profiles","Generation\\","Prosumer "  * jb * ".png")),pGen)

    end
end
