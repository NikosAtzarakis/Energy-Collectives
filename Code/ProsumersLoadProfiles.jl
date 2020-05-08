function ProsumersLoadProfile(Lower_bounds,Upper_bounds)
    pload = Array{Plot}(undef, number_of_prosumers)
    load_upper_diag= Array{Float64}(undef, time_lapses)
    load_lower_diag= Array{Float64}(undef, time_lapses)
    for j in 1:number_of_prosumers
        for i in 1:2:(time_lapses*2)
            load_upper_diag[div(i,2)+1]=Upper_bounds[i,j]
        end
        for i in 1:2:(time_lapses*2)
            load_lower_diag[div(i,2)+1]=Lower_bounds[i,j]
        end
        jb=string(j)

        pload=Gadfly.plot(set_default_plot_size(8inch,5inch),Theme(major_label_font_size=30pt),layer(x=1:time_lapses,y=load_upper_diag,Geom.point,Geom.line,Theme(default_color=color("blue"))), layer(x=1:time_lapses,y=load_lower_diag,Geom.point,Geom.line,Theme(default_color=color("red"))),Guide.XLabel("Ώρες (h)"),Guide.YLabel("Πρόβλεψη Ζήτησης του μέλους " * jb * " (kW)"),Theme(line_width=5pt,point_size=5pt), Guide.xticks(ticks=[1:1:24;]),Guide.yticks(ticks=[0:0.5:8;]) , Guide.manual_color_key("Όρια", ["Άνω", "Κάτω"], ["red", "blue"]))

        draw(PNG(joinpath(pwd(),datapath,"Prosumers Profiles","Load\\","Prosumer "  * jb * ".png")),pload)

    end
end
