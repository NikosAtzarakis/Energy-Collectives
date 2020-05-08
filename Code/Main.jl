using JuMP
using Ipopt
using Gadfly
using DataFrames
using LinearAlgebra
using CSV
using Cairo

##############################################################################
#PATHS
functionspath="Energy_C.optimization\\Energy Collectives 7_11_2019\\Code\\"
datapath= "Energy_C.Optimization\\Energy Collectives 7_11_2019\\Data\\"
##############################################################################

##############################################################################
include(joinpath(pwd() ,functionspath, "ProsumersGenProfile.jl"))
include(joinpath(pwd(),functionspath,"ProsumersLoadProfiles.jl"))
include(joinpath(pwd(),functionspath,"ResultsPlots.jl"))
include(joinpath(pwd(),functionspath,"SetPointsPlots.jl"))



##########################################################################################
m= Model(with_optimizer(Ipopt.Optimizer))

prosumers=DataFrame(CSV.read(joinpath(pwd(),datapath,"15prosumersapril.csv");datarow=3)) #διαβαζει δεδομενα για παραγωγη και φορτιο
prices=DataFrame(CSV.read(joinpath(pwd(),datapath,"Prices.csv"))) #διαβαζει τιμες αγορας
number_of_prosumers=size(prosumers,1)  #αριθμος prosumers
time_lapses=(size(prosumers,2)-1)/4    #χρονικες στιγμες
time_lapses=floor(Int,time_lapses)     #το κανω floor για να γινει int απο float
γcom=prices[1,1] #τιμη συναλλαγης μεσα στην κοινοτητα
γexp=prices[1,2] #τιμη πωλησης στο δικτυο
γimp=prices[1,3] #τιμη αγορας απο το δικτυο
γmax=prices[1,4] #προστιμο σε αυτον που αγοραζει την περισσοτερη ενεργεια απο το δικτυο
a1=prices[1,5] #μεταβλητη που χρειαζονται για την ευρεση των c,d
b1=prices[1,6] #μεταβλητη που χρειαζονται για την ευρεση των c,d



######################################################################
            #ΠΑΝΩ ΚΑΙ ΚΑΤΩ ΟΡΙΑ ΤΩΝ SET-POINTS

Lower_bounds = Array{Float64}(undef, (time_lapses*2) , number_of_prosumers)
for i in 1:(number_of_prosumers)
    for j in 1:(time_lapses*2)
        if j % 2 != 0
            Lower_bounds[j,i] = prosumers[i,j+j+2]
        else
            Lower_bounds[j,i]= prosumers[i,j+(j-2)]
        end

    end
end

Upper_bounds = Array{Float64}(undef, (time_lapses*2) , number_of_prosumers)
for i in 1:(number_of_prosumers)
    for j in 1:(time_lapses*2)
        if j % 2 != 0
            Upper_bounds[j,i] = prosumers[i,j+(j+3)]
        else
            Upper_bounds[j,i]= prosumers[i,j+(j-1)]
        end

    end
end


######################################################################
#ΚΑΜΠΥΛΕΣ ΦΟΡΤΙΟΥ ΚΑΙ ΠΑΡΑΓΩΓΗΣ ΟΛΩΝ ΤΩΝ PROSUMERS ΠΡΙΝ ΓΙΝΕΙ Η ΒΕΛΤΙΣΤΟΠΟΙΗΣΗ
#προφιλ καταναλωτων/παραγωγων

ProsumersLoadProfile(Lower_bounds,Upper_bounds)
ProsumersGenProfile(Lower_bounds,Upper_bounds)


######################################################################
             #ΣΤΑΘΕΡΕΣ C,D ΓΙΑ ΤΗΝ ΣΥΝΑΡΤΗΣΗ ΚΟΣΤΟΥΣ Ψ
             #ΧΡΗΣΙΜΟΠΟΙΗΘΗΚΑΝ ΟΙ ΤΥΠΟΙ ΤΟΥ PAPER

c = Array{Float64}(undef, (time_lapses*2) , number_of_prosumers)  #σταθερα c για την συναρτηση κοστους ψ
for i in 1:number_of_prosumers
  for j in 1:(time_lapses*2)
      if j % 2 !=0    # ΑΝ j%2 != 0  ΤΟΤΕ ΠΡΟΚΕΙΤΑΙ ΓΙΑ ΦΟΡΤΙΟ
          if Lower_bounds[j,i] == Upper_bounds[j,i]
              c[j,i]= (2*b1)/(-Lower_bounds[j,i])
          else
              c[j,i]=(2*b1)/(Upper_bounds[j,i]-Lower_bounds[j,i])
          end

     else           #ΑΝ j%2==0 ΤΟΤΕ ΠΡΟΚΕΙΤΑΙ ΓΙΑ ΠΑΡΑΓΩΓΗ
         if Lower_bounds[j,i] == Upper_bounds[j,i]
             c[j,i]=0
         else
             c[j,i]=-(2*b1)/(Upper_bounds[j,i]-Lower_bounds[j,i])
         end
     end
end
end

d = Array{Float64}(undef, (time_lapses*2) , number_of_prosumers)  #σταθερα d για την συναρτηση κοστους ψ
for i in 1:number_of_prosumers
  for j in 1:(time_lapses*2)
      if j % 2 !=0
          if Lower_bounds[j,i] == Upper_bounds[j,i]
              d[j,i]=a1
          else
             d[j,i]=a1 + (b1)*(-(Upper_bounds[j,i]+Lower_bounds[j,i])/(Upper_bounds[j,i]-Lower_bounds[j,i]))
         end
     else
         if Lower_bounds[j,i] == Upper_bounds[j,i]
             d[j,i]=0
         else
             d[j,i]=-(a1 + (b1)*((Upper_bounds[j,i]+Lower_bounds[j,i])/(Upper_bounds[j,i]-Lower_bounds[j,i])))
         end
     end
end
end


####################################################################
                        #VARIABLES

@variables m begin
    u[1:(time_lapses)*2,1:number_of_prosumers]        #set-points
    a[1:time_lapses,1:number_of_prosumers]            #energy imported from the grid
    b[1:time_lapses,1:number_of_prosumers]            #energy exported to the grid
    q[1:time_lapses,1:number_of_prosumers]            #energy traded within the community
    qb[1:time_lapses,1:number_of_prosumers]           #βοηθητικη μεταβλητη για την απολυτη τιμη στην εκφραση h
    # ab                                                #βοηθητικη μεταβλητη για την απειρη νορμα στην εκφραση g
end


###################################################################
                        #EXPRESSIONS

@expression(m,qimp[j=1:time_lapses] ,sum(a[j,i] for i in 1:number_of_prosumers))    #total energy imported from the grid
@expression(m,qexp[j=1:time_lapses] ,sum(b[j,i] for i in 1:number_of_prosumers))     #total energy exported to the grid
@expression(m, g[i=1:time_lapses] , γimp*qimp[i] + γexp*qexp[i])        #penalty for the maximum importer  #the role of the Community Manager

coefhcom=1    #συντελεστης επηρεαζει την ποσοτητα ενεργειας που ανταλλασεται εντος κοινοτητας
coefhimp=1     #συντελεστης που επηρεαζει την ενεργεια που εισερχεται στην κοινοτητα
coefhexp=1     #συντελεστης που επηρεαζει την ενεργεθα που εξερχεται απο την κοινοτητα
h=@expression(m,[1:time_lapses,1:number_of_prosumers],0)     #homogeneous costs previously agreed within the community
for i in 1:number_of_prosumers
  for j in 1:time_lapses
    h[j,i]=coefhcom*γcom*qb[j,i]+coefhimp*γimp*a[j,i]+coefhexp*γexp*b[j,i]
  end
end


ψ = Matrix{JuMP.QuadExpr}(undef, 2*time_lapses, number_of_prosumers) #συναρτηση κοστους ψ συμφωνα με το paper
for i in 1:2time_lapses
    for j in 1:number_of_prosumers
        ψ[i, j] = zero(JuMP.QuadExpr)
    end
end

for i in 1:number_of_prosumers
   for j in 1:(time_lapses*2)
       JuMP.add_to_expression!(ψ[j,i], c[j,i]u[j,i]+d[j,i], u[j,i])
   end
end


revenues=@expression(m,[1:(time_lapses),1:number_of_prosumers],0)   #εσοδα απο την πωληση ενεργειας στο δικτυο
for i in 1:number_of_prosumers
    for j in 1:(time_lapses)
        revenues[j,i]=γexp*u[j+j,i]
    end
end


Pnet= @expression(m, [1:time_lapses,1:number_of_prosumers], 0 )     #Net Power of each prosumer #καθαρη εξερχομενη ενεργεια
for i in 1:number_of_prosumers
  global t=1
  for j in 1:time_lapses
    for k in t:(t+1)
     Pnet[j,i]= Pnet[j,i]+u[k,i]
    end
    global t+=2
  end
end




################################################################
                        #CONSTRAINTS

@constraint(m ,[i=1:number_of_prosumers,j=1:(time_lapses*2)], Lower_bounds[j,i] <= u[j,i])
@constraint(m ,[i=1:number_of_prosumers,j=1:(time_lapses*2)], Upper_bounds[j,i] >= u[j,i])
# @constraint(m,a[1:end,1:end] .>= -ab)  #περιορισμος για την μεταβλητη a, ωστε να παρω την απειρη νορμα της στην εκφραση g
# @constraint(m,a[1:end,1:end] .<= ab)   #περιορισμος για την μεταβλητη a, ωστε να παρω την απειρη νορμα της στην εκφραση g
@constraint(m,a[1:end,1:end] .>= 0)
@constraint(m,b[1:end,1:end] .>= 0)
@constraint(m,sum(q,dims=2) .== 0)     #το αθροισμα των ενεργειων που ανταλλασονται μεσα στην κοινοτητα πρεπει να ειναι ισο με το μηδεν
@constraint(m, [j=1:time_lapses,i=1:number_of_prosumers] , Pnet[j,i] == b[j,i] - q[j,i]-a[j,i])   #ισοζυγιο ισχυος
@constraint(m,[j=1:time_lapses,i=1:number_of_prosumers] ,q[j,i] <= qb[j,i])     #περιορισμος για την απολυτη τιμη
@constraint(m,[j=1:time_lapses,i=1:number_of_prosumers] ,q[j,i] >= -qb[j,i])    #περιορισμος για την μεταβλητη q ωστε να παρω την απολυτη τιμη της στην εκφραση h



###################################################################
                            #OBJECTIVE
@objective(m,Min, sum(ψ[1:end,1:end]) +sum(h[1:end,1:end]) + sum(g[1:end]) - sum(revenues[1:end,1:end]))
optimize!(m)




##########################################################################
##########################################################################
a_val = JuMP.value.(a)
b_val=JuMP.value.(b)
q_val=JuMP.value.(q)
u_val=JuMP.value.(u)
qimp_val=JuMP.value.(qimp)
qexp_val=JuMP.value.(qexp)
for i in 1:number_of_prosumers #οι τιμες οι οποιες ειναι πολυ κοντα στο 0 γινονται 0
    for j in 1:time_lapses*2
        if abs(u_val[j,i]-0) <= 0.001
            u_val[j,i] = 0
        end
    end
end
for i in 1:number_of_prosumers
    for j in 1:time_lapses
        if abs(a_val[j,i]-0) <= 0.00001
            a_val[j,i]=0
        end
        if abs(b_val[j,i]-0) <= 0.00001
            b_val[j,i]=0
        end
        if abs(q_val[j,i]-0) <= 0.00001
            q_val[j,i]=0
        end
    end
end
for i in 1:time_lapses
    if abs(qimp_val[i] -0) <= 0.001
        qimp_val[i]=0
    end
    if abs(qexp_val[i]-0)<=0.001
        qexp_val[i]=0
    end
end
##########################################################################
##########################################################################



#########################################################################
                                #PLOTTING


ResultsPlots(-q_val,"Energy traded within the community")
ResultsPlots(a_val,"Energy imported from the grid")
ResultsPlots(b_val,"Energy exported to the grid")
SetPointsPlots(u_val,"Set-Points for generation")
SetPointsPlots(u_val,"Set-Points for load")
pqimp=Gadfly.plot((12cm,8cm),x=1:time_lapses,y=qimp_val[1:time_lapses],Geom.point,Geom.line,Guide.xlabel("Time(h)"),Guide.ylabel("Qimp(kW)"))
draw(PDF(joinpath(pwd(),datapath,"results","Total energy imported from the grid", "Qimp.pdf"),20cm,20cm),pqimp)
pqimp=Gadfly.plot((12cm,8cm),x=1:time_lapses,y=qexp_val[1:time_lapses],Geom.point,Geom.line,Guide.xlabel("Time(h)"),Guide.ylabel("Qexp(kW)"))
draw(PDF(joinpath(pwd(),datapath,"results","Total energy exported to the grid", "Qexp.pdf"),20cm,20cm),pqimp)
#######################################################################
