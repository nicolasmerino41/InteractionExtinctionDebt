using CSV,DataFrames,CairoMakie,Statistics
root=abspath(ARGS[1]);out=joinpath(root,"All/outputs/29_degree_distribution_audit");mkpath(out)
mapping=[("Quercus","Quercus","quercus"),("Nahuel","Nahuel","nahuel"),("Salix_Galpar","Salix","galpar"),("Gottin_HP","Gottin","Gottin-HP"),("Gottin_PP","Gottin","Gottin-PP"),("Garraf_HP","Garraf-Montseny-Olot","garraf-HP"),("Garraf_PP","Garraf-Montseny-Olot","garraf-PP"),("Garraf_PP2","Garraf-Montseny-Olot","garraf-pp2"),("Olot","Garraf-Montseny-Olot","olot"),("Montseny","Garraf-Montseny-Olot","montseny")]
current=CSV.read(joinpath(root,"All/CombinedOutputs/29_cumulative_degree_distribution_summary.csv"),DataFrame;delim=';',decimal=',')
rows=NamedTuple[]; differences=NamedTuple[]
for (guild,tag,col) in [("Consumer","pred",:consumer),("Resource","prey",:resource)]
 fig=Figure(size=(1800,900),fontsize=14)
 Label(fig[0,1:5],"Zero-removal $guild degree distributions: script 29 versus original Galiana tables",fontsize=23)
 for (i,(ds,folder,suffix)) in enumerate(mapping)
  original=CSV.read(joinpath(root,folder,"network_real_$(tag)_$(suffix).csv"),DataFrame)
  original=unique(original[:,2:end]);od=Int.(original.interactions)
  inc=CSV.read(joinpath(root,"All/outputs/42",ds,"42_incidence.csv"),DataFrame;types=Dict(:consumer=>String,:resource=>String))
  regional=unique(inc[:,[:consumer,:resource]]);cd=combine(groupby(regional,col),nrow=>:degree)
  x=collect(1:maximum(cd.degree));y=[mean(cd.degree.>=k) for k in x]
  ox=sort(unique(od));oy=[mean(od.>=k) for k in ox]
  oldmap=Dict(String(r.species)=>Int(r.interactions) for r in eachrow(original));newmap=Dict(String(r[col])=>Int(r.degree) for r in eachrow(cd))
  keysall=union(keys(oldmap),keys(newmap)); ndiff=0
  for node in keysall
   a=get(oldmap,node,-1);b=get(newmap,node,-1)
   if a!=b;ndiff+=1;push!(differences,(dataset=ds,guild=guild,species=node,original_degree=a,current_degree=b));end
  end
  exported=current[(current.dataset.==ds).&(current.guild.==guild).&(current.removal_fraction.==0).&(current.distribution_type.=="All original nodes"),:]
  err=maximum(abs(r.median_cumulative_probability-mean(cd.degree.>=r.degree_threshold)) for r in eachrow(exported))
  @assert err<1e-12
  grid=1:max(maximum(od),maximum(cd.degree));delta=maximum(abs(mean(od.>=k)-mean(cd.degree.>=k)) for k in grid)
  push!(rows,(dataset=ds,guild=guild,original_nodes=length(od),current_nodes=nrow(cd),original_links=sum(od),current_links=sum(cd.degree),identical_degree_multiset=sort(od)==sort(cd.degree),species_degree_differences=ndiff,max_cdf_difference=delta,script29_export_error=err))
  ax=Axis(fig[(i-1)÷5+1,(i-1)%5+1],title=ds,xscale=log10,yscale=log10,xlabel="Number of distinct partners",ylabel="P(degree >= x)")
  lines!(ax,x,y,color=:steelblue,linewidth=2);scatter!(ax,x,y,color=:steelblue,markersize=4)
  lines!(ax,ox,oy,color=:darkorange,linestyle=:dash,linewidth=2);scatter!(ax,ox,oy,color=:darkorange,marker=:rect,markersize=5)
 end
 Legend(fig[3,1:5],[LineElement(color=:steelblue),LineElement(color=:darkorange,linestyle=:dash)], ["Script 29: every integer threshold","Original tables: observed degree values only"],orientation=:horizontal)
 save(joinpath(out,"zero_removal_comparison_$(lowercase(guild)).png"),fig)
end
CSV.write(joinpath(out,"degree_audit.csv"),DataFrame(rows))
if !isempty(differences);CSV.write(joinpath(out,"species_degree_differences.csv"),DataFrame(differences));end
println(DataFrame(rows))
