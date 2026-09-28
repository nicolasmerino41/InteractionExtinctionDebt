using CSV, DataFrames, CairoMakie, Statistics, LinearAlgebra
root=abspath(ARGS[1]); out=joinpath(root,"All/outputs/31_hidden_support_main_figures")
pairs=CSV.read(joinpath(out,"31_3B_fresh_pair_support.csv"),DataFrame)
meta=CSV.read(joinpath(out,"31_3B_fresh_dataset_sizes.csv"),DataFrame)
loss(k,N,m)=k>m ? 0.0 : prod((m-j)/(N-j) for j in 0:k-1;init=1.0)
@assert loss(1,10,5)==0.5
@assert isapprox(loss(2,10,5),2/9)
rankmean(x)=[count(<(v),x)+ (count(==(v),x)+1)/2 for v in x]
states=NamedTuple[]; summaries=NamedTuple[]
for ds in meta.dataset
 d=pairs[pairs.dataset.==ds,:]; real=d[d.full_K.>0,:]; N=only(meta.n_sites[meta.dataset.==ds]); xs=[0.,.1,.2,.4,.6,.8]
 pink=Float64[]
 for x in xs
  keep=max(1,round(Int,N*(1-x))); m=N-keep
  grey=mean(loss(k,N,m) for k in real.full_n); extinct=mean(loss(k,N,m) for k in real.full_K)
  green=1-extinct; pp=extinct-grey;push!(pink,pp)
  @assert pp>=-1e-12 && isapprox(green+pp+grey,1)
  push!(states,(dataset=ds,removal_fraction=x,actual_removed=m,green=green,pink=pp,grey=grey))
 end
 auc=sum(diff(xs).*(pink[1:end-1].+pink[2:end])./2)
 push!(summaries,(dataset=ds,p_emp=sum(d.full_K)/sum(d.full_n),f_regional=nrow(real)/nrow(d),
  mean_pair_conversion=mean(d.full_K./d.full_n),mean_realised_pair_conversion=mean(real.full_K./real.full_n),
  pink_auc_0_08=auc,mean_pink_0_08=auc/.8))
end
s=DataFrame(summaries);sort!(s,:pink_auc_0_08); st=DataFrame(states)
CSV.write(joinpath(out,"31_Figure3B_pink_area_and_conversion.csv"),s)
CSV.write(joinpath(out,"31_Figure3B_exact_state_means.csv"),st)
corr=DataFrame(metric=String[],pearson=Float64[],spearman=Float64[])
for col in [:p_emp,:f_regional,:mean_pair_conversion,:mean_realised_pair_conversion]
 push!(corr,(String(col),cor(s[!,col],s.pink_auc_0_08),cor(rankmean(s[!,col]),rankmean(s.pink_auc_0_08))))
end
CSV.write(joinpath(out,"31_Figure3B_area_correlations.csv"),corr)
fig=Figure(size=(1800,930),fontsize=15)
Label(fig[0,1:5],"Original interaction states after site removal",fontsize=25,font=:bold)
colors=["#999999","#CC79A7","#009E73"]
for (i,r) in enumerate(eachrow(s))
 d=st[st.dataset.==r.dataset,:]; x=d.removal_fraction
 title="$(r.dataset)\np_emp = $(round(100r.p_emp,digits=1))% | f = $(round(100r.f_regional,digits=1))%"
 ax=Axis(fig[(i-1)÷5+1,(i-1)%5+1],title=title,titlesize=17,xlabel="Proportion of sites removed",ylabel="Share of original interactions",xticks=[0,.2,.4,.6,.8],yticks=([0,.5,1],["0%","50%","100%"]))
 band!(ax,x,zeros(nrow(d)),d.grey,color=colors[1]);band!(ax,x,d.grey,d.grey+d.pink,color=colors[2]);band!(ax,x,d.grey+d.pink,ones(nrow(d)),color=colors[3]);xlims!(ax,0,.8);ylims!(ax,0,1)
end
Legend(fig[3,1:5],[PolyElement(color=c) for c in reverse(colors)],
 ["Interaction present","Recorded together; interaction absent","No longer recorded together"],orientation=:horizontal,framevisible=false)
Label(fig[4,1:5],"",fontsize=14)
for ext in ["png","pdf"];save(joinpath(out,"31_Figure3B_link_states_by_dataset."*ext),fig);end
write(joinpath(out,"31_Figure3B_method_note.txt"),"Updated rendering uses Julia/CairoMakie. Expected state fractions are calculated exactly for the same uniform-without-replacement removal design, instead of estimating the means with 500 random subsets. For link support K and co-occurrence support n, pink = P(all K sites removed) - P(all n sites removed). All originally realised regional links receive equal weight. p_emp pools all local pair-site opportunities, including never-realised pairs; f is the regional pair fraction. Alternative unweighted pair means and descriptive correlations across ten datasets are in the CSV files. Area is measured only over the displayed 0-0.8 range. No causal interpretation or fitted homogeneous-model p is implied.")
println(s);println(corr)
