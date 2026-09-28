using CSV,DataFrames,CairoMakie,LinearAlgebra,Statistics
root=abspath(ARGS[1]);input=joinpath(root,"All/outputs/42");out=joinpath(root,"All/outputs/42.1 nestedness/global_nestedness");mkpath(out)
function nodf_rows(M)
 n=size(M,1); totals=vec(sum(M,dims=2));shared=M*transpose(M);v=0.0
 for i in 1:n-1,j in i+1:n
  lo=min(totals[i],totals[j]);hi=max(totals[i],totals[j])
  if 0<lo<hi;v+=shared[i,j]/lo;end
 end
 100v/binomial(n,2)
end
@assert nodf_rows([1 0 0;1 1 0;1 1 1])==100
@assert nodf_rows([1 0 0;0 1 1])==0
@assert nodf_rows([1 1;1 1])==0 # standard NODF ties count as zero
meta=CSV.read(joinpath(input,"42_dataset_summary.csv"),DataFrame)
sites=CSV.read(joinpath(input,"42_site_importance.csv"),DataFrame;types=Dict(:site=>String))
rows=NamedTuple[]
for ds in meta.dataset
 d=CSV.read(joinpath(input,ds,"42_incidence.csv"),DataFrame;types=Dict(:site=>String,:consumer=>String,:resource=>String))
 sn=unique(sites.site[sites.dataset.==ds]);links=unique(collect(zip(d.consumer,d.resource)))
 ri=Dict(s=>i for (i,s) in enumerate(sn));ci=Dict(s=>i for (i,s) in enumerate(links));M=zeros(Int,length(sn),length(links))
 for r in eachrow(d);M[ri[r.site],ci[(r.consumer,r.resource)]]=1;end
 a=nodf_rows(M);b=nodf_rows(Matrix(transpose(M)));nr=binomial(size(M,1),2);nc=binomial(size(M,2),2)
 push!(rows,(dataset=ds,sites=size(M,1),links=size(M,2),site_NODF=a,link_NODF=b,whole_matrix_NODF=(nr*a+nc*b)/(nr+nc)))
end
s=DataFrame(rows);sort!(s,:site_NODF);CSV.write(joinpath(out,"global_nestedness.csv"),s)
fig=Figure(size=(1100,650),fontsize=17);ax=Axis(fig[1,1],ylabel="Overall site nestedness (NODF, 0-100)",xticks=(1:nrow(s),s.dataset),xticklabelrotation=pi/4)
barplot!(ax,1:nrow(s),s.site_NODF,color=:steelblue);ylims!(ax,0,100)
text!(ax,1:nrow(s),s.site_NODF.+2;text=string.(round.(s.site_NODF,digits=1)),align=(:center,:bottom),fontsize=15)
Label(fig[2,1],"One score per dataset, across the complete site-by-interaction matrix. Standard NODF gives equal-richness pairs zero contribution.\nDescriptive scores, not a null-model test. The CSV also reports interaction-column and combined whole-matrix NODF.",fontsize=13)
for ext in ["png","pdf"];save(joinpath(out,"global_site_nestedness."*ext),fig);end
write(joinpath(out,"README.md"),"# Overall site nestedness\nRows are sites and columns are realised regional interactions. Site NODF is 100 times the sum of poorer-row overlap fractions for strictly unequal, nonempty row pairs, divided by ALL site pairs. It provides one dataset-level site-nestedness score. Link-column NODF and full matrix NODF are supplied separately; the latter combines both axes and may be dominated by the more numerous interaction columns. This is the standard overlap-and-decreasing-fill metric, not Baselga nestedness-resultant dissimilarity. Unlike the earlier mean containment, its denominator includes tied-richness pairs (zero contribution). NODF still aggregates pair overlaps internally; a global score does not erase their contribution. No raw-score threshold establishes statistical significance. Reference: Almeida-Neto et al. 2008, doi:10.1111/j.0030-1299.2008.16644.x.\n")
println(s)
