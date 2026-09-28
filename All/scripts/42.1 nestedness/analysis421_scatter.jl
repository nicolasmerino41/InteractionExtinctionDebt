using CSV, DataFrames, CairoMakie, Statistics

# Run from any directory: julia analysis421_scatter.jl [repository root]
root = isempty(ARGS) ? normpath(joinpath(@__DIR__, "..", "..", "..")) : abspath(ARGS[1])
input = joinpath(root, "All", "outputs", "42")
output = joinpath(root, "All", "outputs", "42.1 nestedness", "scatter")
mkpath(output)
datasets = CSV.read(joinpath(input,"42_dataset_summary.csv"), DataFrame).dataset
all_sites = CSV.read(joinpath(input,"42_site_importance.csv"), DataFrame; types=Dict(:site=>String))

# Fragility uses K+1 as the denominator after removal:
# contribution = 1/K - 1/(K+1).
# K=1 interactions are therefore retained in the fragility index.
impact(k) = (
    fragility=sum((1.0/x - 1.0/(x+1) for x in k); init=0.0),
    lost=count(==(1),k),
    singleton=count(==(2),k)
)

@assert impact([1,2,3,15]).lost == 1
@assert impact([1,2,3,15]).singleton == 1
@assert isapprox(
    impact([1,2,3,15]).fragility,
    1/2 + 1/6 + 1/12 + 1/240
)
@assert impact(Int[]).fragility == 0
@assert impact([1]).fragility > impact([2]).fragility
@assert impact([3]).fragility > impact([15]).fragility

function ranks(x)
    order=sortperm(x); r=zeros(length(x)); i=1
    while i<=length(x)
        j=i
        while j<length(x) && x[order[j+1]]==x[order[i]]; j+=1; end
        r[order[i:j]].=(i+j)/2; i=j+1
    end
    r
end

safe_cor(x,y) = length(x)>1 && std(x)>0 && std(y)>0 ? cor(x,y) : NaN

rows=NamedTuple[]; summaries=NamedTuple[]

for ds in datasets
    inc=CSV.read(
        joinpath(input,ds,"42_incidence.csv"),
        DataFrame;
        types=Dict(:site=>String,:consumer=>String,:resource=>String)
    )
    
    @assert nrow(unique(inc,[:site,:consumer,:resource]))==nrow(inc)
    
    counts=combine(
        groupby(inc,[:consumer,:resource]),
        nrow=>:calculated_support
    )
    
    inc=leftjoin(inc,counts,on=[:consumer,:resource])
    @assert all(inc.support .== inc.calculated_support)
    
    sites=unique(all_sites.site[all_sites.dataset .== ds])
    
    for site in sites
        k=Int.(inc.support[inc.site .== site])
        v=impact(k)
        
        push!(rows,(
            dataset=ds,
            site=site,
            richness=length(k),
            fragility=v.fragility,
            full_loss=v.lost,
            new_singletons=v.singleton,
            fragility_per_interaction=isempty(k) ? NaN : v.fragility/length(k)
        ))
    end
    
    d=DataFrame(filter(r->r.dataset==ds,rows))
    sort!(d,[:richness,:site])
    d.richness_rank=1:nrow(d)
    
    # Each singleton disappears in one removal, each doubleton becomes singleton in two removals.
    @assert sum(d.full_loss)==count(==(1),counts.calculated_support)
    @assert sum(d.new_singletons)==2count(==(2),counts.calculated_support)
    
    # Every interaction is included in fragility, including K=1.
    # An interaction with support K occurs at K sites, so its total contribution
    # across site removals is K * (1/K - 1/(K+1)) = 1/(K+1).
    @assert isapprox(
        sum(d.fragility),
        sum((1.0/(k+1) for k in counts.calculated_support);init=0.0)
    )
    
    folder=joinpath(output,ds)
    mkpath(folder)
    CSV.write(joinpath(folder,"site_impacts.csv"),d)
    
    for metric in [:fragility,:full_loss,:new_singletons]
        push!(summaries,(
            dataset=ds,
            metric=String(metric),
            sites=nrow(d),
            pearson=safe_cor(d.richness,d[!,metric]),
            spearman=safe_cor(ranks(d.richness),ranks(d[!,metric])),
            minimum=minimum(d[!,metric]),
            maximum=maximum(d[!,metric])
        ))
    end
end

CSV.write(joinpath(output,"all_site_impacts.csv"),DataFrame(rows))
CSV.write(joinpath(output,"richness_associations.csv"),DataFrame(summaries))

set_theme!(Theme(
    fontsize=15,
    Axis=(xgridvisible=false,ygridcolor=(:gray,.15),)
))

specs=[
    (:fragility,"421A_fragility","","",""),
    (:full_loss,"421B_full_loss","","",""),
    (:new_singletons,"421C_new_singletons","","","")
]

for (metric,stem,title,ylabel,definition) in specs
    fig=Figure(size=(1800,950))
    Label(fig[0,1:5],title,fontsize=26,font=:bold)
    
    for (i,ds) in enumerate(datasets)
        d=CSV.read(
            joinpath(output,ds,"site_impacts.csv"),
            DataFrame;
            types=Dict(:site=>String)
        )
        
        ax=Axis(
            fig[(i-1)÷5+1,(i-1)%5+1],
            title=ds,
            xlabel="Site interaction richness",
            ylabel=ylabel,
            titlesize=17,
            xlabelsize=13,
            ylabelsize=13,
            xticklabelsize=12,
            yticklabelsize=12
        )
        
        lim=max(1.0,maximum(d.richness)*1.05)
        
        single=Figure(size=(1000,750))
        Label(single[0,1],"$ds - $title",fontsize=23,font=:bold)
        
        a=Axis(
            single[1,1],
            xlabel="Site interaction richness",
            ylabel=ylabel
        )
        
        for axis in [ax,a]
            # lines!(axis,[0,lim],[0,lim],color=:gray45,linestyle=:dash,linewidth=1.5)
            
            scatter!(
                axis,
                d.richness,
                d[!,metric],
                color=(:steelblue,.65),
                markersize=axis===ax ? 13 : 17
            )
            
            # Linear regression (LM) line
            x=Float64.(d.richness)
            y=Float64.(d[!,metric])
            
            if length(x)>1 && std(x)>0
                slope=cor(x,y)*std(y)/std(x)
                intercept=mean(y)-slope*mean(x)
                xfit=range(minimum(x),maximum(x),length=200)
                yfit=intercept .+ slope .* xfit
                
                lines!(
                    axis,
                    xfit,
                    yfit,
                    color=:black,
                    linewidth=2.5
                )
            end
            
            # xlims!(axis,0,lim);ylims!(axis,0,lim)
        end
        
        for ext in ["png"]
            save(joinpath(output,ds,stem*"."*ext),single)
        end
    end
    
    for ext in ["png"]
        save(joinpath(output,stem*"."*ext),fig)
    end
end

open(joinpath(output,"README.md"),"w") do io
    write(io,"""# Single-site removal: support fragility

All ten datasets use the validated site-interaction incidence exported by analysis 42. Support k is the number of distinct sites recording an interaction, not observation frequency. No 50% removal scenario or null model is used.

- A: sum of 1/k - 1/(k+1) over interactions at the removed site. This reciprocal-support index weights scarce support more heavily. Interactions with k = 1 are included: their contribution is 1 - 1/2 = 0.5, so interaction extinctions remain represented in the fragility index.
- B: number of regional interactions lost immediately (k = 1).
- C: number of interactions newly left at one supporting site (k = 2). Existing singleton interactions elsewhere are not included.

All interactions follow the same weighting rule. Raw sums measure total site impact; fragility_per_interaction in the table separates average contribution from richness. Richness associations are descriptive Pearson/Spearman correlations, not proof that richness explains all variation. Constant outcomes have undefined correlations (NaN).

The output root contains three faceted scatter figures; each dataset folder contains corresponding individual figures and site-level tables. Actual richness is used on x, with no jitter or aggregation; exact duplicate points overlap. Each scatter plot includes an ordinary least-squares linear regression line as a descriptive summary.

These figures measure removal effects, not nestedness or containment.

Validation passed: toy examples, distinct incidences, support recalculated from site identities, and identities summing full loss, new singleton counts and fragility across sites.
""")
end

println("Completed all ten datasets: ",output)
