module RevisedFigure3

using CSV, DataFrames, CairoMakie, Statistics, Random

const SPECIALIST = "#4477AA"
const GENERALIST = "#CC6677"

const GROUPS = ["Specialists", "Generalists"]

const DATA = [
    "Quercus",
    "Nahuel",
    "Salix_Galpar",
    "Gottin_HP",
    "Gottin_PP",
    "Garraf_HP",
    "Garraf_PP",
    "Garraf_PP2",
    "Olot",
    "Montseny"
]

const LEVELS = collect(0.0:0.1:0.8)


# ============================================================
# CALCULATE
# ============================================================

function calculate(root, out)

    baseline = NamedTuple[]
    repout = NamedTuple[]
    groupinfo = NamedTuple[]

    for (di, ds) in enumerate(DATA)

        ints = CSV.read(
            joinpath(out, "inputs", ds * "_interactions.csv"),
            DataFrame;
            types=String
        )

        occ = CSV.read(
            joinpath(out, "inputs", ds * "_occupancy.csv"),
            DataFrame;
            types=String
        )

        sites = sort(unique(vcat(ints.site, occ.site)))
        si = Dict(s => i for (i, s) in enumerate(sites))
        N = length(sites)

        links = unique(
            collect(
                zip(
                    ints.consumer,
                    ints.resource
                )
            )
        )

        li = Dict(
            p => i
            for (i, p) in enumerate(links)
        )

        L = length(links)

        nodes = vcat(
            [
                ("Consumer", s)
                for s in sort(unique(ints.consumer))
            ],
            [
                ("Resource", s)
                for s in sort(unique(ints.resource))
            ]
        )

        ni = Dict(
            s => i
            for (i, s) in enumerate(nodes)
        )

        ends = [
            (
                ni[("Consumer", p[1])],
                ni[("Resource", p[2])]
            )
            for p in links
        ]

        site_links = [
            Int[]
            for _ in sites
        ]

        K = zeros(Int, L)
        degree = zeros(Int, length(nodes))
        support = zeros(length(nodes))

        for r in eachrow(ints)

            j = li[
                (
                    r.consumer,
                    r.resource
                )
            ]

            push!(
                site_links[
                    si[r.site]
                ],
                j
            )

            K[j] += 1
        end

        for (j, (a, b)) in enumerate(ends)

            degree[a] += 1
            degree[b] += 1

            support[a] += K[j]
            support[b] += K[j]
        end

        support ./= degree

        memberships =
            Dict{
                Tuple{String,String},
                Vector{Int}
            }()

        for guild in ["Consumer", "Resource"]

            ids = findall(
                x -> x[1] == guild,
                nodes
            )

            values = sort(
                unique(
                    degree[ids]
                )
            )

            @assert length(values) > 1

            # Same grouping as scripts 36 and 36a:
            # lower/upper halves of distinct degree values.

            cut =
                values[
                    fld(
                        length(values),
                        2
                    )
                ]

            for (group, predicate) in [
                ("Specialists", k -> k <= cut),
                ("Generalists", k -> k > cut)
            ]

                idx = filter(
                    i -> predicate(degree[i]),
                    ids
                )

                @assert !isempty(idx)

                memberships[
                    (
                        guild,
                        group
                    )
                ] = idx

                push!(
                    baseline,
                    (
                        dataset=ds,
                        guild=guild,
                        group=group,
                        support=mean(support[idx])
                    )
                )

                push!(
                    groupinfo,
                    (
                        dataset=ds,
                        guild=guild,
                        group=group,
                        n_species=length(idx),
                        minimum_degree=minimum(degree[idx]),
                        maximum_degree=maximum(degree[idx])
                    )
                )
            end
        end

        rng = MersenneTwister(
            4300 + di
        )

        for rep in 1:500

            order = randperm(rng, N)
            remaining = copy(K)
            last = 0
            retained = copy(degree)

            for f in LEVELS

                m =
                    N -
                    max(
                        1,
                        round(
                            Int,
                            N * (1-f)
                        )
                    )

                for pos in last+1:m

                    for j in site_links[
                        order[pos]
                    ]

                        remaining[j] -= 1

                        if remaining[j] == 0

                            a, b = ends[j]

                            retained[a] -= 1
                            retained[b] -= 1
                        end
                    end
                end

                last = m

                for ((guild, group), ids) in memberships

                    lost = mean(
                        degree[ids] .-
                        retained[ids]
                    )

                    active = mean(
                        retained[ids] .> 0
                    )

                    @assert lost >= 0
                    @assert 0 <= active <= 1

                    if f == 0
                        @assert lost == 0
                        @assert active == 1
                    end

                    push!(
                        repout,
                        (
                            dataset=ds,
                            guild=guild,
                            group=group,
                            removal=f,
                            replicate=rep,
                            links_lost=lost,
                            active_fraction=active
                        )
                    )
                end
            end
        end

        println(
            "Figure 3 group analysis: ",
            ds
        )
    end

    a = DataFrame(baseline)
    r = DataFrame(repout)

    # Average species within each guild, then weight
    # the two guilds equally and datasets equally.

    ag = combine(
        groupby(
            a,
            [:dataset, :group]
        ),
        :support => mean => :support
    )

    rg = combine(
        groupby(
            r,
            [
                :dataset,
                :group,
                :removal,
                :replicate
            ]
        ),
        :links_lost => mean => :links_lost,
        :active_fraction => mean => :active_fraction
    )

    rs = combine(
        groupby(
            rg,
            [
                :dataset,
                :group,
                :removal
            ]
        ),
        :links_lost => mean => :links_lost,
        :active_fraction => mean => :active_fraction
    )

    for (name, t) in [
        (
            "figure3_support_by_guild",
            a
        ),
        (
            "figure3_support_merged",
            ag
        ),
        (
            "figure3_group_definitions",
            DataFrame(groupinfo)
        ),
        (
            "figure3_group_replicates",
            r
        ),
        (
            "figure3_group_curves",
            rs
        )
    ]

        CSV.write(
            joinpath(
                out,
                "tables",
                name * ".csv"
            ),
            t
        )
    end

    ag, rs
end


# ============================================================
# RENDER
# ============================================================

function render(root; plot_only=false)

    out = joinpath(
        root,
        "All/outputs/43 main figures"
    )

    if plot_only

        ag = CSV.read(
            joinpath(
                out,
                "tables/figure3_support_merged.csv"
            ),
            DataFrame
        )

        rs = CSV.read(
            joinpath(
                out,
                "tables/figure3_group_curves.csv"
            ),
            DataFrame
        )

    else

        ag, rs = calculate(
            root,
            out
        )
    end

    dd = CSV.read(
        joinpath(
            out,
            "tables/degree_example.csv"
        ),
        DataFrame
    )


    # ========================================================
    # THEME
    # ========================================================

    set_theme!(
        Theme(
            font="Arial",
            fontsize=22,
            Axis=(
                xgridvisible=false,
                ygridvisible=false,
                topspinevisible=false,
                rightspinevisible=false,
                xticklabelsize=19,
                yticklabelsize=19,
                xlabelsize=22,
                ylabelsize=22,
                spinewidth=1.1
            )
        )
    )


    # ========================================================
    # FIGURE
    #
    # TRUE 1 × 3 STRUCTURE:
    #
    # A | B | C
    #
    # B internally contains two stacked axes.
    # C occupies one normal figure column.
    # ========================================================

    fig = Figure(
        size=(2050, 850),
        figure_padding=(75,35,35,35)
    )

    left = fig[1,1] = GridLayout()
    middle = fig[1,2] = GridLayout()
    right = fig[1,3] = GridLayout()


    # Wider C than before.
    #
    # Old:
    # A = .29
    # B = .44
    # C = .27
    #
    # New:
    # A = .25
    # B = .40
    # C = .35
    #
    # C therefore gets substantially more horizontal space.

    colsize!(
        fig.layout,
        1,
        Relative(.25)
    )

    colsize!(
        fig.layout,
        2,
        Relative(.40)
    )

    colsize!(
        fig.layout,
        3,
        Relative(.35)
    )


    # A-B separation can remain comfortable.
    colgap!(
        fig.layout,
        1,
        55
    )

    # Bring C much closer to B.
    colgap!(
        fig.layout,
        2,
        28
    )


    # ========================================================
    # PANEL A
    # ========================================================

    a = Axis(
        left[1,1],
        ylabel=
            "Mean supporting sites per realised partner",
        xticks=(
            [1,2],
            GROUPS
        ),
        limits=(
            .65,
            2.35,
            0,
            nothing
        )
    )

    for ds in DATA

        d = ag[
            ag.dataset .== ds,
            :
        ]

        v = [
            only(
                d.support[
                    d.group .== g
                ]
            )
            for g in GROUPS
        ]

        lines!(
            a,
            [1,2],
            v,
            color=(:gray50,.28),
            linewidth=1.5
        )

        for i in 1:2

            scatter!(
                a,
                [i],
                [v[i]],
                color=(
                    [
                        SPECIALIST,
                        GENERALIST
                    ][i],
                    .5
                ),
                markersize=10
            )
        end
    end

    means = [
        mean(
            ag.support[
                ag.group .== g
            ]
        )
        for g in GROUPS
    ]

    lines!(
        a,
        [1,2],
        means,
        color=:gray25,
        linewidth=2.5
    )

    for i in 1:2

        v = ag.support[
            ag.group .== GROUPS[i]
        ]

        q = quantile(
            v,
            [.25,.75]
        )

        color = [
            SPECIALIST,
            GENERALIST
        ][i]

        rangebars!(
            a,
            [i],
            [q[1]],
            [q[2]],
            color=color,
            linewidth=4,
            whiskerwidth=16
        )

        scatter!(
            a,
            [i],
            [means[i]],
            color=color,
            markersize=19,
            strokecolor=:white,
            strokewidth=1.5
        )
    end


    # ========================================================
    # PANEL B
    # ========================================================

    ticks = (
        [0,.2,.4,.6,.8],
        ["0","20","40","60","80"]
    )

    b1 = Axis(
        middle[1,1],
        ylabel="Links lost per original species",
        xticks=ticks
    )

    b2 = Axis(
        middle[2,1],
        xlabel="Sites removed (%)",
        ylabel="Species retaining interactions (%)",
        xticks=ticks,
        yticks=(
            [0,.25,.5,.75,1],
            ["0","25","50","75","100"]
        )
    )

    hidexdecorations!(
        b1;
        grid=false
    )

    linkxaxes!(
        b1,
        b2
    )

    for (group, color) in zip(
        GROUPS,
        [
            SPECIALIST,
            GENERALIST
        ]
    )

        d = rs[
            rs.group .== group,
            :
        ]

        for ds in DATA

            v = d[
                d.dataset .== ds,
                :
            ]

            sort!(
                v,
                :removal
            )

            lines!(
                b1,
                v.removal,
                v.links_lost,
                color=(color,.18),
                linewidth=1.3
            )

            lines!(
                b2,
                v.removal,
                v.active_fraction,
                color=(color,.18),
                linewidth=1.3
            )
        end

        p = combine(
            groupby(
                d,
                :removal
            ),
            :links_lost => mean => :links_lost,
            :active_fraction => mean => :active_fraction
        )

        sort!(
            p,
            :removal
        )

        lines!(
            b1,
            p.removal,
            p.links_lost,
            color=color,
            linewidth=3.7
        )

        lines!(
            b2,
            p.removal,
            p.active_fraction,
            color=color,
            linewidth=3.7
        )
    end

    xlims!(
        b1,
        0,
        .8
    )

    xlims!(
        b2,
        0,
        .8
    )

    ylims!(
        b1,
        0,
        nothing
    )

    ylims!(
        b2,
        0,
        1.02
    )

    rowgap!(
        middle,
        25
    )


    # ========================================================
    # PANEL C
    #
    # Wide landscape panel occupying the FULL third column.
    #
    # No extra empty rows.
    # No external vertical legend.
    # ========================================================

    c = Axis(
        right[1,1],

        xlabel="Degree",
        ylabel="Cumulative probability",

        xscale=log10,
        yscale=log10,

        xlabelsize=19,
        ylabelsize=19,

        xticklabelsize=17,
        yticklabelsize=17,

        xticks=(
            [1,2,5,10,20,50],
            ["1","2","5","10","20","50"]
        ),

        yticks=(
            [.01,.1,1],
            ["0.01","0.1","1"]
        )
    )

    colors = [
        "#88CCEE",
        "#4477AA",
        "#332288"
    ]

    for (j, f) in enumerate(
        [0.,.4,.8]
    )

        v = dd[
            dd.removal .== f,
            :
        ]

        sort!(
            v,
            :degree
        )

        y = v.probability

        keep = [
            i
            for i in eachindex(y)
            if (
                y[i] > 0 &&
                (
                    i == length(y) ||
                    y[i] > y[i+1] + 1e-12
                )
            )
        ]

        lines!(
            c,
            v.degree[keep],
            y[keep],
            color=colors[j],
            linewidth=2.5
        )
    end


    # ========================================================
    # C LEGEND
    #
    # Horizontal and INSIDE C.
    #
    # This avoids consuming any extra vertical layout space.
    # ========================================================

    axislegend(
        c,

        [
            LineElement(
                color=x,
                linewidth=3
            )
            for x in colors
        ],

        [
            "0% removed",
            "40% removed",
            "80% removed"
        ];

        position=:lb,
        orientation=:horizontal,
        framevisible=false,
        labelsize=16,
        patchsize=(20,12)
    )


    # ========================================================
    # PANEL LABELS
    # ========================================================

    for (slot, letter) in [

        (
            left[
                1,
                1,
                TopLeft()
            ],
            "a"
        ),

        (
            middle[
                1,
                1,
                TopLeft()
            ],
            "b"
        ),

        (
            right[
                1,
                1,
                TopLeft()
            ],
            "c"
        )

    ]

        Label(
            slot,
            letter,
            font=:bold,
            fontsize=28,
            halign=:left,
            valign=:top,
            padding=(-42,0,10,0),
            tellwidth=false,
            tellheight=false
        )
    end


    # ========================================================
    # SPECIALIST / GENERALIST LEGEND
    # ========================================================

    Legend(
        fig[2,1:3],

        [
            LineElement(
                color=SPECIALIST,
                linewidth=4
            ),

            LineElement(
                color=GENERALIST,
                linewidth=4
            )
        ],

        GROUPS,

        orientation=:horizontal,
        framevisible=false,
        labelsize=21
    )


    # Keep the bottom legend compact.

    rowsize!(
        fig.layout,
        2,
        Auto(0.08)
    )

    rowgap!(
        fig.layout,
        15
    )


    # ========================================================
    # SAVE — PNG ONLY
    # ========================================================

    save(
        joinpath(
            out,
            "Figure3_remnant_network.png"
        ),
        fig;
        px_per_unit=2
    )


    # ========================================================
    # CAPTION
    # ========================================================

    write(
        joinpath(
            out,
            "FIGURE3_REVISED_CAPTION.md"
        ),

"""Figure 3. Generalist–specialist differences in local support and response to site removal. (a) Mean supporting-site count per realised partner, first averaged across species within each dataset, guild and degree group. Consumer and resource means receive equal weight within datasets. Faint points and connecting lines represent ten datasets; large points are dataset-balanced means, and vertical bars show the across-dataset interquartile range. (b, upper) Absolute number of regional partners lost per original species, including species that lose all partners. (b, lower) Fraction of original species retaining at least one recorded interaction. Thin lines show dataset means, thick lines equally weighted means across datasets, after averaging the two guilds equally. These are network-active species, not independent demographic survival observations. 500 uniform-random site-removal permutations per dataset; seed 4300 + dataset index. Groups are fixed at zero removal. As in scripts 36/36a, specialists occupy the lower half of distinct initial degree values within each dataset and guild, and generalists the upper half; this is not a median split of species counts. Definitions and sample sizes are saved in tables/figure3_group_definitions.csv. Guilds are merged for presentation, not species identities; each link contributes to its consumer and resource endpoint. (c) Smaller illustrative consumer degree-distribution example from Salix–Galpar, among active consumers, retaining the existing corrected 43 input and Galiana-style change-point rendering. No claim of a guild comparison is made.
"""
    )


    println(
        "Updated Figure 3 only; Figure 2 untouched"
    )
end


end


# ============================================================
# RUN
# ============================================================

RevisedFigure3.render(
    isempty(ARGS) ?
        normpath(
            joinpath(
                @__DIR__,
                "..",
                "..",
                ".."
            )
        ) :
        abspath(
            ARGS[1]
        );

    plot_only=
        ("--plot-only" in ARGS)
)