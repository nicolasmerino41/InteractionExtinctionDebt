using CSV, DataFrames, CairoMakie, Statistics, Random


# ============================================================
# RELATIVE VERSION OF FIGURE 3b1 ONLY
#
# Response:
#
#   fraction of original partners lost per species
#
#       (initial degree - retained degree) / initial degree
#
# This is calculated species-by-species BEFORE averaging.
#
# Run:
#
# julia relative_b1.jl [repository root]
#
# ============================================================


const SPECIALIST = "#4477AA"
const GENERALIST = "#CC6677"

const GROUPS = [
    "Specialists",
    "Generalists"
]

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

const NREPS = 500


# ============================================================
# PATHS
# ============================================================

ROOT = pwd()

OUT = joinpath(
    ROOT,
    "All",
    "outputs",
    "43 main figures"
)

INPUT = joinpath(
    OUT,
    "inputs"
)

# ============================================================
# CALCULATE RELATIVE PARTNER LOSS
# ============================================================

rows = NamedTuple[]


for (di, ds) in enumerate(DATA)

    ints =
        CSV.read(
            joinpath(
                INPUT,
                ds * "_interactions.csv"
            ),
            DataFrame;
            types=String
        )


    occ =
        CSV.read(
            joinpath(
                INPUT,
                ds * "_occupancy.csv"
            ),
            DataFrame;
            types=String
        )


    # --------------------------------------------------------
    # Sites
    # --------------------------------------------------------

    sites =
        sort(
            unique(
                vcat(
                    ints.site,
                    occ.site
                )
            )
        )


    si =
        Dict(
            s => i
            for (i, s) in enumerate(sites)
        )


    N =
        length(sites)


    # --------------------------------------------------------
    # Regional interactions
    # --------------------------------------------------------

    links =
        unique(
            collect(
                zip(
                    ints.consumer,
                    ints.resource
                )
            )
        )


    li =
        Dict(
            p => i
            for (i, p) in enumerate(links)
        )


    L =
        length(links)


    # --------------------------------------------------------
    # Nodes
    # --------------------------------------------------------

    nodes =
        vcat(

            [
                ("Consumer", s)
                for s in sort(
                    unique(
                        ints.consumer
                    )
                )
            ],

            [
                ("Resource", s)
                for s in sort(
                    unique(
                        ints.resource
                    )
                )
            ]
        )


    ni =
        Dict(
            s => i
            for (i, s) in enumerate(nodes)
        )


    ends = [
        (
            ni[
                (
                    "Consumer",
                    p[1]
                )
            ],

            ni[
                (
                    "Resource",
                    p[2]
                )
            ]
        )

        for p in links
    ]


    # --------------------------------------------------------
    # Initial link support and degree
    # --------------------------------------------------------

    site_links = [
        Int[]
        for _ in sites
    ]


    K =
        zeros(
            Int,
            L
        )


    degree =
        zeros(
            Int,
            length(nodes)
        )


    for r in eachrow(ints)

        j =
            li[
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
    end


    @assert all(
        degree .> 0
    )


    # ========================================================
    # SPECIALIST / GENERALIST MEMBERSHIPS
    #
    # Exactly the same definition as the current Figure 3:
    # lower/upper halves of DISTINCT degree values,
    # separately within consumer and resource guilds.
    # ========================================================

    memberships =
        Dict{
            Tuple{String,String},
            Vector{Int}
        }()


    for guild in [
        "Consumer",
        "Resource"
    ]

        ids =
            findall(
                x -> x[1] == guild,
                nodes
            )


        values =
            sort(
                unique(
                    degree[ids]
                )
            )


        @assert length(values) > 1


        cut =
            values[
                fld(
                    length(values),
                    2
                )
            ]


        for (group, predicate) in [

            (
                "Specialists",
                k -> k <= cut
            ),

            (
                "Generalists",
                k -> k > cut
            )

        ]

            idx =
                filter(
                    i -> predicate(
                        degree[i]
                    ),
                    ids
                )


            @assert !isempty(idx)


            memberships[
                (
                    guild,
                    group
                )
            ] = idx
        end
    end


    # ========================================================
    # RANDOM SITE REMOVAL
    # ========================================================

    rng =
        MersenneTwister(
            4300 + di
        )


    for rep in 1:NREPS

        order =
            randperm(
                rng,
                N
            )


        remaining =
            copy(K)


        retained =
            copy(degree)


        last =
            0


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


            # ------------------------------------------------
            # Remove sites incrementally
            # ------------------------------------------------

            for pos in last+1:m

                for j in site_links[
                    order[pos]
                ]

                    remaining[j] -= 1


                    # Regional interaction disappears only
                    # when its final supporting site is lost.

                    if remaining[j] == 0

                        a, b =
                            ends[j]

                        retained[a] -= 1
                        retained[b] -= 1
                    end
                end
            end


            last =
                m


            # =================================================
            # RELATIVE LOSS
            #
            # Important:
            #
            # calculate the fraction separately for every
            # species first:
            #
            #   (k0 - kf) / k0
            #
            # and THEN average.
            # =================================================

            for ((guild, group), ids) in memberships

                relative_loss =
                    mean(
                        (
                            degree[ids] .-
                            retained[ids]
                        ) ./
                        degree[ids]
                    )


                @assert (
                    -1e-12 <=
                    relative_loss <=
                    1 + 1e-12
                )


                if f == 0

                    @assert isapprox(
                        relative_loss,
                        0.0
                    )
                end


                push!(
                    rows,
                    (
                        dataset=ds,
                        guild=guild,
                        group=group,
                        removal=f,
                        replicate=rep,
                        relative_links_lost=
                            relative_loss
                    )
                )
            end
        end
    end


    println(
        "Computed relative partner loss: ",
        ds
    )
end


# ============================================================
# DATAFRAME
# ============================================================

r =
    DataFrame(rows)


# ============================================================
# MERGE CONSUMER + RESOURCE GUILDS
#
# First average the two guilds equally within:
#
# dataset × group × removal × replicate
#
# This follows the weighting logic of the existing Figure 3.
# ============================================================

rg =
    combine(
        groupby(
            r,
            [
                :dataset,
                :group,
                :removal,
                :replicate
            ]
        ),

        :relative_links_lost =>
            mean =>
            :relative_links_lost
    )


# ============================================================
# AVERAGE REPLICATES WITHIN DATASET
# ============================================================

rs =
    combine(
        groupby(
            rg,
            [
                :dataset,
                :group,
                :removal
            ]
        ),

        :relative_links_lost =>
            mean =>
            :relative_links_lost
    )


# ============================================================
# SAVE THE VALUES TOO
# ============================================================

CSV.write(
    joinpath(
        OUT,
        "tables",
        "figure3_relative_partner_loss.csv"
    ),
    rs
)


# ============================================================
# THEME
# ============================================================

set_theme!(
    Theme(
        font="Arial",
        fontsize=22,

        Axis=(
            backgroundcolor=:white,

            xgridvisible=false,
            ygridvisible=false,

            topspinevisible=false,
            rightspinevisible=false,

            spinewidth=1.1,

            xtickwidth=1.1,
            ytickwidth=1.1,

            xticklabelsize=19,
            yticklabelsize=19,

            xlabelsize=22,
            ylabelsize=22
        ),

        Legend=(
            framevisible=false,
            labelsize=19
        )
    )
)


# ============================================================
# PLOT ONLY THE RELATIVE B1
# ============================================================

fig =
    Figure(
        size=(900,650),
        figure_padding=(
            75,
            35,
            35,
            35
        )
    )


ticks = (
    [
        0,
        .2,
        .4,
        .6,
        .8
    ],

    [
        "0",
        "20",
        "40",
        "60",
        "80"
    ]
)


yticks = (
    [
        0,
        .2,
        .4,
        .6,
        .8,
        1
    ],

    [
        "0",
        "20",
        "40",
        "60",
        "80",
        "100"
    ]
)


ax =
    Axis(
        fig[1,1],

        xlabel=
            "Sites removed (%)",

        ylabel=
            "Original partners lost (%)",

        xticks=
            ticks,

        yticks=
            yticks
    )


# ============================================================
# DRAW GROUPS
# ============================================================

for (group, color) in zip(
    GROUPS,
    [
        SPECIALIST,
        GENERALIST
    ]
)

    d =
        rs[
            rs.group .== group,
            :
        ]


    # --------------------------------------------------------
    # Individual datasets
    # --------------------------------------------------------

    for ds in DATA

        v =
            d[
                d.dataset .== ds,
                :
            ]


        sort!(
            v,
            :removal
        )


        lines!(
            ax,
            v.removal,
            v.relative_links_lost,

            color=(
                color,
                .18
            ),

            linewidth=
                1.3
        )
    end


    # --------------------------------------------------------
    # Dataset-balanced mean
    # --------------------------------------------------------

    p =
        combine(
            groupby(
                d,
                :removal
            ),

            :relative_links_lost =>
                mean =>
                :relative_links_lost
        )


    sort!(
        p,
        :removal
    )


    lines!(
        ax,
        p.removal,
        p.relative_links_lost,

        color=
            color,

        linewidth=
            3.7
    )
end


# ============================================================
# AXIS LIMITS
# ============================================================

xlims!(
    ax,
    0,
    .8
)


ylims!(
    ax,
    0,
    1.02
)


# ============================================================
# LEGEND
# ============================================================

axislegend(
    ax,

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

    GROUPS;

    position=:lt,
    framevisible=false,
    labelsize=19
)


# ============================================================
# SAVE — PNG ONLY
# ============================================================

save(
    joinpath(
        OUT,
        "Figure3_b1_relative_partner_loss.png"
    ),
    fig;
    px_per_unit=2
)


println(
    "Saved relative B1 test plot to ",
    joinpath(
        OUT,
        "Figure3_b1_relative_partner_loss.png"
    )
)