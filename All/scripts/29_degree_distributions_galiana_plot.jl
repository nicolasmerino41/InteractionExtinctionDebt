using CSV, DataFrames, CairoMakie, Statistics


# ============================================================
# PATHS
# ============================================================

root = isempty(ARGS) ?
    normpath(joinpath(@__DIR__, "..", "..")) :
    abspath(ARGS[1])

out = joinpath(root, "All", "CombinedOutputs")

oldfile = joinpath(
    out,
    "29_cumulative_degree_distribution_summary.before_grid_correction.csv"
)

nodefile = joinpath(
    out,
    "29_node_degree_transitions.csv"
)

@assert isfile(oldfile) "Missing old summary file: $oldfile"
@assert isfile(nodefile) "Missing node transition file: $nodefile"


# ============================================================
# READ DATA
# ============================================================

old = CSV.read(
    oldfile,
    DataFrame;
    delim=';',
    decimal=','
)

nodes = CSV.read(
    nodefile,
    DataFrame;
    delim=';',
    decimal=',',
    select=[
        :dataset,
        :guild,
        :removal_fraction,
        :replicate,
        :initial_degree,
        :retained_degree
    ]
)


# ============================================================
# REBUILD CUMULATIVE DEGREE DISTRIBUTIONS
# ============================================================

summaryrows = NamedTuple[]

for group in groupby(
    nodes,
    [:dataset, :guild, :removal_fraction]
)

    maxdeg = maximum(group.initial_degree)
    reps = groupby(group, :replicate)

    allprob = zeros(
        length(reps),
        maxdeg
    )

    activeprob = fill(
        NaN,
        length(reps),
        maxdeg
    )

    active = Int[]


    for (i, r) in enumerate(reps)

        hist = zeros(
            Int,
            maxdeg
        )

        for k in r.retained_degree

            if k > 0
                hist[k] += 1
            end
        end


        tail = reverse(
            cumsum(
                reverse(hist)
            )
        )


        nactive = sum(hist)

        push!(
            active,
            nactive
        )


        # All original nodes
        allprob[i, :] =
            tail ./ nrow(r)


        # Active nodes only
        if nactive > 0

            activeprob[i, :] =
                tail ./ nactive
        end
    end


    for (kind, prob) in [

        (
            "All original nodes",
            allprob
        ),

        (
            "Active nodes only",
            activeprob
        )

    ]

        for k in 1:maxdeg

            v = filter(
                isfinite,
                prob[:, k]
            )

            isempty(v) && continue


            push!(
                summaryrows,
                (
                    dataset =
                        first(group.dataset),

                    guild =
                        first(group.guild),

                    removal_fraction =
                        first(group.removal_fraction),

                    distribution_type =
                        kind,

                    degree_threshold =
                        k,

                    median_cumulative_probability =
                        median(v),

                    q025 =
                        quantile(v, 0.025),

                    q975 =
                        quantile(v, 0.975),

                    number_of_original_nodes =
                        nrow(first(reps)),

                    median_number_of_active_nodes =
                        median(active)
                )
            )
        end
    end
end


d = DataFrame(summaryrows)


# ============================================================
# VERIFY ZERO-REMOVAL VALUES
# ============================================================

for r in eachrow(
    old[old.removal_fraction .== 0, :]
)

    v = d[
        (d.dataset .== r.dataset) .&
        (d.guild .== r.guild) .&
        (d.removal_fraction .== 0) .&
        (d.distribution_type .== r.distribution_type) .&
        (d.degree_threshold .== r.degree_threshold),
        :
    ]


    @assert isapprox(
        only(
            v.median_cumulative_probability
        ),
        r.median_cumulative_probability;
        atol=1e-12
    )
end


# ============================================================
# SAVE CORRECTED SUMMARY
# ============================================================

CSV.write(
    joinpath(
        out,
        "29_cumulative_degree_distribution_summary.csv"
    ),
    d;
    delim=';',
    decimal=','
)


# ============================================================
# DATASET ORDER
#
# Galiana order:
#
# PP:
# Garraf PP
# Garraf PP2
# Montseny
# Gottin PP
# Nahuel
#
# HP:
# Garraf HP
# Quercus
# Olot
# Gottin HP
# Galpar
# ============================================================

datasets = [

    "Garraf_PP",
    "Garraf_PP2",
    "Montseny",
    "Gottin_PP",
    "Nahuel",

    "Garraf_HP",
    "Quercus",
    "Olot",
    "Gottin_HP",
    "Salix_Galpar"
]


# ============================================================
# DISPLAY NAMES
# ============================================================

function display_name(ds)

    if ds == "Salix_Galpar"
        return "Galpar"
    end

    return replace(
        ds,
        "_" => " "
    )
end


# ============================================================
# REMOVAL LEVELS
# ============================================================

levels = [
    0.0,
    0.4,
    0.8
]


colors = [
    "#1b9e77",
    "#d95f02",
    "#7570b3"
]


# ============================================================
# Y TICKS
#
# These control ONLY which labels are printed.
#
# They do NOT control the physical limits of the axis.
# ============================================================


# ------------------------------------------------------------
# CONSUMERS
# ------------------------------------------------------------

consumer_yticks = Dict(

    "Garraf_PP" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Garraf_PP2" =>
        (
            [0.03, 0.10, 0.30, 1.00],
            ["0.03", "0.10", "0.30", "1.00"]
        ),

    "Montseny" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Gottin_PP" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Nahuel" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Garraf_HP" =>
        (
            [0.03, 0.10, 0.30, 1.00],
            ["0.03", "0.10", "0.30", "1.00"]
        ),

    "Quercus" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        ),

    "Olot" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        ),

    "Gottin_HP" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Salix_Galpar" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        )
)


# ------------------------------------------------------------
# RESOURCES
# ------------------------------------------------------------

resource_yticks = Dict(

    "Garraf_PP" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Garraf_PP2" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Montseny" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        ),

    "Gottin_PP" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Nahuel" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        ),

    "Garraf_HP" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        ),

    "Quercus" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        ),

    "Olot" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Gottin_HP" =>
        (
            [0.01, 0.10, 1.00],
            ["0.01", "0.10", "1.00"]
        ),

    "Salix_Galpar" =>
        (
            [0.10, 0.30, 1.00],
            ["0.1", "0.3", "1.0"]
        )
)


# ============================================================
# PANEL-SPECIFIC Y LIMITS
#
# IMPORTANT:
#
# These are independent from the tick positions.
#
# They are deliberately panel-specific rather than calculated
# from the lowest tick.
#
# This avoids forcing every dataset with the same tick pattern
# into exactly the same amount of empty log-space.
#
# Upper limits are slightly above 1 because 1.0 is a tick
# inside the plotting area, not the upper boundary.
# ============================================================


# ------------------------------------------------------------
# CONSUMERS
# ------------------------------------------------------------

consumer_ylims = Dict(

    # Plant-pollinator

    "Garraf_PP" =>
        (0.0045, 1.28),

    "Garraf_PP2" =>
        (0.008, 1.28),

    "Montseny" =>
        (0.0045, 1.28),

    "Gottin_PP" =>
        (0.0055, 1.28),

    "Nahuel" =>
        (0.007, 1.28),


    # Host-parasite

    "Garraf_HP" =>
        (0.0190, 1.28),

    "Quercus" =>
        (0.0200, 1.28),

    "Olot" =>
        (0.0250, 1.28),

    "Gottin_HP" =>
        (0.035, 1.28),

    "Salix_Galpar" =>
        (0.0055, 1.28)
)


# ------------------------------------------------------------
# RESOURCES
# ------------------------------------------------------------

resource_ylims = Dict(

    # Plant-pollinator

    "Garraf_PP" =>
        (0.0050, 1.28),

    "Garraf_PP2" =>
        (0.0045, 1.28),

    "Montseny" =>
        (0.0550, 1.28),

    "Gottin_PP" =>
        (0.0055, 1.28),

    "Nahuel" =>
        (0.0550, 1.28),


    # Host-parasite

    "Garraf_HP" =>
        (0.0550, 1.28),

    "Quercus" =>
        (0.0500, 1.28),

    "Olot" =>
        (0.0050, 1.28),

    "Gottin_HP" =>
        (0.0045, 1.28),

    "Salix_Galpar" =>
        (0.0550, 1.28)
)


# ============================================================
# GET PANEL TICKS
# ============================================================

function panel_yticks(
    ds,
    guild
)

    if guild == "Consumer"

        return consumer_yticks[ds]

    elseif guild == "Resource"

        return resource_yticks[ds]

    else

        error(
            "Unknown guild: $guild"
        )
    end
end


# ============================================================
# GET PANEL Y LIMITS
# ============================================================

function panel_ylims(
    ds,
    guild
)

    if guild == "Consumer"

        return consumer_ylims[ds]

    elseif guild == "Resource"

        return resource_ylims[ds]

    else

        error(
            "Unknown guild: $guild"
        )
    end
end


# ============================================================
# SET PANEL Y LIMITS
# ============================================================

function set_panel_ylims!(
    ax,
    ds,
    guild
)

    yl = panel_ylims(
        ds,
        guild
    )

    ylims!(
        ax,
        yl[1],
        yl[2]
    )
end


# ============================================================
# CCDF KNOTS
# ============================================================

function knots(y)

    [
        i for i in eachindex(y)

        if y[i] > 0 &&

        (
            i == length(y) ||

            y[i] >
            y[i+1] + 1e-12
        )
    ]
end


@assert knots(
    [
        1.0,
        0.75,
        0.25,
        0.25,
        0.25
    ]
) == [
    1,
    2,
    5
]


# ============================================================
# DRAW ONE PANEL
# ============================================================

function draw_panel!(
    ax,
    d,
    ds,
    guild,
    kind,
    levels,
    colors;
    export_to=nothing
)

    for (j, level) in enumerate(
        levels
    )

        a = d[
            (d.dataset .== ds) .&
            (d.guild .== guild) .&
            (d.distribution_type .== kind) .&
            (d.removal_fraction .== level),
            :
        ]


        sort!(
            a,
            :degree_threshold
        )


        isempty(a) &&
            continue


        @assert all(
            diff(
                a.median_cumulative_probability
            ) .<= 1e-10
        )


        a = a[
            knots(
                a.median_cumulative_probability
            ),
            :
        ]


        isempty(a) &&
            continue


        if export_to !== nothing

            append!(
                export_to,
                a
            )
        end


        # Zero cannot be displayed on a log axis.

        lo = map(
            x ->
                x > 0 ?
                x :
                NaN,

            a.q025
        )


        hi = map(
            x ->
                x > 0 ?
                x :
                NaN,

            a.q975
        )


        band!(
            ax,

            a.degree_threshold,

            lo,
            hi,

            color=(
                colors[j],
                0.13
            )
        )


        lines!(
            ax,

            a.degree_threshold,

            a.median_cumulative_probability,

            color=
                colors[j],

            linewidth=2
        )


        scatter!(
            ax,

            a.degree_threshold,

            a.median_cumulative_probability,

            color=
                colors[j],

            markersize=5
        )
    end
end


# ============================================================
# EXPORT TABLE
# ============================================================

exported = DataFrame()


# ============================================================
# FIGURES 1 AND 2
#
# COMBINED FIGURES
#
# Row 1 = Consumers
# Row 2 = Resources
#
# 10 columns = datasets in Galiana order
# ============================================================

for (kind, stem) in [

    (
        "All original nodes",
        "all_nodes"
    ),

    (
        "Active nodes only",
        "active_nodes"
    )

]

    fig = Figure(
        size=(2300, 850),
        fontsize=13
    )


    Label(
        fig[0, 1:10],

        "Cumulative degree distributions - $kind",

        fontsize=23
    )


    for (i, ds) in enumerate(
        datasets
    )

        for (row, guild) in enumerate(
            [
                "Consumer",
                "Resource"
            ]
        )

            yt = panel_yticks(
                ds,
                guild
            )


            ax = Axis(
                fig[row, i],

                title =
                    "$(display_name(ds))\n$guild",

                titlesize=13,

                xscale=log10,

                yscale=log10,

                yticks=yt,

                xlabel=
                    "Degree",

                ylabel=
                    i == 1 ?
                    "P(degree ≥ x)" :
                    ""
            )


            # -----------------------------------------------
            # DATASET + GUILD SPECIFIC LIMITS
            # -----------------------------------------------

            set_panel_ylims!(
                ax,
                ds,
                guild
            )


            draw_panel!(
                ax,
                d,
                ds,
                guild,
                kind,
                levels,
                colors;
                export_to=exported
            )
        end
    end


    Legend(
        fig[3, 1:10],

        [
            LineElement(
                color=c,
                linewidth=2
            )

            for c in colors
        ],

        [
            "0% removed",
            "40% removed",
            "80% removed"
        ],

        orientation=:horizontal,

        framevisible=false
    )


    save(
        joinpath(
            out,

            "29_cumulative_degree_distributions_$(stem).png"
        ),

        fig
    )
end


# ============================================================
# FIGURE 3
#
# ACTIVE CONSUMERS ONLY
#
# 2 × 5
# ============================================================

kind =
    "Active nodes only"

guild =
    "Consumer"


fig_consumers = Figure(
    size=(1500, 800),
    fontsize=15
)


Label(
    fig_consumers[0, 1:5],

    "Cumulative degree distributions - Active consumers",

    fontsize=23
)


for (i, ds) in enumerate(
    datasets
)

    plotrow =
        (i - 1) ÷ 5 + 1

    plotcol =
        (i - 1) % 5 + 1


    yt = panel_yticks(
        ds,
        guild
    )


    ax = Axis(
        fig_consumers[
            plotrow,
            plotcol
        ],

        title=
            display_name(ds),

        titlesize=15,

        xscale=log10,

        yscale=log10,

        yticks=yt,

        xlabel=
            plotrow == 2 ?
            "Degree" :
            "",

        ylabel=
            plotcol == 1 ?
            "P(degree ≥ x)" :
            "",

        xlabelsize=14,

        ylabelsize=14,

        xticklabelsize=12,

        yticklabelsize=12
    )


    # --------------------------------------------------------
    # EACH DATASET HAS ITS OWN Y RANGE
    # --------------------------------------------------------

    set_panel_ylims!(
        ax,
        ds,
        guild
    )


    draw_panel!(
        ax,
        d,
        ds,
        guild,
        kind,
        levels,
        colors
    )
end


Legend(
    fig_consumers[3, 1:5],

    [
        LineElement(
            color=c,
            linewidth=2
        )

        for c in colors
    ],

    [
        "0% removed",
        "40% removed",
        "80% removed"
    ],

    orientation=:horizontal,

    framevisible=false
)


save(
    joinpath(
        out,

        "29_cumulative_degree_distributions_active_nodes_consumers.png"
    ),

    fig_consumers
)


# ============================================================
# FIGURE 4
#
# ACTIVE RESOURCES ONLY
#
# 2 × 5
# ============================================================

kind =
    "Active nodes only"

guild =
    "Resource"


fig_resources = Figure(
    size=(1500, 800),
    fontsize=15
)


Label(
    fig_resources[0, 1:5],

    "Cumulative degree distributions - Active resources",

    fontsize=23
)


for (i, ds) in enumerate(
    datasets
)

    plotrow =
        (i - 1) ÷ 5 + 1

    plotcol =
        (i - 1) % 5 + 1


    yt = panel_yticks(
        ds,
        guild
    )


    ax = Axis(
        fig_resources[
            plotrow,
            plotcol
        ],

        title=
            display_name(ds),

        titlesize=15,

        xscale=log10,

        yscale=log10,

        yticks=yt,

        xlabel=
            plotrow == 2 ?
            "Degree" :
            "",

        ylabel=
            plotcol == 1 ?
            "P(degree ≥ x)" :
            "",

        xlabelsize=14,

        ylabelsize=14,

        xticklabelsize=12,

        yticklabelsize=12
    )


    # --------------------------------------------------------
    # EACH DATASET HAS ITS OWN Y RANGE
    # --------------------------------------------------------

    set_panel_ylims!(
        ax,
        ds,
        guild
    )


    draw_panel!(
        ax,
        d,
        ds,
        guild,
        kind,
        levels,
        colors
    )
end


Legend(
    fig_resources[3, 1:5],

    [
        LineElement(
            color=c,
            linewidth=2
        )

        for c in colors
    ],

    [
        "0% removed",
        "40% removed",
        "80% removed"
    ],

    orientation=:horizontal,

    framevisible=false
)


save(
    joinpath(
        out,

        "29_cumulative_degree_distributions_active_nodes_resources.png"
    ),

    fig_resources
)


# ============================================================
# SAVE PLOTTED KNOTS
# ============================================================

CSV.write(
    joinpath(
        out,
        "29_galiana_plot_knots.csv"
    ),
    exported
)


# ============================================================
# SAVE NOTE
# ============================================================

write(
    joinpath(
        out,
        "29_galiana_plot_note.txt"
    ),

    """
Cumulative degree distributions are plotted on logarithmic x and y axes.

The y-axis labels use ordinary decimal probabilities rather than scientific 10^n notation.

The y-axis tick positions and y-axis limits are independent.

Each dataset and guild has its own y-axis limits. This avoids imposing a common multiplicative padding on datasets with very different probability ranges.

The 1.0 probability tick lies inside the plotting region rather than exactly at the upper panel boundary.

The lower labelled tick also lies inside the plotting region rather than defining the lower panel boundary.

Datasets are displayed in Galiana et al. order:

Plant-pollinator:
Garraf PP
Garraf PP2
Montseny
Gottin PP
Nahuel

Host-parasite:
Garraf HP
Quercus
Olot
Gottin HP
Galpar

Four figures are produced:

1. All original nodes: consumers and resources
2. Active nodes only: consumers and resources
3. Active consumers only: 2 x 5
4. Active resources only: 2 x 5

The cumulative distributions themselves are unchanged by these plotting modifications.
"""
)


# ============================================================
# FINISHED
# ============================================================

println(
    "Finished script 29."
)

println(
    "Figures produced:"
)

println(
    "  1. all original nodes"
)

println(
    "  2. active nodes"
)

println(
    "  3. active consumers, 2 x 5"
)

println(
    "  4. active resources, 2 x 5"
)

println(
    "Y-axis ticks and limits are now controlled independently."
)

println(
    "Y-axis limits are dataset- and guild-specific."
)

println(
    "Output directory: ",
    out
)