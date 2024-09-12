function flexible_smooth_indicator(x, intervals, ε; mode="bounded_intervals")
    if mode == "bounded_intervals"
        return sum(
            smooth_indicator_tanh.(x, interval[1], interval[2], ε) for interval in intervals
        )
    elseif mode == "unbounded_intervals"
        return 0.5 * (tanh((x - intervals[1][1]) / ε) + 1)
    end
end

function smooth_indicator_tanh(x, a, b, ε)
    return 0.5 * (tanh((x - a) / ε) - tanh((x - b) / ε))
end

# For bounded intervals
fNC_bounded = (x, I, ε) -> flexible_smooth_indicator(x, I, ε; mode="bounded_intervals")

# For [a, +∞)
fNC_unboundedplus =
    (x, a, ε) -> flexible_smooth_indicator(x, [(a, 0)], ε; mode="unbounded_intervals")

# For (-∞,a]
fNC_unboundedminus =
    (x, a, ε) -> 1 - flexible_smooth_indicator(x, [(a, 0)], ε; mode="unbounded_intervals")
