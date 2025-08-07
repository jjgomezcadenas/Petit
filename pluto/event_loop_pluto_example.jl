### A Pluto.jl notebook ###
# Example of using event_loop_pluto

using Markdown
using InteractiveUtils

# ╔═╡ 1
begin
    cmdir = joinpath(ENV["DATA"], "HD5t")
    pdir = joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 2
begin
    using Pkg
    Pkg.activate(pdir)
    Pkg.instantiate()
end

# ╔═╡ 3
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ 4
jn = ingredients(string(pdir,"/src/Petit.jl"))

# ╔═╡ 5
md"""
# Example: Using event_loop_pluto

This notebook demonstrates the correct usage of `event_loop_pluto` which returns both results and a progress display function.
"""

# ╔═╡ 6
# Example usage - CORRECT WAY
begin
    # The function returns a tuple: (results, progress_display)
    results, progress_display = jn.Petit.event_loop_pluto(cmdir; 
        input_file="0nubb.next.h5",
        events_to_run=100, 
        voxel_size_mm=5.0,
        max_distance_mm=10.0, 
        energy_threshold_kev=10.0)
    
    md"""
    ## Analysis Complete!
    
    - Processed $(results.n_events_processed) events
    - Single track events: $(results.n_single_track)
    - Two track events: $(results.n_two_track)
    - Three+ track events: $(results.n_three_plus_track)
    - Failed events: $(results.n_failed)
    """
end

# ╔═╡ 7
# Display the progress bar (shows final state)
progress_display()

# ╔═╡ 8
md"""
## Alternative: If you only need the results

If you don't care about the progress display, you can ignore it:
"""

# ╔═╡ 9
begin
    # Just get the results, ignore progress display
    results_only, _ = jn.Petit.event_loop_pluto(cmdir; 
        input_file="0nubb.next.h5",
        events_to_run=50)
    
    md"Processed $(results_only.n_events_processed) events (ignoring progress display)"
end

# ╔═╡ Cell order:
# ╠═1
# ╠═2
# ╠═3
# ╠═4
# ╠═5
# ╠═6
# ╠═7
# ╠═8
# ╠═9