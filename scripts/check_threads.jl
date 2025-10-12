#!/usr/bin/env julia

"""
Diagnostic script to check thread configuration and core affinity on your Mac.
"""

using Base.Threads

println("="^60)
println("JULIA THREADING DIAGNOSTICS")
println("="^60)

# 1. Number of threads Julia is using
println("\n1. Julia Thread Configuration:")
println("   Threads available:     $(Threads.nthreads())")
println("   Thread pool size:      $(Threads.nthreads(:default))")
println("   Interactive threads:   $(Threads.nthreads(:interactive))")

# 2. System CPU information
println("\n2. System CPU Information:")
println("   Physical CPU cores:    $(Sys.CPU_THREADS)")

# Get detailed CPU info on macOS
if Sys.isapple()
    println("\n3. macOS CPU Details:")

    # Get CPU brand
    cpu_brand = read(`sysctl -n machdep.cpu.brand_string`, String) |> strip
    println("   CPU:                   $cpu_brand")

    # Get core counts
    physical_cores = parse(Int, read(`sysctl -n hw.physicalcpu`, String) |> strip)
    logical_cores = parse(Int, read(`sysctl -n hw.logicalcpu`, String) |> strip)

    println("   Physical cores:        $physical_cores")
    println("   Logical cores (HT):    $logical_cores")

    # Performance vs Efficiency cores (Apple Silicon)
    try
        perf_cores = parse(Int, read(`sysctl -n hw.perflevel0.physicalcpu`, String) |> strip)
        eff_cores = parse(Int, read(`sysctl -n hw.perflevel1.physicalcpu`, String) |> strip)
        println("   Performance cores:     $perf_cores")
        println("   Efficiency cores:      $eff_cores")
    catch
        println("   (No performance/efficiency core split detected)")
    end
end

# 4. Thread affinity test
println("\n4. Thread Distribution Test:")
println("   Testing which threads run on which cores...")

# Simple test to see thread IDs
function test_thread_distribution()
    results = Vector{Tuple{Int, Int}}(undef, Threads.nthreads())

    Threads.@threads for i in 1:Threads.nthreads()
        tid = Threads.threadid()
        results[i] = (i, tid)
    end

    return results
end

results = test_thread_distribution()
println("\n   Iteration -> Thread ID mapping:")
for (iter, tid) in results
    println("     Iteration $iter -> Thread $tid")
end

# 5. Compute-intensive test to stress cores
println("\n5. CPU Load Test (5 seconds):")
println("   Running parallel computation to stress cores...")
println("   Monitor CPU usage with: Activity Monitor (GUI) or 'top' command")

function stress_test(duration_sec=5)
    start_time = time()
    counters = zeros(Int, Threads.nthreads())

    Threads.@threads for tid in 1:Threads.nthreads()
        my_tid = Threads.threadid()
        while time() - start_time < duration_sec
            # Compute-intensive work
            x = 0.0
            for i in 1:10000
                x += sin(i) * cos(i)
            end
            counters[my_tid] += 1
        end
    end

    return counters
end

println("   Starting stress test...")
counters = stress_test(5)

println("\n   Work distribution (iterations per thread):")
for (tid, count) in enumerate(counters)
    if count > 0
        println("     Thread $tid: $(count) iterations")
    end
end

# 6. Recommendations
println("\n" * "="^60)
println("RECOMMENDATIONS")
println("="^60)

if Sys.isapple()
    physical_cores = parse(Int, read(`sysctl -n hw.physicalcpu`, String) |> strip)

    if Threads.nthreads() == 1
        println("⚠️  Julia is running with only 1 thread!")
        println("   To use multiple threads, start Julia with:")
        println("   julia -t auto              # Use all logical cores")
        println("   julia -t $physical_cores              # Use all physical cores")
    elseif Threads.nthreads() < physical_cores
        println("ℹ️  Julia is using $(Threads.nthreads()) threads")
        println("   You could use up to $physical_cores physical cores")
    elseif Threads.nthreads() == physical_cores
        println("✓ Julia is using all physical cores ($physical_cores threads)")
    else
        println("✓ Julia is using hyperthreading ($(Threads.nthreads()) threads)")
    end

    println("\nTo monitor core usage in real-time:")
    println("   1. GUI: Open Activity Monitor -> Window -> CPU History")
    println("   2. Terminal: Run 'top' and press '1' to see per-core usage")
    println("   3. Terminal: Run 'htop' (if installed via Homebrew)")
end

println("\n" * "="^60)
println("HOW TO VERIFY SEPARATE CORE USAGE")
println("="^60)
println("""
While running a multi-threaded Julia script:

Method 1: Activity Monitor (macOS GUI)
   1. Open Activity Monitor (/Applications/Utilities/)
   2. Go to Window -> CPU History
   3. You should see multiple CPU cores at ~100% usage

Method 2: Terminal with 'top'
   1. Open Terminal
   2. Run: top
   3. Press '1' to show individual CPU cores
   4. Look for Julia process using multiple cores

Method 3: Terminal with 'htop' (better visualization)
   1. Install: brew install htop
   2. Run: sudo htop
   3. See visual bars for each core

Method 4: Check with powermetrics (requires sudo)
   1. Run: sudo powermetrics --samplers cpu_power -n 1
   2. Shows per-cluster CPU usage (P-cores vs E-cores)
""")

println("="^60)
