#!/usr/bin/env julia

"""
Simple test to verify threads are using separate cores.
Run this while monitoring Activity Monitor or htop.
"""

using Base.Threads

println("Testing thread distribution across CPU cores")
println("Threads available: $(Threads.nthreads())")
println("\nStarting 30-second CPU stress test...")
println("Open Activity Monitor and watch CPU History to see core usage!\n")

# Stress test that should use all threads
function parallel_work(duration=30)
    println("Running for $duration seconds...")
    start_time = time()

    Threads.@threads for tid in 1:Threads.nthreads()
        my_tid = Threads.threadid()
        println("Thread $my_tid started on iteration $tid")

        iteration_count = 0
        while time() - start_time < duration
            # CPU-intensive work
            result = 0.0
            for i in 1:100000
                result += sin(i) * cos(i) * tan(i/1000)
            end
            iteration_count += 1

            # Print progress every 5 seconds
            if iteration_count % 1000 == 0
                elapsed = round(time() - start_time, digits=1)
                println("Thread $my_tid: $(elapsed)s elapsed, $iteration_count iterations")
            end
        end

        println("Thread $my_tid finished with $iteration_count total iterations")
    end
end

parallel_work(30)

println("\n✓ Test complete!")
println("\nIf each thread ran on a separate core, you should have seen:")
println("  - Multiple CPU cores at ~100% in Activity Monitor")
println("  - All threads completed similar number of iterations")
