using Test, Sockets, Reactive
include("NMTA.jl")
clearconsole()

@testset "Nelder-Mead" begin
    # Time limit test
    result1 = Optimize(NonDifferentiable(x -> (1.0 / 2.0) * (x[1]^2 + 0.9 * x[2]^2), [127.0, 921.0]), [127.0, 921.0], Options(time_limit = 0.0))
    @test initial_state(result1) == [127.0, 921.0]
    @test !converged(result1)
    @test time_limit(result1) < time_run(result1)
    @test (length(unique(g_norm_trace(result1))) != 1 || length(unique(f_trace(result1))) != 1 ) && issorted(f_trace(result1)[end:1])
    # Iteration limit test
    result2 = Optimize(NonDifferentiable(x -> (1.0 / 2.0) * (x[1]^2 + 0.9 * x[2]^2), [127.0, 921.0]), [127.0, 921.0], Options(show_trace = true))
    @test !converged(result2)
    @test iteration_limit_reached(result2)
    # Threshold accepting test
    result3 = Optimize(NonDifferentiable(x -> (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2, [0.0, 0.0]), [0.0, 0.0], Options(show_trace = true, ta_rounds = 5, f_abstol = 5.0, ξ = 0.2))
    # Visualisation
    # result4 = Optimize(NonDifferentiable(x -> (1.0 / 2.0) * (x[1]^2 + 0.9 * x[2]^2), [127.0, 921.0]), [127.0, 921.0], Options(trace_simplex = true, extended_trace = true, store_trace = true))
    # trace = trace(result4)
    # f = zeros(Float64, length(trace)); g_norm = zeros(Float64, length(trace)); f_simplex = fill(0.0, length(trace), 3); centroid = fill(0.0, length(trace), 2)
    # for tr in 1:length(trace)
    #     f[i] = tr[i].value
    #     g_norm[i] = tr[i].g_norm
    #     metadata = tr[i].metadata
    #     f_simplex[i, :] = transpose(metadata["simplex_values"])
    #     centroid[i, :] = transpose(metadata["centroid"])
    # end
    # # Objectives
    # plot(1:length(NMtrace), f, seriestype = :line, linecolor = :blue, label = "Objective", xlabel = "Iteration", ylabel = "Value", legend = :topleft)
    # plot!(twinx(), 1:length(NMtrace), g_norm, seriestype = :line, linecolor = :red, label = "NM Objective", legend = :topright)
    # # Simplex values
    # plot(1:length(NMtrace), f_simplex, seriestype = :line, linecolor = [:blue :red :green], xlabel = "Iteration", ylabel = "Simplex objective values")
    # # Centroids
    # plot(1:length(NMtrace), centroid, seriestype = :line, linecolor = [:blue :red], xlabel = "Iteration", ylabel = "Simplex objective values")
end

@testset "Sockets" begin
    receiver = UDPSocket()
    bind(receiver, ip"127.0.0.1", 2000)
    @async message = String(recv(receiver))
    sender = UDPSocket()
    bind(sender, ip"127.0.0.1", 1000)
    send(sender, ip"127.0.0.1", 2000, "Hello World from the UDP")
    @test message == "Hello World from the UDP"
    close(receiver); close(sender)
end

@testset "Reactive Programming" begin
    x = Signal(0)
    preserve(map(println, x))
    @async for i in 1:10
        push!(x, i)
        sleep(3)
    end
end

#=
udpsock = UDPSocket()
z = bind(udpsock,ip"127.0.0.1",40456)
join_multicast_group(udpsock, "224.0.1.1")
close(udpsock)

@async begin
    server = listen(2000)
    while true
       sock = accept(server)
       println("Hello World\n")
    end
end
CancelOrder(client, 2, "John", "Buy", 49)
=#
