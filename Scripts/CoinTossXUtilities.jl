#=
CoinTossXUtilities:
- Julia version: 1.5.3
- Authors: Ivan Jericevich, Patrick Chang, Tim Gebbie
- Function: Provide the necessary functions for running simulations with CoinTossX
- Structure:
    1. Build, deploy, start CoinTossX and initialise Java Virtual Machine with required byte code paths
    2. Order creation structure
    3. Agent structure
    4. Initialize client by logging in to the trading gateway and starting the trading session
    5. Submit an order to CoinTossX
    6. Modify an existing order
    7. Cancel an existing order
    8. Receive updates to the best bid/ask
    9. Receive snapshot of LOB
    10. Receive static price reference at the end of each auction
    11. Receive string of orders submitted/executed
    12. Market data timeouts
    13. Destroy client by logging out and ending the trading session
    14. Shutdown all components of CoinTossX
=#
using JavaCall, Sockets, Dates, Reactive
import JavaCall: iterator
#---------------------------------------------------------------------------------------------------
# Set the directory to the folder where CoinTossX is
directory = "/home/ivanjericevich"
directory = "/Users/patrickchang1/Exchange"
#----- Build, deploy, start CoinTossX and initialise Java Virtual Machine with required byte code paths -----#
# function StartJVM()
#     JavaCall.addClassPath(directory * "/CoinTossX/ClientSimulator/build/classes/main")
#     JavaCall.addClassPath(directory * "/CoinTossX/ClientSimulator/build/install/ClientSimulator/lib/*.jar")
#     JavaCall.init()
#     Juno.notification("JVM started"; kind = :Info, options = Dict(:dismissable => false))
# end
function StartCoinTossX(; build::Bool = true, deploy::Bool = true)
    # run(`sudo sysctl net.core.rmem_max=2097152`)
    # run(`sudo sysctl net.core.wmem_max=2097152`)
    # run(`sudo sysctl vm.nr_hugepages=10000`)
    cd(directory * "/CoinTossX")
	if build
		run(`./gradlew -Penv=local build -x test`)
	end
    if deploy
		run(`./gradlew -Penv=local clean installDist bootWar copyResourcesToInstallDir copyToDeploy deployLocal`)
	end
    cd(directory * "/run/scripts")
    run(`./startAll.sh`)
	Juno.notification("CoinTossX started"; kind = :Info, options = Dict(:dismissable => false))
    # cd(directory)
    JavaCall.addClassPath(directory * "/CoinTossX/ClientSimulator/build/install/ClientSimulator/lib/*.jar")
    JavaCall.init()
    Juno.notification("JVM started"; kind = :Info, options = Dict(:dismissable => false))
end
#---------------------------------------------------------------------------------------------------

#----- Order creation structure -----#
mutable struct Order
    orderId::Int64
    traderMnemonic::String
    volume::Int64
    price::Int64
    side::String
    type::String
    tif::String
    displayVolume::Int64
    mes::Int64
    stopPrice::Int64
    expireTime::String
    function Order(orderId::Int64, traderMnemonic::String, side::String, type::String, volume::Int64, price::Int64 = 0; tif::String = "Day", displayVolume = missing, mes::Int64 = 0, stopPrice::Int64 = 0, expireTime::String = Dates.format(now(), "YYYYmmdd-HH:MM:SS"))
        if ismissing(displayVolume)
            new(orderId, traderMnemonic, volume, price, side, type, tif, volume, mes, stopPrice, expireTime)
        else
            new(orderId, traderMnemonic, volume, price, side, type, tif, displayVolume, mes, stopPrice, expireTime)
        end
    end
    Order() = new(0, "", 0, 0, "", "", "", 0, 0, 0, "")
end
#---------------------------------------------------------------------------------------------------

#----- Agent structure -----#
struct Client
    id::Int64
    securityId::Int64
    javaObject::JavaObject{Symbol("client.Client")}
end
#---------------------------------------------------------------------------------------------------

#----- Initialize client by logging in to the trading gateway and starting the trading session -----#
function Login(clientId::Int64, securityId::Int64)
    utilities = @jimport example.Utilities
    javaObject = jcall(utilities, "loadClientData", JavaObject{Symbol("client.Client")}, (jint, jint), clientId, securityId)
    jcall(javaObject, "sendStartMessage", Nothing, ())
    Juno.notification("Logged in and trading session started"; kind = :Info, options = Dict(:dismissable => false))
    return Client(clientId, securityId, javaObject)
end
#---------------------------------------------------------------------------------------------------

#----- Submit an order to CoinTossX -----#
function SubmitOrder(client::Client, order::Order)
    jcall(client.javaObject, "submitOrder", Nothing, (jint, JString, jlong, jlong, JString, JString, JString, JString, jlong, jlong, jlong), Int32(order.orderId), order.traderMnemonic, order.volume, order.price, order.side, order.type, order.tif, order.expireTime, order.displayVolume, order.mes, order.stopPrice)
end
#---------------------------------------------------------------------------------------------------

#----- Modify an existing order -----#
function ModifyOrder(client::Client, order::Order)
    jcall(client.javaObject, "replaceOrder", Nothing, (JString, jlong, jlong, JString, JString, JString, jlong, jlong, jlong), order.orderId, order.volume, order.price, order.side, order.type, order.tif, order.displayVolume, order.mes, order.stopPrice)
end
#---------------------------------------------------------------------------------------------------

#----- Cancel an existing order -----#
function CancelOrder(client::Client, orderId::Int64, traderMnemonic::String, side::String, price::Int64)
    jcall(client.javaObject, "cancelOrder", Nothing, (jint, JString, JString, jlong), Int32(orderId), traderMnemonic, side, price)
end
#---------------------------------------------------------------------------------------------------

#----- Receive updates to the best bid/ask -----#
function UpdateMarketData!(client::Client; bestBid::Union{Vector{Int64}, Missing} = missing, bestAsk::Union{Vector{Int64}, Missing} = missing)
    if !ismissing(bestBid) && !ismissing(bestAsk)
        bestBid[1] = jcall(client.javaObject, "getBid", jlong, ()); bestBid[2] = jcall(client.javaObject, "getBidQuantity", jlong, ())
        bestAsk[1] = jcall(client.javaObject, "getOffer", jlong, ()); bestAsk[2] = jcall(client.javaObject, "getOfferQuantity", jlong, ())
    elseif !ismissing(bestBid)
        bestBid[1] = jcall(client.javaObject, "getBid", jlong, ()); bestBid[2] = jcall(client.javaObject, "getBidQuantity", jlong, ())
    elseif !ismissing(bestAsk)
        bestAsk[1] = jcall(client.javaObject, "getOffer", jlong, ()); bestAsk[2] = jcall(client.javaObject, "getOfferQuantity", jlong, ())
    end
end
function ReceiveMarketData(client::Client, side::Symbol, type::Symbol)
    if side == :Bid
        data = type == :Price ? jcall(client.javaObject, "getBid", jlong, ()) : jcall(client.javaObject, "getBidQuantity", jlong, ())
    elseif side == :Ask
        data = type == :Price ? jcall(client.javaObject, "getOffer", jlong, ()) : jcall(client.javaObject, "getOfferQuantity", jlong, ())
    end
    return data
end
function ReceiveMarketData(client::Client, side::Symbol)
    best = NamedTuple()
    if side == :Bid
        return (Price = jcall(client.javaObject, "getBid", jlong, ()), Volume = jcall(client.javaObject, "getBidQuantity", jlong, ()))
    else
        return (Price = jcall(client.javaObject, "getOffer", jlong, ()), Volume = jcall(client.javaObject, "getOfferQuantity", jlong, ()))
    end
end
#---------------------------------------------------------------------------------------------------

#----- Receive snapshot of LOB -----#
function ReceiveLOBSnapshot(client::Client, side::String)
    orders = jcall(client.javaObject, "lobSnapshot", JavaObject{Symbol("java.util.ArrayList")}, ())
    LOB = Dict{String, NamedTuple}()
    for order in iterator(orders)
        fields = split(unsafe_string(JavaCall.JNI.GetStringUTFChars(Ptr(order), Ref{JavaCall.JNI.jboolean}())), ",")
        if fields[2] == side
            push!(LOB, fields[1] => (Price = parse(Int, fields[4]), Volume = parse(Int, fields[3])))
        end
    end
    return LOB
end
function ReceiveLOBSnapshot(client::Client)
    orders = jcall(client.javaObject, "lobSnapshot", JavaObject{Symbol("java.util.ArrayList")}, ())
    LOB = Dict{String, NamedTuple}()
    for order in iterator(orders)
        fields = split(unsafe_string(JavaCall.JNI.GetStringUTFChars(Ptr(order), Ref{JavaCall.JNI.jboolean}())), ",")
        push!(LOB, fields[1] => (Side = fields[2], Price = parse(Int, fields[4]), Volume = parse(Int, fields[3])))
    end
    return LOB
end
#---------------------------------------------------------------------------------------------------

#----- Receive static price reference at the end of each auction -----#
function StaticPriceReference(client::Client)
    return jcall(client.javaObject, "getStaticPriceReference", jlong, ())
end
#---------------------------------------------------------------------------------------------------

#----- Receive string of orders submitted/executed -----#
function ExecutionReport(client::Client)
    return jcall(client.javaObject, "getExecutionReport", JString, ())
end
#---------------------------------------------------------------------------------------------------

#----- Market data timeouts -----#
function WaitForMarketData(client::Client)
    jcall(client.javaObject, "waitForMarketDataUpdate", Nothing, ())
end
#---------------------------------------------------------------------------------------------------

#----- Destroy client by logging out and ending the trading session -----#
function Logout(client::Client)
    jcall(client.javaObject, "sendEndMessage", Nothing, ())
    jcall(client.javaObject, "close", Nothing, ())
    Juno.notification("Logged out and trading session ended"; kind = :Info, options = Dict(:dismissable => false))
end
#---------------------------------------------------------------------------------------------------

#----- Shutdown all components of CoinTossX -----#
function StopCoinTossX()
    cd(directory * "/run/scripts")
    run(`./stopAll.sh`)
    Juno.notification("CoinTossX stopped"; kind = :Info, options = Dict(:dismissable => false))
    exit()
end
#---------------------------------------------------------------------------------------------------

#----- Market data listener -----#
function Listen(signal)
    receiver = UDPSocket()
    connected = bind(receiver, ip"127.0.0.1", 1234)
    if connected
		try
			while true
	            message = (String(recv(receiver)))
				push!(signal, message)
	        end
		finally
			close(receiver)
		end
    end
end
function OnEvent(signal)
	println(signal)
end
#---------------------------------------------------------------------------------------------------
using Reactive

# Start and login
StartCoinTossX(build = false, deploy = false)
client = Login(1, 1)
# Reactive signal
signal = Signal("")
preserve(map(OnEvent, signal))
# Start listener
status = @async Listen(signal)
# Submit orders
SubmitOrder(client, Order(5, "John", "Buy", "Limit", 100, 50))
@time SubmitOrder(client, Order(2, "John", "Buy", "Limit", 100, 49))
SubmitOrder(client, Order(4, "John", "Sell", "Market", 50))
SubmitOrder(client, Order(4, "John", "Buy", "Limit", 100, 50; tif = "GTT", expireTime = Dates.format(now() + Second(20), "YYYYmmdd-HH:MM:SS")))
# Logout and stop
Logout(client)
StopCoinTossX()

CancelOrder(client, 5, "John", "Buy", 50)

# Market data requests
ReceiveMarketData(client, :Bid, :Price)
ReceiveLOBSnapshot(client)


























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
