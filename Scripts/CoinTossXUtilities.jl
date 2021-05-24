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
    7. Cancel an existing order
    8. Receive updates to the best bid/ask
    13. Destroy client by logging out and ending the trading session
=#
using JavaCall, Sockets, Dates
directory = "/home/ivanjericevich"
directory = "/Users/patrickchang1/Exchange"
#---------------------------------------------------------------------------------------------------

#----- Build, deploy, start CoinTossX and initialise Java Virtual Machine with required byte code paths -----#
function StartCoinTossX(; build::Bool = true, deploy::Bool = true)
    run(`sudo sysctl net.core.rmem_max=2097152`)
    run(`sudo sysctl net.core.wmem_max=2097152`)
    run(`sudo sysctl vm.nr_hugepages=10000`)
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
    cd(directory)
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
    cd(directory)
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

#----- Cancel an existing order -----#
function CancelOrder(client::Client, orderId::Int64, traderMnemonic::String, side::String, price::Int64)
    jcall(client.javaObject, "cancelOrder", Nothing, (jint, JString, JString, jlong), Int32(orderId), traderMnemonic, side, price)
end
#---------------------------------------------------------------------------------------------------

#----- Receive updates to the best bid/ask -----#
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
end
#---------------------------------------------------------------------------------------------------
