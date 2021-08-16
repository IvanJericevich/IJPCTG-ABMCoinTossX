# Simulation and estimation of an agent-based market-model with a matching engine

## Authors
* Ivan Jericevich
* Patrick Chang
* Tim Gebbie


## Description
We present an agent-based model in an artificial exchange framework that is able to recover price-impact while replicating the conventional stylised facts by directly interacting with a decoupled Matching Engine (ME) that manages orders away from the software relating to agent interactions and rules. This repo contains scripts for the replication of results in the pubished article. The scripts in this repository provide a number of different functions:
* A flexible asynchronous ABM for interacting with CoinTossX
* A Julia API for CoinTossX along with the functionality for polling the market-data feed
* Code for the recovery of sensitivity analysis and calibration results
* Methodology for the replication of many stylized facts (including price impact)
* A simple implementation of the Nelder-Mead with Threshold Accepting optimization heuristic
* Construction of the method of moments objective and implementation of the simulated minimum distance method in the context of financial market ABMs
* JSE and CoinTossX data cleaning functions


## Prerequisites
* [Julia](https://julialang.org) programming langauge
* A text editor or IDE such as [Atom](https://flight-manual.atom.io/getting-started/sections/installing-atom/), [VS Code](https://code.visualstudio.com/download) or [Jupyter](https://jupyter.org/install)
* CoinTossX matching engine. The version of CoinTossX used is a fork from the official release of version [v1.1.0](https://github.com/dharmeshsing/CoinTossX/tree/v1.1.0) and contains added functionalty necessary for ABM simulation and calibration. This modified version can be found [here](https://github.com/IvanJericevich/CoinTossX/tree/abm-experiment). For instructions on how to install, start, and use CoinTossX refer to [CoinTossX](https://github.com/dharmeshsing/CoinTossX)

## Usage
Clone the repository
```sh
git clone https://github.com/IvanJericevich/IJPCTG-ABMCoinTossX.git
```
Packages can be installed using
```julia
Pkg.add("...")
```
Each script contains examples of how to run the relevant functions. [CoinTossXUtilities.jl](Scripts/CoinTossXUtilities.jl) contains the important functions for interacting with CoinTossX. [ABMVolatilityAuctionProxy.jl](Scripts/ABMVolatilityAuctionProxy.jl) is the script used to simulate using the ABM and submit generated orders via the model to CoinTossX. [Calibration.jl](Scripts/Calibration.jl) calibrates parameters to empirical JSE TAQ data. [SensitivityAnalysis.jl](Scripts/SensitivityAnalysis.jl) provides code for the recovery of sensitivity analysis results and figures. Important to the calibration of parameters and sensitivity analysis is the algorithm for cleaning and re-formatting raw data from both COinTossX and JSE. This is demonstrated in [DataCleaning](DataCleaning).
