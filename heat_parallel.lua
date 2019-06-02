
ug_load_script("xbraid_util.lua")

numSpatialProcs = util.GetParamNumber("-numSpatialProcs", 2, "simulated time frame in seconds")


numWorldRanks = NumProcs()
xCommunicator = XCommunicator()

if numWorldRanks % numSpatialProcs == 0 then -- make sure numSpatialProcs*numTemporal == numWorldRanks
    xCommunicator:split(numSpatialProcs)
    print("Using: ".. numSpatialProcs .." of " ..numWorldRanks .." for spatial")
else
    xCommunicator:split(1)
    print("Using: ".. 1 .." of " ..numWorldRanks .." for spatial")
end

