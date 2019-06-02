xbraid_util = xbraid_util or {}

function brutil.SolveLinearTimeProblem(comm,
                                       u,
                                       domainDisc,
                                       linSolver,
                                       out,
                                       filename,
                                       timeScheme, -- implEuler
                                       orderOrTheta,
                                       startTime,
                                       endTime,
                                       n,
                                       maxStepSize,
                                       minStepSize,
                                       reductionFactor)
    if u == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No grid function for the solution specified.")
        util.PrintUsageOfSolveTimeProblem()
        exit()
    end

    if domainDisc == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No domain discretization specified.")
        util.PrintUsageOfSolveTimeProblem()
        exit()
    end

    if linSolver == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No lin. solver specified.")
        util.PrintUsageOfSolveTimeProblem()
        exit()
    end

    if timeScheme == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No time scheme specified.")
        util.PrintUsageOfSolveTimeProblem()
        exit()
    end

    if startTime == nil or endTime == nil then
        print("SolveLinearTimeProblem: Illegal parameters: Start or end time not specified.")
        util.PrintUsageOfSolveTimeProblem()
        exit()
    end

    if maxStepSize == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No max. time step specified.")
        util.PrintUsageOfSolveTimeProblem()
        exit()
    end

    -- check parameters
    if minStepSize == nil then
        minStepSize = maxStepSize
    end
    if reductionFactor == nil then
        reductionFactor = 0.5
    end

    -- create time disc
    local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)


    local ux = u:clone()
    ux:set(-10)

    domainDisc:adjust_solution(ux)
    -- set order for bdf to 1 (initially)
    if timeScheme:lower() == "bdf" then
        timeDisc:set_order(1)
    end

    local assembled_dt = nil
    local braid = Braid2dCPU1()

    braid:setCommunicator(comm);

    braid:setStartTime(startTime)
    braid:setEndTime(endTime)
    braid:setNumberOfTimesteps(n)

    braid:setStartVector(u)
    braid:setRemainingVector(ux)

    braid:setTimeDisc(timeDisc)
    braid:setLinearSolver(linSolver)

    braid:setOutput(out);
    braid:setFilename(filename);
    braid:setVerbose(true);

    braid:run()

end