xbraid_util = xbraid_util or {}


-- todo renew function
function xbraid_util.SolveLinearTimeProblem(comm,
                                            u, generator, cmp,
                                            domainDisc,
                                            linSolver,
                                            out,
                                            filename,
                                            timeScheme, -- implEuler
                                            orderOrTheta,
                                            startTime,
                                            endTime,
                                            n)
    if u == nil or generator == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No grid function for the solution specified.")
        exit()
    end

    if domainDisc == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No domain discretization specified.")
        exit()
    end

    if linSolver == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No lin. solver specified.")
        exit()
    end

    if timeScheme == nil then
        print("SolveLinearTimeProblem: Illegal parameters: No time scheme specified.")
        exit()
    end

    if startTime == nil or endTime == nil then
        print("SolveLinearTimeProblem: Illegal parameters: Start or end time not specified.")
        exit()
    end


    -- create time disc
    local timeDisc = util.CreateTimeDisc(domainDisc, timeScheme, orderOrTheta)

    local braid = Braid2dCPU1()

    braid:setCommunicator(comm);

    braid:setStartTime(startTime)
    braid:setEndTime(endTime)
    braid:setNumberOfTimesteps(n)

    braid:setStartVector(u)
    braid:setVectorGenerator(generator)
    braid:setGeneratorComponent(cmp)

    braid:setTimeDisc(timeDisc)
    braid:setLinearSolver(linSolver)

    braid:setOutput(out);
    braid:setFilename(filename);
    braid:setVerbose(true);
    braid:setResidual(true);
    braid:setMaxIterations(4);

    braid:run()

end

defaultBraidSettings = {
    type = "uniform", -- uniform, leveldependend
    time = { t0 = 0, tn = 10, n = 100 },
    maxLevels = 16,

    -- sequential = false,
    -- timeRefinement = false,
    -- maxRefinements = "",
    -- ----- level dependend

    -- cycle = "",
    -- vcycle = "",
    -- FMG, n*FMG, FMGV
    -- nrelax = "", -- number of relax
    solver = "",
    timeDisc = "",
    cfactor = 2, -- coarsening factor
    -- ---- level independend
    storeMatrix = 8, -- save Assembled Operator for level <= x
    storeValues = 0, -- store points for level >= x

    skipDowncycleWork = true,
    spatialCoarsenAndRefine = false, -- todo methods!

    minCoarsening = "",
    -- temporalNorm = 2, -- possible values 1 (one norm), 2 (two norm), 3 (infinity norm)
    exactSolver = false, -- using residual
    verbose = false, -- still printing braid output

    convCheck = { -- convCheck for braid
        -- type = "standard", -- todo ask
        iterations = 100,
        absolute = 1e-9,
        -- reduction = 1e-12, -- relTol
    }
}

function xbraid_util.SetConvCheck(desc, braid)
    if (desc.iterations ~= nil) then
        braid:setMaxIterations(desc.iterations)
    else
        braid:setMaxIterations(defaultBraidSettings.convCheck.iterations)
    end

    if (desc.absolute == nil and desc.reduction == nil) then
        braid.setAbsoluteTol(defaultBraidSettings.convCheck.absolute)
    else
        if (desc.absolute ~= nil) then
            braid:setAbsoluteTol(desc.absolute)
        end
        if (desc.reduction ~= nil) then
            braid:setReduction(desc.reduction)
        end
    end
end

function xbraid_util.ConfigurateUniform()
end

function xbraid_util.LinearLevelDependend()
end

function xbraid_util.NonLinearLevelDependend(desc, braid, app)

end

function xbraid_util.CreateBraidSolver(desc, XComm, domainDisc)
    if type(desc) == "table" then
        app = nil
        braid = nil
        if (desc.type == 'uniform') then
            print("using uniform linear level configuration")
            app = RGFBraidApp()
            -- app:setRecurringStrongIteration(3)
            if (desc.strongfirst ~= nil) then
                if (desc.strongfirst == 1) then
                    app:setStrongFirstIteration(true)
                end
            end
            if (desc.adaptiveSolver ~= nil) then

                if (desc.adaptiveSolver == 1) then
                    print("use adaptive tol")
                    app:setAdaptConv(true)
                    app:setLooseTol(desc.looseTol)
                    app:setTightTol(desc.tightTol)
                end
            end
        elseif (desc.type == 'multistageleveldependend') then
            print("using multi stage level dependend configuration")
            app = ITSGFBraidApp()
            if desc.exactSolver ~= nil then
                if desc.exactSolver == 0 then
                    app:setResidual(false,0)
                else
                    print("XBraid Residual set")
                    app:setResidual(true,0)
                end
            end
        elseif (desc.type == 'leveldependend') then
            print("using linear level dependend configuration")
            app = MGFBraidApp()
        else
            print("Error Unknown type for configuration")
            exit()
        end
        app:setDomainDisc(domainDisc)

        if (desc.time ~= nil) then
            app:setTimeValues(desc.time.t0, desc.time.tn, desc.time.n)
        else
            app:setTimeValues(defaultBraidSettings.time.t0, defaultBraidSettings.time.tn, defaultBraidSettings.time.n)
        end

        braid = Braid(XComm, app)

        maxLevels = 0
        if desc.maxLevels == nil then
            maxLevels = defaultBraidSettings.maxLevels
        else
            maxLevels = desc.maxLevels
        end
        print("using a maximum of " .. maxLevels .. " level")
        braid:setMaxLevels(maxLevels)

        if (desc.type == 'multistageleveldependend') then
            app:setMaxLevels(maxLevels)

            if type(desc.level) == "table" then
                leveliter = 0
                leveldesc = desc.level[1] -- idx == 0!

                idx = leveldesc.from

                solverdesc = leveldesc.solver
                solver = util.solver.CreateSolver(solverdesc)

                timedesc = leveldesc.timeDisc

                orderOrTheta = 1
                if leveldesc.orderOrTheta ~= nil then
                    orderOrTheta = leveldesc.orderOrTheta
                end

                storeOperator = true
                if leveldesc.storeOperator ~= nil then
                    storeOperator = leveldesc.storeOperator
                end

                forceConvergence = false
                if leveldesc.forceConvergence ~= nil then
                    forceConvergence = leveldesc.forceConvergence
                end
                timedisc = util.CreateTimeDisc(domainDisc, timedesc, orderOrTheta)

                if timedisc:num_stages() > 1 then
                    if solverdesc.type ~= "newton" then
                        -- todo other nonlinear solver?
                        print(solverdesc.type)
                        print("multistage needs a non linear solver (newton) @level 0")
                        exit()
                    end
                end

                linear = true
                if solverdesc.type == "newton" then
                    -- todo other nonlinear solver?
                    linear = false
                    else
                    linear = true
                end

                solvername = ""
                if solverdesc.name ~= nil then
                    solvername = solverdesc.name
                end

                for i = 1, table.getn(desc.level) - 1 do
                    leveldesc = desc.level[i + 1]
                    nextidx = leveldesc.from
                    while leveliter < nextidx do
                        print(leveliter .. " : " .. timedesc .. " (" .. orderOrTheta .. ") " .. solverdesc.type .. " - " .. solvername.. "\tstore?" ..tostring(storeOperator) )
                        app:setTimeDisc(timedisc, leveliter)
                        if linear then
                            app:setLinearSolver(solver, leveliter)
                        else
                            app:setSolver(solver, leveliter)
                        end
                        app:setStoreOperator(storeOperator, leveliter)
                        app:setForceConv(forceConvergence, leveliter)

                        leveliter = leveliter + 1
                    end


                    if leveldesc.timeDisc ~= nil then
                        timedesc = leveldesc.timeDisc
                    end

                    if leveldesc.orderOrTheta ~= nil then
                        orderOrTheta = leveldesc.orderOrTheta
                    end


                    if leveldesc.solver ~= nil then
                        solverdesc = leveldesc.solver
                    end

                    solvername = ""
                    if solverdesc.name ~= nil then
                        solvername = solverdesc.name
                    end

                    if solverdesc.type == "newton" then
                        -- todo other nonlinear solver?
                        linear = false
                    else
                        linear = true
                    end
                    solver = util.solver.CreateSolver(solverdesc)
                    timedisc = util.CreateTimeDisc(domainDisc, timedesc, orderOrTheta) -- create an object for every level
                     if timedisc:num_stages() > 1 then
                        if solverdesc.type ~= "newton" then
                            -- todo other nonlinear solver? && lowercase compare?
                            print("multistage needs a non linear solver (newton) @level: " .. leveliter)
                            exit()
                        end
                    end

                    if leveldesc.forceConvergence ~= nil then
                        forceConvergence = leveldesc.forceConvergence
                    end

                    if leveldesc.storeOperator ~= nil then
                        storeOperator = leveldesc.storeOperator
                    end
                end

                while leveliter < maxLevels do
                    -- fill up the remaining levels
                    print(leveliter .. " : " .. timedesc .. " (" .. orderOrTheta .. ") " .. solverdesc.type .. " - " .. solvername.. "\tstore?" ..tostring(storeOperator) )
                    app:setTimeDisc(timedisc, leveliter)
                    if linear then
                        app:setLinearSolver(solver, leveliter)
                    else
                        app:setSolver(solver, leveliter)
                    end
                    app:setStoreOperator(storeOperator, leveliter)
                    app:setForceConv(forceConvergence, leveliter)

                    leveliter = leveliter + 1
                end
            else
                print("for this configuration a level dependend description is neccessary")
            end
            print("finished multistage leveldependend configuration")



        elseif (desc.type == 'leveldependend') then
            app:setMaxLevels(maxLevels)
            if desc.storeOperator ~= nil then
                print("store operator for level 0 to " .. desc.storeOperator)
                app:setStoreOperator(desc.storeOperator)
            end

            if desc.forceConvergence ~= nil then
                if type(desc.forceConvergence) == "table" then
                    -- todo
                else
                    for i = 0, maxLevels - 1 do
                        app:setForceConv(desc.forceConvergence, i)
                    end
                end
            end

            if type(desc.timeDisc) == "table" then
                leveliter = 0
                idx = desc.timeDisc[1][1] -- idx == 0!
                timeDesc = desc.timeDisc[1][2]
                timeDisc = util.CreateTimeDisc(domainDisc, timeDesc, 1)
                for i = 1, table.getn(desc.timeDisc) do
                    nextidx = desc.timeDisc[i][1]
                    while leveliter < nextidx do
                        app:setTimeDisc(timeDisc, leveliter)
                        leveliter = leveliter + 1
                    end
                    timeDesc = desc.timeDisc[i][2]
                    timeDisc = util.CreateTimeDisc(domainDisc, timeDesc, 1)
                end
                while leveliter < maxLevels do
                    app:setTimeDisc(timeDisc, leveliter)
                    leveliter = leveliter + 1
                end
            else
                timeDisc = util.CreateTimeDisc(domainDisc, desc.timeDisc, 1)
                for i = 0, maxLevels - 1 do
                    app:setTimeDisc(timeDisc, i)
                end
            end

            if desc.solver.type == nil then
                print("solver NIL")
                leveliter = 0
                idx = desc.solver[1][1] -- idx == 0!
                solverDesc = desc.solver[1][2]
                -- solver = util.solver.CreateSolver(solverDesc)
                for i = 1, table.getn(desc.solver) do
                    nextidx = desc.solver[i][1]

                    while leveliter < nextidx do
                        solver = util.solver.CreateSolver(solverDesc) -- sds
                        app:setLinearSolver(solver, leveliter)
                        leveliter = leveliter + 1
                    end
                    solverDesc = desc.solver[i][2]
                    -- solver = util.solver.CreateSolver(solverDesc)
                end
                while leveliter < maxLevels do
                    solver = util.solver.CreateSolver(solverDesc)
                    app:setLinearSolver(solver, leveliter)
                    leveliter = leveliter + 1
                end
            else
                solverDesc = desc.solver
                -- solver = util.solver.CreateSolver(desc.solver)

                for i = 0, maxLevels - 1 do
                    solver = util.solver.CreateSolver(solverDesc)
                    app:setLinearSolver(solver, i)
                end


            end
        else
            -- initial RGFBraidApp
            solver = util.solver.CreateSolver(desc.solver)
            timeDisc = util.CreateTimeDisc(domainDisc, desc.timeDisc, 1)
            app:setLinearSolver(solver, 0)
            app:setTimeDisc(timeDisc, 0)
            if desc.forceConvergence ~= nil then
                app:setForceConv(desc.forceConvergence)
            end

        end

        -- Braid related configuration
        -- todo defaultCFactor = x,
        if desc.convCheck ~= nil then
            xbraid_util.SetConvCheck(desc.convCheck, braid)
        else
            xbraid_util.SetConvCheck(defaultBraidSettings.convCheck, braid)
        end

        if desc.temporalNorm ~= nil then
            norm = desc.temporalNorm
            braid:setTemporalNorm(norm)
        end

        if desc.verbose ~= nil then
            app:setVerbose(desc.verbose)
        end

        if desc.exactSolver ~= nil then
            if desc.exactSolver == 0 then
                print("User defined residual set")
                braid:setResidual()
            else
                print("XBraid Residual set")
            end
        end

        if desc.nrelax ~= nil then
            braid:setNRelax(-1, desc.nrelax)
        end

        if type(desc.cfactor) == "table" then
            leveliter = 0
            idx = desc.cfactor[1][1] -- idx == 0!
            cfactor = desc.cfactor[1][2]

            for i = 1, table.getn(desc.cfactor) do
                nextidx = desc.cfactor[i][1]
                while leveliter < nextidx do

                    braid:setCFactor(leveliter, cfactor);
                    leveliter = leveliter + 1
                end
                cfactor = desc.cfactor[i][2]

            end
            while leveliter < maxLevels do

                braid:setCFactor(leveliter, cfactor);
                leveliter = leveliter + 1
            end

        else
            if cfactor ~= nil then
                cfactor = desc.cfactor
                braid:setCFactor(-1, cfactor)
            end

        end

        if desc.skipDowncycleWork ~= nil then
            braid:setSkip(desc.skipDowncycleWork)
        end

        if desc.storeValues ~= nil then
            print("store values from level " .. desc.storeValues .. " up to max Level")
            braid:setStoreValues(desc.storeValues)
        end

        if desc.storeLevel ~= nil then
            -- alternative for storeValues
            print("store values from level " .. desc.storeValues .. " up to max Level")
            braid:setStoreValues(desc.storeLevel)
        end

        if desc.sequential ~= nil then
            print("Using sequential time stepping")
            braid:setSequential()
        end

        if desc.fmg ~= nil then
            print("set FMG to " .. desc.fmg)
            braid:setCycleNFMG(desc.fmg)
        end
        -- todo mincoarsening
        -- todo spatial coarse and refine
        -- todo cycle type
        -- todo time refinement
        -- todo max refinements
        -- todo access level
        -- todo print level
        -- todo print filename
        print("XBraid Object created")
        return braid
    else
        return xbraid_util.CreateBraidSolver(defaultBraidSettings)
    end
end
