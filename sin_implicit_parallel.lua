walltime = Talasma() -- get Time of whole execution
walltime:start()

ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")
ug_load_script("xbraid_util.lua") -- load neccessary XBraid lua interfaces

p_numSpatialProcs = util.GetParamNumber("-npX", 1, "number of spatial procs (must divide totalproc number)")  -- numSpatialProcs * numTimeProcs = numGlobalProcs
p_coarseningStrategy = util.GetParamNumber("-strategy", 1, " see list below for specific coarsening strategies")
p_max_iter = util.GetParamNumber("-maxiter", 100, " maximum number of iterations")
p_max_level = util.GetParamNumber("-maxlevel", 15, " maximum number of levels")
p_cycle = util.GetParamNumber("-cycle", 0, " cycletype 1-3 V-Cycle, 4-6 F-Cycle with relaxation type FCF, F-FCF or F-relaxation")
p_useResidual = util.GetParamNumber("-useResidual", 0, " 0 use residual, 1 xbraid residual")
p_adaptiveSolver = util.GetParamNumber("-adaptiveSolver", 0, " use a tight and a loose tol for solving")
p_solver = util.GetParamNumber("-solver", 0, "see list below; 0 for GMG(V,1,1,jac(0.66))[reduction], 1 for GMG(V,1,1, ilu)[absolute]")
p_strongfirst = util.GetParamNumber("-first", 0, " Use a high accuracy in the first iteration")
alpha = util.GetParamNumber("-alpha", 0.1, " diffusion constant (default: 0.1)")

numWorldRanks = NumProcs()
xCommunicator = SpaceTimeCommunicator()
numTemporalProcs = numWorldRanks;
if numWorldRanks % p_numSpatialProcs == 0 then
    -- make sure numSpatialProcs*numTemporal == numWorldRanks
    xCommunicator:split(p_numSpatialProcs)
    numTemporalProcs = numWorldRanks / p_numSpatialProcs;
    print("Using: " .. p_numSpatialProcs .. " of " .. numWorldRanks .. " for spatial")
else
    xCommunicator:split(1)
    print("Using: " .. 1 .. " of " .. numWorldRanks .. " for spatial")
end

-- Parse parameters and print help
dim = util.GetParamNumber("-dim", 2, "simulated time frame in seconds")
p_gridName = util.GetParam("-grid", "grids/cube_"..dim.."d.ugx","filename of underlying grid")
p_numRefs = util.GetParamNumber("-numRefs", 5, "number of refinements") -- grid with 129x129 grid points
p_N = util.GetParamNumber("-N", 16384, "simulated time frame in seconds")
p_startTime = util.GetParamNumber("-startTime", 0, "simulated time frame in seconds")
p_endTime = util.GetParamNumber("-endTime", 4 * math.pi, "simulated time frame in")
p_dt = util.GetParamNumber("-dt", p_endTime / p_N, "time step size")
p_modal = util.GetParamNumber("-Mod", 512, " divisor to save every Mod'th result")
solFileName = "p" .. numWorldRanks .. "s" .. p_numSpatialProcs .. "_sin2"..dim
numSteps = p_N
InitUG(dim, AlgebraType("CPU", 1));

-- Lua problem definition ----------------------------------------------------------------------------------------------
------- 3d -------------------------------------------------------------------------------------------------------------
function sinSource3d(x, y, z, t)
    return -math.sin(math.pi * x) * math.sin(math.pi * y) * math.sin(math.pi * z) * (math.sin(t) - 3 * alpha * math.pi * math.pi * math.cos(t))
end

function sinBoundary3d(x, y, z, t)
    return true, 0
end

function sinAnalyticSolution3d(x, y, z, t)
    return math.sin(math.pi * x) * math.sin(math.pi * y) * math.sin(z) * math.cos(t)
end
------- 2d -------------------------------------------------------------------------------------------------------------
function sinSource2d(x, y, t)
    return -math.sin(math.pi * x) * math.sin(math.pi * y) * (math.sin(t) - 2 * alpha * math.pi * math.pi * math.cos(t))
end

function sinBoundary2d(x, y, t)
    return true, 0
end

function sinAnalyticSolution2d(x, y, t)
    return math.sin(math.pi * x) * math.sin(math.pi * y) * math.cos(t)
end
------- 1d -------------------------------------------------------------------------------------------------------------
function sinSource1d(x, t)
    return -math.sin(math.pi * x) * (math.sin(t) - 1 * alpha * math.pi * math.pi * math.cos(t))
end

function sinBoundary1d(x, t)
    return true, 0
end

function sinAnalyticSolution1d(x, t)
    return math.sin(math.pi * x) * math.cos(t)
end
-- C++ problem definition ----------------------------------------------------------------------------------------------
source = SinSourceOneCube()
source:setAlpha(alpha)
boundary = 0
analyticsolution = SinAnalyticSolutionOneCube()
analyticsolution:setAlpha(alpha)

-- guess function ------------------------------------------------------------------------------------------------------
function originator(x, y, z, t, si)
    -- guess function to initiate other time steps of vector v ( v[0] = u0)
    return 0
end
-- Prepare Domain ------------------------------------------------------------------------------------------------------

requiredSubsets = { "Inner", "Boundary" }
dom = util.CreateDomain(p_gridName, 0, requiredSubsets)
-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, p_numRefs, true)


-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("t", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


-- set up discretization
convection = ConvectionDiffusion("t", "Inner", "fe")
convection:set_diffusion(alpha)
convection:set_source(source)

dboundary = DirichletBoundary()
dboundary:add(boundary, "t", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(convection)
domainDisc:add(dboundary)

-- -------------------------------------- Solver -----------------------------------

VSolver = {
    type = "bicgstab",
    name = "Linear Weak V Solver",
    precond = {
        type = "gmg",
        approxSpace = approxSpace,
        smoother = "jac", -- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
        baseSolver = "lu",
        cycle = "V", -- gmg-cycle ["V", "F", "W"]
        preSmooth = 1, -- number presmoothing steps
        postSmooth = 1, -- number postsmoothing steps
    },
    convCheck = {
        type = "standard",
        iterations = 2,
        absolute = 1e-64,
        reduction = 1e-6,
        verbose = false
    }
}

exactSolver = {
    type = "bicgstab",
    name = "Linear Exact Solver",
    precond = {
        type = "gmg",
        approxSpace = approxSpace,
        smoother = "ilu",
        cycle = "V",
        preSmooth = 1,
        postSmooth = 1,
        baseSolver = "lu"
    },
    convCheck = {
        type = "standard",
        iterations = 100,
       reduction = 1e-9, -- relTol
--         absolute = 1e-9,
        absolute = 1e-32,
        verbose = false
    }
}

newtonSolver = {
    type = "newton",
    name = "Newton Solver",
    convCheck = {
        type		= "standard",
        iterations 	= 30,		-- maximum number of iterations
        absolute	= 5e-5,		-- absolut value of defect to be reached; usually 1e-7 - 1e-9
        reduction	= 1e-3,	-- reduction factor of defect to be reached; usually 1e-6 - 1e-8
        verbose		= false			-- print convergence rates if true
    },
    linSolver =
    {
        type = "bicgstab",
        precond = {
            type            = "gmg",        -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
            adaptive        = false,
            smoother		= "ilu",
            cycle           = "V",          -- gmg-cycle ["V", "F", "W"]
            preSmooth       = 3,            -- number presmoothing steps
            postSmooth      = 3,            -- number postsmoothing steps
            rap             = false,        -- comutes RAP-product instead of assembling if true
            rim				= false,			-- smooth on surface rim
            emulateFullRefined	= false,	-- emulate full grid (works with rap=true only)
            baseLevel       = 0,               -- gmg - baselevel
            gatheredBaseSolverIfAmbiguous = true,
            baseSolver = "lu",
            approxSpace = approxSpace
        },
        convCheck = {
            type            = "standard",
            iterations      = 100,          -- number of iterations
            absolute        = 1e-8,        -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be larger than in newton section)
            reduction       = 1e-6,         -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be larger than in newton section)
            verbose         = false,         -- print convergence rates if true
        }
    }
}

if p_solver == 1 then
    argsolver = VSolver
elseif p_solver == 2 then
    argsolver = newtonSolver
else
    argsolver = exactSolver
end

-- -------------------------------------- Coarsening ---------------------------------
coarsening = 2
function estimatecoarsening(numTemporalProcs, N, cfactor_proc, cfactor_comm)
    numPerProc = N / numTemporalProcs;
    coarse = {}
    cN = N
    level = 0
    while numPerProc / cfactor_proc >= 1 do
        table.insert(coarse, { level, cfactor_proc })
        level = level + 1
        numPerProc = numPerProc / cfactor_proc
        cN = cN / cfactor_proc
    end
    while cN > 1 do
        table.insert(coarse, { level, cfactor_comm })
        level = level + 1
        cN = cN / cfactor_comm
    end
    while cfactor_comm > 1 do
        table.insert(coarse, { level, 2 })
        cfactor_comm = cfactor_comm / 2
    end
    return coarse
end
print("")
print("")
print("using coarsening strategy nr " .. p_coarseningStrategy)
-- processor dependend coarsening

if p_coarseningStrategy == 1 then
    coarsening = { { 0, 2}, { 1, 2}, { 2, 2}, { 3, 2}, { 4, 2}, { 5, 2}, { 6,2}, { 7, 2}, { 8,2}, { 9, 2}, { 10, 2}, { 11, 2} }
elseif p_coarseningStrategy == 51299 then
    coarsening = { { 0, 64}, { 1, 2 }, { 2, 4 }, { 3, 4 }, { 4, 8}, { 5, 4}, { 6, 4}, { 7, 2}, { 8,2}, { 9, 2}, { 10, 2}, { 11, 2} }
elseif p_coarseningStrategy == 77777 then -- 193.340141964
    coarsening = { { 0, 16}, { 1, 8}, { 2, 8}, { 3, 16 }, { 4, 8}, { 5, 4}, { 6, 4}, { 7, 2}, { 8,2}, { 9, 2}, { 10, 2}, { 11, 2} }
end
-- -------------------------------------- Braid object definition ---------------------------------
defaultBraidSettingsA = {
    type = "uniform",
    time = { t0 = p_startTime, tn = p_endTime, n = numSteps },
    maxLevels = p_max_level,
    level = { { from = 0,
                solver = exactSolver, --newtonSolver,
                timeDisc = "impleuler",
                storeOperator = true,
                orderOrTheta = 1 } ,
              --[[{ from = 0,
                solver = newtonSolver,
               timeDisc = "sdirk",
                orderOrTheta = 4
              }, {
            from = 1,
            solver = exactSolver,
            timeDisc = "crank-nicolson",
            orderOrTheta = 0.5
        }--]]
    },
    timeDisc = "ImplEuler",
    cfactor = coarsening,
    --fmg = 2,
    storeOperator = 7,

    solver = argsolver, -- VSolver,

    exactSolver = p_useResidual,
    forceConvergence = false,
    adaptiveSolver = p_adaptiveSolver,
    looseTol = 5e-6,
    tightTol = 5e-8,
    strongfirst = p_strongfirst,

    temporalNorm = 3, -- {1,2,3}

    convCheck = { -- convCheck for braid
        iterations = p_max_iter,
        -- reduction = 1e-9,
        absolute = 5e-7
    },
    verbose = true
}

u = GridFunction(approxSpace)
u:set(0.0)
Interpolate(analyticsolution, u, "t", "Inner", 0)
domainDisc:adjust_solution(u)


braid = xbraid_util.CreateBraidSolver(defaultBraidSettingsA, xCommunicator, domainDisc)

braid:setStoreValues(0)
-- -------------------------------------- predefined cycle and relaxation types ---------------------------------
if p_cycle == 1 then
    print("V FCF Cycle")
    braid:setNRelax(-1, 1)
    braid:setNRelax(0, 1)
elseif p_cycle == 2 then
    print("V F-FCF Cycle")
    braid:setNRelax(-1, 1)
    braid:setNRelax(0, 0)
elseif p_cycle == 3 then
    print("V F Cycle")
    braid:setNRelax(-1, 0)
    braid:setNRelax(0, 0)
elseif p_cycle == 4 then
    print("F FCF Cycle")
    braid:setNRelax(-1, 1)
    braid:setNRelax(0, 1)
    braid:setCycleFMG()
elseif p_cycle == 5 then
    print("F F-FCF Cycle")
    braid:setNRelax(-1, 1)
    braid:setNRelax(0, 0)
    braid:setCycleFMG()
elseif p_cycle == 6 then
    print("F F Cycle")
    braid:setNRelax(-1, 0)
    braid:setNRelax(0, 0)
    braid:setCycleFMG()
end


-- -------------------------------------- Access definition ---------------------------------
out = MultiScriptor()

eval_out = EvalScriptor()
eval_out:setRelative(true)
eval_out:setFile("difference"..xCommunicator:getTemporalRank())
eval_out:setGeneratorComponent("t")
eval_out:setVectorGenerator(analyticsolution)
eval_out:setDomain(domainDisc)

out:addScriptor(eval_out)

-- if p_modal >= 0 then
-- modal_out = VTKModScriptor(VTKOutput(), solFileName)
-- modal_out:setModal(p_modal)
-- out:addScriptor(modal_out)
-- end

-- -------------------------------------- Access and information settins --------------------------
braid:setPrintLevel(3)
braid:setAccessLevel(3)

-- -------------------------------------- Braid object definition ---------------------------------
time = Talasma()
time:start()


braid:run(u, "originator", "t", out)

time:stop()
print(time:get() .. " seconds for parallel time stepping")

walltime:stop()
print(walltime:get() .. " seconds wall time")
