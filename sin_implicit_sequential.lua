walltime = Talasma()
walltime:start()

ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

-- Parse parameters and print help
dim = util.GetParamNumber("-dim", 2, "simulated time frame in seconds")
p_gridName = util.GetParam("-grid", "grids/cube_" .. dim .. "d.ugx", "filename of underlying grid")
p_numRefs = util.GetParamNumber("-numRefs", 3, "number of refinements")
endTime = util.GetParamNumber("-endTime", 6, "simulated time frame in seconds")
p_N = util.GetParamNumber("-N", 16384, "simulated time frame in seconds")
alpha = util.GetParamNumber("-alpha", 0.1, "simulated time frame in seconds")

dt = util.GetParamNumber("-dt", endTime / p_N, "time step size")
-- modal = util.GetParamNumber("-Mod", 512, "simulated time frame in seconds")

startTime = 0

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


-- set up solver (using 'util/solver_util.lua')
solverDesc = {
    type = "bicgstab",
    precond = {
        type = "gmg",
        approxSpace = approxSpace,
        smoother = "jac",
        damping = 1,
        cycle = "V",
        preSmooth = 1,
        postSmooth = 1,
        baseSolver = "lu"
    },
    convCheck = {
        type = "standard",
        iterations = 100,
	reduction = 1e-9, -- relTol
        absolute = 1e-32,
        verbose = true
    }
}

solver = util.solver.CreateSolver(solverDesc)

print("\nsolving...")
u = GridFunction(approxSpace)
u:set(0.0)
Interpolate(analyticsolution, u, "t", "Inner", 0)
domainDisc:adjust_solution(u)

solFileName = "sequential"

out = MultiScriptor()

eval_out = EvalScriptor()
eval_out:setFile("difference")
eval_out:setGeneratorComponent("t")
eval_out:setVectorGenerator(analyticsolution)
eval_out:setDomain(domainDisc)
out:addScriptor(eval_out)

-- modal_out = VTKModScriptor(VTKOutput(), solFileName);
-- modal_out:setModal(modal)
-- out:addScriptor(modal_out)

time = Talasma()
time:start()

util.SolveLinearTimeProblem(u, domainDisc, solver, out, solFileName, "ImplEuler", 1, startTime, endTime, dt);

out:write_time_pvd(solFileName, u)
time:stop()
print(time:get() .. " seconds for time stepping")

walltime:stop()
print(walltime:get() .. " seconds wall time")
