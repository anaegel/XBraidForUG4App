
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

function DirichletValue2d(x, y, t)
    return true, 0
end

-- Parse parameters and print help
gridName	= util.GetParam("-grid", "grids/simple_grid_2d.ugx",
        "filename of underlying grid")
numRefs		= util.GetParamNumber("-numRefs", 0, "number of refinements")

steadyState	= util.HasParamOption("-steadyState", "If specified, the steady state of the problem is computed. Else a time-dependent problem is computed.")

endTime 	= util.GetParamNumber("-endTime", 1.6, "simulated time frame in seconds")
dt			= util.GetParamNumber("-dt", 0.1, "time step size")


InitUG(2, AlgebraType("CPU", 1));


-- Load a domain without initial refinements.
requiredSubsets = {"Inner", "Boundary"}
dom = util.CreateDomain(gridName, 0, requiredSubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, numRefs, true)


-- set up approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("t", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


-- set up discretization
coolerDisc = ConvectionDiffusion("t", "Inner", "fv1")
coolerDisc:set_diffusion(-2.3/100)

flowBnd = DirichletBoundary()
flowBnd:add("DirichletValue2d", "t", "Boundary")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(coolerDisc)
domainDisc:add(flowBnd)


-- set up solver (using 'util/solver_util.lua')
solverDesc = {
    type = "bicgstab",
    precond = {
        type		= "gmg",
        approxSpace	= approxSpace,
        smoother	= "ilu",
        baseSolver	= "lu"
    }
}

solver = util.solver.CreateSolver(solverDesc)


print("\nsolving...")
u = GridFunction(approxSpace)
u:set(35.0)
domainDisc:adjust_solution(u)


if steadyState then
    local A = AssembledLinearOperator(domainDisc)
    local b = GridFunction(approxSpace)
    domainDisc:adjust_solution(u)
    domainDisc:assemble_linear(A, b)

    solver:init(A, u)
    solver:apply(u, b)

    solFileName = "sol_cooler"
    print("writing solution to '" .. solFileName .. "'...")
    WriteGridFunctionToVTK(u, solFileName)
    SaveVectorForConnectionViewer(u, solFileName .. ".vec")
else
    local startTime = 0
    util.SolveLinearTimeProblem(u, domainDisc, solver, VTKOutput(), "sequential",
            "ImplEuler", 1, startTime, endTime, dt);
end

print("done")
