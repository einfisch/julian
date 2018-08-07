-- Copyright (c) 2010-2016:  G-CSC, Goethe University Frankfurt
-- Authors: Andreas Vogel, Sebastian Reiter
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

-- Parse parameters and print help
local gridName	= util.GetParam("--grid", "laplace_sample_grid_2d.ugx", --WAS 2D
							"filename of underlying grid")
local numRefs		= util.GetParamNumber("--numRefs", 3, "number of refinements")

util.CheckAndPrintHelp("Waermeleitungsgleichung mit nicht-linear(in 2D)");


-- initialize ug with the world dimension dim=2 and the algebra type
local blockSize = 2 -- select 2 for point-block
InitUG(2, AlgebraType("CPU", blockSize));  


-- Load a domain without initial refinements.
local mandatorySubsets = {"Inner", "Boundary"}
dom = util.CreateDomain(gridName, 0, mandatorySubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, numRefs, true)


-- set up approximation space: linear functions
local approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("v", "Lagrange", 1)
approxSpace:add_fct("w", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("approximation space:")
approxSpace:print_statistic()


-- use finite-volume scheme
local elemDisc = {}
elemDisc["v"] = ConvectionDiffusion("v", "Inner", "fv1")  -- potential
elemDisc["w"] = ConvectionDiffusion("w", "Inner", "fv1")  -- gating variable


--[[ set up discretization for 
  $$ \frac{\partial c}{\partial t} + \triangle c + r(c)*c = 0$$
  where  r(c) = alpha*(1-c/c_mx)*c 
--]]

local alpha = -0.8
local beta = 0.0
local tau = 10.0
local eps = 0.01

-- Initial values ("Anfangswerte")
-- start with c=1 in circle around center (and 0 elsewhere) 
function MyInitialValueV(x, y)
 -- if ((x-0.5)*(x-0.5)+(y)*(y)<0.0625) then return 1.0 
 -- else return 0.0 end
 return alpha
end


function MyInitialValueW(x, y)
   local wref =  alpha - (alpha^3)/3.0
    if ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<0.1) then 
      return wref-0.1
    else
      return wref 
    end
end


function DiffusionV(x,y)
  return eps, 0.0, 0.0, eps*0.1
end

-- \dot v = v - v^3 - w
function ReactionV(v,w)
  return (v - (1.0/3.0)*v*v*v - w) * (-1.0) --(v - (1.0/3.0)*v*v*v - w) * (-1.0)
end

function ReactionV_dV(v,w) return (1.0 - v*v) * (-1.0) end
function ReactionV_dW(v,w) return -1.0 * (-1.0) end

-- external current
function ISourceV(x,y,t,si) if (t<10) then return 10.0 else return 0.0 end end

-- \tau \dot w = v - a - b*w
function ReactionW(v,w) 
  return (v - alpha - beta*w)/tau * (-1.0)
end

function ReactionW_dV(v,w) return (1.0/tau)* (-1.0) end
function ReactionW_dW(v,w) return (-beta/tau)* (-1.0) end


local nonlinearGrowth = {}
nonlinearGrowth["v"] = LuaUserFunctionNumber("ReactionV", 2)
nonlinearGrowth["v"]:set_input(0, elemDisc["v"]:value())
nonlinearGrowth["v"]:set_input(1, elemDisc["w"]:value())
nonlinearGrowth["v"]:set_deriv(0, "ReactionV_dV")
nonlinearGrowth["v"]:set_deriv(1, "ReactionV_dW")

nonlinearGrowth["w"] = LuaUserFunctionNumber("ReactionW", 2)
nonlinearGrowth["w"]:set_input(0, elemDisc["v"]:value())
nonlinearGrowth["w"]:set_input(1, elemDisc["w"]:value())
nonlinearGrowth["w"]:set_deriv(0, "ReactionW_dV")
nonlinearGrowth["w"]:set_deriv(1, "ReactionW_dW")


-- elemDisc["v"]:set_mass_scale(1.0)  -- default
elemDisc["v"]:set_diffusion("DiffusionV")
elemDisc["v"]:set_reaction(nonlinearGrowth["v"])          
--   elemDisc["v"]:set_source("ISourceV") 

--elemDisc["w"]:set_mass_scale(tau)  
elemDisc["w"]:set_diffusion(0.0)
elemDisc["w"]:set_reaction(nonlinearGrowth["w"])  


local dirichletBND = DirichletBoundary()
dirichletBND:add(alpha, "v", "Boundary")
dirichletBND:add(alpha-alpha^3/3.0, "w", "Boundary")

local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc["v"])
domainDisc:add(elemDisc["w"])
domainDisc:add(dirichletBND)


-- set up solver (using 'util/solver_util.lua')
local solverDesc = {
	
	type = "newton",
	
	linSolver = {
	 type = "bicgstab",
	
	 precond = {
		  type		= "gmg",
		  approxSpace	= approxSpace,
		  smoother	= "sgs",
		  baseSolver	= "lu"
	 },
	
  },
	
}

local nlsolver = util.solver.CreateSolver(solverDesc)
print (nlsolver)

print("\nsolving...")
local A = AssembledLinearOperator(domainDisc)
local u = GridFunction(approxSpace)
local b = GridFunction(approxSpace)
u:set(0.0)
domainDisc:adjust_solution(u)
domainDisc:assemble_linear(A, b)


Interpolate("MyInitialValueV", u, "v")
Interpolate("MyInitialValueW", u, "w")

local startTime = 0
local endTime = 100.0
local dt = (endTime-startTime)/200.0
util.SolveNonlinearTimeProblem(u, domainDisc, nlsolver, VTKOutput(),
"FitzHughNagumo", "ImplEuler", 1, startTime, endTime, dt);


print("done")
