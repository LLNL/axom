-- Basic example of a simulation Input File
thermal_solver={}
thermal_solver.mesh = {
 filename = "/data/star.mesh",
 serial = 1,
 parallel = 1
}
thermal_solver.order = 2
thermal_solver.timestepper = "quasistatic"
-- define initial conditions
thermal_solver.u0 = { type = "function", func = "BoundaryTemperature"}
-- define conducitivity
thermal_solver.kappa = { type = "constant", constant = 0.5}
-- linear solver parameters
thermal_solver.solver = {
 rel_tol = 1.e-6,
 abs_tol = 1.e-12,
 print_level = 0,
 max_iter = 100,
 dt = 1.0,
 steps = 1 
}
-- boundary conditions
thermal_solver.bcs = {
  ["temperature_1"] = {
    attrs = {3, 4, 7},
    coef = function (v)
      -- Constant is defined as a function
      return 12.55
    end
  },
  ["temperature_2"] = {
    attrs = {4, 6, 1},
    coef = function (v)
      return v.x * 0.12
    end
  },
  ["flux"] = {
    attrs = {14, 62, 11},
    vec_coef = function (v)
      scale = 0.12
      return v * scale
    end
  }
}
