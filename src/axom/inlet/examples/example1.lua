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
  [1] = {
    attrs = {3, 4, 7},
    coef = function (x, y, z)
      -- Constant is defined as a function
      return 12.55
    end
  },
  [4] = {
    attrs = {4, 6, 1},
    coef = function (x, y, z)
      return x * 0.12
    end
  },
  [8] = {
    attrs = {14, 62, 11},
    vec_coef = function (x, y, z)
      scale = 0.12
      return x * scale, y * scale, z * scale
    end
  }
}
