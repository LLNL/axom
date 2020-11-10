-- boundary conditions
-- names are used to disambiguate between functions with different signatures
-- as lua provides no way of inspecting the argument/return types of a function
bcs = {
  [1] = {
    attrs = { [1] = 3, [2] = 4, [3] = 7},
    coef = function (x, y, z)
      -- Constant is defined as a function
      return 12.55
    end
  },
  -- time-independent scalar coefficient
  [4] = {
    attrs = { [7] = 4, [12] = 6, [9] = 1},
    coef = function (x, y, z)
      return x * 0.12
    end
  },
  -- time-independent vector coefficient
  [8] = {
    attrs = { [4] = 14, [8] = 62, [6] = 11},
    vec_coef = function (x, y, z)
      scale = 0.12
      return x * scale, y * scale, z * scale
    end
  },
  -- time-dependent scalar coefficient
  [11] = {
    attrs = { [4] = 14, [8] = 62, [6] = 11},
    coef_t = function (x, y, z, t)
      return y * t * 0.12
    end
  },
  -- time-dependent vector coefficient
  [6] = {
    attrs = { [4] = 14, [8] = 62, [6] = 11},
    vec_coef_t = function (x, y, z, t)
      return x * t, y * t, z * t
    end
  },
}
