-- boundary conditions
-- names are used to disambiguate between functions with different signatures
-- as lua provides no way of inspecting the argument/return types of a function
-- all coefficients are defined as time-dependent, though t can be left unused
bcs = {
  [1] = {
    attrs = {3, 4, 7},
    coef = function (v, t)
      -- Constant is defined as a function
      return 12.55
    end
  },
  -- time-independent scalar coefficient
  [4] = {
    attrs = {4, 6, 1},
    coef = function (v, t)
      return v.x * 0.12
    end
  },
  -- time-independent vector coefficient
  [8] = {
    attrs = {14, 62, 11},
    vec_coef = function (v, t)
      x = v.x
      s = 0.1 / 64
      first = -s * x * x
      last = s * x * x * (8.0 - x)
      if v.dim == 2 then
        return Vector.new(first, last)
      else
        return Vector.new(first, 0, last)
      end
    end
  },
  -- time-dependent scalar coefficient
  [11] = {
    attrs = {14, 62, 11},
    coef = function (v, t)
      return v.y * t * 0.12
    end
  },
  -- time-dependent vector coefficient
  [6] = {
    attrs = {14, 62, 11},
    vec_coef = function (v, t)
      return v * t
    end
  },
}
