--
-- Basic tests
--
nx = tonumber(args[2]) or 200
vskip = math.floor(nx/200)
batch = tonumber(args[3]) or 1 -- batch time steps (batch 1 equals 2 steps)
block_n = tonumber(args[4]) or 64 -- block size of x or y axis

pond = {
  init = function(x,y) return 1, 0, 0 end,
  out = "pond.out",
  nx = nx,
  vskip = vskip,
  batch = batch,
  block_n = block_n
}

river = {
  init = function(x,y) return 1, 1, 0 end,
  out = "river.out",
  nx = nx,
  vskip = vskip,
  batch = batch,
  block_n = block_n
}

dam = {
  init = function(x,y)
    if (x-1)*(x-1) + (y-1)*(y-1) < 0.25 then
      return 1.5, 0, 0
    else
      return 1, 0, 0
    end
  end,
  out = "dam_break.out",
  nx = nx,
  vskip = vskip,
  batch = batch,
  block_n = block_n
}

wave = {
  init = function(x,y)
    return 1.0 + 0.2 * math.sin(math.pi * x), 1, 0
  end,
  out = "wave.out",
  frames = 100,
  nx = nx,
  vskip = vskip,
  batch = batch,
  block_n = block_n
}

simulate(_G[args[1]])
