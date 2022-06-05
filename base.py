def cartesianProd(pools):
  result = [[]]
  for pool in pools:
    result = [x + [y] for x in result for y in pool]
  return result