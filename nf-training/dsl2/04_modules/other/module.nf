/*
 * module parameters are overridden by the invoking context
 */

params.alpha = 'x'
params.delta = 'x'

process foo {
  echo true
  """
  echo alpha=${params.alpha} and delta=${params.delta}
  """
}



