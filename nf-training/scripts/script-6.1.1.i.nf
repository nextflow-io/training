params.data = 'World'

process foo {
  script:
  """
  echo Hello $params.data
  """
}