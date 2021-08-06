process foo {
  cpus 2
  memory 8.GB
  container 'image/name'

  script:
  """
  your_command --this --that
  """
}
