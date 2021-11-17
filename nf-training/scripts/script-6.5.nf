process foo {
  cpus 2
  memory 1.GB
  container 'image/name'

  script:
  """
  echo your_command --this --that
  """
}