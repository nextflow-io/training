process foo {
  script:
  """
  echo "The current directory is \$PWD"
  """
}