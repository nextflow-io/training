process example {
  script:
  """
  echo 'Hello world!\nHola mundo!\nCiao mondo!\nHallo Welt!' > file
  cat file | head -n 1 | head -c 5 > chunk_1.txt
  gzip -c chunk_1.txt  > chunk_archive.gz
  """
}
