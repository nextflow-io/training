/*
 * The `splitText` operator allow the parsing of text file one or more lines at time
 */

count=0

Channel
   .fromPath('data/meta/random.txt')
   .splitText()
   .view { "${count++}: ${it.toUpperCase().trim()}" }

 