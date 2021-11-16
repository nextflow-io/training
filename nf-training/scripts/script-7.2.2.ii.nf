Channel
    .from( 'hello', 'world' )
    .map { word -> [word, word.size()] }
    .view { word, len -> "$word contains $len letters" }