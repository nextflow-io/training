foo = [1,2,3]
bar = [4,5,6]

Channel
    .from(foo, bar)
    .flatten()
    .view()