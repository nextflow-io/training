Channel
     .from( 'a', 'b', 'c' )
     .into{ foo; bar }

foo.view{ "Foo emits: " + it }
bar.view{ "Bar emits: " + it }