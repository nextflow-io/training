Channel
      .from('foo', 'bar', 'baz')
      .view { "- $it" }