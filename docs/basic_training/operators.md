# Operators

Operators are methods that allow you to connect channels, transform values emitted by a channel, or apply some user-provided rules.

There are seven main groups of operators are described in greater detail within the Nextflow Reference Documentation, linked below:

1. [Filtering operators](https://www.nextflow.io/docs/latest/operator.html#filtering-operators)
2. [Transforming operators](https://www.nextflow.io/docs/latest/operator.html#transforming-operators)
3. [Splitting operators](https://www.nextflow.io/docs/latest/operator.html#splitting-operators)
4. [Combining operators](https://www.nextflow.io/docs/latest/operator.html#combining-operators)
5. [Forking operators](https://www.nextflow.io/docs/latest/operator.html#forking-operators)
6. [Maths operators](https://www.nextflow.io/docs/latest/operator.html#maths-operators)
7. [Other operators](https://www.nextflow.io/docs/latest/operator.html#other-operators)

## Basic example

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1"
nums = Channel.of(1,2,3,4) // (1)!
square = nums.map { it -> it * it } // (2)!
square.view() // (3)!
```

1. Creates a queue channel emitting four values
2. Creates a new channel, transforming each number into its square
3. Prints the channel content

![Channel map](img/channel-map.png)

Operators can also be chained to implement custom behaviors, so the previous snippet can also be written as:

```groovy linenums="1"
Channel
    .of(1,2,3,4)
    .map { it -> it * it }
    .view()
```

## Basic operators

Here we explore some of the most commonly used operators.

### `view()`

The `view` operator prints the items emitted by a channel to the console standard output, appending a _new line_ character to each item. For example:

```groovy linenums="1"
Channel
    .of('foo', 'bar', 'baz')
    .view()
```

```console title="Output"
foo
bar
baz
```

An optional _closure_ parameter can be specified to customize how items are printed. For example:

```groovy linenums="1"
Channel
    .of('foo', 'bar', 'baz')
    .view { "- $it" }
```

```console title="Output"
- foo
- bar
- baz
```

### `map()`

The `map` operator applies a function of your choosing to every item emitted by a channel and returns the items obtained as a new channel. The function applied is called the _mapping_ function and is expressed with a _closure_ as shown in the example below:

```groovy linenums="1"
Channel
    .of( 'hello', 'world' )
    .map { it -> it.reverse() }
    .view()
```

A `map` can associate a generic _tuple_ to each element and can contain any data.

```groovy linenums="1"
Channel
    .of( 'hello', 'world' )
    .map { word -> [word, word.size()] }
    .view { word, len -> "$word contains $len letters" }
```

!!! exercise

    Use `fromPath` to create a channel emitting the _fastq_ files matching the pattern `data/ggal/*.fq`, then use `map` to return a pair containing the file name and the path itself, and finally, use `view` to print the resulting channel.

    ??? result "Solution"

        ```groovy linenums="1"
        Channel
            .fromPath('data/ggal/*.fq')
            .map { file -> [ file.name, file ] }
            .view { name, file -> "> $name : $file" }
        ```

### `mix()`

The `mix` operator combines the items emitted by two (or more) channels into a single channel.

```groovy linenums="1"
c1 = Channel.of( 1,2,3 )
c2 = Channel.of( 'a','b' )
c3 = Channel.of( 'z' )

c1 .mix(c2,c3).view()
```

```console title="Output"
1
2
a
3
b
z
```

!!! warning

    The items in the resulting channel have the same order as in the respective original channels. However, there is no guarantee that the element of the second channel are appended after the elements of the first. Indeed, in the example above, the element `a` has been printed before `3`.

### `flatten()`

The `flatten` operator transforms a channel in such a way that every _tuple_ is flattened so that each entry is emitted as a sole element by the resulting channel.

```groovy linenums="1"
foo = [1,2,3]
bar = [4,5,6]

Channel
    .of(foo, bar)
    .flatten()
    .view()
```

```console title="Output"
1
2
3
4
5
6
```

### `collect()`

The `collect` operator collects all of the items emitted by a channel in a list and returns the object as a sole emission.

```groovy linenums="1"
Channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
```

It prints a single value:

```console title="Output"
[1,2,3,4]
```

!!! info

    The result of the `collect` operator is a **value** channel.

### `groupTuple()`

The `groupTuple` operator collects tuples (or lists) of values emitted by the source channel, grouping the elements that share the same key. Finally, it emits a new tuple object for each distinct key collected.

Try the following example:

```groovy linenums="1"
Channel
    .of( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
    .groupTuple()
    .view()
```

```console title="Output"
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

This operator is useful to process a group together with all the elements that share a common property or grouping key.

!!! exercise

    Use `fromPath` to create a channel emitting all of the files in the folder `data/meta/`, then use a `map` to associate the `baseName` prefix to each file. Finally, group all files that have the same common prefix.

    ??? result "Solution

        ```groovy linenums="1"
        Channel.fromPath('data/meta/*')
            .map { file -> tuple(file.baseName, file) }
            .groupTuple()
            .view { baseName, file -> "> $baseName : $file" }
        ```

### `join()`

The `join` operator creates a channel that joins together the items emitted by two channels with a matching key. The key is defined, by default, as the first element in each item emitted.

```groovy linenums="1"
left = Channel.of(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.of(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).view()
```

```console title="Output"
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
```

!!! note

    Notice _P_ is missing in the final result.

### `branch()`

The `branch` operator allows you to forward the items emitted by a source channel to one or more output channels.

The selection criterion is defined by specifying a closure that provides one or more boolean expressions, each of which is identified by a unique label. For the first expression that evaluates to a true value, the item is bound to a named channel as the label identifier. For example:

```groovy linenums="1"
Channel
    .of(1,2,3,40,50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }

result.small.view { "$it is small" }
result.large.view { "$it is large" }
```

!!! info

    The `branch` operator returns a multi-channel object (i.e., a variable that holds more than one channel object).

!!! note

    In the above example, what would happen to a value of 10? To deal with this, you can also use `>=`.

## More resources

Check the [operators documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site.
