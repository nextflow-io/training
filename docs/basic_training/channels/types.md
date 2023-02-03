# Channel types

Nextflow distinguishes two different kinds of channels: **queue** channels and **value** channels.

## Queue channel

A _queue_ channel is an _asynchronous_ unidirectional _FIFO_ queue that connects two processes or operators.

-   _asynchronous_ means that operations are non-blocking.
-   _unidirectional_ means that data flows from a producer to a consumer.
-   _FIFO_ means that the data is guaranteed to be delivered in the same order as it is produced. First In, First Out.

A queue channel is implicitly created by process output definitions or using channel factories such as [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).

Try the following snippets:

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1"
ch = Channel.of(1,2,3)
println(ch) // (1)!
ch.view() // (2)!
```

1. Use the built-in print line function `println` to print the `ch` channel
2. Apply the `view` method to the `ch` channel prints each item emitted by the channels

!!! exercise

    Try to execute this snippet. You can do that by creating a new `.nf` file or by editing an already existing `.nf` file.

    ```groovy linenums="1"
    ch = Channel.of(1,2,3)
    ch.view()
    ```

## Value channels

A **value** channel (a.k.a. singleton channel) by definition is bound to a single value and it can be read unlimited times without consuming its contents. A `value` channel is created using the [value](https://www.nextflow.io/docs/latest/channel.html#value) factory method or by operators returning a single value, such as [first](https://www.nextflow.io/docs/latest/operator.html#first), [last](https://www.nextflow.io/docs/latest/operator.html#last), [collect](https://www.nextflow.io/docs/latest/operator.html#operator-collect), [count](https://www.nextflow.io/docs/latest/operator.html#operator-count), [min](https://www.nextflow.io/docs/latest/operator.html#operator-min), [max](https://www.nextflow.io/docs/latest/operator.html#operator-max), [reduce](https://www.nextflow.io/docs/latest/operator.html#operator-reduce), and [sum](https://www.nextflow.io/docs/latest/operator.html#operator-sum).

To better understand the difference between value and queue channels, save the snippet below as `example.nf`.

```groovy linenums="1" title="example.nf" linenums="1"
ch1 = Channel.of(1,2,3)
ch2 = Channel.of(1)

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    SUM(ch1,ch2).view()
}
```

When you run the script, it prints only 2, as you can see below:

```console
2
```

To understand why, we can inspect the queue channel and running Nextflow with DSL1 gives us a more explicit comprehension of what is behind the curtains.

```groovy linenums="1"
ch1 = Channel.of(1)
println ch1
```

```console
$ nextflow run example.nf -dsl1
...
DataflowQueue(queue=[DataflowVariable(value=1), DataflowVariable(value=groovyx.gpars.dataflow.operator.PoisonPill@34be065a)])
```

We have the value 1 as the single element of our queue channel and a poison pill, which will tell the process that there’s nothing left to be consumed. That’s why we only have one output for the example above, which is 2. Let’s inspect a value channel now.

```groovy linenums="1"
ch1 = Channel.value(1)
println ch1
```

```console
$ nextflow run example.nf -dsl1
...
DataflowVariable(value=1)
```

There is no poison pill, and that’s why we get a different output with the code below, where `ch2` is turned into a value channel through the `first` operator.

```groovy linenums="1"
ch1 = Channel.of(1,2,3)
ch2 = Channel.of(1)

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    SUM(ch1,ch2.first()).view()
}
```

```console title="Output"
4

3

2
```

Besides, in many situations, Nextflow will implicitly convert variables to value channels when they are used in a process invocation. For example, when you invoke a process with a pipeline parameter (`params.example`) which has a string value, it is automatically cast into a value channel.
