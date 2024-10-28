---
title: Groovy introduction
description: A quick intro to Groovy basic structures and idioms
---

# Groovy basic structures and idioms

Nextflow is a domain specific language (DSL) implemented on top of the Groovy programming language, which in turn is a super-set of the Java programming language. This means that Nextflow can run any Groovy or Java code.

You have already been using some Groovy code in the previous sections, but now it's time to learn more about it.

## Printing values

To print something is as easy as using one of the `print` or `println` methods.

```groovy linenums="1" title="snippet.nf"
println("Hello, World!")
```

The only difference between the two is that the `println` method implicitly appends a new line character to the printed string.

!!! tip

    Parentheses for function invocations are optional. Therefore, the following syntax is also valid:

    ```groovy linenums="1" title="snippet.nf"
    println "Hello, World!"
    ```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to print a string to the console

## Comments

Comments use the same syntax as C-family programming languages:

```groovy linenums="1" title="snippet.nf"
// comment a single line

/*
    a comment spanning
    multiple lines
*/
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to write comments in Groovy

## Variables

To define a variable, simply assign a value to it:

```groovy linenums="1" title="snippet.nf"
x = 1
println x

x = new java.util.Date()
println x

x = -3.1499392
println x

x = false
println x

x = "Hi"
println x
```

Local variables are defined using the `def` keyword:

```groovy linenums="1" title="snippet.nf"
def x = 'foo'
```

The `def` should be always used when defining variables local to a function or a closure.

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define variables in Groovy

## Lists

A List object can be defined by placing the list items in square brackets:

```groovy linenums="1" title="snippet.nf"
list = [10, 20, 30, 40]

println list
```

You can access a given item in the list with square-bracket notation (indexes start at `0`) or using the `get` method:

```groovy linenums="1" title="snippet.nf"
list = [10, 20, 30, 40]

println list[0]
println list.get(0)
```

In order to get the length of a list you can use the `size` method:

```groovy linenums="1" title="snippet.nf"
list = [10, 20, 30, 40]

println list.size()
```

You can use the `assert` keyword to test if a condition is true (similar to an `if` function).

Here, Groovy will print nothing if it is correct, else it will raise an AssertionError message.

```groovy linenums="1" title="snippet.nf"
list = [10, 20, 30, 40]

assert list[0] == 10
```

!!! question "Exercise"

    This assertion should be correct, try changing it to an incorrect one.

    ```groovy linenums="1" title="snippet.nf"
    list = [10, 20, 30, 40]

    assert list[0] == 10
    ```

    ??? Solution

        Your solution could look something similar to this:

        ```groovy linenums="1" title="snippet.nf"
        list = [10, 20, 30, 40]

        assert list[0] == 11
        ```

        You should see an error message similar to this:

        ```console title="Output"
        ERROR ~ assert list[0] == 11
           |   |
           |   10
           [10, 20, 30, 40]
        ```

Lists can also be indexed with negative indexes and reversed ranges.

```groovy linenums="1" title="snippet.nf"
list = [0, 1, 2]
assert list[-1] == 2
assert list[-1..0] == list.reverse()
```

!!! info

    In the last assert line you are referencing the initial list and converting this with a "shorthand" range (`..`), to run from the -1th element (2) to the 0th element (0).

List objects implement all methods provided by the [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html) interface, plus the extension methods provided by [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html).

```groovy linenums="1" title="snippet.nf"
assert [1, 2, 3] << 1 == [1, 2, 3, 1]
assert [1, 2, 3] + [1] == [1, 2, 3, 1]
assert [1, 2, 3, 1] - [1] == [2, 3]
assert [1, 2, 3] * 2 == [1, 2, 3, 1, 2, 3]
assert [1, [2, 3]].flatten() == [1, 2, 3]
assert [1, 2, 3].reverse() == [3, 2, 1]
assert [1, 2, 3].collect { it + 3 } == [4, 5, 6]
assert [1, 2, 3, 1].unique().size() == 3
assert [1, 2, 3, 1].count(1) == 2
assert [1, 2, 3, 4].min() == 1
assert [1, 2, 3, 4].max() == 4
assert [1, 2, 3, 4].sum() == 10
assert [4, 2, 1, 3].sort() == [1, 2, 3, 4]
assert [4, 2, 1, 3].find { it % 2 == 0 } == 4
assert [4, 2, 1, 3].findAll { it % 2 == 0 } == [4, 2]
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define a list in Groovy
    2. How to access a list item
    3. How to apply methods
    4. How to use the `assert` keyword

## Maps

Maps are like lists that have an arbitrary key instead of an integer. Therefore, the syntax is very much aligned.

```groovy linenums="1" title="snippet.nf"
map = [a: 0, b: 1, c: 2]
```

Maps can be accessed in a conventional square-bracket syntax or as if the key was a property of the map.

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```groovy linenums="1" title="snippet.nf"
map = [a: 0, b: 1, c: 2]

assert map['a'] == 0 // (1)!
assert map.b == 1 // (2)!
assert map.get('c') == 2 // (3)!
```

1. Using square brackets.
2. Using dot notation.
3. Using the `get` method.

To add data or to modify a map, the syntax is similar to adding values to a list:

```groovy linenums="1" title="snippet.nf"
map = [a: 0, b: 1, c: 2]

map['a'] = 'x' // (1)!
map.b = 'y' // (2)!
map.put('c', 'z') // (3)!
assert map == [a: 'x', b: 'y', c: 'z']
```

1. Using square brackets.
2. Using dot notation.
3. Using the put method.

Map objects implement all methods provided by the [java.util.Map](https://docs.oracle.com/javase/8/docs/api/java/util/Map.html) interface, plus the extension methods provided by [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html).

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define a map in Groovy
    2. How to access and update maps

## String interpolation

String literals can be defined by enclosing them with either _single-_ ('') or _double-_ ("") quotation marks.

```groovy linenums="1" title="snippet.nf"
foxtype = 'quick'
foxcolor = ['b', 'r', 'o', 'w', 'n']
println "The $foxtype ${foxcolor.join()} fox"

x = 'Hello'
y = 'World'
println '$x $y'
```

```console title="Output"
The quick brown fox
$x $y
```

!!! info

    Note the different use of `$` and `${..}` syntax to interpolate value expressions in a string literal.
    The `$x` variable was _not_ expanded, as it was enclosed by single quotes.

!!! question "Exercise"

    Modify the script above to print `Hello World` instead of `$x $y`.

    ??? Solution

        Modify `println '$x $y'` to `println "$x $y"`.

        ```groovy linenums="1" title="snippet.nf"
        foxtype = 'quick'
        foxcolor = ['b', 'r', 'o', 'w', 'n']
        println "The $foxtype ${foxcolor.join()} fox"

        x = 'Hello'
        y = 'World'
        println "$x $y"
        ```

Finally, string literals can also be defined using the `/` character as a delimiter. They are known as **slashy** strings and are useful for defining regular expressions and patterns, as there is no need to escape backslashes. As with double-quote strings they allow to interpolate variables prefixed with a `$` character.

See the difference below:

```groovy linenums="1" title="snippet.nf"
x = /tic\tac\toe/
y = 'tic\tac\toe'
z = "tic\tac\toe"

println x
println y
```

```console title="Output"
tic\tac\toe
tic    ac    oe
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define string literals in Groovy
    2. How to interpolate variables in string literals

## Multi-line strings

A block of text that spans multiple lines can be defined by delimiting it with triple single or double quotes:

```groovy linenums="1" title="snippet.nf"
text = """
    Hello there James.
    How are you today?
    """
println text
```

```console title="Output"
Hello there James.
How are you today?
```

Finally, multi-line strings can also be defined with slashy strings. For example:

```groovy linenums="1" title="snippet.nf"
text = /
    This is a multi-line
    slashy string!
    It's cool, isn't it?!
    /
println text
```

```console title="Output"
This is a multi-line
slashy string!
It's cool, isn't it?!
```

!!! info

    Like before, multi-line strings inside double quotes and slash characters support variable interpolation, while single-quoted multi-line strings do not.

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define multi-line strings in Groovy

## If statement

The `if` statement uses the same syntax common in other programming languages, such as Java, C, and JavaScript.

```groovy linenums="1"
if (< boolean expression >) {
    // true branch
}
else {
    // false branch
}
```

The `else` branch is optional. Also, the curly brackets are optional when the branch defines just a single statement.

```groovy linenums="1" title="snippet.nf"
x = 11
if (x > 10)
    println 'Hello'
```

```console title="Output"
Hello
```

!!! tip

    `null`, empty strings, and empty collections are evaluated to `false`.

    Therefore a statement like:

    ```groovy linenums="1" title="snippet.nf"
    list = [1, 2, 3]
    if (list != null && list.size() > 0) {
        println list
    }
    else {
        println 'The list is empty'
    }
    ```

    Can be written as:

    ```groovy linenums="1" title="snippet.nf"
    list = [1, 2, 3]
    if (list)
        println list
    else
        println 'The list is empty'
    ```

    See the [Groovy-Truth](http://groovy-lang.org/semantics.html#Groovy-Truth) for further details.

!!! tip

    In some cases it can be useful to replace the `if` statement with a ternary expression (aka a conditional expression):

    ```groovy linenums="1" title="snippet.nf"
    println list ? list : 'The list is empty'
    ```

    The previous statement can be further simplified using the [Elvis operator](http://groovy-lang.org/operators.html#_elvis_operator):

    ```groovy linenums="1" title="snippet.nf"
    println list ?: 'The list is empty'
    ```

!!! Question "Exercise"

    Write an if statement that prints `Hello` if the variable `x` is greater than 10 and `Goodbye` if it is less than 10.

    ??? Solution

        Your solution could look something similar to this:

        ```groovy linenums="1" title="snippet.nf"
        x = 11

        if (x > 10)
            println 'Hello'
        else
            println 'Goodbye'
        ```

        Or this:

        ```groovy linenums="1" title="snippet.nf"
        x = 11

        println x > 10 ? 'Hello' : 'Goodbye'
        ```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define an `if` statement in Groovy
    2. How to use the ternary operator
    3. How to use the Elvis operator

## For statement

The classical `for` loop syntax is supported:

```groovy linenums="1" title="snippet.nf"
for (int i = 0; i < 3; i++) {
    println("Hello World $i")
}
```

Iteration over list objects is also possible using the syntax below:

```groovy linenums="1" title="snippet.nf"
list = ['a', 'b', 'c']

for (String elem : list) {
    println elem
}
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define a `for` loop in Groovy

## Functions

It is possible to define a custom function into a script:

```groovy linenums="1" title="snippet.nf"
def fib(int n) {
    return n < 2 ? 1 : fib(n - 1) + fib(n - 2)
}

assert fib(10)==89
```

A function can take multiple arguments separating them with a comma.

The `return` keyword can be omitted and the function implicitly returns the value of the last evaluated expression. Also, explicit types can be omitted, though not recommended:

```groovy linenums="1" title="snippet.nf"
def fact(n) {
    n > 1 ? n * fact(n - 1) : 1
}

assert fact(5) == 120
```

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define a function in Groovy

## Closures

Closures are the Swiss army knife of Nextflow/Groovy programming. In a nutshell, a closure is a block of code that can be passed as an argument to a function.

A closure can also be used to define an anonymous function.

More formally, a closure allows the definition of functions as first-class objects.

```groovy linenums="1" title="snippet.nf"
square = { it * it }
```

The curly brackets around the expression `it * it` tells the script interpreter to treat this expression as code. The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.

Once compiled, the function object is assigned to the variable `square` as any other variable assignment shown previously.

To invoke the closure execution use the special method `call` or just use the round parentheses to specify the closure parameter(s):

```groovy linenums="1" title="snippet.nf"
assert square.call(5) == 25
assert square(9) == 81
```

As is, this may not seem interesting, but you can now pass the `square` function as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists:

```groovy linenums="1" title="snippet.nf"
x = [1, 2, 3, 4].collect(square)
println x
```

```console linenums="1" title="Output"
[1, 4, 9, 16]
```

By default, closures take a single parameter called `it`.

To give it a different name use the `->` syntax. For example:

```groovy linenums="1" title="snippet.nf"
square = { num -> num * num }
```

It’s also possible to define closures with multiple, custom-named parameters.

For example, when the method `each()` is applied to a map it can take a closure with two arguments, to which it passes the _key-value_ pair for each entry in the `map` object. For example:

```groovy linenums="1" title="snippet.nf"
printMap = { a, b -> println "$a with value $b" } // (1)!
values = ["Yue": "Wu", "Mark": "Williams", "Sudha": "Kumari"] // (2)!
values.each(printMap) // (3)!
```

1. Closure object that prints the key-value pair.
2. Defines a map object with three entries.
3. Invokes the `each` method passing the closure object.

```console title="Output"
Yue with value Wu
Mark with value Williams
Sudha with value Kumari
```

A closure has two other important features.

First, it can access and _modify_ variables in the scope where it is defined.

Second, a closure can be defined in an _anonymous_ manner, meaning that it is not given a name, and is only defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment:

```groovy linenums="1" title="snippet.nf"
result = 0 // (1)!
values = ["China": 1, "India": 2, "USA": 3] // (2)!
values.keySet().each { result += values[it] } // (3)!
println result
```

1. Defines a global variable.
2. Defines a map object.
3. Invokes the `each` method passing the closure object which modifies the `result` variable.

```console title="Output"
6
```

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html).

!!! cboard-list-2 "Summary"

    In this step you have learned:

    1. How to define a closure in Groovy
    2. How to invoke a closure

## More resources

The complete Groovy language documentation is available at [this link](http://groovy-lang.org/documentation.html#languagespecification).

A great resource to master Apache Groovy syntax is the book: [Groovy in Action](https://www.manning.com/books/groovy-in-action-second-edition).
