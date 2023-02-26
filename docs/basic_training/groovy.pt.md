---
title: Introdução ao Groovy
description: Uma introdução rápida às estruturas básicas e expressões do Groovy
---

# Estruturas básicas e expressões do Groovy

Nextflow é uma linguagem específica de domínio (DSL) implementada sobre a linguagem de programação Groovy, que por sua vez é um superconjunto da linguagem de programação Java. Isso significa que o Nextflow pode executar qualquer código Groovy ou Java.

Aqui estão algumas sintaxes Groovy importantes que são comumente usadas no Nextflow.

## Imprimindo valores

Imprimir algo é tão fácil quanto usar um dos métodos `print` ou `println`.

```groovy linenums="1"
println("Olá, mundo!")
```

A única diferença entre os dois é que o método `println` anexa implicitamente um caractere de nova linha à string impressa.

!!! tip

    Parênteses para invocações de função são opcionais. Portanto, a seguinte sintaxe também é válida:

    ```groovy linenums="1"
    println "Olá, mundo!"
    ```

## Comentários

Os comentários usam a mesma sintaxe das linguagens de programação da família C:

```groovy linenums="1"
// comente um único arquivo de configuração

/*
    um comentário abrangendo
    várias linhas
*/
```

## Variáveis

Para definir uma variável, basta atribuir um valor a ela:

```groovy linenums="1"
x = 1
println x

x = new java.util.Date()
println x

x = -3.1499392
println x

x = false
println x

x = "Oi"
println x
```

As variáveis locais são definidas usando a palavra-chave `def`:

```groovy linenums="1"
def x = 'foo'
```

O `def` deve ser sempre usado ao definir variáveis locais para uma função ou clausura.

## Listas

Um objeto List pode ser definido colocando os itens da lista entre colchetes:

```groovy linenums="1"
lista = [10,20,30,40]
```

Você pode acessar um determinado item na lista com a notação de colchetes (índices começam em `0`) ou usando o método get:

```groovy linenums="1"
println lista[0]
println lista.get(0)
```

Para obter o comprimento de uma lista, você pode usar o método size:

```groovy linenums="1"
println lista.size()
```

Usamos a palavra-chave `assert` para testar se uma condição é verdadeira (semelhante a uma função `if`). Aqui, o Groovy não imprimirá nada se estiver correto, caso contrário, gerará uma mensagem AssertionError.

```groovy linenums="1"
assert lista[0] == 10
```

!!! note

    Esta afirmação deve estar correta, tente alterá-la para uma incorreta.

As listas também podem ser indexadas com índices negativos e intervalos invertidos.

```groovy linenums="1"
lista = [0,1,2]
assert lista[-1] == 2
assert lista[-1..0] == lista.reverse()
```

!!! info

    Na afirmação da última linha, estamos referenciando a lista inicial e convertendo-a com um intervalo "abreviado" (`..`), para executar do -1º elemento (2), o último, ao 0º elemento (0), o primeiro.

Objetos List implementam todos os métodos fornecidos pela interface [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html), mais os métodos de extensão fornecidos pelo [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html).

```groovy linenums="1"
assert [1,2,3] << 1 == [1,2,3,1]
assert [1,2,3] + [1] == [1,2,3,1]
assert [1,2,3,1] - [1] == [2,3]
assert [1,2,3] * 2 == [1,2,3,1,2,3]
assert [1,[2,3]].flatten() == [1,2,3]
assert [1,2,3].reverse() == [3,2,1]
assert [1,2,3].collect{ it+3 } == [4,5,6]
assert [1,2,3,1].unique().size() == 3
assert [1,2,3,1].count(1) == 2
assert [1,2,3,4].min() == 1
assert [1,2,3,4].max() == 4
assert [1,2,3,4].sum() == 10
assert [4,2,1,3].sort() == [1,2,3,4]
assert [4,2,1,3].find{it%2 == 0} == 4
assert [4,2,1,3].findAll{it%2 == 0} == [4,2]
```

## Mapas

Os mapas são como listas que possuem uma chave arbitrária em vez de um número inteiro. Portanto, a sintaxe é bem parecida.

```groovy linenums="1"
mapa = [a:0, b:1, c:2]
```

Os mapas podem ser acessados em uma sintaxe convencional de colchetes ou como se a chave fosse uma propriedade do mapa.

!!! info ""

    Clique no ícone :material-plus-circle: para ver explicações no código.

```groovy linenums="1"
assert mapa['a'] == 0 // (1)!
assert mapa.b == 1 // (2)!
assert mapa.get('c') == 2 // (3)!
```

1. Usando colchetes.
2. Usando a notação de ponto.
3. Usando o método get.

Para adicionar dados ou modificar um mapa, a sintaxe é semelhante à adição de valores a uma lista:

```groovy linenums="1"
mapa['a'] = 'x' // (1)!
mapa.b = 'y' // (2)!
mapa.put('c', 'z') // (3)!
assert mapa == [a:'x', b:'y', c:'z']
```

1. Usando colchetes.
2. Usando a notação de ponto.
3. Usando o método put.

Objetos Map implementam todos os métodos fornecidos pela interface [java.util.Map](https://docs.oracle.com/javase/8/docs/api/java/util/Map.html), mais os métodos de extensão fornecidos pelo [Groovy](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html).

## Interpolação de Strings

Strings podem ser definidas colocando-as entre aspas _simples_ ('') ou _duplas_ ("").

```groovy linenums="1"
tipoderaposa = 'rápida'
cordaraposa = ['m', 'a', 'r', 'r', 'o', 'm']
println "A $tipoderaposa raposa ${cordaraposa.join()}"

x = 'Olá'
println '$x + $y'
```

```console title="Output"
A rápida raposa marrom
$x + $y
```

!!! info

    Observe o uso diferente das sintaxes `$` e `${..}` para interpolar expressões de valor em uma string.

Por fim, strings também podem ser definidas usando o caractere `/` como delimitador. Elas são conhecidas como strings **com barras** e são úteis para definir expressões e padrões regulares, pois não há necessidade de escapar das barras invertidas. Assim como as strings de aspas duplas, elas permitem interpolar variáveis prefixadas com um caractere `$`.

Tente o seguinte para ver a diferença:

```groovy linenums="1"
x = /tic\tac\toe/
y = 'tic\tac\toe'

println x
println y
```

```console title="Output"
tic\tac\toe
tic    ac    oe
```

## Strings de várias linhas

Um bloco de texto que abrange várias linhas pode ser definido delimitando-o com aspas simples ou duplas triplas:

```groovy linenums="1"
text = """
    Hello there James.
    How are you today?
    """
println text
```

Por fim, strings de várias linhas também podem ser definidas com strings com barras. Por exemplo:

```groovy linenums="1"
text = /
    Esta é uma string abrangendo
    várias linhas com barras!
    Super legal, né?!
    /
println text
```

!!! info

    Como antes, strings de várias linhas dentro de aspas duplas e caracteres de barra suportam interpolação de variável, enquanto strings de várias linhas com aspas simples não.

## Declarações condicionais com if

The `if` statement uses the same syntax common in other programming languages, such as Java, C, JavaScript, etc.

```groovy linenums="1"
if( < boolean expression > ) {
    // true branch
}
else {
    // false branch
}
```

The `else` branch is optional. Also, the curly brackets are optional when the branch defines just a single statement.

```groovy linenums="1"
x = 1
if( x > 10 )
    println 'Hello'
```

!!! tip

    `null`, empty strings, and empty collections are evaluated to `false`.

    Therefore a statement like:

    ```groovy linenums="1"
    list = [1,2,3]
    if( list != null && list.size() > 0 ) {
        println list
    }
    else {
        println 'The list is empty'
    }
    ```

    Can be written as:

    ```groovy linenums="1"
    list = [1,2,3]
    if( list )
        println list
    else
        println 'The list is empty'
    ```

    See the [Groovy-Truth](http://groovy-lang.org/semantics.html#Groovy-Truth) for further details.

!!! tip

    In some cases it can be useful to replace the `if` statement with a ternary expression (aka conditional expression). For example:

    ```groovy linenums="1"
    println list ? list : 'The list is empty'
    ```

    The previous statement can be further simplified using the [Elvis operator](http://groovy-lang.org/operators.html#_elvis_operator), as shown below:

    ```groovy linenums="1"
    println list ?: 'The list is empty'
    ```

## Declarações de loop com for

The classical `for` loop syntax is supported as shown here:

```groovy linenums="1"
for (int i = 0; i <3; i++) {
    println("Hello World $i")
}
```

Iteration over list objects is also possible using the syntax below:

```groovy linenums="1"
list = ['a','b','c']

for( String elem : list ) {
    println elem
}
```

## Funções

It is possible to define a custom function into a script, as shown here:

```groovy linenums="1"
int fib(int n) {
    return n < 2 ? 1 : fib(n-1) + fib(n-2)
}

assert fib(10)==89
```

A function can take multiple arguments separating them with a comma. The `return` keyword can be omitted and the function implicitly returns the value of the last evaluated expression. Also, explicit types can be omitted, though not recommended:

```groovy linenums="1"
def fact( n ) {
    n > 1 ? n * fact(n-1) : 1
}

assert fact(5) == 120
```

## Clausuras

Closures are the Swiss army knife of Nextflow/Groovy programming. In a nutshell, a closure is a block of code that can be passed as an argument to a function. A closure can also be used to define an anonymous function.

More formally, a closure allows the definition of functions as first-class objects.

```groovy linenums="1"
square = { it * it }
```

The curly brackets around the expression `it * it` tells the script interpreter to treat this expression as code. The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.

Once compiled, the function object is assigned to the variable `square` as any other variable assignment shown previously. To invoke the closure execution use the special method `call` or just use the round parentheses to specify the closure parameter(s). For example:

```groovy linenums="1"
assert square.call(5) == 25
assert square(9) == 81
```

As is, this may not seem interesting, but we can now pass the `square` function as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists:

```groovy linenums="1"
x = [ 1, 2, 3, 4 ].collect(square)
println x
```

```console title="Output"
[ 1, 4, 9, 16 ]
```

By default, closures take a single parameter called `it`, to give it a different name use the `\->` syntax. For example:

```groovy linenums="1"
square = { num -> num * num }
```

It’s also possible to define closures with multiple, custom-named parameters.

For example, when the method `each()` is applied to a map it can take a closure with two arguments, to which it passes the _key-value_ pair for each entry in the `map` object. For example:

```groovy linenums="1"
printMap = { a, b -> println "$a with value $b" }
values = [ "Yue" : "Wu", "Mark" : "Williams", "Sudha" : "Kumari" ]
values.each(printMap)
```

```console title="Output"
Yue with value Wu
Mark with value Williams
Sudha with value Kumari
```

A closure has two other important features.

First, it can access and _modify_ variables in the scope where it is defined.

Second, a closure can be defined in an _anonymous_ manner, meaning that it is not given a name, and is defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment:

```groovy linenums="1"
result = 0 // (1)!
values = ["China": 1 , "India" : 2, "USA" : 3] // (2)!
values.keySet().each { result += values[it] } // (3)!
println result
```

1. Defines a global variable.
2. Defines a map object.
3. Invokes the `each` method passing the closure object which modifies the `result` variable.

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html).

## Mais recursos

The complete Groovy language documentation is available at [this link](http://groovy-lang.org/documentation.html#languagespecification).

A great resource to master Apache Groovy syntax is the book: [Groovy in Action](https://www.manning.com/books/groovy-in-action-second-edition).
