# Parte 4: Pruebas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Los plugins son software independiente en el que los desarrolladores de pipelines necesitan confiar.
Probar cada funcionalidad de forma independiente, fuera de un pipeline, garantiza que el plugin funcione correctamente antes de que alguien lo integre en un workflow.
En esta sección, escribirá y ejecutará pruebas usando el framework de pruebas Spock.

!!! tip "¿Empezando desde aquí?"

    Si se une en esta parte, copie la solución de la Parte 3 para usarla como punto de partida:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    Luego cambie al directorio del plugin:

    ```bash
    cd nf-greeting
    ```

Asegúrese de estar en el directorio del plugin:

```bash
cd nf-greeting
```

---

## 1. ¿Por qué hacer pruebas?

Una compilación exitosa significa que el código compila, pero no verifica que funcione como se espera.
Las pruebas unitarias son pequeñas piezas de código que verifican automáticamente si sus funciones producen la salida correcta para una entrada dada.
Por ejemplo, una prueba podría verificar que `#!groovy reverseGreeting("Hello")` retorna `"olleH"`.

Las pruebas son valiosas porque:

- Detectan errores antes de que los usuarios los encuentren
- Le dan confianza para hacer cambios sin romper nada
- Sirven como documentación que muestra cómo deben usarse las funciones

---

## 2. Entendiendo las pruebas de Spock

La plantilla del plugin usa [Spock](https://spockframework.org/), un framework de pruebas para Groovy.
Spock ya está configurado en el proyecto (a través de `build.gradle`), por lo que no necesita agregar nada.

Si ha usado herramientas de prueba antes (como `pytest` en Python o `testthat` en R), Spock cumple el mismo rol: usted escribe pequeñas funciones que llaman a su código con entradas conocidas y verifican las salidas.
La diferencia es que Spock usa bloques etiquetados (`given:`, `expect:`, `when:`, `then:`) que son similares a un proceso o workflow de Nextflow.

Esta es la estructura básica:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **Nombre de la prueba entre comillas**: Describe qué verifica la prueba. Use lenguaje simple.
2. **Bloque `given:`**: Configure lo que necesita para la prueba (crear objetos, preparar datos)
3. **Bloque `expect:`**: Las verificaciones reales. Cada línea debe ser `true` para que la prueba pase

Esta estructura hace que las pruebas sean legibles: "Dado un objeto de extensión, se espera que `reverseGreeting('Hello')` sea igual a `'olleH'`."

---

## 3. Escribir las pruebas

Escriba pruebas para las dos funciones que creó en la Parte 3: `reverseGreeting` y `decorateGreeting`.

### 3.1. Crear la clase de prueba

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

Ábralo en su editor y agregue el esqueleto vacío de la clase de prueba:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Pruebas para las funciones de extensión de saludo
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. Todas las clases de prueba de Spock extienden `Specification`. Este es el punto de partida para cualquier archivo de prueba de Spock.

### 3.2. Probar reverseGreeting

Agregue un método de prueba dentro del cuerpo de la clase.
El bloque `given:` crea una instancia de `GreetingExtension`, y el bloque `expect:` verifica que `reverseGreeting` invierta correctamente dos entradas diferentes.
Esto prueba la función directamente, sin ejecutar un pipeline.

=== "Después"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * Pruebas para las funciones de extensión de saludo
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()            // (1)!

            expect:
            ext.reverseGreeting('Hello') == 'olleH'     // (2)!
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

    1. Cree una instancia de su extensión para probarla directamente, sin ejecutar un pipeline
    2. Cada línea en `expect:` es una aserción; la prueba pasa solo si todas son `true`

=== "Antes"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * Pruebas para las funciones de extensión de saludo
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. Probar decorateGreeting

Agregue un segundo método de prueba después del primero.
Este verifica que `decorateGreeting` envuelva la cadena de entrada con `***` en cada lado.

=== "Después"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * Pruebas para las funciones de extensión de saludo
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }

        def 'should decorate a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.decorateGreeting('Hello') == '*** Hello ***'
        }
    }
    ```

=== "Antes"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * Pruebas para las funciones de extensión de saludo
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

---

## 4. Ejecutar las pruebas

```bash
make test
```

??? example "Salida de las pruebas"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **¿Dónde están los resultados de las pruebas?** Gradle oculta la salida detallada cuando todas las pruebas pasan.
    "BUILD SUCCESSFUL" significa que todo funcionó correctamente.
    Si alguna prueba falla, verá mensajes de error detallados.

??? exercise "Agregar una prueba de caso límite"

    Agregue una prueba que verifique que `reverseGreeting` maneje una cadena vacía.
    ¿Qué debería retornar `reverseGreeting('')`?
    Agregue la prueba, ejecute `make test` y verifique que pase.

    ??? solution "Solución"

        Agregue este método de prueba a `GreetingExtensionTest.groovy`:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        Una cadena vacía invertida sigue siendo una cadena vacía.

---

## 5. Ver el reporte de pruebas

Gradle genera un reporte de pruebas en HTML con resultados detallados para cada prueba.
Inicie un servidor web en el directorio del reporte:

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code le pedirá que abra la aplicación en su navegador.
Navegue hasta su clase de prueba para ver los resultados individuales de cada prueba:

![Reporte de pruebas mostrando que todas las pruebas pasaron](./img/test_report.png)

El reporte muestra cada método de prueba y si pasó o falló.

Presione ++ctrl+c++ para detener el servidor, luego regrese al directorio anterior:

```bash
popd
```

Regrese al directorio principal del proyecto:

```bash
cd ..
```

---

## Conclusión

Aprendió que:

- Las pruebas de Spock usan una estructura legible `given:`/`expect:`
- Use `make test` para ejecutar las pruebas y `build/reports/tests/test/` para el reporte HTML
- Las pruebas verifican el comportamiento y sirven como documentación sobre cómo deben usarse las funciones

---

## ¿Qué sigue?

Hasta ahora, su plugin agrega funciones personalizadas que los pipelines pueden llamar.
Los plugins también pueden reaccionar a eventos del workflow (una tarea completándose, un archivo siendo publicado, el pipeline finalizando) usando observadores de traza.
En la siguiente sección, construirá un observador que cuenta las tareas completadas e imprime un resumen cuando el pipeline finaliza.

[Continuar a la Parte 5 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
