# Parte 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Hacer que la vaca cite científicos famosos

Esta sección contiene algunos ejercicios adicionales para practicar lo que ha aprendido hasta ahora.
Realizar estos ejercicios _no es obligatorio_ para comprender las partes posteriores del entrenamiento, pero proporcionan una forma divertida de reforzar sus aprendizajes al descubrir cómo hacer que la vaca cite científicos famosos.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 1.1. Modificar el script `hello-containers.nf` para usar un proceso getQuote

Tenemos una lista de pioneros de la informática y la biología en el archivo `containers/data/pioneers.csv`.
A nivel general, para completar este ejercicio necesitará:

- Modificar el `params.input_file` predeterminado para que apunte al archivo `pioneers.csv`.
- Crear un proceso `getQuote` que use el contenedor `quote` para obtener una cita para cada entrada.
- Conectar la salida del proceso `getQuote` al proceso `cowsay` para mostrar la cita.

Para la imagen del contenedor `quote`, puede usar la que construyó usted mismo en el ejercicio adicional anterior o usar la que obtuvo de Seqera Containers.

!!! Hint

    Una buena opción para el bloque `script` de su proceso getQuote podría ser:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Puede encontrar una solución a este ejercicio en `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modificar su pipeline de Nextflow para permitir que se ejecute en los modos `quote` y `sayHello`.

Agregue algo de lógica de ramificación a su pipeline para permitir que acepte entradas destinadas tanto a `quote` como a `sayHello`.
Aquí hay un ejemplo de cómo usar una sentencia `if` en un workflow de Nextflow:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint

    Puede usar `new_ch = processName.out` para asignar un nombre al canal de salida de un proceso.

Puede encontrar una solución a este ejercicio en `containers/solutions/hello-containers-4.2.nf`.

### Conclusión

¡Ahora sabe cómo usar contenedores en Nextflow para ejecutar procesos y cómo construir algo de lógica de ramificación en sus pipelines!

### ¿Qué sigue?

¡Celebre, tome un descanso para estirarse y beba algo de agua!

Cuando esté listo, continúe con la Parte 3 de esta serie de entrenamiento para aprender cómo aplicar lo que ha aprendido hasta ahora a un caso de uso de análisis de datos más realista.
