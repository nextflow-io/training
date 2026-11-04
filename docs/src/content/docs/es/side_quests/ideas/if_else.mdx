# Parte 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Haz que la vaca cite científicos famosos

Esta sección contiene algunos ejercicios adicionales para practicar lo que has aprendido hasta ahora.
Realizar estos ejercicios _no es obligatorio_ para comprender las partes posteriores de la capacitación, pero proporcionan una forma divertida de reforzar tu aprendizaje al descubrir cómo hacer que la vaca cite científicos famosos.

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

### 1.1. Modifica el script `hello-containers.nf` para usar un proceso getQuote

Tenemos una lista de pioneros de la computación y la biología en el archivo `containers/data/pioneers.csv`.
A nivel general, para completar este ejercicio necesitarás:

- Modificar el `params.input_file` predeterminado para que apunte al archivo `pioneers.csv`.
- Crear un proceso `getQuote` que use el contenedor `quote` para obtener una cita para cada entrada.
- Conectar la salida del proceso `getQuote` al proceso `cowsay` para mostrar la cita.

Para la imagen de contenedor `quote`, puedes usar la que construiste tú mismo en el ejercicio adicional anterior o usar la que obtuviste de Seqera Containers.

!!! Hint "Pista"

    Una buena opción para el bloque `script` de tu proceso getQuote podría ser:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Puedes encontrar una solución a este ejercicio en `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modifica tu pipeline de Nextflow para permitir que se ejecute en los modos `quote` y `sayHello`.

Agrega algo de lógica de ramificación a tu pipeline para permitir que acepte entradas destinadas tanto a `quote` como a `sayHello`.
Aquí hay un ejemplo de cómo usar una declaración `if` en un workflow de Nextflow:

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

!!! Hint "Pista"

    Puedes usar `new_ch = processName.out` para asignar un nombre al canal de salida de un proceso.

Puedes encontrar una solución a este ejercicio en `containers/solutions/hello-containers-4.2.nf`.

### Conclusión

¡Sabes cómo usar contenedores en Nextflow para ejecutar procesos y cómo construir algo de lógica de ramificación en tus pipelines!

### ¿Qué sigue?

¡Celebra, toma un descanso para estirarte y bebe algo de agua!

Cuando estés listo/a, continúa con la Parte 3 de esta serie de capacitación para aprender cómo aplicar lo que has aprendido hasta ahora a un caso de uso de análisis de datos más realista.
