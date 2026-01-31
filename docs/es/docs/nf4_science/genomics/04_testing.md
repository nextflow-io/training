# Parte 4: Agregando pruebas

En la primera parte de este curso, construiste un pipeline de llamado de variantes que fue completamente lineal y procesó los datos de cada muestra de forma independiente de las otras.

En la segunda parte, te mostramos cómo usar canales y operadores de canal para implementar el llamado conjunto de variantes con GATK.

En la tercera parte, modularizamos el pipeline.

En esta parte del entrenamiento, vamos a mostrarte cómo usar [**nf-test**](https://www.nf-test.com/), un marco de pruebas que se integra bien con Nextflow y facilita agregar pruebas tanto a nivel de módulo como a nivel de flujo de trabajo a tu pipeline. Para seguir esta parte del entrenamiento, deberías haber completado la Parte 1, Parte 2 y Parte 3, así como la [misión secundaria de nf-test](../../side_quests/nf-test.md), que cubre los conceptos básicos de nf-test y por qué las pruebas son importantes.

---

## 0. Calentamiento

!!! note "Nota"

    Asegúrate de estar en el directorio de trabajo correcto:
    `cd /workspaces/training/nf4-science/genomics`

Si trabajaste en las partes anteriores de este curso de entrenamiento, deberías tener una versión funcional del pipeline de genómica con la estructura de directorios de módulos apropiada.

??? abstract "Contenido del directorio"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

Este directorio de módulos se puede encontrar en el directorio `solutions` si lo necesitas.

Vamos a comenzar con el mismo flujo de trabajo que en la Parte 3, que hemos proporcionado en el archivo `genomics-4.nf`. Exactamente como en la [misión secundaria de nf-test](../../side_quests/nf-test.md), vamos a agregar algunos tipos diferentes de pruebas a los tres procesos en este pipeline, así como una prueba a nivel de flujo de trabajo.

### 0.1. Verificar que el flujo de trabajo se ejecute

Antes de comenzar a agregar pruebas, asegúrate de que el flujo de trabajo se ejecute como se espera.

```bash
nextflow run genomics-4.nf -resume
```

Esto debería verse muy familiar a estas alturas si has estado trabajando en este curso de entrenamiento desde el principio.

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Como anteriormente, ahora habrá un directorio `work` y un directorio `results_genomics` dentro de tu directorio de proyecto. Realmente haremos uso de estos resultados más adelante en nuestras pruebas. Pero a partir de ahora vamos a usar el paquete `nf-test` para probar el pipeline.

### 0.2. Inicializar `nf-test`

Como en la [misión secundaria de nf-test](../../side_quests/nf-test.md), necesitamos inicializar el paquete `nf-test`.

```bash
nf-test init
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "Contenido de nf-test.config"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

También crea un directorio `tests` que contiene un esqueleto de archivo de configuración.

### Conclusión

Ahora estamos listos para comenzar a escribir pruebas para nuestro pipeline de genómica.

### ¿Qué sigue?

Escribir pruebas básicas que evalúen si las llamadas a los procesos fueron exitosas y produjeron las salidas correctas.

---

## 1. Probar un proceso para éxito y salidas coincidentes

Comenzaremos probando el proceso `SAMTOOLS_INDEX`, que crea archivos de índice para archivos BAM para permitir un acceso aleatorio eficiente. Este es un buen primer caso de prueba porque:

1. Tiene una única entrada bien definida (un archivo BAM)
2. Produce una salida predecible (un archivo de índice BAI)
3. La salida debería ser idéntica para entradas idénticas

### 1.1. Generar un esqueleto de archivo de prueba

Primero, genera un esqueleto de archivo de prueba:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Salida del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Esto crea un archivo en el mismo directorio que `main.nf`.
Puedes navegar al directorio en el explorador de archivos y abrir el archivo, que debería contener el siguiente código:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Las afirmaciones iniciales deberían ser familiares de la [misión secundaria de nf-test](../../side_quests/nf-test.md):

- `assert process.success` establece que esperamos que el proceso se ejecute exitosamente y se complete sin fallas.
- `snapshot(process.out).match()` establece que esperamos que el resultado de la ejecución sea idéntico al resultado obtenido en una ejecución anterior (si aplica).
  Discutimos esto con más detalle más adelante.

Usando esto como punto de partida, necesitamos agregar las entradas de prueba correctas para el proceso de índice de samtools, y cualquier parámetro si aplica.

### 1.2. Mover el archivo de prueba y actualizar la ruta del script

Antes de comenzar a trabajar en completar la prueba, necesitamos mover el archivo a su ubicación definitiva. Parte de la razón por la que agregamos un directorio para cada módulo es que ahora podemos enviar pruebas en un directorio `tests` co-ubicado con el archivo `main.nf` de cada módulo. Crea ese directorio y mueve el archivo de prueba allí.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Ahora podemos simplificar la sección `script` del archivo de prueba a una ruta relativa:

=== "Después"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Antes"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

Esto le dice a la prueba dónde encontrar el archivo `main.nf` del módulo, sin tener que especificar la ruta completa.

### 1.3. Proporcionar entradas de prueba para SAMTOOLS_INDEX

El archivo esqueleto incluye un marcador de posición que necesitamos reemplazar con una entrada de prueba real, apropiada para la entrada de `samtools index`. La entrada apropiada es un archivo BAM, que tenemos disponible en el directorio `data/bam`.

=== "Después"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Nombrar la prueba basándose en la funcionalidad

Como aprendimos antes, es una buena práctica renombrar la prueba a algo que tenga sentido en el contexto de la prueba.

=== "Después"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    Esto toma una cadena arbitraria, así que podríamos poner lo que queramos.
    Aquí elegimos referirnos al nombre del archivo y su formato.

=== "Antes"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Ejecutar la prueba y examinar la salida

Ejecuta la prueba:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Como aprendimos anteriormente, esto verificó la afirmación básica sobre el éxito del proceso y creó un archivo de snapshot basado en la salida del proceso. Podemos ver el contenido del archivo snapshot en el archivo `tests/modules/samtools/index/tests/main.nf.test.snap`:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

También podemos ejecutar la prueba nuevamente y ver que pasa, porque la salida es idéntica al snapshot:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. Agregar más pruebas a `SAMTOOLS_INDEX`

A veces es útil probar un rango de diferentes archivos de entrada para asegurarnos de que estamos probando una variedad de problemas potenciales. Agrega pruebas para los archivos BAM de la madre y el padre en el trío de nuestros datos de prueba.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Luego puedes ejecutar la prueba nuevamente:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

Observa la advertencia, refiriéndose al efecto del parámetro `--update-snapshot`.

!!! note "Nota"

    Aquí estamos usando datos de prueba que usamos anteriormente para demostrar las salidas científicas del pipeline.
    Si hubiéramos estado planeando operar estas pruebas en un entorno de producción, habríamos generado entradas más pequeñas para propósitos de prueba.

    En general es importante mantener las pruebas unitarias lo más ligeras posible usando las piezas de datos más pequeñas necesarias y suficientes para evaluar la funcionalidad del proceso, de lo contrario el tiempo de ejecución total puede acumularse bastante seriamente.
    Una suite de pruebas que tarda demasiado en ejecutarse regularmente es una suite de pruebas que probablemente se omitirá en interés de la conveniencia.

### Conclusión

Has escrito tu primera prueba de módulo para un proceso de genómica, verificando que `SAMTOOLS_INDEX` crea correctamente archivos de índice para diferentes archivos BAM. La suite de pruebas asegura que:

1. El proceso se ejecuta exitosamente
2. Se crean archivos de índice
3. Las salidas son consistentes a través de las ejecuciones
4. El proceso funciona para todos los archivos BAM de muestra

### ¿Qué sigue?

Aprender cómo escribir pruebas para otros procesos en nuestro flujo de trabajo de genómica, usando el método setup para manejar procesos encadenados. También evaluaremos si las salidas, específicamente nuestros archivos VCF, contienen las llamadas de variantes esperadas.

---

## 2. Agregar pruebas a un proceso encadenado y probar contenido

Para probar `GATK_HAPLOTYPECALLER`, necesitamos proporcionar al proceso la salida de `SAMTOOLS_INDEX` como entrada. Podríamos hacer eso ejecutando `SAMTOOLS_INDEX`, recuperando sus salidas, y almacenándolas con los datos de prueba para el flujo de trabajo. Ese es en realidad el enfoque recomendado para un pipeline pulido, pero nf-test proporciona un enfoque alternativo, usando el método `setup`.

Con el método setup, podemos activar el proceso `SAMTOOLS_INDEX` como parte de la configuración de la prueba, y luego usar su salida como entrada para `GATK_HAPLOTYPECALLER`. Esto tiene un costo: vamos a tener que ejecutar el proceso `SAMTOOLS_INDEX` cada vez que ejecutemos la prueba para `GATK_HAPLOTYPECALLER`. Sin embargo, tal vez todavía estamos desarrollando el flujo de trabajo y no queremos pre-generar datos de prueba que podríamos tener que cambiar más tarde. El proceso `SAMTOOLS_INDEX` también es muy rápido, así que tal vez los beneficios de pre-generar y almacenar sus salidas son insignificantes. Así es como funciona el método setup.

### 2.1. Generar y ubicar el archivo de prueba

Como anteriormente, primero generamos el esqueleto de archivo:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Salida del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Esto produce el siguiente esqueleto de prueba:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. Mover el archivo de prueba y actualizar la ruta del script

Creamos un directorio para el archivo de prueba co-ubicado con el archivo `main.nf` del módulo:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

Y movemos el archivo esqueleto de prueba allí:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Finalmente, no olvides actualizar la ruta del script:

=== "Después"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Antes"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Proporcionar entradas usando el método setup

Insertamos un bloque `setup` antes del bloque `when`, donde podemos activar una ejecución del proceso `SAMTOOLS_INDEX` en uno de nuestros archivos de entrada originales. Además, recuerda como antes cambiar el nombre de la prueba a algo significativo.

=== "Después"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Antes"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Luego podemos referirnos a la salida de ese proceso en el bloque `when` donde especificamos las entradas de prueba:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Realiza ese cambio y ejecuta la prueba nuevamente:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

También produce un archivo snapshot como anteriormente.

### 2.4. Ejecutar nuevamente y observar falla

Curiosamente, si ejecutas exactamente el mismo comando nuevamente, esta vez la prueba fallará.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

El mensaje de error te dice que hubo diferencias entre los snapshots para las dos ejecuciones; específicamente, los valores md5sum son diferentes para los archivos VCF.

¿Por qué? Para hacer una larga historia corta, la herramienta HaplotypeCaller incluye una marca de tiempo en el encabezado del VCF que es diferente cada vez (por definición).
Como resultado, no podemos simplemente esperar que los archivos tengan md5sums idénticos incluso si tienen contenido idéntico en términos de las llamadas de variantes en sí.

¿Cómo lidiamos con eso?

### 2.5. Usar un método de afirmación de contenido para verificar una variante específica

Una forma de resolver el problema es usar un [tipo diferente de afirmación](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).
En este caso, vamos a verificar contenido específico en lugar de afirmar identidad.
Más exactamente, haremos que la herramienta lea las líneas del archivo VCF y verifique la existencia de líneas específicas.

En la práctica, reemplazamos la segunda afirmación en el bloque `then` de la siguiente manera:

=== "Después"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Antes"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Aquí estamos leyendo el contenido completo del archivo de salida VCF y buscando una coincidencia de contenido, lo cual está bien hacer en un archivo de prueba pequeño, pero no querrías hacer eso en un archivo más grande.
En su lugar, podrías elegir leer líneas específicas.

Este enfoque requiere elegir más cuidadosamente qué queremos usar como 'señal' para probar.
En el lado positivo, se puede usar para probar con gran precisión si una herramienta de análisis puede identificar consistentemente características 'difíciles' (como variantes raras) a medida que experimenta un mayor desarrollo.

### 2.6. Ejecutar nuevamente y observar éxito

Una vez que hemos modificado la prueba de esta manera, podemos ejecutar la prueba múltiples veces, y pasará consistentemente.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Agregar más pruebas

Agrega pruebas similares para las muestras de la madre y el padre:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Ejecutar el comando de prueba

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

Eso completa el plan de prueba básico para este segundo paso en el pipeline. ¡Hacia la tercera y última prueba a nivel de módulo!

### Conclusión

Has aprendido cómo:

1. Probar procesos que dependen de salidas de otros procesos
2. Verificar variantes genómicas específicas en archivos de salida VCF
3. Manejar salidas no determinísticas verificando contenido específico
4. Probar llamado de variantes a través de múltiples muestras

### ¿Qué sigue?

Aprender cómo escribir pruebas que usen datos de prueba pre-generados para el paso de genotipado conjunto.

---

## 3. Usar datos de prueba pre-generados

Para el paso de genotipado conjunto, usaremos un enfoque diferente - usar datos de prueba pre-generados. Esto suele ser preferible para:

1. Procesos complejos con múltiples dependencias
2. Procesos que toman mucho tiempo en ejecutarse
3. Procesos que son parte de un pipeline estable, de producción

### 3.1. Generar datos de prueba

Inspecciona los resultados que generamos al inicio de esta sección:

```bash
tree results_genomics/
```

```console title="Contenido del directorio de resultados"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

El paso de genotipado conjunto necesita los archivos VCF producidos por los pasos del llamador de haplotipo como entradas, junto con los índices. Así que copiemos los resultados que tenemos al directorio de prueba del módulo `jointgenotyping`.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Ahora podemos usar estos archivos como entradas a la prueba que vamos a escribir para el paso de genotipado conjunto.

### 3.2. Generar el esqueleto de archivo de prueba

Como anteriormente, primero generamos el esqueleto de archivo:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Salida del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Esto produce el siguiente esqueleto de prueba:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. Mover el archivo de prueba y actualizar la ruta del script

Esta vez ya tenemos un directorio para pruebas co-ubicado con el archivo `main.nf` del módulo, así que podemos mover el archivo esqueleto de prueba allí:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

Y no olvides actualizar la ruta del script:

=== "Después"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Antes"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Proporcionar entradas

Completa las entradas basándose en las definiciones de entrada del proceso y renombra la prueba en consecuencia:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Usar afirmaciones de contenido

La salida del paso de genotipado conjunto es otro archivo VCF, así que vamos a usar una afirmación de contenido nuevamente.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

Al verificar el contenido de una variante específica en el archivo de salida, esta prueba verifica que:

1. El proceso de genotipado conjunto se ejecuta exitosamente
2. El VCF de salida contiene las tres muestras en el orden correcto
3. Una variante específica es llamada correctamente con:
   - Genotipos precisos para cada muestra (0/1 para el padre, 1/1 para la madre y el hijo)
   - Profundidades de lectura y calidades de genotipo correctas
   - Estadísticas a nivel de población como la frecuencia alélica (AF=0.833)

No hemos tomado un snapshot del archivo completo, pero al verificar una variante específica, podemos estar seguros de que el proceso de genotipado conjunto está funcionando como se espera.

### 3.6. Ejecutar la prueba

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

La prueba pasa, verificando que nuestro proceso de genotipado conjunto correctamente:

1. Combina VCFs de muestras individuales
2. Realiza llamado de variantes conjunto
3. Produce un VCF de múltiples muestras con llamados de genotipo consistentes a través de las ejecuciones

### Conclusión

Sabes cómo:

- Usar resultados generados previamente como entradas para pruebas
- Escribir pruebas usando datos de prueba pre-generados

### ¿Qué sigue?

Agregar una prueba a nivel de flujo de trabajo para verificar que todo el pipeline de llamado de variantes funcione de principio a fin.

---

## 4. Agregar una prueba a nivel de flujo de trabajo

Ahora probaremos el pipeline completo de llamado de variantes, desde archivos BAM hasta genotipos conjuntos. Esto verifica que:

1. Todos los procesos funcionan juntos correctamente
2. Los datos fluyen apropiadamente entre pasos
3. Las llamadas de variantes finales son consistentes

### 4.1. Generar la prueba de flujo de trabajo

Genera un archivo de prueba para el pipeline completo:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Salida del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Esto crea un esqueleto de prueba básico:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Solo corrige el nombre a algo significativo (verás por qué esto es útil en breve).

=== "Después"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Antes"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note "Nota"

    En este caso el archivo de prueba puede permanecer donde `nf-test` lo creó.

### 4.2. Especificar parámetros de entrada

Todavía necesitamos especificar entradas, lo cual se hace ligeramente diferente a nivel de flujo de trabajo comparado con las pruebas a nivel de módulo.
Hay varias formas de hacer esto, incluyendo especificando un perfil.
Sin embargo, una forma más simple es configurar un bloque `params {}` en el archivo `nextflow.config` que `nf-test init` creó originalmente en el directorio `tests`.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

// Directorio de salida para salidas del flujo de trabajo
outputDir = 'results_genomics'

/*
 * Parámetros del pipeline
 */

params {
    // Entrada primaria (archivo de archivos de entrada, uno por línea)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Archivos accesorios
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Nombre base para el archivo de salida final
    cohort_name = "family_trio"
}
```

Cuando ejecutemos la prueba, `nf-test` recogerá este archivo de configuración y extraerá las entradas en consecuencia.

### 4.3. Ejecutar la prueba de flujo de trabajo

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

La prueba pasa, confirmando que nuestro pipeline completo de llamado de variantes:

1. Procesa exitosamente todas las muestras
2. Encadena correctamente todos los pasos

### 4.4. Ejecutar TODAS las pruebas

nf-test tiene un truco más bajo la manga. ¡Podemos ejecutar todas las pruebas a la vez! Modifica el archivo `nf-test.config` para que nf-test busque en cada directorio archivos de nf-test. Puedes hacer esto modificando el parámetro `testsDir`:

=== "Después"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Antes"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Ahora, simplemente podemos ejecutar nf-test y ejecutará _cada prueba individual_ en nuestro repositorio:

```bash
nf-test test
```

??? success "Salida del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

¡8 pruebas en 1 comando! Pasamos mucho tiempo configurando muchas pruebas, pero cuando llegó el momento de ejecutarlas fue muy rápido y fácil. Puedes ver lo útil que es esto cuando se mantiene un pipeline grande, que podría incluir cientos de elementos diferentes. Pasamos tiempo escribiendo pruebas una vez para poder ahorrar tiempo ejecutándolas muchas veces.

Además, ¡podemos automatizar esto! Imagina pruebas ejecutándose cada vez que tú o un colega intenta agregar código nuevo. Así es como aseguramos que nuestros pipelines mantengan un alto estándar.

## Conclusión

Ahora sabes cómo escribir y ejecutar varios tipos de pruebas para tu pipeline de genómica usando nf-test. Este marco de pruebas ayuda a asegurar que tu flujo de trabajo de llamado de variantes produzca resultados consistentes y confiables en diferentes entornos y a medida que realizas cambios en el código.

Has aprendido a probar componentes críticos como:

- El proceso `SAMTOOLS_INDEX` que prepara archivos BAM para el llamado de variantes
- El proceso `GATK_HAPLOTYPECALLER` que identifica variantes en muestras individuales
- El proceso `GATK_JOINTGENOTYPING` que combina llamadas de variantes a través de una cohorte

También has implementado diferentes estrategias de prueba específicas para datos genómicos:

- Verificar que los archivos VCF contengan las llamadas de variantes esperadas a pesar de elementos no determinísticos como marcas de tiempo
- Probar con un conjunto de datos de trío familiar para asegurar la identificación apropiada de variantes a través de muestras relacionadas
- Verificar coordenadas genómicas específicas e información de variantes en tus archivos de salida

Estas habilidades de prueba son esenciales para desarrollar pipelines bioinformáticos robustos que puedan procesar datos genómicos de manera confiable y producir llamadas de variantes precisas. A medida que continúes trabajando con Nextflow para análisis genómico, esta base de pruebas te ayudará a mantener código de alta calidad que produzca resultados científicos confiables.
