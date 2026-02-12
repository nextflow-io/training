# Plantilla de Curso nf4_science

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Esta es una plantilla genérica para crear un nuevo curso específico de dominio bajo `nf4_science/`.
Está basada en los patrones comunes encontrados en los cursos de Genómica y RNAseq.

## Cómo usar

1. Copie este directorio y el directorio de scripts correspondiente para crear su curso:

   ```bash
   # Documentación
   cp -r docs/en/docs/nf4_science/_template docs/en/docs/nf4_science/{your_domain}

   # Scripts
   cp -r nf4-science/_template nf4-science/{your_domain}
   ```

2. Renombre el script del workflow:

   ```bash
   mv nf4-science/{your_domain}/_template.nf nf4-science/{your_domain}/{your_domain}.nf
   ```

3. Renombre los archivos de módulos para que coincidan con sus herramientas:

   ```bash
   mv nf4-science/{your_domain}/modules/_tool_a.nf nf4-science/{your_domain}/modules/{tool_a}.nf
   mv nf4-science/{your_domain}/modules/_tool_b.nf nf4-science/{your_domain}/modules/{tool_b}.nf
   ```

4. Busque todos los valores `{PLACEHOLDER}` en la documentación y scripts y reemplácelos con su contenido.

5. Agregue al nav de `mkdocs.yml` (vea el fragmento a continuación).

6. Llene el directorio `data/` con conjuntos de datos de prueba.

7. Complete el directorio `solutions/` con código funcional para cada parte.

## Fragmento de nav para mkdocs.yml

```yaml
- { Your Domain }:
    - nf4_science/{your_domain}/index.md
    - nf4_science/{your_domain}/00_orientation.md
    - nf4_science/{your_domain}/01_method.md
    - nf4_science/{your_domain}/02_single_sample.md
    - nf4_science/{your_domain}/03_multi_sample.md
    - nf4_science/{your_domain}/survey.md
    - nf4_science/{your_domain}/next_steps.md
```

## Referencia de marcadores de posición

### Marcadores de posición a nivel de curso

| Marcador de posición         | Descripción                                | Ejemplo (Genómica)                                                |
| ---------------------------- | ------------------------------------------ | ----------------------------------------------------------------- |
| `{DOMAIN}`                   | Nombre del dominio (mayúsculas y minúsculas) | Genomics                                                          |
| `{DOMAIN_DIR}`               | Nombre del directorio (minúsculas)         | genomics                                                          |
| `{METHOD}`                   | Nombre del método de análisis              | variant calling                                                   |
| `{METHOD_SHORT_DESCRIPTION}` | Descripción del método en una línea        | variant calling with GATK                                         |
| `{ACCESSORY_FILES}`          | Tipos de archivos accesorios utilizados    | index files and reference genome resources                        |

### Marcadores de posición de herramientas

| Marcador de posición                  | Descripción                         | Ejemplo (Genómica)                                  |
| ------------------------------------- | ----------------------------------- | --------------------------------------------------- |
| `{TOOL_A}` / `{TOOL_B}`               | Nombres de visualización de herramientas | Samtools / GATK                                     |
| `{TOOL_A_URL}` / `{TOOL_B_URL}`       | URLs de documentación de herramientas | https://www.htslib.org/                             |
| `{TOOL_A_CONTAINER_URI}`              | URI completo del contenedor         | community.wave.seqera.io/library/samtools:1.20--... |
| `{TOOL_A_MODULE}` / `{TOOL_B_MODULE}` | Nombres de archivos de módulos (sin .nf) | samtools_index / gatk_haplotypecaller               |
| `{TOOL_A_PROCESS_NAME}`               | Nombre del proceso en MAYÚSCULAS    | SAMTOOLS_INDEX                                      |
| `{TOOL_A_COMMAND}`                    | El comando shell a ejecutar         | samtools index '<input_bam>'                        |

### Marcadores de posición de entrada/salida

| Marcador de posición   | Descripción                          | Ejemplo (Genómica)   |
| ---------------------- | ------------------------------------ | -------------------- |
| `{PRIMARY_INPUT_TYPE}` | Tipo de archivo de entrada principal | BAM file             |
| `{PRIMARY_PARAM_NAME}` | Nombre del parámetro de Nextflow     | reads_bam            |
| `{TEST_INPUT_PATH}`    | Ruta a la entrada de prueba en data/ | bam/reads_mother.bam |
| `{FINAL_OUTPUT_TYPE}`  | Tipo de salida final del pipeline    | joint-called VCFs    |

### Marcadores de posición de contenido

| Marcador de posición    | Descripción                                      |
| ----------------------- | ------------------------------------------------ |
| `{PART2_SUMMARY}`       | Resumen de una línea del alcance de la Parte 2  |
| `{AGGREGATION_SUMMARY}` | Resumen de una línea del paso de agregación     |
| `{PROCESS_LIST}`        | Lista separada por comas de nombres de procesos |
| `{TYPEFORM_ID}`         | ID de inserción de Typeform para la encuesta    |

## Estructura pedagógica

La plantilla sigue un arco de tres partes:

1. **Parte 1 (01_method.md)**: Pruebas manuales en contenedores Docker para comprender la metodología
2. **Parte 2 (02_single_sample.md)**: Envolver comandos en Nextflow; muestra única, luego por lotes
3. **Parte 3 (03_multi_sample.md)**: Agregar agregación de múltiples muestras usando operadores de canal

### Convenciones clave

- Use **bloques de código con pestañas Antes/Después** con `hl_lines` para mostrar cambios de código
- Use `??? success "Command output"` para salida esperada colapsable
- Use `??? abstract "Directory contents"` para árboles de directorios colapsables
- Termine cada sección principal con subsecciones **Conclusión** y **¿Qué sigue?**
- Use reglas horizontales `---` para separar secciones numeradas principales
- Proporcione archivos esqueleto (workflow + módulos) para que los estudiantes los completen
- Organice las soluciones por parte (`solutions/part2/`, `solutions/part3/`)

## Estructura de directorios

```
docs/en/docs/nf4_science/{domain}/
├── index.md                    # Course overview with frontmatter
├── 00_orientation.md           # Environment setup
├── 01_method.md                # Manual testing in containers
├── 02_single_sample.md         # Single-sample Nextflow implementation
├── 03_multi_sample.md          # Multi-sample aggregation
├── survey.md                   # Typeform feedback survey
├── next_steps.md               # Course summary and suggestions
└── img/                        # Diagrams (.excalidraw.svg, .png)

nf4-science/{domain}/
├── {domain}.nf                 # Skeleton workflow file
├── nextflow.config             # Minimal config (docker.enabled = true)
├── data/                       # Test datasets and resources
│   └── samplesheet.csv         # Sample metadata
├── modules/                    # Skeleton module files
│   ├── {tool_a}.nf
│   └── {tool_b}.nf
└── solutions/
    ├── part2/                  # Complete Part 2 solution
    │   ├── {domain}-2.nf
    │   ├── nextflow.config
    │   └── modules/
    └── part3/                  # Complete Part 3 solution
        ├── {domain}-3.nf
        ├── nextflow.config
        └── modules/
```
