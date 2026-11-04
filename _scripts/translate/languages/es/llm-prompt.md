# Translation Rules for Spanish

The target language for this translation is **Spanish** (`es`).

## 1. Grammar & Tone

- Use Latin American Spanish spelling conventions
- For instructions, use imperative forms that work with both tú and usted (e.g., "Ejecute", "Consulte")
- Use proper Spanish punctuation: ¡...! and ¿...?
- Prefer active voice when possible
- Consider gender-inclusive greetings: "¡Te damos la bienvenida!" or "¡Bienvenido/a!" rather than just "¡Bienvenido!"

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Spanish (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "El canal de entrada recibe los archivos..." (translate "channel" to "canal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// emite un saludo`

## 3. Code Comments

**Always translate code comments to Spanish.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Spanish translation
params.greeting = "Hello" // define el saludo predeterminado
```

## 4. Common Mistakes

Avoid these translation errors specific to Spanish:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Canal.fromPath('*.fastq')
proceso FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  versión 24.04.0
executor >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Missing Spanish punctuation

```markdown
// Wrong - missing opening punctuation marks
Que sigue?
Felicidades!

// Correct - proper Spanish punctuation
¿Qué sigue?
¡Felicidades!
```

### ❌ Using Spain Spanish for Latin American audience

```markdown
// Wrong - Spain Spanish vocabulary
El ordenador ejecuta el script...

// Correct - Latin American Spanish
La computadora ejecuta el script...
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code):

| English     | Spanish              | Notes                          |
| ----------- | -------------------- | ------------------------------ |
| channel     | canal                | In prose; keep English in code |
| process     | proceso              | In prose; keep English in code |
| workflow    | workflow             | Keep in English (common usage) |
| pipeline    | pipeline             | Keep in English (common usage) |
| directive   | directiva            |                                |
| container   | contenedor           |                                |
| input       | entrada              |                                |
| output      | salida               |                                |
| task        | tarea                |                                |
| tuple       | tupla                |                                |
| operator    | operador             |                                |
| parameter   | parámetro            |                                |
| environment | entorno              |                                |
| directory   | directorio           |                                |
| file        | archivo              |                                |
| sample      | muestra              |                                |
| alignment   | alineamiento         |                                |
| reference   | referencia           |                                |
| training    | capacitación         | Latin American preference      |
| workshop    | taller               |                                |
| module      | módulo               |                                |
| command     | comando              |                                |
| index       | índice               |                                |
| run         | ejecutar / ejecución |                                |
| open source | código abierto       |                                |
| community   | comunidad            |                                |

## 6. Admonition Titles

| English  | Spanish     |
| -------- | ----------- |
| Note     | Nota        |
| Tip      | Consejo     |
| Warning  | Advertencia |
| Exercise | Ejercicio   |
| Solution | Solución    |
| Example  | Ejemplo     |

## 7. Section Headers

These recurring section headers should be translated consistently:

| English           | Spanish                   |
| ----------------- | ------------------------- |
| Takeaway          | Conclusión                |
| What's next?      | ¿Qué sigue?               |
| Warmup            | Calentamiento             |
| Environment Setup | Configuración del entorno |
| Getting Started   | Primeros pasos            |

## 8. Tab Labels

| English | Spanish |
| ------- | ------- |
| After   | Después |
| Before  | Antes   |
| Gitpod  | Gitpod  |
| Local   | Local   |
