# Translation Rules for Italian

The target language for this translation is **Italian** (`it`).

## 1. Grammar & Tone

- Use first person plural inclusive (noi/voi) rather than formal Lei - this creates a more collaborative, educational tone
- Example: "ci addentriamo", "apriamo", "eseguiamo", "vediamo insieme"
- Follow standard Italian spelling conventions
- Prefer active voice when possible
- Use natural Italian expressions (e.g., "Voilà!" is acceptable as it's commonly used in Italian)

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Italian (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Il canale di input riceve i file..." (translate "channel" to "canale")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// emette un saluto`

## 3. Code Comments

**Always translate code comments to Italian.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Italian translation
params.greeting = "Hello" // imposta il saluto predefinito
```

## 4. Common Mistakes

Avoid these translation errors specific to Italian:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Canale.fromPath('*.fastq')
processo FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  versione 24.04.0
esecutore >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Using formal Lei when voi/noi is more appropriate

```markdown
// Wrong - too formal for tutorial
Lei dovrà eseguire il comando seguente...

// Correct - collaborative, educational tone
Eseguiamo il comando seguente...
// or
Eseguite il comando seguente...
```

### ❌ Translating English loanwords common in Italian tech

Italian tech writing commonly keeps these in English:

```markdown
// Acceptable - these are standard in Italian tech docs
Il file di output...
La directory di lavoro...
Il container Docker...

// "Archivio" and "cartella" are NOT preferred
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code).

Note: Italian tech writing commonly keeps "file", "directory", "output", "input", and "container" in English.

| English         | Italian                          | Notes                                 |
| --------------- | -------------------------------- | ------------------------------------- |
| channel         | canale / canali                  | Translate in prose                    |
| process         | processo                         | Translate in prose                    |
| workflow        | flusso di lavoro                 | Always translate in prose             |
| pipeline        | pipeline                         | Keep in English                       |
| channel factory | fabbrica di canali               |                                       |
| queue channel   | canale di coda                   |                                       |
| work directory  | directory di lavoro              | Note: "directory" stays English       |
| directive       | direttiva                        |                                       |
| container       | container                        | Keep in English (standard in Italian) |
| input           | input                            | Keep in English (standard in Italian) |
| output          | output                           | Keep in English (standard in Italian) |
| file            | file                             | Keep in English (standard in Italian) |
| directory       | directory                        | Keep in English (standard in Italian) |
| task            | attività                         |                                       |
| tuple           | tupla                            |                                       |
| operator        | operatore                        |                                       |
| parameter       | parametro                        |                                       |
| qualifier       | qualificatore                    |                                       |
| environment     | ambiente                         |                                       |
| sample          | campione                         |                                       |
| alignment       | allineamento                     |                                       |
| reference       | riferimento                      |                                       |
| training        | formazione                       |                                       |
| module          | modulo                           |                                       |
| command         | comando                          |                                       |
| index           | indice                           |                                       |
| run             | eseguire / esecuzione / lanciare |                                       |
| greeting        | saluto                           | In context of Hello World examples    |
| log             | registri / log                   | Both acceptable                       |

## 6. Admonition Titles

| English  | Italian             |
| -------- | ------------------- |
| Note     | Nota                |
| Tip      | Suggerimento        |
| Warning  | Avviso / Attenzione |
| Exercise | Esercizio           |
| Solution | Soluzione           |
| Example  | Esempio             |

## 7. Section Headers

These recurring section headers should be translated consistently:

| English            | Italian            | Notes                                                   |
| ------------------ | ------------------ | ------------------------------------------------------- |
| Takeaway           | Takeaway           | Keep in English - common in Italian educational content |
| What's next?       | Cosa c'è dopo?     |                                                         |
| Warmup             | Riscaldamento      |                                                         |
| Directory contents | Directory contents | Keep in code block titles                               |
| Output             | Output             | Keep in code block titles                               |

## 8. Tab Labels

| English | Italian |
| ------- | ------- |
| After   | Dopo    |
| Before  | Prima   |
| Gitpod  | Gitpod  |
| Local   | Locale  |

## 9. Common Expressions

| English            | Italian                        |
| ------------------ | ------------------------------ |
| Congratulations!   | Congratulazioni!               |
| Good news          | Buone notizie                  |
| Voilà!             | Voilà! (acceptable in Italian) |
| Take a short break | Prendetevi una piccola pausa   |
| you've earned it   | ve la siete meritata           |
