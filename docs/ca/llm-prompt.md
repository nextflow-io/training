# Translation Rules for Catalan

The target language for this translation is **Catalan** (`ca`).

## 1. Grammar & Tone

- Use Catalan standard (central) spelling conventions
- For instructions, use formal imperative where appropriate (e.g., "Executeu", "Consulteu")
- Use proper Catalan punctuation (e.g., do not use Spanish opening punctuation; use standard Catalan punctuation)
- Prefer active voice when possible
- Use inclusive greetings when appropriate: "Et donem la benvinguda!" or "Benvingut/da!"

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Catalan (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "El canal d'entrada rep els fitxers..." (translate "channel" to "canal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` -> `// emet una salutacio`

## 3. Code Comments

**Always translate code comments to Catalan.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Catalan translation
params.greeting = "Hello" // defineix la salutacio per defecte
```

## 4. Common Mistakes

Avoid these translation errors specific to Catalan:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Canal.fromPath('*.fastq')
proces FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  versió 24.04.0
executor >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Using non-standard punctuation

```markdown
// Wrong - using Spanish opening punctuation
¿Què segueix?
¡Felicitats!

// Correct - standard Catalan punctuation
Què segueix?
Felicitats!
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code):

| English     | Catalan             | Notes                          |
| ----------- | ------------------- | ------------------------------ |
| channel     | canal               | In prose; keep English in code |
| process     | procés              | In prose; keep English in code |
| workflow    | workflow            | Keep in English (common usage) |
| pipeline    | pipeline            | Keep in English (common usage) |
| directive   | directiva           |                                |
| container   | contenidor          |                                |
| input       | entrada             |                                |
| output      | sortida             |                                |
| task        | tasca               |                                |
| tuple       | tupla               |                                |
| operator    | operador            |                                |
| parameter   | parametre           |                                |
| environment | entorn              |                                |
| directory   | directori           |                                |
| file        | fitxer              |                                |
| sample      | mostra              |                                |
| alignment   | alineament          |                                |
| reference   | referència          |                                |
| training    | formació            |                                |
| workshop    | taller              |                                |
| module      | mòdul               |                                |
| command     | comanda             |                                |
| index       | índex               |                                |
| run         | executar / execució |                                |
| open source | codi obert          |                                |
| community   | comunitat           |                                |

## 6. Admonition Titles

| English  | Catalan     |
| -------- | ----------- |
| Note     | Nota        |
| Tip      | Consell     |
| Warning  | Advertència |
| Exercise | Exercici    |
| Solution | Solució     |
| Example  | Exemple     |

## 7. Section Headers

These recurring section headers should be translated consistently:

| English           | Catalan                  |
| ----------------- | ------------------------ |
| Takeaway          | Conclusio                |
| What's next?      | Què segueix?             |
| Warmup            | Escalfament              |
| Environment Setup | Configuració de l'entorn |
| Getting Started   | Primers passos           |

## 8. Tab Labels

| English | Catalan |
| ------- | ------- |
| After   | Després |
| Before  | Abans   |
| Gitpod  | Gitpod  |
| Local   | Local   |
