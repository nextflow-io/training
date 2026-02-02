# Translation Rules for Polish

The target language for this translation is **Polish** (`pl`).

## 1. Grammar & Tone

- Use informal tone (Ty for one person, Wy for multiple, not Pan/Pani)
- Always capitalize second-person pronouns, including possessives (Twój, Wasz); do not capitalize other pronouns
- Follow standard Polish spelling conventions
- Prefer active voice when possible
- Use natural Polish syntax and sentence structure; feel free to rephrase the English text wherever it improves clarity and naturalness
- Always use correct declension, even with English terms. Use apostrophe endings (e.g. "workflow'u") where necessary
- Use inanimate genders when referring to programming terms. Exception: use virile (animate masculine) to refer to Nextflow (i.e. Nextflow'a)
- Avoid repeating nouns, verbs, adjectives, and adverbs within a sentence or in two consecutive sentences wherever possible.
- Avoid sentences with more than 3 verbs. Split long sentences as required.

### Terms with apostrophe endings

Some English terms require Polish declension with apostrophe:

- pipeline (pipeline'u, pipeline'owi, etc.)
- workflow (workflow'u, workflow'owi, etc.)

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Polish (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations; use Polish words where available unless referring to code keywords verbatim; always use Markdown inline code for keywords

For example:

- In prose: "Kanał wejściowy otrzymuje pliki..." (translate "channel" to "kanał")
- In prose with verbatim syntax: "kanały można tworzyć przy pomocy przestrzeni nazw `channel`" (keep `channel` in code formatting)
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// Wyemituj powitanie`

## 3. Code Comments

**Always translate code comments to Polish.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Polish translation
params.greeting = "Hello" // ustaw domyślne powitanie
```

## 4. Common Mistakes

Avoid these translation errors specific to Polish:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Kanał.fromPath('*.fastq')
proces FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  wersja 24.04.0
wykonawca >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Missing declension on English terms

```markdown
// Wrong - no declension
Zainstaluj workflow na swoim komputerze.
Pracujemy z pipeline.

// Correct - with proper Polish declension
Zainstaluj workflow'a na swoim komputerze.
Pracujemy z pipeline'em.
```

### ❌ Using formal Pan/Pani instead of Ty

```markdown
// Wrong - too formal
Proszę uruchomić workflow...
Pan/Pani zobaczy wyniki...

// Correct - informal Ty form
Uruchom workflow'a...
Zobaczysz wyniki...
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code).

| English         | Polish                   | Notes                        |
| --------------- | ------------------------ | ---------------------------- |
| channel         | kanał / kanały           | In prose                     |
| process         | proces / procesy         | In prose                     |
| workflow        | workflow / workflow'y    | Keep English, add declension |
| pipeline        | pipeline / pipeline'y    | Keep English, add declension |
| queue channel   | kanał kolejki            |                              |
| value channel   | kanał wartości           |                              |
| script          | skrypt                   |                              |
| shell           | powłoka                  |                              |
| params          | parametry                | In prose                     |
| directive       | dyrektywa                |                              |
| container       | kontener                 |                              |
| input           | wejście                  |                              |
| output          | wyjście                  |                              |
| task            | zadanie                  |                              |
| tuple           | krotka                   |                              |
| operator        | operator                 |                              |
| parameter       | parametr                 |                              |
| environment     | środowisko               |                              |
| directory       | katalog                  |                              |
| file            | plik                     |                              |
| sample          | próbka                   |                              |
| alignment       | dopasowanie              |                              |
| reference       | referencja               |                              |
| training        | szkolenie                |                              |
| module          | moduł                    |                              |
| command         | polecenie                |                              |
| index           | indeks                   |                              |
| run             | uruchomić / uruchomienie |                              |
| parallelization | paralelizacja            |                              |
| parallelize     | paralelizować            | Never "równoleglić"          |
| domain          | dziedzina                | domain of knowledge          |
| journey (fig.)  | droga / ścieżka          |                              |

## 6. Admonition Titles

| English  | Polish      |
| -------- | ----------- |
| Note     | Uwaga       |
| Tip      | Wskazówka   |
| Warning  | Ostrzeżenie |
| Exercise | Ćwiczenie   |
| Solution | Rozwiązanie |
| Example  | Przykład    |

## 7. Section Headers

| English           | Polish                  |
| ----------------- | ----------------------- |
| Takeaway          | Podsumowanie            |
| What's next?      | Co dalej?               |
| Warmup            | Rozgrzewka              |
| Environment Setup | Konfiguracja środowiska |
| Getting Started   | Pierwsze kroki          |

## 8. Tab Labels

| English | Polish   |
| ------- | -------- |
| After   | Po       |
| Before  | Przed    |
| Gitpod  | Gitpod   |
| Local   | Lokalnie |
