# Translation Rules for German

The target language for this translation is **German** (`de`).

**Language name**: Use lowercase "deutsch" (not "Deutsch") in the language selector to match the style of other Romance languages (español, français, italiano, português).

## 1. Grammar & Tone

- Use informal tone (du instead of Sie)
- Use standard German (Hochdeutsch) spelling conventions
- Prefer active voice when possible
- Use gender-neutral language where appropriate, including gender stars for gender-specific suffixes, e.g. "Entwickler\*innen"

## 2. Writing Style

- Keep sentences short and direct; this is tutorial text, not academic writing
- Avoid overly formal or dated phrasing (e.g., "wissenschaftliche Rechenanforderungen" sounds archaic)
- Prefer modern, conversational German that feels natural for technical documentation
- Break complex sentences into simpler ones rather than using long subordinate clause chains

Example:

- Avoid: "Am Ende dieses Kurses wirst du gut vorbereitet sein, um die nächsten Schritte auf deiner Reise zur Entwicklung reproduzierbarer Workflows für deine wissenschaftlichen Rechenanforderungen anzugehen."
- Better: "Nach Abschluss dieses Kurses bist du bereit, reproduzierbare Workflows für deine eigenen Projekte zu entwickeln."

## 3. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to German (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Der Eingabekanal empfängt die Dateien..." (translate "channel" to "Kanal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// Eine Begrüßung ausgeben`

## 4. Code Comments

**Always translate code comments to German.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// German translation
params.greeting = "Hello" // Standard-Begrüßung festlegen
```

## 5. Common Mistakes

Avoid these translation errors specific to German:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Kanal.fromPath('*.fastq')
Prozess FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  Version 24.04.0
Ausführer >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Overly formal or academic language

```markdown
// Wrong - too formal/academic
Sie werden die Möglichkeit haben, wissenschaftliche Rechenanforderungen zu bewältigen...

// Correct - conversational tutorial style
Du lernst, wie du Workflows für deine Projekte entwickelst...
```

### ❌ Long subordinate clause chains

```markdown
// Wrong - German "Schachtelsätze"
Du wirst sehen, dass, wenn du den Workflow ausführst, der, wie wir bereits besprochen haben, mehrere Prozesse enthält, die Ergebnisse im Ausgabeverzeichnis erscheinen.

// Correct - break into shorter sentences
Führe den Workflow aus. Er enthält mehrere Prozesse, wie wir bereits besprochen haben. Die Ergebnisse erscheinen im Ausgabeverzeichnis.
```

## 6. Terms to Translate

These terms should be translated in prose (but kept in English in code):

| English     | German                 |
| ----------- | ---------------------- |
| channel     | Kanal                  |
| process     | Prozess                |
| workflow    | Workflow               |
| pipeline    | Pipeline               |
| directive   | Direktive              |
| container   | Container              |
| input       | Eingabe                |
| output      | Ausgabe                |
| task        | Aufgabe                |
| tuple       | Tupel                  |
| operator    | Operator               |
| parameter   | Parameter              |
| environment | Umgebung               |
| directory   | Verzeichnis            |
| file        | Datei                  |
| sample      | Probe                  |
| alignment   | Alignment              |
| reference   | Referenz               |
| training    | Training               |
| module      | Modul                  |
| command     | Befehl                 |
| index       | Index                  |
| run         | ausführen / Ausführung |

## 7. Admonition Titles

| English  | German   |
| -------- | -------- |
| Note     | Hinweis  |
| Tip      | Tipp     |
| Warning  | Warnung  |
| Exercise | Übung    |
| Solution | Lösung   |
| Example  | Beispiel |

## 8. Section Headers

These recurring section headers should be translated consistently:

| English           | German              |
| ----------------- | ------------------- |
| Takeaway          | Fazit               |
| What's next?      | Wie geht es weiter? |
| Warmup            | Aufwärmen           |
| Environment Setup | Umgebung einrichten |
| Getting Started   | Erste Schritte      |

## 9. Tab Labels

| English | German |
| ------- | ------ |
| After   | Danach |
| Before  | Vorher |
| Gitpod  | Gitpod |
| Local   | Lokal  |
