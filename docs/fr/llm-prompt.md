# Translation Rules for French

The target language for this translation is **French** (`fr`).

## 1. Grammar & Tone

- Use formal tone (vous instead of tu)
- Follow standard French spelling conventions
- Prefer active voice when possible
- Standard quotation marks ("") are acceptable; French quotation marks (« ») are optional

### Inclusive Writing (Écriture Inclusive)

Consider using inclusive writing with middle dots for gender-neutral terms:

- "débutant·es" (beginners, both genders)
- "prêt·e" (ready, both genders)
- "chercheur·euses" (researchers, both genders)
- "formateur·trice" (trainer, both genders)

This is optional but aligns with modern French writing practices.

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to French (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "Le canal d'entrée reçoit les fichiers..." (translate "channel" to "canal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// émet un message de bienvenue`

## 3. Code Comments

**Always translate code comments to French.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// French translation
params.greeting = "Hello" // définit le message d'accueil par défaut
```

## 4. Common Mistakes

Avoid these translation errors specific to French:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Canal.fromPath('*.fastq')
processus FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  version 24.04.0
exécuteur >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Translating terms commonly kept in English

French tech writing commonly keeps these terms in English:

```markdown
// Acceptable - these terms are standard in French tech docs
Le workflow utilise plusieurs pipelines...
Configurez le container Docker...

// Also acceptable if translating
Le flux de travail utilise plusieurs pipelines...
Configurez le conteneur Docker...
```

### ❌ Missing spaces before punctuation

French requires a non-breaking space before `:`, `;`, `!`, and `?`:

```markdown
// Wrong
Attention: ceci est important!

// Correct (with non-breaking space)
Attention : ceci est important !
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code).

Note: French technical writing commonly keeps "workflow" and "pipeline" in English.

| English         | French               | Notes                              |
| --------------- | -------------------- | ---------------------------------- |
| channel         | canal                | In prose; keep English in code     |
| process         | processus            | In prose; keep English in code     |
| workflow        | workflow             | Keep in English (common in French) |
| pipeline        | pipeline             | Keep in English (common in French) |
| directive       | directive            |                                    |
| container       | conteneur            |                                    |
| input           | entrée               |                                    |
| output          | sortie               |                                    |
| task            | tâche                |                                    |
| tuple           | tuple                |                                    |
| operator        | opérateur            |                                    |
| parameter       | paramètre            |                                    |
| environment     | environnement        |                                    |
| directory       | répertoire           |                                    |
| file            | fichier              |                                    |
| sample          | échantillon          |                                    |
| alignment       | alignement           |                                    |
| reference       | référence            |                                    |
| training        | formation            |                                    |
| workshop        | atelier              |                                    |
| module          | module               |                                    |
| command         | commande             |                                    |
| index           | index                |                                    |
| run             | exécuter / exécution |                                    |
| Side Quests     | Quêtes secondaires   | Creative translation, keep theme   |
| quick reference | référence rapide     |                                    |

## 6. Admonition Titles

| English  | French        |
| -------- | ------------- |
| Note     | Note          |
| Tip      | Astuce        |
| Warning  | Avertissement |
| Exercise | Exercice      |
| Solution | Solution      |
| Example  | Exemple       |

## 7. Section Headers

These recurring section headers should be translated consistently:

| English           | French                           |
| ----------------- | -------------------------------- |
| Takeaway          | À retenir                        |
| What's next?      | Et ensuite ?                     |
| Warmup            | Échauffement                     |
| Environment Setup | Configuration de l'environnement |
| Getting Started   | Premiers pas                     |

## 8. Tab Labels

| English | French |
| ------- | ------ |
| After   | Après  |
| Before  | Avant  |
| Gitpod  | Gitpod |
| Local   | Local  |
