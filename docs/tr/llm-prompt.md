# Translation Rules for Turkish

The target language for this translation is **Turkish** (`tr`).

## 1. Grammar & Tone

- Use formal tone (`siz` instead of `sen`) throughout. Use imperative with `-in`/`-ın`/`-un`/`-ün` suffix, not bare verb stem.
- Follow Turkish Language Association (TDK) spelling conventions (kurallar).
- Prefer active voice when possible.
- Pay strict attention to vowel harmony in all suffixes and borrowed-word inflections.

### 1.1. Sentence Structure (Word Order)

Turkish is a **verb-final, left-branching** language. Subordinate and participial clauses that appear at the **end** of an English sentence must be moved to **before** the main clause in Turkish.

**Pattern**: English `[main clause], [participial/subordinate phrase]` → Turkish `[participial/subordinate phrase], [main clause]`

```
Wrong (English word order preserved):
  Mevcut bir iş akışını nf-core şablon iskeletine uyarlayın,
  Hello Nextflow kursunda üretilen basit iş akışından başlayarak.

Correct (Turkish word order):
  Hello Nextflow kursunda üretilen basit iş akışından başlayarak,
  mevcut bir iş akışını nf-core şablon iskeletine uyarlayın.
```

This applies to all trailing English phrases such as:
- "starting from..." → "...dan/den başlayarak, [ana cümle]"
- "using..." → "...kullanarak, [ana cümle]"
- "by running..." → "...çalıştırarak, [ana cümle]"
- "before doing X..." → "X yapmadan önce, [ana cümle]"

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run).
2. **In code comments**: TRANSLATE comments to Turkish (they are not executable).
3. **In prose/explanatory text**: Follow the glossary below for translations.

For example:

- In prose: "Giriş kanalı dosyaları alır..." (translate "channel" to "kanal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// bir selamlama yayınla`

## 3. Code Comments

**Always translate code comments to Turkish.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Turkish translation
params.greeting = "Hello" // varsayılan selamlamayı ayarla
```

## 4. Apostrophes with English Words

When an English technical term is used in a Turkish sentence, attach Turkish suffixes with an apostrophe. Apply vowel harmony based on the **last vowel sound** of the English word as it is pronounced in Turkish.

| English term | Last vowel (sound) | Suffix example       | Turkish form        |
| ------------ | ------------------ | -------------------- | ------------------- |
| workflow     | o → back-rounded   | locative: -'da/-'de  | workflow'da         |
| channel      | e → front          | ablative: -'den/-'dan| channel'dan         |
| pipeline     | a → back           | dative: -'a/-'e      | pipeline'a          |
| process      | e → front          | genitive: -'in/-'ın  | process'in          |
| script       | i → front          | locative: -'de/-'da  | script'te           |
| container    | e → front          | plural: -'ler/-'lar  | container'lar       |

**Rule**: Never add a suffix directly without an apostrophe (e.g., `workflowda` ✗ → `workflow'da` ✓).

## 5. Common Mistakes

Avoid these translation errors specific to Turkish:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Kanal.fromPath('*.fastq')
süreç FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```
Wrong:
  N E X T F L O W  ~  sürüm 24.04.0
  yürütücü >  local (3)

Correct (leave exactly as-is):
  N E X T F L O W  ~  version 24.04.0
  executor >  local (3)
```

### ❌ Incorrect vowel harmony

```
Wrong (vowel harmony violation):
  workflow'a koşalım
  container'lar

Correct (proper vowel harmony):
  workflow'u çalıştıralım
  container'lar
```

### ❌ Using informal sen instead of siz

```
Wrong (too informal):
  Workflow'u çalıştır. Sonuçları göreceksin.

Correct (formal siz form):
  Workflow'u çalıştırın. Sonuçları göreceksiniz.
```

### ❌ Compound nouns without possessive suffix

In Turkish, noun + verbal noun compounds require the third-person possessive suffix (`-sı/-si/-su/-sü`). Without it, the phrase sounds unnatural.

```
Wrong (missing possessive suffix):
  girdi doğrulama       (input validation)
  hata yönetim          (error management)
  süreç çalıştırma      (process execution)

Correct (with possessive suffix):
  girdi doğrulaması
  hata yönetimi
  süreç çalıştırması
```

Prefer verbal rephrasing over literal noun compounds. Expand vague nouns ("input") into what they actually refer to in context ("parameters", "data files", etc.):

```
Instead of: "girdi doğrulaması uygulayın"
Prefer:      "komut satırı parametrelerini ve veri dosyalarını doğrulayın"
```

### ❌ Translating terms that should stay in English

File paths, flag names, and CLI options must never be translated:

```
Wrong:
  `--çıktı-dizini` bayrağını kullanın

Correct:
  `--outdir` bayrağını kullanın
```

## 6. Terms to Translate

These terms should be translated in prose (but kept in English in code):

| English      | Turkish                  |
| ------------ | ------------------------ |
| channel      | kanal                    |
| process      | süreç                    |
| workflow     | iş akışı                 |
| pipeline     | boru hattı / pipeline    |
| directive    | yönerge                  |
| container    | konteyner                |
| input        | girdi                    |
| output       | çıktı                    |
| task         | görev                    |
| tuple        | demet                    |
| operator     | operatör                 |
| parameter    | parametre                |
| environment  | ortam                    |
| directory    | dizin                    |
| file         | dosya                    |
| sample       | örnek                    |
| alignment    | hizalama                 |
| reference    | referans                 |
| training     | eğitim                   |
| module       | modül                    |
| command      | komut                    |
| index        | dizin (index for files)  |
| run          | çalıştırmak / çalıştırma |
| conventions  | kurallar                 |
| script       | betik / script           |
| executor     | yürütücü                 |
| configuration| yapılandırma             |
| resume       | devam ettirme            |
| publish      | yayımlamak               |
| emit         | yayınlamak               |
| collect      | toplamak                 |
| wrap         | kapsamak                 |
| overview     | giriş                    | 


## 7. Admonition Titles

| English  | Turkish   |
| -------- | --------- |
| Note     | Not       |
| Tip      | İpucu     |
| Warning  | Uyarı     |
| Exercise | Alıştırma |
| Solution | Çözüm     |
| Example  | Örnek     |

## 8. Section Headers

| English           | Turkish        |
| ----------------- | -------------- |
| Takeaway          | Özet           |
| What's next?      | Sırada ne var? |
| Warmup            | Isınma         |
| Environment Setup | Ortam Kurulumu |
| Getting Started   | Başlarken      |

## 9. Tab Labels

| English | Turkish |
| ------- | ------- |
| After   | Sonra   |
| Before  | Önce    |
| Gitpod  | Gitpod  |
| Local   | Yerel   |