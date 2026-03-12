# Translation Rules for Turkish

The target language for this translation is **Turkish** (`tr`).

## 1. Grammar & Tone

- Use formal tone (`siz` instead of `sen`) throughout. Use imperative with `-in`/`-ın`/`-un`/`-ün` suffix, not bare verb stem.
- Follow Turkish Language Association (TDK) spelling conventions (kurallar).
- Prefer active voice when possible.
- Pay strict attention to vowel harmony in all suffixes and borrowed-word inflections.
- **Tone**: Maintain the educational future tense (`öğreneceksiniz`, `yapacaksınız`) found in the source, but ensure the phrasing is natural and free of "translation smell" (anlamsız kelime dizilimleri).

### 1.1. Punctuation: Comma, Semicolon, and Period

Turkish uses commas, semicolons, and periods differently from English. Do not mirror English punctuation mechanically.

**Comma (virgül)** — use sparingly:

- No comma before "ve" (and) in lists:
  ```
  Wrong:  modüler, ölçeklenebilir, ve taşınabilir
  Correct: modüler, ölçeklenebilir ve taşınabilir
  ```
- No comma before participial (sıfat-fiil) phrases — these are embedded in Turkish with no comma:
  ```
  Wrong:  Nextflow kullanılarak oluşturulmuş, küratörlüğü yapılmış pipeline'lar
  Correct: Nextflow kullanılarak oluşturulmuş pipeline'lar
  ```
- Comma before clause connectors like "bu da", "bu nedenle", "bu sayede" is correct:
  ```
  Correct: Pipeline'lar taşınabilir olacak şekilde tasarlanmıştır, bu da
           araştırmacıların kendi verileriyle kolayca çalıştırmasını sağlar.
  ```

**Semicolon (noktalı virgül)** — use to join closely related independent clauses, or to separate list items that already contain commas:

```
Correct: Pipeline'lar modüler ve ölçeklenebilirdir; araştırmacılar bunları
         kendi hesaplama kaynaklarıyla kolayca çalıştırabilir.
```

Do not use a semicolon where a period would be cleaner.

**Period (nokta)** — prefer over comma chains:

When an English sentence is long and comma-heavy, split it into shorter Turkish sentences with periods rather than preserving the comma chain.

```
Wrong (comma chain):
  nf-core, açık geliştirmeyi, test etmeyi ve akran değerlendirmesini teşvik eden,
  topluluk tarafından geliştirilen pipeline'lar sunan bir projedir.

Correct (split into sentences):
  nf-core, topluluk tarafından geliştirilen pipeline'lar sunan bir projedir.
  Proje; açık geliştirmeyi, test etmeyi ve akran değerlendirmesini teşvik eder.
```

### 1.2. Sentence Structure (Word Order)

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

### 1.3. Handling Long Introductory Phrases and Complex Sentences

When translating structures like "You will be introduced to [X]...":

- If [X] is very long, start with "X ile tanışacaksınız" or move the main subject to the beginning as an independent sentence.
- If the source sentence has too many subordinate clauses, split it into two meaningful parts instead of forcing a single complex inverted sentence.
  - **Preferred structure**: "X ile tanışacaksınız. Bu oluşum aynı zamanda Y ve Z'yi de içerir."

### 1.4. Sentence-Initial '-ing' Phrases

When translating sentences starting with long '-ing' phrases (e.g., "Working through..."), do not detach them from the main clause. Use Turkish adverbial participles (zarf-fiil: `-arak`, `-erek`) to connect them fluently to the main verb.

- **Example**: "Working through practical examples..." → "Pratik örnekler üzerinden çalışarak..."

### 1.5. Handling "Because" in Long Sentences

If a sentence is connected by "because" and is very long, split it into two sentences. End the first sentence with a period and start the new sentence with "Zira," or "Bunun nedeni," to improve readability.

### 1.6. Actions in Summaries/Lists

When translating 'Summary/Özet' sections in tables or lists, use a consistent mood. If an action is described (Run, Explore, Create), use either the infinitive (`-mak`/`-mek`) or simple present tense.

- **Example**: "Mevcut bir pipeline'ı çalıştırmak" or "...çalıştırılır".

### 1.7. Specific Idioms and Definitions

- **Definitions**: Avoid clunky passive structures like "X olarak adlandırılır". Instead, use more natural phrasing like "X adı verilen...".
  - *Example*: Instead of "Metro haritası olarak adlandırılır", use "metro haritası (subway map) adı verilen görsel temsil".
- **Direct Narrative**: Instead of structures like "Buna şu denir ve bu şudur", use fluid engineering language like "X olarak adlandırılan bu yapı...".
- **"And so forth"**: Translate as "ve benzerleri" to complete the sentence naturally.
- **"Out of the way"**: Translate as "ayrı bir yerde" or "erişimi doğrudan olmayan" depending on context.
- **"Poke around"**: Translate as "kodu incelemek" or "kurcalamak".
- **Fluency**: Use natural equivalents like "Sizin belirleyeceğiniz bir isim" instead of translation-heavy phrases like "Bunun için uydurulabilecek bir ad".
- **Placeholders**: When describing symbols like `< >` in command line examples, refer to them as "yer tutucular".
- **Engineering Terminology**: If a term refers to both a biological process and a software component, prioritize engineering terminology (e.g., "Read trimming" → "Okuma kırpma işlemi").

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run).
2. **In code comments**: TRANSLATE comments to Turkish (they are not executable).
3. **In prose/explanatory text**: Follow the glossary below for translations.
4. **Filenames**: Never translate technical filenames (e.g., `main.nf`, `demo.nf`).

For example:

- In prose: "Giriş kanalı dosyaları alır..." (translate "channel" to "kanal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// bir selamlama yayınla`

### Specific Nextflow Keywords

- **`take`, `main`, `topic`**: Keep these keywords in English within code blocks.
- **Prose Translation**: When explaining these structures in text, use Turkish equivalents:
  - `take` / `main` → "iş akışı girişi" or "kanal tanımı".
  - `emit` vs `publish`: Use "yayınlamak" for internal data flow (`emit`) and "yayımlamak" for writing files externally (`publish`).

### Groovy Data Types

- Keep specific Groovy/Nextflow types in English to maintain technical accuracy: **Closure**, **Map**, **List**, **Channel**, **Value**, **Queue**.
- Translate general references: "list of files" → "dosya listesi".

### Code Block Titles

Translate descriptive titles in code blocks, but keep filenames intact.
- `title="Syntax"` → `title="Sözdizimi"`
- `title="main.nf"` → `title="main.nf"`

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

| English term | Last vowel (sound) | Suffix example        | Turkish form  |
| ------------ | ------------------ | --------------------- | ------------- |
| workflow     | o → back-rounded   | locative: -'da/-'de   | workflow'da   |
| channel      | e → front          | ablative: -'den/-'dan | channel'dan   |
| pipeline     | a → back           | dative: -'a/-'e       | pipeline'a    |
| process      | e → front          | genitive: -'in/-'ın   | process'in    |
| script       | i → front          | locative: -'de/-'da   | script'te     |
| container    | e → front          | plural: -'ler/-'lar   | container'lar |

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

| English       | Turkish                  |
| ------------- | ------------------------ |
| channel       | kanal                    |
| process       | süreç                    |
| workflow      | iş akışı                 |
| pipeline      | pipeline / iş hattı      |
| directive     | yönerge                  |
| container     | konteyner                |
| input         | girdi                    |
| output        | çıktı                    |
| task          | görev                    |
| tuple         | demet                    |
| operator      | operatör                 |
| parameter     | parametre                |
| environment   | ortam                    |
| directory     | dizin                    |
| file          | dosya                    |
| sample        | örnek                    |
| alignment     | hizalama                 |
| reference     | referans                 |
| training      | eğitim                   |
| module        | modül                    |
| command       | komut                    |
| index         | indeks (BWA/Samtools) / indis (dizi) |
| run           | çalıştırmak / çalıştırma |
| conventions   | kurallar / gelenekler    |
| best practices| en iyi uygulamalar       |
| portable      | taşınabilir              |
| robust        | dayanıklı / güçlü        |
| validation    | doğrulama                |
| reproducibility | tekrarlanabilirlik     |
| script        | betik / script           |
| executor      | yürütücü                 |
| configuration | yapılandırma             |
| resume        | devam ettirme            |
| publish       | yayımlamak (dosya çıkışı)|
| emit          | yayınlamak (veri akışı)  |
| collect       | toplamak                 |
| wrap          | kapsamak                 |
| overview      | giriş                    |
| peer review   | hakem değerlendirmesi    |
| curated set   | özenle seçilmiş / denetlenmiş koleksiyon |
| community effort | topluluk girişimi / topluluk çabası |
| retrieve      | çekmek / almak / indirmek|
| pull          | çekmek / indirmek        |
| plain code    | yalın kod / standart kod |

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
| ----------------- | ------------------------ |
| Takeaway          | Özetle / Öğrendiklerimiz |
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
