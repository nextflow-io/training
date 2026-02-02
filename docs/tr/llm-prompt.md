# Translation Rules for Turkish

The target language for this translation is **Turkish** (`tr`).

## 1. Grammar & Tone

- Use formal tone (siz instead of sen)
- Follow Turkish Language Association (TDK) spelling conventions
- Prefer active voice when possible
- Pay attention to vowel harmony in suffixes

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Turkish (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

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

## 4. Common Mistakes

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

```console
// Wrong
N E X T F L O W  ~  sürüm 24.04.0
yürütücü >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Incorrect vowel harmony

```markdown
// Wrong - vowel harmony violation
workflow'a koşalım
container'lar

// Correct - proper vowel harmony
workflow'u çalıştıralım
container'lar (or konteynırlar)
```

### ❌ Using informal sen instead of siz

```markdown
// Wrong - too informal
Workflow'u çalıştır. Sonuçları göreceksin.

// Correct - formal siz form
Workflow'u çalıştırın. Sonuçları göreceksiniz.
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code):

| English     | Turkish                  |
| ----------- | ------------------------ |
| channel     | kanal                    |
| process     | süreç                    |
| workflow    | iş akışı                 |
| pipeline    | boru hattı / pipeline    |
| directive   | yönerge                  |
| container   | konteyner                |
| input       | girdi                    |
| output      | çıktı                    |
| task        | görev                    |
| tuple       | demet                    |
| operator    | operatör                 |
| parameter   | parametre                |
| environment | ortam                    |
| directory   | dizin                    |
| file        | dosya                    |
| sample      | örnek                    |
| alignment   | hizalama                 |
| reference   | referans                 |
| training    | eğitim                   |
| module      | modül                    |
| command     | komut                    |
| index       | dizin (index for files)  |
| run         | çalıştırmak / çalıştırma |

## 6. Admonition Titles

| English  | Turkish   |
| -------- | --------- |
| Note     | Not       |
| Tip      | İpucu     |
| Warning  | Uyarı     |
| Exercise | Alıştırma |
| Solution | Çözüm     |
| Example  | Örnek     |

## 7. Section Headers

| English           | Turkish        |
| ----------------- | -------------- |
| Takeaway          | Özet           |
| What's next?      | Sırada ne var? |
| Warmup            | Isınma         |
| Environment Setup | Ortam Kurulumu |
| Getting Started   | Başlarken      |

## 8. Tab Labels

| English | Turkish |
| ------- | ------- |
| After   | Sonra   |
| Before  | Önce    |
| Gitpod  | Gitpod  |
| Local   | Yerel   |
