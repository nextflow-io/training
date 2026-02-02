# Translation Rules for Hindi

The target language for this translation is **Hindi** (`hi`).

## 1. Grammar & Tone

- Use informal tone (तुम instead of आप where appropriate for technical tutorials)
- Use standard Devanagari script
- Prefer active voice when possible
- Use Hindi numerals or Arabic numerals as appropriate for context
- Maintain a conversational yet professional tone

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Hindi (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "इनपुट चैनल फ़ाइलें प्राप्त करता है..." (translate "channel" to "चैनल" or keep in English with Hindi text)
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// एक अभिवादन emit करें`

## 3. Code Comments

**Always translate code comments to Hindi.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Hindi translation
params.greeting = "Hello" // डिफ़ॉल्ट अभिवादन सेट करें
```

## 4. Common Mistakes

Avoid these translation errors specific to Hindi:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
चैनल.fromPath('*.fastq')
प्रोसेस FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  संस्करण 24.04.0
निष्पादक >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Using overly formal आप when तुम is more appropriate

```markdown
// Wrong - too formal for tutorial
आपको workflow चलाना होगा...
आप परिणाम देखेंगे...

// Correct - appropriate informal tone
तुम्हें workflow चलाना होगा...
तुम परिणाम देखोगे...
```

### ❌ Mixing scripts inconsistently

```markdown
// Wrong - inconsistent script usage
यह workflow बहुत powerful है और तुम इसे आसानी से use कर सकते हो।

// Correct - consistent approach (either keep English terms or transliterate consistently)
यह workflow बहुत शक्तिशाली है और तुम इसे आसानी से उपयोग कर सकते हो।
// or
यह वर्कफ़्लो बहुत पावरफुल है और तुम इसे आसानी से यूज़ कर सकते हो।
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code).

Note: Hindi tech writing often keeps English terms with Devanagari transliteration. Both approaches are acceptable - be consistent within a document.

| English     | Hindi      | Notes              |
| ----------- | ---------- | ------------------ |
| channel     | चैनल       | Transliteration    |
| process     | प्रोसेस    | Transliteration    |
| workflow    | वर्कफ़्लो  | Transliteration    |
| pipeline    | पाइपलाइन   | Transliteration    |
| directive   | निर्देश    | Actual translation |
| container   | कंटेनर     | Transliteration    |
| input       | इनपुट      | Transliteration    |
| output      | आउटपुट     | Transliteration    |
| task        | कार्य      | Actual translation |
| tuple       | टपल        | Transliteration    |
| operator    | ऑपरेटर     | Transliteration    |
| parameter   | पैरामीटर   | Transliteration    |
| environment | वातावरण    | Actual translation |
| directory   | डायरेक्टरी | Transliteration    |
| file        | फ़ाइल      | Transliteration    |
| sample      | नमूना      | Actual translation |
| alignment   | संरेखण     | Actual translation |
| reference   | संदर्भ     | Actual translation |
| training    | प्रशिक्षण  | Actual translation |
| module      | मॉड्यूल    | Transliteration    |
| command     | कमांड      | Transliteration    |
| index       | इंडेक्स    | Transliteration    |
| run         | चलाना / रन | Both acceptable    |

## 6. Admonition Titles

| English  | Hindi   |
| -------- | ------- |
| Note     | नोट     |
| Tip      | सुझाव   |
| Warning  | चेतावनी |
| Exercise | अभ्यास  |
| Solution | समाधान  |
| Example  | उदाहरण  |

## 7. Section Headers

| English           | Hindi          |
| ----------------- | -------------- |
| Takeaway          | सारांश         |
| What's next?      | आगे क्या है?   |
| Warmup            | वार्मअप        |
| Environment Setup | पर्यावरण सेटअप |
| Getting Started   | शुरू करना      |

## 8. Tab Labels

| English | Hindi   |
| ------- | ------- |
| After   | बाद में |
| Before  | पहले    |
| Gitpod  | Gitpod  |
| Local   | लोकल    |
