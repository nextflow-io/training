# Translation Rules for Korean

The target language for this translation is **Korean** (`ko`).

## 1. Grammar & Tone

- Use formal polite speech level (합쇼체/하십시오체)
- Sentence endings: "~입니다", "~습니다", "~하세요", "~수 있습니다"
- Maintain consistent honorific usage throughout
- Follow standard Korean spelling conventions (한글 맞춤법)
- Preserve personality and humor where present (e.g., "그동안 차 한 잔을 준비하거나" for suggesting a tea break)
- **Avoid run-on sentences**: Split long, complex sentences into shorter ones or use commas (,) appropriately to improve readability
- **Natural phrasing**: When translating "introduce information", prefer "학습하다" (learn) or "다루다" (cover) instead of direct translation like "소개하다"

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Korean (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "입력 채널이 파일을 받습니다..." (translate "channel" to "채널")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// 인사말을 내보냅니다`

## 3. Code Comments

**Always translate code comments to Korean.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Korean translation
params.greeting = "Hello" // 기본 인사말 설정
```

## 4. English Terms with Korean Particles

Korean uses particles attached to words. When English terms appear in Korean text, attach Korean particles directly:

- "README.md 파일" (README.md file)
- "VSCode IDE에서" (in VSCode IDE)
- "Docker를" (Docker + object marker)

## 5. UI Element Pattern

For UI elements, consider including English in parentheses to help users match translations to what they see in the (often English) interface:

- "**사이드바(Sidebar)**"
- "**파일 탐색기(File Explorer)**"

## 6. Image Alt Text

Translate image alt text to Korean for accessibility:

- "GitHub Codespace 세션 리스트"
- "GitHub Codespaces 환영 화면"

## 7. Common Mistakes

Avoid these translation errors specific to Korean:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
채널.fromPath('*.fastq')
프로세스 FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  버전 24.04.0
실행자 >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Inconsistent honorific levels

```markdown
// Wrong - mixing speech levels
워크플로우를 실행하세요. 결과가 나온다.

// Correct - consistent formal polite level
워크플로우를 실행하세요. 결과가 나옵니다.
```

### ❌ Using overly academic translations

```markdown
// Wrong - sounds like political independence
독립 실행형 스크립트...

// Correct - appropriate technical meaning
단독 실행형 스크립트...
```

## 8. Terms to Translate

These terms should be translated in prose (but kept in English in code).

Note: Some terms use transliteration (음차) while others use actual Korean translations. Follow these patterns:

### Transliterations (음차)

| English   | Korean     | Notes                      |
| --------- | ---------- | -------------------------- |
| workflow  | 워크플로우 | Common transliteration     |
| pipeline  | 파이프라인 | Common transliteration     |
| terminal  | 터미널     | Common transliteration     |
| container | 컨테이너   | Common transliteration     |
| sidebar   | 사이드바   | Common transliteration     |
| module    | 모듈       | Common transliteration     |
| sample    | 샘플       | Common transliteration     |
| tuple     | 튜플       | Common transliteration     |
| index     | 인덱스     | Common transliteration     |
| core      | 코어       | Common transliteration     |
| directory | 디렉토리   | Transliteration acceptable |
| file      | 파일       | Transliteration            |

### Actual Translations

| English       | Korean      | Notes                                                                   |
| ------------- | ----------- | ----------------------------------------------------------------------- |
| archived      | 아카이브된  | Transliteration, not 보관된 (standard for tech contexts)                |
| channel       | 채널        | Transliteration but treated as translation                              |
| dependency    | 의존성      | Not 종속성 (의존성 is standard in Korean tech community)                |
| gradually     | 단계적으로  | Not 점진적으로 (단계적으로 = step-by-step, better for tutorials)        |
| parser        | 분석기      | Not 파서 (more formal, e.g., AWS uses 로그 분석기 for "Log parser")     |
| process       | 프로세스    | Transliteration but treated as translation                              |
| sleep         | 중지(sleep) | Not 절전 (power saving). Context: "Codespace goes to sleep" = 중지 상태 |
| standalone    | 단독        | Not 독립 (political independence)                                       |
| mini-course   | 단기 과정   | Not 미니 과정 (too casual)                                              |
| wrap          | 적용하다    | Not 래핑하다 (sounds like physical packaging)                           |
| training      | 교육        | Actual translation                                                      |
| environment   | 환경        | Actual translation                                                      |
| repository    | 저장소      | Actual translation                                                      |
| file explorer | 파일 탐색기 | Mixed                                                                   |
| main editor   | 메인 편집기 | Mixed                                                                   |
| alignment     | 정렬        | Actual translation                                                      |
| command       | 명령        | Actual translation                                                      |
| directive     | 지시문      | Actual translation                                                      |
| input         | 입력        | Actual translation                                                      |
| output        | 출력        | Actual translation                                                      |
| operator      | 연산자      | Actual translation                                                      |
| parameter     | 매개변수    | Actual translation                                                      |
| reference     | 참조        | Actual translation                                                      |
| run           | 실행        | Actual translation                                                      |
| task          | 작업        | Actual translation                                                      |
| lowercase     | 소문자      | Actual translation                                                      |
| materials     | 자료        | Actual translation                                                      |
| course        | 과정        | Actual translation                                                      |

### Parsing Note

Note: "parser" is translated as 분석기, but "parsing" is typically transliterated as 파싱.

## 9. Admonition Titles

| English  | Korean |
| -------- | ------ |
| Note     | 참고   |
| Tip      | 팁     |
| Warning  | 경고   |
| Exercise | 연습   |
| Solution | 해결책 |
| Example  | 예제   |

## 10. Section Headers

| English           | Korean    |
| ----------------- | --------- |
| Takeaway          | 핵심 정리 |
| What's next?      | 다음 단계 |
| Warmup            | 준비 운동 |
| Environment Setup | 환경 설정 |
| Getting Started   | 시작하기  |

## 11. Tab Labels

| English | Korean |
| ------- | ------ |
| After   | 후     |
| Before  | 전     |
| Gitpod  | Gitpod |
| Local   | 로컬   |
