# Translation Rules for Portuguese

The target language for this translation is **Brazilian Portuguese** (`pt`).

## 1. Grammar & Tone

- Use informal tone (você instead of o senhor/a senhora)
- Use Brazilian Portuguese spelling conventions (e.g., "arquivo" not "ficheiro", "diretório" not "pasta")
- Prefer active voice when possible
- Use a friendly, encouraging tone where appropriate (e.g., "Parabéns!" for congratulations, "Boas notícias" for good news)

## 2. Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Portuguese (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "O canal de entrada recebe os arquivos..." (translate "channel" to "canal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// emite uma saudação`

## 3. Code Comments

**Always translate code comments to Portuguese.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Portuguese translation
params.greeting = "Hello" // define a saudação padrão
```

## 4. Common Mistakes

Avoid these translation errors specific to Portuguese:

### ❌ Translating code syntax

```groovy
// Wrong - translating Nextflow keywords
Canal.fromPath('*.fastq')
processo FOO { }

// Correct - keep Nextflow keywords in English
Channel.fromPath('*.fastq')
process FOO { }
```

### ❌ Translating console output

Console output shows exactly what users will see and must not be translated:

```console
// Wrong
N E X T F L O W  ~  versão 24.04.0
executor >  local (3)

// Correct - leave exactly as-is
N E X T F L O W  ~  version 24.04.0
executor >  local (3)
```

### ❌ Mixing English keywords in prose

```markdown
// Wrong - keeping "channel" in English when discussing concepts
O channel de entrada recebe os arquivos...

// Correct - translate "channel" to "canal" in prose
O canal de entrada recebe os arquivos...
```

### ❌ Using European Portuguese

```markdown
// Wrong - European Portuguese spelling
O ficheiro está na pasta...

// Correct - Brazilian Portuguese
O arquivo está no diretório...
```

## 5. Terms to Translate

These terms should be translated in prose (but kept in English in code):

| English       | Portuguese              |
| ------------- | ----------------------- |
| channel       | canal / canais          |
| process       | processo / processos    |
| workflow      | fluxo de trabalho       |
| directive     | diretiva / diretivas    |
| container     | contêiner / contêineres |
| input         | entrada                 |
| output        | saída                   |
| task          | tarefa                  |
| tuple         | tupla                   |
| queue channel | canal de fila           |
| value channel | canal de valor          |
| operator      | operador                |
| parameter     | parâmetro               |
| environment   | ambiente                |
| directory     | diretório               |
| file          | arquivo                 |
| sample        | amostra                 |
| alignment     | alinhamento             |
| reference     | referência              |
| training      | treinamento             |
| module        | módulo                  |
| command       | comando                 |
| index         | índice                  |
| run           | executar / execução     |

## 6. Admonition Titles

| English  | Portuguese |
| -------- | ---------- |
| Note     | Nota       |
| Tip      | Dica       |
| Warning  | Aviso      |
| Exercise | Exercício  |
| Solution | Solução    |
| Example  | Exemplo    |

## 7. Section Headers

These recurring section headers should be translated consistently:

| English            | Portuguese                                       |
| ------------------ | ------------------------------------------------ |
| Takeaway           | Conclusão                                        |
| What's next?       | O que vem a seguir?                              |
| Warmup             | Aquecimento                                      |
| Directory contents | Conteúdo do diretório                            |
| Output             | Saída (in prose) / Output (in code block titles) |

## 8. Tab Labels

| English | Portuguese |
| ------- | ---------- |
| After   | Depois     |
| Before  | Antes      |
| Gitpod  | Gitpod     |
| Local   | Local      |

## 9. Common Expressions

Brazilian Portuguese translations for common expressions:

| English            | Portuguese                      |
| ------------------ | ------------------------------- |
| Congratulations!   | Parabéns!                       |
| Good news          | Boas notícias                   |
| That said          | Dito isso                       |
| free bonus         | bônus gratuito                  |
| symbolic link      | link simbólico                  |
| Ta-da! / Voilà!    | Tcharam! (Brazilian colloquial) |
| Take a short break | Faça uma pequena pausa          |
| you've earned it   | você mereceu                    |

## 10. UI Elements

| English       | Portuguese             |
| ------------- | ---------------------- |
| file explorer | explorador de arquivos |
| sidebar       | barra lateral          |
| editor pane   | painel do editor       |
| terminal      | terminal               |
