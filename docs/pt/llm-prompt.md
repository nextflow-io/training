# Translation Rules for Portuguese

The target language for this translation is **Brazilian Portuguese** (`pt`).

## Grammar Preferences

- Use informal tone (você instead of o senhor/a senhora)
- Use Brazilian Portuguese spelling conventions (e.g., "arquivo" not "ficheiro", "diretório" not "pasta")
- Prefer active voice when possible
- Use a friendly, encouraging tone where appropriate (e.g., "Parabéns!" for congratulations, "Boas notícias" for good news)

## Translation Context Rules

**Important distinction**: Some technical terms have different translation rules depending on context:

1. **In code blocks**: Keep ALL Nextflow syntax in English (the code must run)
2. **In code comments**: TRANSLATE comments to Portuguese (they are not executable)
3. **In prose/explanatory text**: Follow the glossary below for translations

For example:

- In prose: "O canal de entrada recebe os arquivos..." (translate "channel" to "canal")
- In code: `channel.fromPath('*.fastq')` (keep "channel" in English)
- In comments: `// emit a greeting` → `// emite uma saudação`

## Code Comments

**Always translate code comments to Portuguese.** Comments are not executable code and should be in the target language for better comprehension.

```groovy
// English original
params.greeting = "Hello" // set default greeting

// Portuguese translation
params.greeting = "Hello" // define a saudação padrão
```

## Glossary

### Terms to Keep in English (DO NOT TRANSLATE)

These terms should remain in English in both code and prose:

- Nextflow
- DSL2
- pipeline / pipelines
- script
- cache / caching
- bucket
- dataflow
- job
- pipe
- quick start
- resume
- shell
- staging
- framework
- log
- debug (noun)
- nf-core
- Docker
- Singularity
- Apptainer
- Conda
- Mamba
- GitHub
- GitLab
- Bitbucket
- Gitpod
- GitHub Codespaces
- Seqera Platform
- Wave
- Fusion
- Groovy
- closure
- shebang
- glob / globbing
- regex

### Nextflow Keywords (keep in English in code only)

These terms appear as Nextflow syntax and must remain in English in code blocks,
but should be translated when used in explanatory prose:

- channel / channels → canal / canais (in prose)
- process / processes → processo / processos (in prose)
- workflow / workflows → fluxo de trabalho / fluxos de trabalho (in prose)
- emit → emitir (in prose)
- take → receber (in prose)
- main → principal (in prose)
- params → parâmetros (in prose)
- input → entrada (in prose)
- output → saída (in prose)

### Operators (keep in English)

Operator names should always remain in English:

- map, filter, collect, flatten, groupTuple, join, combine, mix, merge
- view, first, last, take, branch, multiMap
- splitCsv, splitFastq, splitFasta, splitText

### Directives (keep in English)

Directive names should always remain in English:

- publishDir, container, conda, memory, cpus, time
- errorStrategy, maxRetries, maxErrors, queue, scratch
- storeDir, tag, label, cache, executor

### Data Types (keep in English in code)

- val, path, env, stdin, stdout, tuple, file

### File Formats (keep in English)

- FASTQ, FASTA, BAM, SAM, VCF, BED, GFF, GTF, CSV, TSV, JSON, YAML
- `.nf`, `nextflow.config`, `main.nf`, `modules.nf`, `workflows.nf`

### Terms to Translate (prose only, keep English in code)

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

### Admonition Titles

| English  | Portuguese |
| -------- | ---------- |
| Note     | Nota       |
| Tip      | Dica       |
| Warning  | Aviso      |
| Exercise | Exercício  |
| Solution | Solução    |
| Example  | Exemplo    |

### Section Headers

These recurring section headers should be translated consistently:

| English            | Portuguese                                       |
| ------------------ | ------------------------------------------------ |
| Takeaway           | Conclusão                                        |
| What's next?       | O que vem a seguir?                              |
| Warmup             | Aquecimento                                      |
| Directory contents | Conteúdo do diretório                            |
| Output             | Saída (in prose) / Output (in code block titles) |

### Tab Labels

| English | Portuguese |
| ------- | ---------- |
| After   | Depois     |
| Before  | Antes      |
| Gitpod  | Gitpod     |
| Local   | Local      |

### UI Elements

| English       | Portuguese             |
| ------------- | ---------------------- |
| file explorer | explorador de arquivos |
| sidebar       | barra lateral          |
| editor pane   | painel do editor       |
| terminal      | terminal               |

### Common Expressions

Brazilian Portuguese translations often use these natural expressions:

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
