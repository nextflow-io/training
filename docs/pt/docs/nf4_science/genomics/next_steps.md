# Resumo do curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Parabéns por concluir o curso de treinamento Nextflow para Genômica! 🎉

## Sua jornada

Você começou executando ferramentas de chamada de variantes manualmente no terminal para entender a metodologia.
Depois, construiu um pipeline Nextflow para uma única amostra para automatizar o processo, escalou-o para lidar com múltiplas amostras em paralelo e adicionou genotipagem conjunta de múltiplas amostras usando operadores de canal.

### O que você construiu

- Um pipeline de chamada de variantes que recebe arquivos BAM como entrada e produz VCFs com chamadas conjuntas como saída.
- Três processos (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` e `GATK_JOINTGENOTYPING`) armazenados em arquivos de módulo separados.
- O pipeline paraleliza automaticamente o processamento de amostras de entrada usando o paradigma de fluxo de dados do Nextflow.
- Os resultados são publicados em um diretório chamado `results/`.

### Habilidades adquiridas

Através deste curso prático, você aprendeu como:

- Escrever um fluxo de trabalho linear para aplicar chamada de variantes a uma única amostra
- Lidar adequadamente com arquivos acessórios, como arquivos de índice e recursos do genoma de referência
- Aproveitar o paradigma de fluxo de dados do Nextflow para paralelizar a chamada de variantes por amostra
- Implementar chamada conjunta de múltiplas amostras usando operadores de canal relevantes

Você está agora preparado para começar a aplicar o Nextflow a fluxos de trabalho de análise genômica em seu próprio trabalho.

## Próximos passos para desenvolver suas habilidades

Aqui estão nossas principais sugestões para o que fazer a seguir:

- Aplique o Nextflow a outros casos de uso de análise científica com [Nextflow for Science](../index.md)
- Comece com o nf-core através do [Hello nf-core](../../hello_nf-core/index.md)
- Explore recursos mais avançados do Nextflow com as [Side Quests](../../side_quests/index.md)

Por fim, recomendamos que você dê uma olhada na [**Seqera Platform**](https://seqera.io/), uma plataforma baseada em nuvem desenvolvida pelos criadores do Nextflow que torna ainda mais fácil lançar e gerenciar seus fluxos de trabalho, além de gerenciar seus dados e executar análises interativamente em qualquer ambiente.

## Obtendo ajuda

Para recursos de ajuda e suporte da comunidade, consulte a [página de Ajuda](../../help.md).

## Pesquisa de feedback

Antes de prosseguir, por favor, reserve um minuto para completar a pesquisa do curso! Seu feedback nos ajuda a melhorar nossos materiais de treinamento para todos.

[Responder à pesquisa :material-arrow-right:](survey.md){ .md-button .md-button--primary }
