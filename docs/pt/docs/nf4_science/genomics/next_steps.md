# Resumo do Curso

Parabéns por concluir o curso de treinamento Nextflow for Genomics! 🎉

## Sua jornada

Você começou executando ferramentas de chamada de variantes manualmente no terminal para entender a metodologia.
Então construiu um pipeline Nextflow de amostra única para automatizar o processo, escalou-o para lidar com múltiplas amostras em paralelo e adicionou genotipagem conjunta de múltiplas amostras usando operadores de canal.

### O que você construiu

- Um pipeline de chamada de variantes que recebe arquivos BAM como entrada e produz VCFs com chamadas conjuntas como saída.
- Três processos (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` e `GATK_JOINTGENOTYPING`) armazenados em arquivos de módulo separados.
- O pipeline escala automaticamente para qualquer número de amostras de entrada usando o paradigma de fluxo de dados do Nextflow.
- Os resultados são publicados em um diretório chamado `results/`.

### Habilidades adquiridas

Através deste curso prático, você aprendeu como:

- Escrever um fluxo de trabalho linear para aplicar chamada de variantes a uma única amostra
- Lidar com arquivos acessórios como arquivos de índice e recursos de genoma de referência de forma apropriada
- Aproveitar o paradigma de fluxo de dados do Nextflow para paralelizar a chamada de variantes por amostra
- Implementar chamada conjunta de múltiplas amostras usando operadores de canal relevantes

Você está agora preparado para começar a aplicar Nextflow a fluxos de trabalho de análise genômica em seu próprio trabalho.

## Próximos passos para aprimorar suas habilidades

Aqui estão nossas principais sugestões sobre o que fazer a seguir:

- Aplicar Nextflow a outros casos de uso de análise científica com [Nextflow for Science](../index.md)
- Começar com nf-core através do [Hello nf-core](../../hello_nf-core/index.md)
- Explorar recursos mais avançados do Nextflow com as [Side Quests](../../side_quests/index.md)

Por fim, recomendamos que você dê uma olhada no [**Seqera Platform**](https://seqera.io/), uma plataforma baseada em nuvem desenvolvida pelos criadores do Nextflow que facilita ainda mais o lançamento e gerenciamento de seus fluxos de trabalho, além de gerenciar seus dados e executar análises de forma interativa em qualquer ambiente.

## Obtendo ajuda

Para recursos de ajuda e suporte da comunidade, consulte a [página de Ajuda](../../help.md).

## Pesquisa de feedback

Antes de prosseguir, por favor reserve um minuto para completar a pesquisa do curso! Seu feedback nos ajuda a melhorar nossos materiais de treinamento para todos.

[Responder à pesquisa :material-arrow-right:](survey.md){ .md-button .md-button--primary }
