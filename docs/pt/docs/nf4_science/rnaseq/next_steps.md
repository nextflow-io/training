# Resumo do Curso

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Parabéns por concluir o curso de treinamento Nextflow para RNAseq!

## Sua jornada

Você começou executando ferramentas de processamento de RNAseq manualmente no terminal para entender a metodologia.
Depois, você construiu um pipeline Nextflow para uma única amostra para automatizar o processo, escalou-o para lidar com múltiplas amostras em paralelo e o estendeu para lidar com dados paired-end e agregar relatórios de QC entre amostras.

### O que você construiu

- Um pipeline de processamento de RNAseq que recebe arquivos FASTQ como entrada e produz reads aparados, alinhamentos e relatórios de QC agregados como saída.
- Processos para aparar (Trim Galore), alinhar (HISAT2), controle de qualidade (FastQC) e agregação de relatórios (MultiQC) armazenados em arquivos de módulo separados.
- O pipeline paraleliza automaticamente o processamento de amostras de entrada usando o paradigma de fluxo de dados do Nextflow.
- O pipeline final lida com dados de sequenciamento paired-end.

### Habilidades adquiridas

Através deste curso prático, você aprendeu como:

- Escrever um fluxo de trabalho linear para aplicar métodos básicos de processamento e QC de RNAseq
- Lidar adequadamente com arquivos específicos do domínio, como FASTQ e recursos de genoma de referência
- Lidar com dados de sequenciamento single-end e paired-end
- Aproveitar o paradigma de fluxo de dados do Nextflow para paralelizar o processamento de RNAseq por amostra
- Agregar relatórios de QC em múltiplas etapas e amostras usando operadores de canal relevantes

Você está agora preparado para começar a aplicar Nextflow a fluxos de trabalho de análise de RNAseq em seu próprio trabalho.

## Próximos passos para aprimorar suas habilidades

Aqui estão nossas principais sugestões sobre o que fazer a seguir:

- Aplicar Nextflow a outros casos de uso de análise científica com [Nextflow para Ciência](../index.md)
- Começar com nf-core com [Hello nf-core](../../hello_nf-core/index.md)
- Explorar recursos mais avançados do Nextflow com as [Missões Secundárias](../../side_quests/index.md)

Por fim, recomendamos que você dê uma olhada no [**Seqera Platform**](https://seqera.io/), uma plataforma baseada em nuvem desenvolvida pelos criadores do Nextflow que facilita ainda mais o lançamento e gerenciamento de seus fluxos de trabalho, bem como o gerenciamento de seus dados e a execução de análises interativamente em qualquer ambiente.

## Obtendo ajuda

Para recursos de ajuda e suporte da comunidade, consulte a [página de Ajuda](../../help.md).

## Pesquisa de feedback

Antes de prosseguir, por favor, reserve um minuto para completar a pesquisa do curso! Seu feedback nos ajuda a melhorar nossos materiais de treinamento para todos.

[Responder à pesquisa :material-arrow-right:](survey.md){ .md-button .md-button--primary }
