---
title: Executar Nextflow para Imagens
hide:
  - toc
---

# Executar Nextflow para Imagens

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Este curso de treinamento é destinado a pesquisadores em imagens e biologia espacial que estão interessados em executar e personalizar pipelines de análise de dados.
Ele ensina conceitos fundamentais do Nextflow relacionados à execução, organização e configuração de workflows usando o [nf-core/molkart](https://nf-co.re/molkart), um pipeline para processamento de dados de transcriptômica espacial Molecular Cartography.
As habilidades que você aprenderá aqui são transferíveis para qualquer pipeline Nextflow ou nf-core.

Vamos começar! Clique no botão "Open in GitHub Codespaces" abaixo para iniciar o ambiente de treinamento (de preferência em uma aba separada), depois continue lendo enquanto ele carrega.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objetivos de aprendizagem

Ao trabalhar neste curso, você aprenderá como aplicar conceitos fundamentais e ferramentas do Nextflow para executar pipelines de análise de imagens.

Ao final deste workshop, você será capaz de:

- Iniciar um workflow Nextflow localmente e monitorar a execução
- Encontrar e interpretar saídas (resultados) e arquivos de log gerados pelo Nextflow
- Executar um pipeline nf-core com dados de teste e entradas personalizadas
- Configurar a execução do pipeline usando perfis e arquivos de parâmetros
- Gerenciar entradas usando samplesheets e parâmetros de linha de comando

## Público e pré-requisitos

Este curso pressupõe alguma familiaridade mínima com o seguinte:

- Experiência com a linha de comando
- Familiaridade básica com formatos de arquivos de imagem (imagens TIFF, dados tabulares)

Para requisitos técnicos e configuração do ambiente, consulte o mini-curso de [Configuração do Ambiente](../../envsetup/).
