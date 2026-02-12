---
title: Imaging
hide:
  - toc
---

# Nextflow per a Imatge

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Aquest curs de formació està destinat a investigadors en imatge i biologia espacial que estiguin interessats en executar i personalitzar pipelines d'anàlisi de dades.
Ensenya conceptes fonamentals de Nextflow relacionats amb l'execució, organització i configuració de workflows utilitzant [nf-core/molkart](https://nf-co.re/molkart), un pipeline per processar dades de transcriptòmica espacial de Cartografia Molecular.
Les habilitats que aprendràs aquí són transferibles a qualsevol pipeline de Nextflow o nf-core.

Comencem! Fes clic al botó "Open in GitHub Codespaces" a continuació per llançar l'entorn de formació (preferiblement en una pestanya separada), i després continua llegint mentre es carrega.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Objectius d'aprenentatge

En treballar a través d'aquest curs, aprendràs a aplicar conceptes i eines fonamentals de Nextflow per executar pipelines d'anàlisi d'imatge.

Al final d'aquest taller seràs capaç de:

- Llançar un workflow de Nextflow localment i monitoritzar l'execució
- Trobar i interpretar sortides (resultats) i fitxers de registre generats per Nextflow
- Executar un pipeline nf-core amb dades de prova i entrades personalitzades
- Configurar l'execució del pipeline utilitzant perfils i fitxers de paràmetres
- Gestionar entrades utilitzant fulls de mostres i paràmetres de línia de comandes

## Audiència i prerequisits

Aquest curs assumeix una familiaritat mínima amb el següent:

- Experiència amb la línia de comandes
- Familiaritat bàsica amb formats de fitxer d'imatge (imatges TIFF, dades tabulars)

Per als requisits tècnics i la configuració de l'entorn, consulteu el mini-curs [Configuració de l'entorn](../../envsetup/).
