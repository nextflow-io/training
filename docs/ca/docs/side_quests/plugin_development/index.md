---
title: Desenvolupament de Plugins
hide:
  - toc
---

# Desenvolupament de Plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

El sistema de plugins de Nextflow us permet estendre el llenguatge amb funcions personalitzades, hooks de monitoratge, backends d'execució i molt més.
Els plugins permeten a la comunitat afegir funcionalitats a Nextflow sense modificar el seu nucli, cosa que els fa ideals per compartir funcionalitats reutilitzables entre pipelines.

Durant aquesta formació, aprendreu a utilitzar plugins existents i, opcionalment, a crear-ne de propis.

## Audiència i prerequisits

La Part 1 cobreix l'ús de plugins existents i és rellevant per a tots els usuaris de Nextflow.
Les Parts 2-6 cobreixen la construcció dels vostres propis plugins i impliquen codi Groovy i eines de construcció.
No cal experiència prèvia amb Java ni Groovy.

**Prerequisits**

- Un compte de GitHub O una instal·lació local tal com es descriu [aquí](../../envsetup/02_local).
- Haver completat el curs [Hello Nextflow](../../hello_nextflow/index.md) o equivalent.
- Java 21 o posterior (inclòs a l'entorn de formació; només necessari per a les Parts 2-6).

**Directori de treball:** `side-quests/plugin_development`

## Objectius d'aprenentatge

Al final d'aquesta formació, sereu capaços de:

**Ús de plugins (Part 1):**

- Instal·lar i configurar plugins existents als vostres workflows
- Importar i utilitzar funcions de plugins

**Desenvolupament de plugins (Parts 2-6):**

- Crear un nou projecte de plugin utilitzant el generador de projectes integrat de Nextflow
- Implementar funcions personalitzades invocables des de workflows
- Construir, provar i instal·lar el vostre plugin localment
- Monitorar esdeveniments del workflow (p. ex., finalització de tasques, inici/fi del pipeline) per a registres o notificacions personalitzades
- Afegir opcions de configuració per fer els plugins personalitzables
- Distribuir el vostre plugin

## Pla de lliçons

#### Part 1: Conceptes bàsics de plugins

Utilitzeu plugins existents en un workflow de Nextflow i configureu el seu comportament.

#### Part 2: Creació d'un projecte de plugin

Genereu un nou projecte de plugin i examineu-ne l'estructura.

#### Part 3: Funcions personalitzades

Implementeu funcions personalitzades, construïu el vostre plugin i executeu-lo en un workflow.

#### Part 4: Proves

Escriviu i executeu proves unitàries utilitzant el framework Spock.

#### Part 5: Monitoratge del workflow

Responeu a esdeveniments com la finalització de tasques per construir un comptador de tasques.

#### Part 6: Configuració i distribució

Llegiu paràmetres de `nextflow.config` per fer el vostre plugin personalitzable i, a continuació, apreneu com compartir-lo.

Esteu preparats per fer el curs?

[Comenceu a aprendre :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
