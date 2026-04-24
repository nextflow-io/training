---
title: Desarrollo de Plugins
hide:
  - toc
---

# Desarrollo de Plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

El sistema de plugins de Nextflow le permite extender el lenguaje con funciones personalizadas, hooks de monitoreo, backends de ejecución y más.
Los plugins permiten a la comunidad agregar funcionalidades a Nextflow sin modificar su núcleo, lo que los hace ideales para compartir funcionalidades reutilizables entre pipelines.

Durante esta capacitación, aprenderá a usar plugins existentes y, opcionalmente, a crear los suyos propios.

## Audiencia y requisitos previos

La Parte 1 cubre el uso de plugins existentes y es relevante para todos los usuarios de Nextflow.
Las Partes 2-6 cubren la creación de sus propios plugins e involucran código Groovy y herramientas de compilación.
No se requiere experiencia previa en Java o Groovy.

**Requisitos previos**

- Una cuenta de GitHub O una instalación local como se describe [aquí](../../envsetup/02_local).
- Haber completado el curso [Hello Nextflow](../../hello_nextflow/index.md) o equivalente.
- Java 21 o posterior (incluido en el entorno de capacitación; solo necesario para las Partes 2-6).

**Directorio de trabajo:** `side-quests/plugin_development`

## Objetivos de aprendizaje

Al finalizar esta capacitación, será capaz de:

**Uso de plugins (Parte 1):**

- Instalar y configurar plugins existentes en sus workflows
- Importar y usar funciones de plugins

**Desarrollo de plugins (Partes 2-6):**

- Crear un nuevo proyecto de plugin usando el generador de proyectos integrado de Nextflow
- Implementar funciones personalizadas que se puedan invocar desde workflows
- Compilar, probar e instalar su plugin localmente
- Monitorear eventos del workflow (por ejemplo, finalización de tareas, inicio/fin del pipeline) para registros personalizados o notificaciones
- Agregar opciones de configuración para hacer los plugins personalizables
- Distribuir su plugin

## Plan de lecciones

#### Parte 1: Conceptos básicos de plugins

Use plugins existentes en un workflow de Nextflow y configure su comportamiento.

#### Parte 2: Crear un proyecto de plugin

Genere un nuevo proyecto de plugin y examine su estructura.

#### Parte 3: Funciones personalizadas

Implemente funciones personalizadas, compile su plugin y ejecútelo en un workflow.

#### Parte 4: Pruebas

Escriba y ejecute pruebas unitarias usando el framework Spock.

#### Parte 5: Monitoreo de workflows

Responda a eventos como la finalización de tareas para construir un contador de tareas.

#### Parte 6: Configuración y distribución

Lea configuraciones desde `nextflow.config` para hacer su plugin personalizable y luego aprenda cómo compartirlo.

¿Listo para tomar el curso?

[Comenzar a aprender :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
