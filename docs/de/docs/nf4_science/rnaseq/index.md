# Nextflow für RNAseq

Dieser Trainingskurs richtet sich an Forscher\*innen im Bereich Transkriptomik und verwandten Feldern, die daran interessiert sind, Datenanalyse-Pipelines zu entwickeln oder anzupassen.
Er baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und zeigt, wie du Nextflow im spezifischen Kontext der Bulk-RNAseq-Analyse einsetzen kannst.

Konkret demonstriert dieser Kurs, wie du eine einfache Bulk-RNAseq-Verarbeitungspipeline implementierst, um Adaptersequenzen zu trimmen, die Reads gegen eine Genomreferenz zu alignieren und Quality Control (QC) in mehreren Schritten durchzuführen.

Los geht's! Klicke auf den Button „Open in GitHub Codespaces" unten, um die Trainingsumgebung zu starten (am besten in einem separaten Tab), und lies dann weiter, während sie lädt.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Lernziele

Durch die Bearbeitung dieses Kurses lernst du, wie du grundlegende Nextflow-Konzepte und -Tools auf einen typischen RNAseq-Anwendungsfall anwendest.

Am Ende dieses Workshops wirst du in der Lage sein:

- Einen linearen Workflow zu schreiben, um grundlegende RNAseq-Verarbeitungs- und QC-Methoden anzuwenden
- Domänenspezifische Dateien wie FASTQ und Genomreferenz-Ressourcen angemessen zu handhaben
- Single-End- und Paired-End-Sequenzierungsdaten zu verarbeiten
- Nextflows Dataflow-Paradigma zu nutzen, um die RNAseq-Verarbeitung pro Probe zu parallelisieren
- QC-Reports über mehrere Schritte und Proben hinweg mithilfe relevanter Kanal-Operatoren zu aggregieren

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Voraussetzungen

Der Kurs setzt minimale Vertrautheit mit Folgendem voraus:

- Tools und Dateiformate, die in diesem wissenschaftlichen Bereich häufig verwendet werden
- Erfahrung mit der Kommandozeile
- Grundlegende Nextflow-Konzepte und -Tools, die im [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training behandelt werden

Für technische Anforderungen und Umgebungseinrichtung siehe den Mini-Kurs [Environment Setup](../../envsetup/).
