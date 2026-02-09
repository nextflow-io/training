# Nextflow für Genomik

**Ein praxisorientierter Kurs zur Anwendung von Nextflow auf einen realen Genomik-Anwendungsfall: Variantenerkennung mit GATK.**

Dieser Kurs baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und zeigt, wie du Nextflow im spezifischen Kontext der Genomik einsetzen kannst.
Du wirst eine Pipeline zur Variantenerkennung mit [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) implementieren, einem weit verbreiteten Softwarepaket zur Analyse von Hochdurchsatz-Sequenzierungsdaten.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die so strukturiert sind, dass Informationen schrittweise eingeführt werden.

Du beginnst damit, die Tools zur Variantenerkennung manuell im Terminal auszuführen, um die Methodik zu verstehen. Anschließend baust du schrittweise eine Nextflow-Pipeline auf, die die Analyse automatisiert und skaliert.

### Lektionsplan

Wir haben den Kurs in drei Teile unterteilt, die sich jeweils auf spezifische Aspekte der Anwendung von Nextflow auf einen Genomik-Anwendungsfall konzentrieren.

| Kurskapitel                                                              | Zusammenfassung                                                                                                      | Geschätzte Dauer |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Methodenübersicht](./01_method.md)                              | Die Methodik zur Variantenerkennung verstehen und die Tools manuell ausführen                                        | 30 Min.          |
| [Teil 2: Variantenerkennung pro Probe](./02_per_sample_variant_calling.md) | Eine Pipeline erstellen, die BAM-Dateien indiziert und Varianten erkennt, dann auf mehrere Proben skalieren         | 60 Min.          |
| [Teil 3: Joint Calling auf einer Kohorte](./03_joint_calling.md)         | Multi-Sample Joint Genotyping hinzufügen, indem Kanal-Operatoren verwendet werden, um Ausgaben pro Probe zu aggregieren | 45 Min.          |

Am Ende dieses Kurses kannst du grundlegende Nextflow-Konzepte und -Tools auf einen typischen Genomik-Anwendungsfall anwenden.

Bereit, den Kurs zu beginnen?

[Los geht's :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
