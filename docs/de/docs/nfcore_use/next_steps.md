# Kursübersicht

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herzlichen Glückwunsch zum Abschluss des Kurses „Use nf-core"! 🎉

<!-- placeholder for video -->

## Deine Reise

Du hast damit begonnen, die `nf-core/demo`-Pipeline zu finden und abzurufen, und dann gelernt, sie mit ihrem Test-Profil auszuführen und ihre Ausgaben zu untersuchen.
Anschließend hast du die Ausführung über Pipeline-Parameter und Konfigurationsdateien angepasst und gesehen, wie nf-core-Pipelines Parameter und Eingabedaten validieren.
Schließlich hast du dieselben Fähigkeiten auf `nf-core/rnaseq` angewendet – eine Pipeline im Produktionsmaßstab – und gelernt, wie du ihre Standard-Ressourcenzuweisungen an die dir verfügbare Hardware anpassen kannst.

### Was du gelernt hast

Du kannst jetzt nf-core-Pipelines finden, abrufen, ausführen und konfigurieren.

- nf-core-Pipelines werden mit `nextflow pull` abgerufen und folgen einer standardisierten Code-Organisation.
- Jede nf-core-Pipeline enthält ein `test`-Profil zur schnellen Validierung mit einem kleinen Datensatz.
- Pipeline-Parameter (gesetzt über `--param_name` oder `-params-file`) und Konfiguration (gesetzt über `-c`) dienen unterschiedlichen Zwecken: Eingaben und Analyseoptionen einerseits, Ausführungsdetails wie Ressourcenzuweisung andererseits.
- nf-core-Pipelines validieren Parameter und Eingabedateien automatisch und erkennen Fehler, bevor irgendeine Arbeit erledigt wird.
- Ressourcen-Standardwerte werden über Labels (`process_low`, `process_medium`, `process_high`) vergeben, die in `conf/base.config` definiert sind und mit einer eigenen Konfigurationsdatei überschrieben werden können.

### Erworbene Fähigkeiten

In diesem praxisorientierten Kurs hast du gelernt, wie du:

- Eine nf-core-Pipeline auf der nf-co.re-Website findest und ihren Quellcode abrufst
- Eine Pipeline mit ihrem integrierten Test-Profil ausführst und ihre Ausgaben untersuchst
- Hilfe bekommst, Parameter setzt und Parameter- sowie Eingabevalidierung verstehst
- Ressourcenzuweisung und Tool-Argumente über Konfigurationsdateien anpasst
- Eine Pipeline im Produktionsmaßstab abrufst und ausführst und ihre Standard-Ressourcen-Labels überschreibst

Du verfügst jetzt über das grundlegende Wissen, um nf-core-Pipelines für deine eigenen Analysen einzusetzen.

## Nächste Schritte zum Ausbau deiner Fähigkeiten

Hier sind unsere wichtigsten Empfehlungen für das weitere Vorgehen:

- Starte und überwache diese Pipelines in großem Maßstab mit [Scale with Seqera](../seqera_scale/index.md)
- Führe nf-core-Pipelines nicht nur aus, sondern entwickle sie! Lerne nf-core Best Practices mit [Build with nf-core](../nfcore_build/index.md)
- Neu bei Nextflow? Starte mit [Nextflow Run](../nextflow_run/index.md)
- Wende Nextflow auf einen wissenschaftlichen Anwendungsfall an mit [Nextflow for Science](../nf4_science/index.md)
- Entdecke fortgeschrittenere Nextflow-Funktionen mit den [Side Quests](../side_quests/index.md)

## Hilfe erhalten

Hilferessourcen und Community-Support findest du auf der [Hilfe-Seite](../help.md).

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
