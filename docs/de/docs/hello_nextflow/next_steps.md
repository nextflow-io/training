# Kurszusammenfassung

Herzlichen Glückwunsch zum Abschluss des Hello Nextflow Trainingskurses! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Sieh dir die [gesamte Playlist auf dem Nextflow YouTube-Kanal](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) an.

:green_book: Du kannst das [Video-Transkript](./transcripts/07_next_steps.md) parallel zum Video lesen.
///

## Deine Lernreise

Du hast mit einem sehr einfachen Workflow begonnen, der einen fest codierten Befehl ausgeführt hat.
Im Verlauf von sechs Teilen hast du diesen einfachen Workflow in eine modulare mehrstufige Pipeline verwandelt, die wichtige Funktionen von Nextflow nutzt, darunter Kanäle, Operatoren, integrierte Container-Unterstützung und Konfigurationsoptionen.

### Was du entwickelt hast

- Die finale Version des Hello-Workflows nimmt als Eingabe eine CSV-Datei mit Textbegrüßungen entgegen.
- Die vier Schritte sind als Nextflow-Prozesse (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`) implementiert, die in separaten Moduldateien gespeichert sind.
- Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht.
- Die finale Ausgabe der Pipeline ist eine Textdatei mit ASCII-Art eines Charakters, der die großgeschriebenen Begrüßungen sagt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Schreibt jede Begrüßung in eine eigene Ausgabedatei (_z.B._ "Hello-output.txt")
2. **`convertToUpper`:** Wandelt jede Begrüßung in Großbuchstaben um (_z.B._ "HELLO")
3. **`collectGreetings`:** Sammelt alle großgeschriebenen Begrüßungen in einer einzigen Batch-Datei
4. **`cowpy`:** Erzeugt ASCII-Art mit dem `cowpy`-Tool

Die Workflow-Konfiguration unterstützt die flexible und reproduzierbare Bereitstellung von Eingaben und Parametern.

### Erworbene Fähigkeiten

In diesem praxisorientierten Kurs hast du gelernt, wie du:

- Zentrale Nextflow-Komponenten beschreibst und nutzt, um einen einfachen mehrstufigen Workflow zu erstellen
- Weiterführende Konzepte wie Operatoren und Channel Factories beschreibst
- Einen Nextflow-Workflow lokal ausführst
- Ausgaben (Ergebnisse) und Log-Dateien, die von Nextflow generiert werden, findest und interpretierst
- Grundlegende Probleme behebst

Du bist jetzt mit dem grundlegenden Wissen ausgestattet, um mit der Entwicklung deiner eigenen Pipelines in Nextflow zu beginnen.

## Nächste Schritte zum Ausbau deiner Fähigkeiten

Hier sind unsere Top-3-Empfehlungen, was du als Nächstes tun kannst:

- Wende Nextflow auf einen wissenschaftlichen Analyse-Anwendungsfall an mit [Nextflow for Science](../nf4_science/index.md)
- Steig ein in nf-core mit [Hello nf-core](../hello_nf-core/index.md)
- Erkunde fortgeschrittenere Nextflow-Funktionen mit den [Side Quests](../side_quests/index.md)

Abschließend empfehlen wir dir, einen Blick auf [**Seqera Platform**](https://seqera.io/) zu werfen, eine cloudbasierte Plattform, die von den Entwickler\*innen von Nextflow entwickelt wurde und es noch einfacher macht, deine Workflows zu starten und zu verwalten sowie deine Daten zu managen und Analysen interaktiv in jeder Umgebung auszuführen.

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
