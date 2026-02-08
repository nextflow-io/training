# Kurszusammenfassung

Herzlichen Glückwunsch zum Abschluss des Hello Nextflow-Trainingskurses! 🎉

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/wnmUAfDL01E?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe die [ganze Playlist auf dem Nextflow YouTube-Kanal](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n).

:green_book: Du kannst das [Video-Transkript](./transcripts/07_next_steps.md) neben dem Video lesen.
///

## Dein Weg

Du hast mit einem einfachen Workflow begonnen, der einen fest codierten Befehl ausführte.
Im Laufe von sechs Teilen hast du diesen einfachen Workflow in eine modulare mehrstufige Pipeline verwandelt, die wichtige Funktionen von Nextflow nutzt, einschließlich channels, Operatoren, integrierter Container-Unterstützung und Konfigurationsoptionen.

### Was du gebaut hast

- Die endgültige Form des Hello-Workflows nimmt als Eingabe eine CSV-Datei mit Textbegrüßungen.
- Die vier Schritte sind als Nextflow-Prozesse (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`) implementiert, die in separaten Moduldateien gespeichert sind.
- Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht.
- Die endgültige Ausgabe der Pipeline ist eine einfache Textdatei mit ASCII-Kunst einer Figur, die die in Großbuchstaben konvertierten Begrüßungen sagt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Schreibt jede Begrüßung in ihre eigene Ausgabedatei (_z.B._ "Hello-output.txt")
2. **`convertToUpper`:** Konvertiert jede Begrüßung in Großbuchstaben (_z.B._ "HELLO")
3. **`collectGreetings`:** Sammelt alle Großbuchstaben-Begrüßungen in einer einzigen Batch-Datei
4. **`cowpy`:** Generiert ASCII-Kunst mit dem `cowpy`-Tool

Die Workflow-Konfiguration unterstützt die flexible und reproduzierbare Bereitstellung von Eingaben und Parametern.

### Erworbene Fähigkeiten

Durch diesen praxisorientierten Kurs hast du gelernt, wie du:

- Nextflow-Kernkomponenten beschreiben und nutzen kannst, um einen einfachen mehrstufigen Workflow zu erstellen
- Weiterführende Konzepte wie Operatoren und channel factories beschreiben kannst
- Einen Nextflow-Workflow lokal starten kannst
- Ausgaben (Ergebnisse) und Log-Dateien, die von Nextflow generiert werden, finden und interpretieren kannst
- Grundlegende Probleme beheben kannst

Du bist jetzt mit dem Grundwissen ausgestattet, um eigene Pipelines in Nextflow zu entwickeln.

## Nächste Schritte zum Aufbau deiner Fähigkeiten

Hier sind unsere Top-3-Vorschläge, was du als Nächstes tun solltest:

- Wende Nextflow auf einen wissenschaftlichen Analyse-Anwendungsfall an mit [Nextflow für die Wissenschaft](../nf4_science/index.md)
- Starte mit nf-core mit [Hello nf-core](../hello_nf-core/index.md)
- Erkunde fortgeschrittenere Nextflow-Funktionen mit den [Side Quests](../side_quests/index.md)

Schließlich empfehlen wir dir, einen Blick auf [**Seqera Platform**](https://seqera.io/) zu werfen, eine Cloud-basierte Plattform, die von den Entwickler\*innen von Nextflow entwickelt wurde und es noch einfacher macht, deine Workflows zu starten und zu verwalten, sowie deine Daten zu verwalten und Analysen interaktiv in jeder Umgebung auszuführen.

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
