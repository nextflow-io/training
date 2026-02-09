# Kurszusammenfassung

Herzlichen Glückwunsch zum Abschluss des Nextflow Run Trainingskurses! 🎉

<!-- placeholder for video -->

## Deine Lernreise

Du hast mit einem sehr einfachen Workflow begonnen und gelernt, ihn auszuführen, die Ausgaben zu finden und seine Ausführung zu verwalten.
Dann hast du dich durch zunehmend komplexere Versionen dieses Workflows gearbeitet und gelernt, die wesentlichen Konzepte und Mechanismen zu erkennen, die Nextflow-Pipelines antreiben, einschließlich Kanäle und Operatoren, Code-Modularisierung und Container.
Schließlich hast du gelernt, wie du die Konfiguration einer Pipeline an deine Präferenzen und deine Recheninfrastruktur anpasst.

### Was du gelernt hast

Du kannst jetzt die Ausführung der Hello-Pipeline verwalten, beschreiben, wie sie strukturiert ist, und die wichtigsten beteiligten Code-Teile identifizieren.

- Die finale Version des Hello-Workflows nimmt als Eingabe eine CSV-Datei mit Textbegrüßungen entgegen.
- Die vier Schritte sind als Nextflow-Prozesse (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`) implementiert, die in separaten Moduldateien gespeichert sind.
- Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht.
- Die finale Ausgabe der Pipeline ist eine Textdatei mit ASCII-Art eines Charakters, der die großgeschriebenen Begrüßungen sagt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Schreibt jede Begrüßung in eine eigene Ausgabedatei (_z.B._ "Hello-output.txt")
2. **`convertToUpper`:** Konvertiert jede Begrüßung in Großbuchstaben (_z.B._ "HELLO")
3. **`collectGreetings`:** Sammelt alle großgeschriebenen Begrüßungen in einer einzigen Batch-Datei
4. **`cowpy`:** Generiert ASCII-Art mit dem `cowpy`-Tool

Die Workflow-Konfiguration unterstützt die flexible und reproduzierbare Bereitstellung von Eingaben und Parametern.

### Erworbene Fähigkeiten

Durch diesen praxisorientierten Kurs hast du gelernt, wie du:

- Einen Nextflow-Workflow lokal startest
- Ausgaben (Ergebnisse) und Log-Dateien, die von Nextflow generiert werden, findest und interpretierst
- Die Kern-Komponenten von Nextflow erkennst, die einen einfachen mehrstufigen Workflow bilden
- Weiterführende Konzepte wie Operatoren und Channel Factories beschreibst
- Pipelines für verschiedene Rechenumgebungen konfigurierst

Du bist jetzt mit dem grundlegenden Wissen ausgestattet, um bestehende Nextflow-Pipelines in deine eigene Arbeit zu integrieren.

## Nächste Schritte zum Ausbau deiner Fähigkeiten

Hier sind unsere Top-Empfehlungen, was du als Nächstes tun kannst:

- Führe nicht nur Nextflow aus, schreibe es! Werde Nextflow-Entwickler\*in mit [Hello Nextflow](../hello_nextflow/index.md)
- Wende Nextflow auf einen wissenschaftlichen Analyse-Anwendungsfall an mit [Nextflow for Science](../nf4_science/index.md)
- Starte mit nf-core durch [Hello nf-core](../hello_nf-core/index.md)
- Lerne Troubleshooting-Techniken mit der [Debugging Side Quest](../side_quests/debugging.md)

Schließlich empfehlen wir dir, einen Blick auf [**Seqera Platform**](https://seqera.io/) zu werfen, eine cloudbasierte Plattform, die von den Entwickler\*innen von Nextflow entwickelt wurde und es noch einfacher macht, deine Workflows zu starten und zu verwalten sowie deine Daten zu verwalten und Analysen interaktiv in jeder Umgebung auszuführen.

## Hilfe erhalten

Für Hilfsressourcen und Community-Support siehe die [Hilfe-Seite](../help.md).

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
