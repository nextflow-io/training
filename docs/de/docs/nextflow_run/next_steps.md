# Kurszusammenfassung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herzlichen Glückwunsch zum Abschluss des Nextflow Run-Trainingskurses! 🎉

<!-- placeholder for video -->

## Dein Fortschritt

Du hast mit einem sehr grundlegenden Workflow begonnen und gelernt, ihn auszuführen, die Ausgaben zu finden und seine Ausführung zu verwalten.
Dann hast du dich durch immer komplexere Versionen dieses Workflows gearbeitet und gelernt, die wesentlichen Konzepte und Mechanismen zu erkennen, die Nextflow-Pipelines antreiben, einschließlich Channels und Operatoren, Code-Modularisierung und Container.
Schließlich hast du gelernt, wie du die Konfiguration einer Pipeline an deine Präferenzen und deine Recheninfrastruktur anpasst.

### Was du gelernt hast

Du bist jetzt in der Lage, die Ausführung der Hello-Pipeline zu verwalten, zu beschreiben, wie sie strukturiert ist, und die wichtigsten Code-Teile zu identifizieren.

- Die endgültige Form des Hello-Workflows nimmt als Eingabe eine CSV-Datei, die Text-Grüße enthält.
- Die vier Schritte sind als Nextflow-Prozesse (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`) implementiert und in separaten Modul-Dateien gespeichert.
- Die Ergebnisse werden in ein Verzeichnis namens `results/` veröffentlicht.
- Die finale Ausgabe der Pipeline ist eine einfache Textdatei, die ASCII-Kunst eines Charakters enthält, der die großgeschriebenen Grüße sagt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Schreibt jeden Gruß in seine eigene Ausgabedatei (_z.B._ "Hello-output.txt")
2. **`convertToUpper`:** Konvertiert jeden Gruß in Großbuchstaben (_z.B._ "HELLO")
3. **`collectGreetings`:** Sammelt alle großgeschriebenen Grüße in einer einzelnen Batch-Datei
4. **`cowpy`:** Generiert ASCII-Kunst mit dem `cowpy`-Tool

Die Workflow-Konfiguration unterstützt die flexible und reproduzierbare Bereitstellung von Eingaben und Parametern.

### Erworbene Fähigkeiten

Durch diesen praxisorientierten Kurs hast du gelernt, wie man:

- Einen Nextflow-Workflow lokal startet
- Ausgaben (Ergebnisse) und Protokolldateien, die von Nextflow generiert werden, findet und interpretiert
- Die Kernkomponenten von Nextflow erkennt, die einen einfachen mehrstufigen Workflow ausmachen
- Fortgeschrittene Konzepte wie Operatoren und Channel-Factories beschreibt
- Pipelines für verschiedene Rechenumgebungen konfiguriert

Du bist jetzt mit dem grundlegenden Wissen ausgestattet, um bestehende Nextflow-Pipelines in deine eigene Arbeit zu integrieren.

## Nächste Schritte zum Aufbau deiner Fähigkeiten

Hier sind unsere Top-Vorschläge, was du als Nächstes tun kannst:

- Führe Nextflow nicht nur aus, schreibe es! Lerne die Nextflow-Entwicklung mit [Hello Nextflow](../hello_nextflow/index.md)
- Wende Nextflow auf einen wissenschaftlichen Analyse-Anwendungsfall an mit [Nextflow for Science](../nf4_science/index.md)
- Starte mit nf-core mit [Hello nf-core](../hello_nf-core/index.md)
- Lerne Troubleshooting-Techniken mit der [Debugging Side Quest](../side_quests/debugging.md)

Schließlich empfehlen wir dir einen Blick auf [**Seqera Platform**](https://seqera.io/), eine Cloud-basierte Plattform, die von den Erstellern von Nextflow entwickelt wurde und es noch einfacher macht, deine Workflows zu starten und zu verwalten, deine Daten zu verwalten und Analysen interaktiv in jeder Umgebung auszuführen.

## Hilfe bekommen

Für Hilfe-Ressourcen und Community-Support siehe die [Hilfe-Seite](../help.md).

## Feedback-Umfrage

Bevor du weitergehst, nimm dir bitte eine Minute Zeit, um die Kurs-Umfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
