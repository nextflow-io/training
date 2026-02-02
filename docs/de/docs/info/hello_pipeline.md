---
title: Die Hello-Pipeline
description: Zusammenfassung dessen, was die Hello-Pipeline tut und wie sie strukturiert ist.
hide:
  - toc
  - footer
---

# Die Hello-Pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Die meisten unserer Trainingskurse verwenden eine einfache domänenunabhängige pipeline, um Nextflow-Konzepte und -Mechanismen zu demonstrieren.
Der Hello Nextflow-Kurs zeigt, wie diese pipeline Schritt für Schritt entwickelt wird, wobei jede Design- und Implementierungsentscheidung erklärt wird.
Andere Trainings verwenden diese pipeline oder Teile davon als Ausgangspunkt.

Diese Seite fasst den Zustand der pipeline nach Abschluss des Hello Nextflow-Kurses zusammen.

### Kurzbeschreibung

Der Hello-workflow nimmt eine CSV-Datei mit Grüßen, schreibt sie in separate Dateien, konvertiert jede in Großbuchstaben, sammelt sie wieder zusammen und gibt eine einzelne Textdatei aus, die ein ASCII-Bild einer lustigen Figur enthält, die die Grüße sagt.

### Workflow-Schritte (Prozesse)

Die vier Schritte sind als Nextflow-processes implementiert (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`), die in separaten Modul-Dateien gespeichert sind.

1. **`sayHello`:** Schreibt jeden Gruß in eine eigene Ausgabedatei (z.B. "Hello-output.txt")
2. **`convertToUpper`:** Konvertiert jeden Gruß in Großbuchstaben (z.B. "HELLO")
3. **`collectGreetings`:** Sammelt alle Großbuchstaben-Grüße in einer einzelnen Batch-Datei
4. **`cowpy`:** Generiert ASCII-Kunst mit dem `cowpy`-Tool

### Diagramm

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Ergebnisse

Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht, und die endgültige Ausgabe der pipeline (bei Ausführung mit Standardparametern) ist eine Klartextdatei mit ASCII-Kunst eines Truthahns, der die Grüße in Großbuchstaben sagt.

```txt title="results/cowpy-COLLECTED-test-batch-output.txt"
  _________
/ BONJOUR \
| HELLO   |
\ HOLà    /
---------
  \                                  ,+*^^*+___+++_
  \                           ,*^^^^              )
    \                       _+*                     ^**+_
    \                    +^       _ _++*+_+++_,         )
              _+^^*+_    (     ,+*^ ^          \+_        )
            {       )  (    ,(    ,_+--+--,      ^)      ^\
            { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
          {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
          ( /  (    (        ,___    ^*+_+* )   <    <      \
          U _/     )    *--<  ) ^\-----++__)   )    )       )
            (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
          (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
        (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
          *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
          \             \_)^)_)) ))^^^^^^^^^^))^^^^)
          (_             ^\__^^^^^^^^^^^^))^^^^^^^)
            ^\___            ^\__^^^^^^))^^^^^^^^)\\
                  ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                    ___) >____) >___   ^\_\_\_\_\_\_\)
                    ^^^//\\_^^//\\_^       ^(\_\_\_\)
                      ^^^ ^^ ^^^ ^
```

Je nach Kurs, in dem die pipeline vorgestellt wird, können kleine Abweichungen in den Details auftreten.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
