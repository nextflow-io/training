---
title: Die Hello-Pipeline
description: Zusammenfassung dessen, was die Hello-Pipeline macht und wie sie strukturiert ist.
hide:
  - toc
  - footer
---

# Die Hello-Pipeline

Die meisten unserer Trainingskurse verwenden eine einfache, domänenunabhängige Pipeline, um Nextflow-Konzepte und -Mechanismen zu demonstrieren.
Der Kurs Hello Nextflow zeigt Schritt für Schritt, wie diese Pipeline entwickelt wird, und erklärt dabei jede Design- und Implementierungsentscheidung.
Andere Trainings verwenden diese Pipeline oder Teile davon als Ausgangspunkt.

Diese Seite fasst den Zustand der Pipeline zusammen, wie er am Ende des Kurses Hello Nextflow vorliegt.

### Zusammenfassung

Der Hello-Workflow nimmt eine CSV-Datei mit Begrüßungen entgegen, schreibt sie in separate Dateien, konvertiert jede in Großbuchstaben, sammelt sie wieder zusammen und gibt eine einzelne Textdatei aus, die ein ASCII-Bild einer lustigen Figur enthält, die die Begrüßungen sagt.

### Workflow-Schritte (Prozesse)

Die vier Schritte sind als Nextflow-Prozesse (`sayHello`, `convertToUpper`, `collectGreetings` und `cowpy`) implementiert und in separaten Moduldateien gespeichert.

1. **`sayHello`:** Schreibt jede Begrüßung in eine eigene Ausgabedatei (z. B. "Hello-output.txt")
2. **`convertToUpper`:** Konvertiert jede Begrüßung in Großbuchstaben (z. B. "HELLO")
3. **`collectGreetings`:** Sammelt alle Begrüßungen in Großbuchstaben in einer einzigen Batch-Datei
4. **`cowpy`:** Generiert ASCII-Art mit dem Tool `cowpy`

### Diagramm

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

### Ergebnisse

Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht, und die finale Ausgabe der Pipeline (bei Ausführung mit Standardparametern) ist eine Klartextdatei mit ASCII-Art eines Truthahns, der die Begrüßungen in Großbuchstaben sagt.

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

Je nach Kurs, in dem die Pipeline vorgestellt wird, können einige Details variieren.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
