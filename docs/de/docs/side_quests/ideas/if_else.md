# Teil 2: If - Else

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Lass die Kuh berühmte Wissenschaftler\*innen zitieren

Dieser Abschnitt enthält einige zusätzliche Übungen, um das bisher Gelernte zu üben.
Diese Übungen sind _nicht erforderlich_, um spätere Teile des Trainings zu verstehen, bieten aber eine unterhaltsame Möglichkeit, dein Wissen zu festigen, indem du herausfindest, wie du die Kuh berühmte Wissenschaftler\*innen zitieren lassen kannst.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 1.1. Modifiziere das `hello-containers.nf`-Skript, um einen getQuote-Prozess zu verwenden

Wir haben eine Liste von Computer- und Biologie-Pionier\*innen in der Datei `containers/data/pioneers.csv`.
Um diese Übung zu lösen, musst du im Wesentlichen:

- Den Standard-`params.input_file` so ändern, dass er auf die Datei `pioneers.csv` zeigt.
- Einen `getQuote`-Prozess erstellen, der den `quote`-Container verwendet, um für jede Eingabe ein Zitat abzurufen.
- Die Ausgabe des `getQuote`-Prozesses mit dem `cowsay`-Prozess verbinden, um das Zitat anzuzeigen.

Für das `quote`-Container-Image kannst du entweder das verwenden, das du selbst in der vorherigen zusätzlichen Übung erstellt hast, oder das von Seqera Containers.

!!! Hint "Hinweis"

    Eine gute Wahl für den `script`-Block deines getQuote-Prozesses könnte sein:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Eine Lösung zu dieser Übung findest du in `containers/solutions/hello-containers-4.1.nf`.

### 1.2. Modifiziere deine Nextflow-Pipeline, um sie in den Modi `quote` und `sayHello` ausführen zu können.

Füge deiner Pipeline eine Verzweigungslogik hinzu, damit sie Eingaben sowohl für `quote` als auch für `sayHello` akzeptieren kann.
Hier ist ein Beispiel, wie du eine `if`-Anweisung in einem Nextflow-Workflow verwendest:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! Hint "Hinweis"

    Du kannst `new_ch = processName.out` verwenden, um dem Ausgabekanal eines Prozesses einen Namen zuzuweisen.

Eine Lösung zu dieser Übung findest du in `containers/solutions/hello-containers-4.2.nf`.

### Fazit

Du weißt jetzt, wie du Container in Nextflow verwendest, um Prozesse auszuführen, und wie du Verzweigungslogik in deine Pipelines einbaust!

### Wie geht es weiter?

Feiere deinen Erfolg, mach eine Pause und trink etwas Wasser!

Wenn du bereit bist, gehe zu Teil 3 dieser Trainingsreihe über, um zu lernen, wie du das bisher Gelernte auf einen realistischeren Datenanalyse-Anwendungsfall anwendest.
