# Teil 3: Ressourcen-Profiling und -Optimierung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

DIES IST EIN PLATZHALTER

!!!note "Hinweis"

    Dieses Trainingsmodul wird derzeit überarbeitet.

---

TODO

### 1.1. Führe den Workflow aus, um einen Bericht zur Ressourcennutzung zu erstellen

Damit Nextflow den Bericht automatisch generiert, füge einfach `-with-report <dateiname>.html` zu deinem Befehl hinzu.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

Der Bericht ist eine HTML-Datei, die du herunterladen und in deinem Browser öffnen kannst. Du kannst auch mit der rechten Maustaste im Datei-Explorer auf der linken Seite darauf klicken und `Show preview` auswählen, um ihn in VS Code anzuzeigen.

Nimm dir ein paar Minuten Zeit, um den Bericht durchzusehen und zu prüfen, ob du einige Möglichkeiten zur Anpassung der Ressourcen identifizieren kannst.
Klicke unbedingt auf die Tabs, die die Nutzungsergebnisse als Prozentsatz der zugewiesenen Ressourcen anzeigen.
Es gibt eine [Dokumentation](https://www.nextflow.io/docs/latest/reports.html), die alle verfügbaren Funktionen beschreibt.

<!-- TODO: insert images -->

Eine Beobachtung ist, dass `GATK_JOINTGENOTYPING` sehr CPU-hungrig zu sein scheint, was Sinn macht, da es viele komplexe Berechnungen durchführt.
Wir könnten also versuchen, diese zu erhöhen und sehen, ob sich die Laufzeit dadurch verkürzt.

Allerdings scheinen wir bei den Speicherzuweisungen über das Ziel hinausgeschossen zu sein; alle Prozesse nutzen nur einen Bruchteil dessen, was wir ihnen geben.
Das sollten wir wieder herunterregeln und Ressourcen sparen.

### 1.2. Passe Ressourcenzuweisungen für einen bestimmten Prozess an

Wir können Ressourcenzuweisungen für einen bestimmten Prozess mit dem `withName`-Prozess-Selektor angeben.
Die Syntax sieht folgendermaßen aus, wenn sie allein in einem process-Block steht:

```groovy title="Syntax"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Fügen wir das zum bestehenden process-Block in der `nextflow.config`-Datei hinzu.

```groovy title="nextflow.config" linenums="11"
process {
    // Standardwerte für alle Prozesse
    cpus = 2
    memory = 2.GB
    // Zuweisungen für einen bestimmten Prozess
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

Mit dieser Angabe gelten die Standardeinstellungen für alle Prozesse **außer** dem `GATK_JOINTGENOTYPING`-Prozess, der eine besondere Ausnahme ist und deutlich mehr CPU bekommt.
Hoffentlich sollte das einen Effekt haben.

### 1.3. Führe den Workflow erneut mit der geänderten Konfiguration aus

Lass uns den Workflow erneut mit der geänderten Konfiguration und mit aktiviertem Reporting-Flag ausführen, aber beachte, dass wir dem Bericht einen anderen Namen geben, damit wir sie unterscheiden können.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

Auch hier wirst du wahrscheinlich keinen wesentlichen Unterschied in der Laufzeit bemerken, da dies eine so kleine Arbeitslast ist und die Tools mehr Zeit mit Nebenaufgaben verbringen als mit der eigentlichen 'echten' Arbeit.

Der zweite Bericht zeigt jedoch, dass unsere Ressourcennutzung jetzt ausgewogener ist.

<!-- **TODO: screenshots?** -->

Wie du sehen kannst, ist dieser Ansatz nützlich, wenn deine Prozesse unterschiedliche Ressourcenanforderungen haben. Er ermöglicht es dir, die Ressourcenzuweisungen, die du für jeden Prozess einrichtest, basierend auf tatsächlichen Daten und nicht auf Vermutungen optimal anzupassen.

!!!note "Hinweis"

    Dies ist nur ein kleiner Vorgeschmack darauf, was du tun kannst, um deine Ressourcennutzung zu optimieren.
    Nextflow selbst hat eine wirklich clevere [dynamische Wiederholungslogik](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) eingebaut, um Jobs, die aufgrund von Ressourcenbeschränkungen fehlschlagen, erneut zu versuchen.
    Darüber hinaus bietet die Seqera Platform KI-gestützte Tools zur automatischen Optimierung deiner Ressourcenzuweisungen.

    Wir werden beide Ansätze in einem kommenden Teil dieses Trainingskurses behandeln.

Allerdings kann es je nach verwendetem Computing-Executor und Compute-Infrastruktur einige Einschränkungen geben, was du zuweisen kannst (oder musst). Dein Cluster könnte beispielsweise verlangen, dass du innerhalb bestimmter Grenzen bleibst, die nicht gelten, wenn du woanders ausführst.
