---
title: Nextflow-Versionen
description: Die Entwicklung der Nextflow-Syntaxversionen verstehen und verwalten
hide:
  - toc
  - footer
---

## Aktuell unterstützte Nextflow-Syntaxversion & Anforderungen

Ab Version 3.0 des Trainingsportals basieren alle unsere Trainingskurse auf der Nextflow-Version 25.10.2, sofern auf der Kursübersichtsseite nicht anders angegeben (mit Ausnahme veralteter oder archivierter Materialien, die möglicherweise keinen Versionshinweis enthalten).

Da die Kurse jetzt typisierte Eingaben auf Workflow-Ebene sowie Workflow-Ausgabedirektiven verwenden, benötigen sie den V2-Syntaxparser.
Wenn du die Umgebung verwendest, die wir über [Github Codespaces](../envsetup/01_setup.md) oder [lokale Devcontainer](../envsetup/03_devcontainer.md) bereitstellen, musst du nichts tun, sofern in den Kursanweisungen nicht ausdrücklich etwas anderes angegeben ist.
Wenn du jedoch planst, die Trainings in deiner eigenen Umgebung zu absolvieren ([Manuelle Installation](../envsetup/02_local.md)), musst du sicherstellen, dass du Nextflow Version 25.10.2 oder neuer mit aktiviertem v2-Syntaxparser verwendest.

## Ältere Versionen der Trainingsmaterialien

Unsere Trainingsmaterialien werden seit Februar 2025 versioniert.

Du kannst auf ältere Versionen der Trainingsmaterialien zugreifen, die mit Nextflow-Versionen **vor 25.10.2** funktionieren, über das Dropdown-Menü oben auf jeder Seite, das die nummerierte Version der Trainingsmaterialien anzeigt.
Wenn du eine ältere Version der Trainingsmaterialien auswählst, geben Links zur Trainingsumgebung automatisch die entsprechende Version der Umgebung an.

## Weitere Informationen zu Nextflow-Syntaxversionen

Nextflow hat zwei unterschiedliche Versionierungskonzepte, die manchmal verwechselt werden: **DSL-Versionen** und **Syntaxparser-Versionen**.

**DSL1 vs DSL2** bezieht sich auf grundlegend unterschiedliche Arten, Nextflow-Pipelines zu schreiben.
DSL1 war die ursprüngliche Syntax, bei der Prozesse implizit über Kanäle verbunden wurden.
DSL2, eingeführt in Nextflow 20.07, fügte Modularitätsfunktionen hinzu: die Möglichkeit, Prozesse und Workflows aus anderen Dateien zu importieren, explizite `workflow`-Blöcke und benannte Prozessausgaben.
DSL1 wurde in Nextflow 22.03 als veraltet markiert und in 22.12 entfernt.
Aller moderner Nextflow-Code verwendet DSL2.

**Syntaxparser v1 vs v2** bezieht sich auf verschiedene Parser, die beide mit DSL2-Code funktionieren.
Der v1-Parser ist der ursprüngliche, tolerantere Parser.
Der v2-Parser ist strenger und ermöglicht neue Sprachfunktionen wie statische Typisierung (typisierte Ein- und Ausgaben) und Ausgabedirektiven auf Workflow-Ebene.
Der v2-Parser liefert auch bessere Fehlermeldungen und erkennt mehr Fehler zur Parse-Zeit statt zur Laufzeit.
Der v2-Parser wird in Nextflow 26.04 zum Standard.

Zusammengefasst: DSL2 ist die Sprache, die du schreibst; die Syntaxparser-Version bestimmt, wie streng diese Sprache interpretiert wird und welche erweiterten Funktionen verfügbar sind.

### Nextflow-Version prüfen und festlegen

Du kannst mit dem Befehl `nextflow --version` prüfen, welche Nextflow-Version auf deinem System installiert ist.

Weitere Informationen zum Aktualisieren deiner Nextflow-Version findest du in der Referenzdokumentation unter [Updating Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Den v2-Syntaxparser aktivieren

Um den v2-Syntaxparser für deine aktuelle Sitzung zu **aktivieren**, führe folgenden Befehl in deinem Terminal aus:

```bash
export NXF_SYNTAX_PARSER=v2
```

Um dies dauerhaft zu machen (bis v2 in Nextflow 26.04 zum Standard wird), füge den export-Befehl zu deinem Shell-Profil hinzu (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export NXF_SYNTAX_PARSER=v2' >> ~/.bashrc
source ~/.bashrc
```

Beachte, dass die Umgebungsvariable `NXF_SYNTAX_PARSER=v2` eine vorübergehende Anforderung ist.
Ab Nextflow 26.04 wird der v2-Parser zum Standard und diese Einstellung wird nicht mehr benötigt.

### Den v2-Syntaxparser deaktivieren

Um den v2-Syntaxparser für deine aktuelle Sitzung zu **deaktivieren**, führe folgenden Befehl in deinem Terminal aus:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Bestehenden Code migrieren

Für Anleitungen zur Migration von bestehendem Code, um mit neueren Nextflow-Versionen kompatibel zu sein, siehe die [Migration Notes](https://www.nextflow.io/docs/latest/migrations/index.html) in der Referenzdokumentation.

Diese beiden Artikel sind besonders hilfreich für die Migration zur aktuellsten Version:

- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrating to static types](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Beide Funktionen werden ab Version 3.0 der Trainingsmaterialien im Anfängertraining behandelt.

Je nachdem, welche Generation von Nextflow-Code du migrieren möchtest, kannst du möglicherweise den Großteil der Arbeit mit dem Nextflow-Linter über den Befehl `nextflow lint -format` erledigen.
Siehe die CLI-Referenz für [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint) für weitere Details.

Wir hoffen, dass dies hilfreich ist.
Wenn du Hilfe benötigst, melde dich auf Slack oder im Forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
