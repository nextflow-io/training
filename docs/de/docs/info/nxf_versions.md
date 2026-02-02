---
title: Nextflow-Versionen
description: Verständnis und Verwaltung der Entwicklung von Nextflow-Syntaxversionen
hide:
  - toc
  - footer
---

## Aktuell unterstützte Nextflow-Syntaxversion & Anforderungen

Ab Version 3.0 des Trainingsportals basieren alle unsere Trainingskurse auf der Nextflow-Version 25.10.2, sofern auf der Kursindexseite nicht anders angegeben (außer veraltete oder anderweitig archivierte Materialien, die möglicherweise keinen Versionshinweis enthalten).

Da die Kurse jetzt typisierte Eingaben auf workflow-Ebene sowie workflow-Ebene-Ausgabe-Direktiven verwenden, erfordern sie die Verwendung des V2-Syntaxparsers.
Wenn du die Umgebung verwendest, die wir über [Github Codespaces](../envsetup/01_setup.md) oder [lokale Devcontainer](../envsetup/03_devcontainer.md) bereitstellen, musst du nichts unternehmen, es sei denn, dies wird in den Kursanweisungen speziell erwähnt.
Wenn du jedoch planst, die Trainings in deiner eigenen Umgebung durchzuarbeiten ([Manuelle Installation](../envsetup/02_local.md)), musst du sicherstellen, dass du Nextflow Version 25.10.2 oder höher mit aktiviertem v2-Syntaxparser verwendest.

## Ältere Versionen der Trainingsmaterialien

Unsere Trainingsmaterialien sind seit Februar 2025 versioniert.

Du kannst auf ältere Versionen der Trainingsmaterialien zugreifen, die mit Nextflow-Versionen **vor 25.10.2** funktionieren, über das Dropdown-Menü oben auf jeder Seite, das die nummerierte Version der Trainingsmaterialien anzeigt.
Wenn du eine ältere Version der Trainingsmaterialien auswählst, geben Links zur Trainingsumgebung automatisch die entsprechende Version der Umgebung an.

## Weitere Informationen zu Nextflow-Syntaxversionen

Nextflow hat zwei unterschiedliche Versionierungskonzepte, die manchmal verwechselt werden: **DSL-Versionen** und **Syntaxparser-Versionen**.

**DSL1 vs DSL2** bezieht sich auf grundlegend unterschiedliche Arten, Nextflow-pipelines zu schreiben.
DSL1 war die ursprüngliche Syntax, bei der processes implizit über channels verbunden wurden.
DSL2, eingeführt in Nextflow 20.07, fügte Modularitätsfunktionen hinzu: die Möglichkeit, processes und workflows aus anderen Dateien zu importieren, explizite `workflow`-Blöcke und benannte process-Ausgaben.
DSL1 wurde in Nextflow 22.03 als veraltet markiert und in 22.12 entfernt.
Aller moderner Nextflow-Code verwendet DSL2.

**Syntaxparser v1 vs v2** bezieht sich auf verschiedene Parser, die beide mit DSL2-Code funktionieren.
Der v1-Parser ist der ursprüngliche, tolerantere Parser.
Der v2-Parser ist strenger und ermöglicht neue Sprachfunktionen wie statische Typisierung (typisierte Eingaben und Ausgaben) und workflow-Ebene-Ausgabe-Direktiven.
Der v2-Parser bietet auch bessere Fehlermeldungen und erkennt mehr Fehler zur Parse-Zeit statt zur Laufzeit.
Der v2-Parser wird in Nextflow 26.04 zum Standard.

Zusammenfassend: DSL2 ist die Sprache, die du schreibst; die Syntaxparser-Version bestimmt, wie streng diese Sprache interpretiert wird und welche erweiterten Funktionen verfügbar sind.

### Überprüfen und Einstellen der Nextflow-Version

Du kannst überprüfen, welche Version von Nextflow auf deinem System installiert ist, indem du den Befehl `nextflow --version` verwendest.

Weitere Informationen zum Aktualisieren deiner Nextflow-Version findest du in der Referenzdokumentation unter [Updating Nextflow](https://www.nextflow.io/docs/latest/updating-nextflow.html).

### Aktivieren des v2-Syntaxparsers

Um den v2-Syntaxparser für deine aktuelle Sitzung zu **aktivieren**, führe den folgenden Befehl in deinem Terminal aus:

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

### Deaktivieren des v2-Syntaxparsers

Um den v2-Syntaxparser für deine aktuelle Sitzung zu **deaktivieren**, führe den folgenden Befehl in deinem Terminal aus:

```bash
export NXF_SYNTAX_PARSER=v1
```

<!-- Will it be possible to disable it in versions after 26.04? -->

### Migrieren von bestehendem Code

Anleitungen zur Migration von bestehendem Code zur Einhaltung neuerer Nextflow-Versionen findest du in den [Migration Notes](https://www.nextflow.io/docs/latest/migrations/index.html) in der Referenzdokumentation.

Diese beiden Artikel sind besonders hilfreich für die Migration zur neuesten Version:

- [Migrating to workflow outputs](https://www.nextflow.io/docs/latest/tutorials/workflow-outputs.html)
- [Migrating to static types](https://www.nextflow.io/docs/latest/tutorials/static-types.html)

Beide Funktionen werden ab Version 3.0 der Trainingsmaterialien im Einsteiger-Training behandelt.

Je nach Generation des Nextflow-Codes, den du migrieren möchtest, kannst du möglicherweise den größten Teil davon mit dem Nextflow-Linter mithilfe des Befehls `nextflow lint -format` erledigen.
Weitere Details findest du in der CLI-Referenz für [`lint`](https://www.nextflow.io/docs/latest/reference/cli.html#lint).

Wir hoffen, dass dies hilfreich sein wird.
Wenn du Hilfe benötigst, melde dich auf Slack oder im Forum.

---

<div markdown class="homepage_logos">

![Seqera](../assets/img/seqera_logo.png#only-light)

![Seqera](../assets/img/seqera_logo_dark.png#only-dark)

</div>
