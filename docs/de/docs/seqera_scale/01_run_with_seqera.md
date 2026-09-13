---

# Teil 1: Pipelines über die Web-Oberfläche starten

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem Teil des Kurses „Scale with Seqera" richtest du den Zugang zur Seqera Platform ein und startest eine produktionsreife Pipeline über die Web-Oberfläche.

Stelle sicher, dass dein Arbeitsverzeichnis auf `seqera-scale/` gesetzt ist, wie auf der Seite [Erste Schritte](./00_orientation.md) beschrieben.

---

## 1. Erste Schritte mit Seqera

Seqera bietet eine umfassende Platform zum Starten, Überwachen und Verwalten von Nextflow-Pipelines.
Dieser Abschnitt führt dich durch die Registrierung und gibt dir einen ersten Überblick, bevor du deine erste Pipeline ausführst.

### 1.1. Kostenloses Konto erstellen

Gehe zu [cloud.seqera.io](https://cloud.seqera.io) und erstelle ein kostenloses Konto.
Du kannst dich mit deiner E-Mail-Adresse, GitHub oder Google anmelden.

Ein kostenloses Konto bietet dir:

- **Persönlicher Workspace**: dein eigener Bereich, um Pipelines hinzuzufügen, Compute-Umgebungen zu konfigurieren und Ausführungen zu verwalten
- **Zugang zur Community Showcase**: eine kuratierte Sammlung von nf-core- und Community-Pipelines mit vorkonfigurierten Einstellungen und Beispieldaten

Eine vollständige Übersicht der Kontostufen und verfügbaren Funktionen findest du in der [Seqera-Dokumentation](https://docs.seqera.io).

### 1.2. Die Community Showcase erkunden

Bevor du eigene Pipelines startest, nimm dir ein paar Minuten, um die Community Showcase zu erkunden.
Sie gibt dir einen realistischen Eindruck davon, wie die Platform mit echten Pipelines und Daten aussieht.

1. Melde dich bei [cloud.seqera.io](https://cloud.seqera.io) an.
2. Klicke in der linken Seitenleiste auf **Showcase**.
3. Stöbere durch die verfügbaren Pipelines — du wirst mehrere nf-core-Pipelines aus dem Kurs „Use nf-core" wiedererkennen.
4. Klicke auf eine Pipeline, um ihre Konfiguration und Starteinstellungen anzusehen.
5. Klicke auf **Runs**, um Beispiel-Ausführungsverläufe zu erkunden, einschließlich Details auf Aufgaben-Ebene und Berichte aus früheren Ausführungen.

Dies ist eine schreibgeschützte Ansicht, aber sie zeigt dir, wie die Oberfläche funktioniert, bevor du selbst etwas ausführst.

### 1.3. Auf einen Workspace mit Compute zugreifen

Zum Starten von Pipelines benötigst du einen Workspace mit einer konfigurierten Compute-Umgebung.

Seqera unterstützt zwei Möglichkeiten, Compute bereitzustellen:

- **Eigene Infrastruktur verbinden**: AWS, Azure, Google Cloud und HPC-Scheduler (SLURM, LSF, PBS und andere).
  Einrichtungsanleitungen findest du in der [Dokumentation zu Compute-Umgebungen](https://docs.seqera.io).
- **Seqera Compute**: ein verwalteter Dienst, der vorkonfigurierte Compute-Umgebungen auf AWS bereitstellt — kostenpflichtig, ohne dass ein eigenes Cloud-Konto eingerichtet werden muss.
  Du kannst ihn direkt über deine Workspace-Einstellungen aktivieren.

**Gruppentraining:**
Wenn du an einer Gruppentraining-Session teilnimmst, wurdest du möglicherweise einer Organisation und einem Workspace hinzugefügt, in dem bereits Compute konfiguriert ist.
Deine\*r Trainer\*in gibt dir den Namen der Organisation, den Namen des Workspace und alle weiteren benötigten Details.

**Selbstständiges Arbeiten:**
Wenn du dieses Training eigenständig durcharbeitest, musst du eine Compute-Umgebung in deinem persönlichen Workspace mit einer der oben genannten Optionen einrichten.
Kostenlose Credits zum Ausprobieren von Seqera Compute sind [auf Anfrage verfügbar](https://seqera.io/platform/compute/).

!!! note "Hinweis"

    Der Rest dieses Kurses setzt voraus, dass du Zugang zu einem Workspace mit einer konfigurierten Compute-Umgebung hast.
    Wenn du an einer Gruppentraining-Session teilnimmst, wird deine\*r Trainer\*in bestätigen, welchen Workspace und welche Compute-Umgebung du verwenden sollst.

### Fazit

Du hast ein Seqera-Konto, hast die Community Showcase erkundet und kannst auf einen Workspace mit Compute zugreifen.

### Wie geht es weiter?

Starte eine produktionsreife RNA-seq-Pipeline über die Seqera Cloud Web-Oberfläche.

---

## 2. nf-core/rnaseq über die Web-Oberfläche starten

Wie im Kurs „Use nf-core" behandelt, ist die nf-core/rnaseq-Pipeline eine von der Community kuratierte Pipeline zur Analyse von Bulk-RNA-Sequenzdaten.

In diesem Abschnitt fügst du die Pipeline zu deinem Workspace hinzu, startest eine Ausführung und überwachst deren Ablauf.

### 2.1. Die Pipeline zum Workspace hinzufügen

Praktischerweise ist nf-core/rnaseq Teil einer kuratierten Pipeline-Sammlung, die über den Seqera Pipelines-Dienst mit wenigen Klicks zu deinem Workspace hinzugefügt werden kann.

_Wie du eigene Pipelines hinzufügst, zeigen wir dir später in diesem Kurs._

1. Navigiere zu [**Seqera Pipelines**](https://seqera.io/pipelines), um die Community-Sammlung zu durchsuchen.
2. Suche nach `rnaseq` und wähle **nf-core/rnaseq** aus.
3. Klicke auf **Launch Pipeline** oder scrolle zum Ende der Seite zum Abschnitt **Launch Pipeline**.
4. Stelle sicher, dass du angemeldet bist, und wähle die passenden Werte aus den Dropdown-Menüs **Organizations**, **Workspace** und **Compute Environment** aus.
   **Tipp für Gruppen:** Wenn du einen gemeinsamen Workspace verwendest, füge dem Pipeline-Namen eine eindeutige Kennung hinzu (z. B. deinen Benutzernamen).
5. Klicke auf **Add pipeline to your Seqera account**.

Es erscheint ein Hinweisfeld mit der Meldung: **Pipeline added: View Pipeline**.
Ein Klick auf den Link bringt dich zum Pipeline-Eintrag in deinem Launchpad.

Die Pipeline ist jetzt im **Launchpad**-Panel deines Workspace aufgelistet und kann gestartet werden.

### 2.2. Die Pipeline starten

Klicke auf die Schaltfläche **Launch** der Pipeline, entweder im **Launchpad**-Panel oder auf der Pipeline-Detailseite.
Dadurch öffnet sich die Konfigurationsoberfläche.

Die Pipeline ist bereits mit dem `test`-Profil konfiguriert, sodass die Eingabedaten, das Ausgabeverzeichnis und die Genomreferenz bereits vorausgefüllt sind.
Die übrigen Parameter und erweiterten Einstellungen kannst du vorerst ignorieren.

Klicke auf die blaue Schaltfläche **Launch**, um die Ausführung tatsächlich zu starten.

### 2.3. Ausführung überwachen

Nach dem Start wirst du zum **Runs**-Panel deiner Pipeline weitergeleitet.

Die Ausführungsansicht zeigt:

- **Status**: aktueller Zustand der Ausführung (submitted, running, succeeded, failed)
- **Command line**: der genaue `nextflow run`-Befehl, den die Platform zusammengestellt und übermittelt hat
- **Parameters**: alle für diese Ausführung verwendeten Parameterwerte
- **Tasks**: eine Tabelle aller Prozessaufrufe mit Status, Dauer und Ressourcennutzung

Klicke auf eine beliebige Aufgaben-Zeile, um die Ausführungsdetails zu prüfen, darunter:

- Das `.command.sh`-Skript, das ausgeführt wurde
- stdout- und stderr-Logs
- CPU-, Arbeitsspeicher- und I/O-Metriken

Der Tab **Reports** zeigt einen MultiQC-Bericht, sobald die Ausführung abgeschlossen ist, der Qualitätskontroll-Metriken über alle Proben hinweg zusammenfasst.

Die Ausführung dauert eine Weile — wir machen daher zunächst weiter und schauen uns die Ausgaben später an.

### Fazit

Du weißt, wie du eine Pipeline zu einem Seqera-Workspace hinzufügst, eine Ausführung konfigurierst und startest und die Ausführung in großem Maßstab überwachst.

### Wie geht es weiter?

Weiter zu [Teil 2](./02_launch_from_cli.md), wo du lernst, all das über die Befehlszeile mit der `tw` CLI zu erledigen.

---

## Zusammenfassung

In diesem Teil hast du gelernt, wie du:

- Ein Seqera-Konto erstellst und die Community Showcase erkundest
- Eine Pipeline aus dem kuratierten Katalog hinzufügst, eine produktionsreife Ausführung startest und die Ausführung überwachst
