# GitHub Codespaces

GitHub Codespaces ist eine webbasierte Plattform, die es uns ermöglicht, eine vorkonfigurierte Umgebung für das Training bereitzustellen, die von virtuellen Maschinen in der Cloud unterstützt wird.
Die Plattform wird von GitHub (das zu Microsoft gehört) betrieben und ist kostenlos (mit Nutzungskontingenten) für alle mit einem GitHub-Account zugänglich.

!!! warning "Warnung"

    Accounts, die mit Organisationen verbunden sind, können bestimmten zusätzlichen Einschränkungen unterliegen.
    Falls das auf dich zutrifft, musst du möglicherweise einen unabhängigen persönlichen Account verwenden oder stattdessen eine lokale Installation nutzen.

## Einen GitHub-Account erstellen

Du kannst einen kostenlosen GitHub-Account auf der [GitHub-Startseite](https://github.com/) erstellen.

## Deinen GitHub Codespace starten

Sobald du bei GitHub angemeldet bist, öffne diesen Link in deinem Browser, um die Nextflow-Trainingsumgebung zu öffnen: <https://codespaces.new/nextflow-io/training?quickstart=1&ref=master>

Alternativ kannst du auf den unten gezeigten Button klicken, der in jedem Trainingskurs wiederholt wird (typischerweise auf der Orientierungsseite).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Du solltest eine Seite sehen, auf der du einen neuen GitHub Codespace erstellen kannst:

![Create a GitHub Codespace](img/codespaces_create.png)

### Konfiguration

Für die allgemeine Nutzung solltest du nichts konfigurieren müssen.
Sofern im Kurs, den du beginnst, nicht anders angegeben, kannst du einfach auf den Hauptbutton klicken, um fortzufahren.

Es ist jedoch möglich, die Umgebung anzupassen, indem du auf den Button "Change options" klickst.

??? info "Konfigurationsoptionen"

    Wenn du auf den Button "Change options" klickst, hast du die Möglichkeit, Folgendes anzupassen:

    #### Branch

    Hier kannst du eine andere Version der Trainingsmaterialien auswählen.
    Der `master`-Branch enthält in der Regel Fehlerbehebungen und Materialien, die kürzlich entwickelt und genehmigt wurden, aber noch nicht auf der Website veröffentlicht sind.
    Andere Branches enthalten laufende Arbeiten, die möglicherweise noch nicht vollständig funktionsfähig sind.

    #### Machine type

    Hier kannst du die virtuelle Maschine anpassen, die du für das Training verwenden wirst.

    Die Verwendung einer Maschine mit mehr Kernen ermöglicht es dir, die Fähigkeit von Nextflow zur Parallelisierung der Workflow-Ausführung besser zu nutzen.
    Allerdings verbraucht dies dein kostenloses Kontingent schneller, daher empfehlen wir, diese Einstellung nicht zu ändern, es sei denn, dies wird in den Anweisungen für den Kurs, den du absolvieren möchtest, empfohlen.

    Siehe 'GitHub Codespaces-Kontingente' weiter unten für weitere Details zu Kontingenten.

### Startzeit

Das erstmalige Öffnen einer neuen GitHub Codespaces-Umgebung kann mehrere Minuten dauern, da das System deine virtuelle Maschine einrichten muss. Mach dir also keine Sorgen, wenn es eine Wartezeit gibt.
Es sollte jedoch nicht länger als fünf Minuten dauern.

## Navigation in der Trainingsumgebung

Sobald dein GitHub Codespace geladen ist, solltest du etwas Ähnliches wie das Folgende sehen (das je nach deinen Account-Einstellungen im Light-Modus geöffnet werden kann):

![GitHub Codespaces welcome](img/codespaces_welcome.png)

Dies ist die Benutzeroberfläche der VSCode IDE, einer beliebten Code-Entwicklungsanwendung, die wir für die Nextflow-Entwicklung empfehlen.

- **Der Haupteditor** ist der Bereich, in dem Nextflow-Code und andere Textdateien geöffnet werden. Hier wirst du Code bearbeiten. Wenn du den Codespace öffnest, wird hier eine Vorschau der `README.md`-Datei angezeigt.
- **Das Terminal** unterhalb des Haupteditors ermöglicht es dir, Befehle auszuführen. Hier wirst du alle Befehlszeilen ausführen, die in den Kursanweisungen angegeben sind.
- **Die Seitenleiste** ermöglicht es dir, deine Umgebung anzupassen und grundlegende Aufgaben auszuführen (kopieren, einfügen, Dateien öffnen, suchen, Git usw.). Standardmäßig ist sie im Datei-Explorer geöffnet, der es dir ermöglicht, den Inhalt des Repositorys zu durchsuchen. Wenn du auf eine Datei im Explorer klickst, wird sie im Haupteditorfenster geöffnet.

Du kannst die relativen Proportionen der Fensterbereiche nach Belieben anpassen.

<!-- TODO (future) Link to development best practices side quest? -->

## Weitere Hinweise zur Verwendung von GitHub Codespaces

### Eine Sitzung fortsetzen

Sobald du eine Umgebung erstellt hast, kannst du sie einfach fortsetzen oder neu starten und dort weitermachen, wo du aufgehört hast.
Deine Umgebung läuft nach 30 Minuten Inaktivität ab und speichert deine Änderungen bis zu 2 Wochen lang.

Du kannst eine Umgebung von <https://github.com/codespaces/> aus erneut öffnen.
Frühere Umgebungen werden aufgelistet.
Klicke auf eine Sitzung, um sie fortzusetzen.

![List GitHub Codespace sessions](img/codespaces_list.png)

Wenn du die URL für deine vorherige GitHub Codespaces-Umgebung gespeichert hast, kannst du sie einfach in deinem Browser öffnen.
Alternativ klicke auf denselben Button, den du verwendet hast, um sie ursprünglich zu erstellen:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

Du solltest die vorherige Sitzung sehen, die Standardoption ist, sie fortzusetzen:

![Resume a GitHub Codespace](img/codespaces_resume.png)

### Dateien auf deinem lokalen Rechner speichern

Um eine Datei aus dem Explorer-Panel zu speichern, klicke mit der rechten Maustaste auf die Datei und wähle `Download`.

### GitHub Codespaces-Kontingente verwalten

GitHub Codespaces gibt dir bis zu 15 GB-Monat Speicher pro Monat und 120 Kern-Stunden pro Monat.
Dies entspricht etwa 60 Stunden Laufzeit der Standardumgebung mit dem Standard-Workspace (2 Kerne, 8 GB RAM und 32 GB Speicher).

Du kannst sie mit mehr Ressourcen erstellen (siehe Erklärung oben), aber dies verbraucht deine kostenlose Nutzung schneller und du hast weniger Stunden Zugriff auf diesen Space.
Wenn du beispielsweise eine 4-Kern-Maschine anstelle der 2-Kern-Standardmaschine auswählst, läuft dein Kontingent in der Hälfte der Zeit ab.

Optional kannst du Zugriff auf mehr Ressourcen kaufen.

Weitere Informationen findest du in der GitHub-Dokumentation:
[About billing for GitHub Codespaces](https://docs.github.com/en/billing/managing-billing-for-your-products/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)
