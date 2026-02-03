# Einführung - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtiger Hinweis"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../00_orientation.md) zurück.

## Willkommen

Hallo, willkommen zu Hello Nextflow. Mein Name ist Phil Ewels. Ich bin Product Manager für Open Source bei Seqera und freue mich sehr, dich heute durch diesen ersten Nextflow-Trainingskurs zu führen.

Wir werden die Grundlagen von Nextflow durchgehen und erklären, wie man Pipelines schreibt, ausführt und konfiguriert.

Und du wirst deine eigene einfache mehrstufige Pipeline erstellen. Wir werden Begriffe wie Operatoren und Channel Factories behandeln, und am Ende des Kurses bist du bereit, mit dem Erstellen deiner eigenen Bioinformatik-Pipelines zu beginnen.

Falls du Fragen hast, melde dich bitte auf community.seqera.io. Wir haben dort eine sehr aktive Nextflow-Community mit einem Bereich, der speziell dem Training gewidmet ist. Lass uns einfach wissen, wo du nicht weiterkommst, und jemand wird dir helfen können.

Gut. Dann lass uns anfangen.

## Trainings-Website

Das gesamte Trainingsmaterial für die Nextflow-Kurse findest du auf training.nextflow.io. Du kannst es in deinem Webbrowser aufrufen. Starte das jetzt und wir können uns gemeinsam umschauen.

Ich werde das mit Version 2.1.1 durchführen. Wir veröffentlichen hier und da kleinere Updates und Korrekturen. Wenn es also ein wenig abweicht, ist das kein Problem. Sollte das Material aber zu weit abgedriftet sein, kannst du immer diese Versionsauswahl oben verwenden, um genau die Version des Materials auszuwählen, über die ich sprechen werde.

Wenn du eher ein Heller-Modus-Mensch bist, kannst du das Theme für die Website hier ändern.

Hier findest du auch Übersetzungen, auch wenn zum Zeitpunkt der Aufnahme wirklich nur Englisch verfügbar ist, das dieses neue Material abdeckt.

Und du kannst auch den gesamten Quellcode für die Trainings-Website und alles, womit wir arbeiten werden, auf GitHub einsehen.

Die Startseite hier listet alle verschiedenen Trainingsmaterial-Kurse auf, die wir haben. Wenn ich nach unten scrolle, sehen wir Nextflow for newcomers mit dem Hello Nextflow-Kurs, den wir hier machen werden. Du siehst auch alle anderen Kurse, die wir haben und die auf ähnliche Weise funktionieren.

## Umgebungseinrichtung

Ich werde tatsächlich mit diesem ersten ganz oben beginnen, der für alle Trainingskurse gemeinsam ist, und es geht speziell um die Einrichtung unserer Umgebung.

Ich klicke darauf, und es bringt mich zu diesem Abschnitt, wo wir Anweisungen für die lokale Entwicklung sehen können. Wenn du deinen eigenen Laptop mit deiner eigenen Kopie von VS Code und deinen eigenen Software-Installationen verwenden möchtest, oder was wir von den meisten Leuten erwarten, nämlich etwas namens GitHub Codespaces zu verwenden.

Codespaces ist ein von GitHub bereitgestellter Service, bei dem ein Webserver in der Cloud ausgeführt wird, mit dem du dich verbinden kannst. Auf diesem Server ist VS Code installiert, das du in deinem Webbrowser ausführen kannst, oder wenn du möchtest, mit deiner lokalen Installation von VS Code verbinden kannst. Die gesamte Berechnung, alle Dateien, alle Bearbeitungen geschehen remote, was bedeutet, dass die gesamte Software, die du benötigst, vorinstalliert ist und für alle gleich ist.

## Erstellen eines GitHub Codespace

Um den Codespace mit allem, was wir brauchen, zu erstellen, suche nach den Buttons im Dokumentationsmaterial, die "Open in GitHub Codespaces" sagen. Ich werde das jetzt anklicken und in einem neuen Tab öffnen. Und mir wird diese Webseite präsentiert. Du siehst, dass dies vorkonfiguriert ist mit nextflow-io training.

Ich kann einfach auf create new codespace klicken. Aber tatsächlich empfehlen wir, eine etwas größere Maschine für das Nextflow-Training mit vier CPUs statt zwei zu verwenden. Du kannst ändern, welche Version des Materials verwendet wird. Das ist standardmäßig auf 2.1.1 eingestellt, weil das die Version der Dokumentation ist, von der ich den Link gefolgt bin. Aber ich könnte es auch auf einen bestimmten Branch des Repositorys setzen, wenn ich möchte.

Jetzt werde ich auf create codespace klicken. Und es wird beginnen, die Umgebung für mich einzurichten.

## Codespace-Erstellung

Beim ersten Mal wird das ziemlich lange dauern, also ist jetzt ein guter Zeitpunkt, um eine Tasse Tee zu holen. Mach es dir bequem, unterhalte dich mit der Person neben dir.

Wenn es dich interessiert, kannst du hier auf building codespace klicken, um die Logs der Einrichtung zu sehen. Und du siehst hier, dass es ein Docker-Image mit allem, was ich brauche, herunterlädt und die Umgebung konfiguriert.

Du musst nur beim ersten Mal so lange warten, wenn du einen Codespace erstellst. Wenn du zu github.com/codespaces gehst, siehst du alle verschiedenen Codespaces, die du geöffnet hast. Hier ist der, den ich gerade erstellt habe. Beim nächsten Mal kannst du hierher gehen und den vorherigen Codespace auswählen und direkt wieder hineinspringen. Und es ist ein viel, viel schnellerer Prozess, diese bestehende Umgebung aufzuwärmen. Das behält auch alle Änderungen bei, die du an VS Code und an den Dateien vorgenommen hast, sodass du deinen Fortschritt nicht verlierst, wenn du gehst und zurückkommst.

Du kannst hier auf die drei Punkte klicken, um andere Aktionen durchzuführen. Zum Beispiel, wenn du ihn mit zwei CPUs konfiguriert hast und jetzt vier möchtest, kannst du den Maschinentyp ändern. Oder wenn du von vorne anfangen und einen frischen Start haben möchtest, kannst du den Codespace löschen.

## Einführung in VS Code

Okay, Codespaces ist mit der Einrichtung meiner Umgebung fertig und präsentiert mir jetzt VS Code im Webbrowser.

Wenn du VS Code kennst, wird sich das sehr vertraut anfühlen. Falls du es noch nicht benutzt hast, ist es ziemlich einfach. Es gibt ein paar verschiedene Teile der Seite, die du kennen solltest.

Hier auf der linken Seite haben wir die Seitenleiste. Du siehst den Explorer mit all den verschiedenen Dateien aus dem GitHub-Repository des Training-Repos.

Mit diesen Buttons unten links können verschiedene Werkzeuge in der Seitenleiste angezeigt werden. Ich kann alle Dateien im gesamten Projekt durchsuchen. Ich kann mit Git arbeiten, kann mit GitHub arbeiten, all solche verschiedenen Dinge.

Oben ist das Hauptmenü. Der Datei-Explorer ist derjenige, den wir hier am meisten offen haben werden, und du kannst mit der rechten Maustaste auf eine dieser Dateien klicken und die normalen Dinge tun, die du erwarten würdest. Möglicherweise musst du einige Warnungen wie diese durchklicken, wo es um Ausschneiden und Kopieren geht, und du kannst auch auf deine lokale Maschine herunterladen.

Wenn der Codespace lädt, zeigt er uns eine Vorschau der Markdown-Datei in diesem Hauptbereich hier. Das ist genau dasselbe wie das, was auf github.com gerendert wird. Ich kann das schließen, und wenn ich doppelt auf diese Readme-Datei klicke, siehst du, dass sie als Code im Code-Editor geöffnet wird, und genau wie bei jeder anderen Datei können wir diesen Code direkt bearbeiten.

Schließlich haben wir hier unten das Terminalfenster. Ich habe mir die Logs beim Erstellen angesehen, also ist das, was es gerade zeigt. Ich kann auch diesen Plus-Button drücken, um eine neue Terminal-Sitzung zu starten. Das läuft nicht auf meiner Maschine. Denk daran, das läuft in der Cloud, und wenn ich tree mit Tiefe zwei mache, siehst du alle dieselben Dateien hier, die auch links waren.

## Nur "hello-nextflow"-Dateien anzeigen

Dieses GitHub-Repository enthält alle verschiedenen Trainingssets, nicht nur das, was wir machen. Wenn du möchtest, kannst du dich nur auf den Hello-Nextflow-Ordner konzentrieren. Eine Möglichkeit, das ein wenig aufzuräumen, ist, zum Menü Datei zu gehen und dann "Ordner zum Arbeitsbereich hinzufügen" auszuwählen.

Wir klicken darauf, gehen zu training, Hello nextflow, und klicken auf Hinzufügen. Es wird deinen Bildschirm aktualisieren. Und dann haben wir im Explorer jetzt zwei verschiedene Arbeitsbereiche, den, den wir vorher für training hatten, und einen mit nur Hello Nextflow.

Wenn du möchtest, kannst du mit der rechten Maustaste auf training klicken und "Ordner aus Arbeitsbereich entfernen" klicken, um ihn vollständig aus der Seitenleiste zu entfernen.

Jetzt haben wir nur die Dateien für diesen bestimmten Trainingskurs in der Seitenleiste. Ich kann diese Warnung ausblenden, und jetzt kann ich dasselbe im Terminal hier machen und cd für Verzeichniswechsel machen. Hello Nextflow. Und wieder haben wir hier dieselben Dateien, die in der Seitenleiste sind.

## Hello Nextflow: Dateien

Schauen wir uns diese Dateien für den Hello-Nextflow-Kurs an.

Wir haben eine Reihe von .nf-Dateien, die für Nextflow sind, und es gibt eine dieser Dateien für jedes Kapitel des Trainingskurses. Wir werden diese Dateien durcharbeiten und sie in den Übungen modifizieren.

Wir haben auch eine nextflow.config-Datei, die nur grundlegende Konfigurationseinstellungen für das Ausführen von Nextflow in dieser Umgebung enthält, um die du dir zu diesem Zeitpunkt keine Sorgen machen musst. Eine greetings.csv-Datei, die wir für die Verarbeitung von Daten verwenden werden und die im nächsten Teil dieses Kurses eingeführt wird, und eine test-params.json-Datei, die in Teil sechs verwendet wird und die du vorerst ignorieren kannst.

Diese Nextflow-Dateien sind nur der Anfang jeder Übung. Wenn du sehen möchtest, wie sie aussehen sollten, wenn sie fertig sind, kannst du in ein solutions-Verzeichnis gehen, und dort sind die Antworten für jeden Teil des Trainingskurses, sodass du eine funktionierende Version dessen sehen kannst, worauf du hinarbeitest.

## Ein Terminal öffnen

Falls du irgendwann das Terminal schließt und dich nicht erinnern kannst, wie du zurückkommst, mach dir keine Sorgen. Diese Buttons oben rechts öffnen und schließen verschiedene Panels im Arbeitsbereich. Klick also auf dieses hier für das untere Panel und es wird wieder erscheinen. Und stelle nur sicher, dass du hier terminal ausgewählt hast. Du kannst auch diesen Button hier klicken, den Pfeil auf der rechten Seite eines Terminals, um es im Vollbildmodus anzuzeigen.

Du wirst sehen, dass ich das ziemlich oft mache, weil ich VS Code vergrößert habe, damit du den Text lesen kannst. Abhängig von deiner Bildschirmgröße musst du das möglicherweise tun oder auch nicht. Dasselbe gilt für das Minimieren des Seitenpanels.

Gut. Das reicht für die Umgebung. Ich denke, wir sind bereit anzufangen. Komm im nächsten Video zu mir zurück für Kapitel eins.

[Nächstes Video-Transkript :octicons-arrow-right-24:](01_hello_world.md)
