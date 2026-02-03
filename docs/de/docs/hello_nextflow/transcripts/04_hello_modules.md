# Teil 4: Hello Modules - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../04_hello_modules.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen

Hi, willkommen zu Teil Vier des Hello Nextflow Trainingskurses.

Dieses Kapitel heißt Hello Modules, und wir werden darüber sprechen, wie man Nextflow-Code modularisiert. Was wir tun werden, ist unser Workflow-Skript zu nehmen und es in separate Dateien aufzuteilen.

Das macht den Code einfacher zu navigieren und zu warten, wenn dein Workflow größer wird, und ermöglicht es auch, Module zwischen Pipelines zu teilen, sodass du, wenn du mehrere Pipelines mit demselben Tool hast, diesen Prozess nur einmal schreiben musst.

Ein klassisches Beispiel dafür ist das nf-core Modules Repository, das Tausende verschiedener Tools in gebrauchsfertigen Modulen hat, die du installieren und in deinem Workflow verwenden kannst.

Nextflow kann auch mit Sub-Workflows arbeiten, die wie Module sind, aber mehrere Prozesse haben. Das liegt außerhalb des Umfangs dieses Trainings, aber es funktioniert im Grunde auf die gleiche Weise.

Okay. Lass uns einen Blick darauf werfen.

Wie üblich beginne mit dem Aufrufen von training.nextflow.io.

Gehe zu "Hello Nextflow" in der Seitenleiste, und wir machen Teil vier: "Hello Modules".

Ich werde jetzt in meine GitHub Codespaces Umgebung springen und mir die Datei "hello-modules" ansehen.

Genau wie zuvor starten wir am Endpunkt des vorherigen Kapitels, also sollte dir dieses Skript bekannt vorkommen. Wir haben unsere drei Prozesse, say hello, convert to upper und collect greetings, und einen einfachen Workflow, der diese drei Befehle ausführt und am Ende eine Nachricht ausgibt. Wir haben zwei Parameter namens greeting und batch, die den Namen angeben, der für die gesammelte Ausgabedatei am Ende verwendet wird.

## 0. Aufwärmen: hello-modules.nf ausführen

Wir können überprüfen, dass dieser Workflow noch wie erwartet funktioniert, indem wir nextflow run hello modules ausführen.

Großartig. Es hat drei Aufgaben mit jedem dieser Prozesse ausgeführt, eine collect-Aufgabe, und es hat uns mitgeteilt, dass es drei Grüße in diesem Batch gibt. Wenn wir in results gehen, haben wir hier unsere verschiedenen Ausgabedateien, einschließlich der gesammelten test batch Ausgabe.

## 1. Ein Verzeichnis zum Speichern von Modulen erstellen

Richtig. Lass uns etwas Modularisierung vornehmen.

Es ist generell eine gute Idee, Module in einen Unterordner in deinem Pipeline-Repository zu legen, einfach um die Dinge aufgeräumt zu halten. Du kannst das nennen, wie du willst, aber per Konvention nennen wir es normalerweise modules.

Also lass uns loslegen, gehe in ein Terminal und mache make dir modules. Du kannst sehen, wie es in der Seitenleiste in VS Code hier auftaucht.

## 2. Ein Modul für sayHello() erstellen

Ich werde dann eine neue Datei für mein erstes Modul erstellen. Du kannst "touch" oder "code" verwenden oder du kannst das in der Seitenleiste machen, es ist wirklich egal. Also werde ich code modules machen und ich werde es nach dem Prozess benennen. Also sayHello.nf. NF ist eine traditionelle Dateierweiterung für Nextflow-Dateien.

Werde hier speichern und wir können sehen, wie unsere neue Moduldatei auftaucht.

## 2.2. Den sayHello Prozess-Code in die Moduldatei verschieben

Richtig, als Nächstes werde ich den Modulcode aus dem Workflow nehmen. Ich werde auch den Shebang hier nehmen und den zuerst kopieren, damit es eindeutig eine Nextflow-Datei ist. Und dann werde ich diesen Prozess nehmen und ich werde ihn ausschneiden. Also werde ich ihn aus meinem Haupt-Workflow-Skript entfernen und ich werde ihn in dieses neue Modul einfügen.

Das ist der gesamte Inhalt, den diese Moduldatei enthalten wird. Nur ein einzelner Prozess, kein Workflow, keine Logik, nur ein Prozess allein.

Ich kann diese Datei jetzt schließen.

## 2.3. Eine Import-Deklaration vor dem Workflow-Block hinzufügen

Jetzt fehlt meinem Workflow dieser erste Prozess, also müssen wir ihn zurückbringen, indem wir ihn importieren. Die Syntax dafür ist sehr ähnlich zu anderen Programmiersprachen, also könnte sie sich vertraut anfühlen. Wir machen include geschweifte Klammern, den Namen des Prozesses, say hello, und dann from dem Dateipfad modules, say hello, nf. Fantastisch.

Ein paar Tricks hier. Die VS Code Erweiterung ist clever in dieser Hinsicht. Sie erkennt diesen Dateipfad und du kannst mit der Maus darüber fahren und auf Link folgen klicken. Oder ich bin auf Mac, ich kann Option-Klick machen und es öffnet diese Datei. So können wir schnell zu ihr springen.

Dieser Prozessname wird jetzt vom Workflow hier unten verwendet, und wir können dasselbe hier machen. Es zeigt uns ein bisschen Information über diesen Prozess, und wieder kann ich Option halten, darauf klicken, und es wird ihn im Editor öffnen.

Es ist also eine wirklich schnelle Möglichkeit, wenn du viele Dateien für deine verschiedenen Prozesse hast, schnell in deiner Codebasis in VS Code zu navigieren.

Okay. Das ist im Grunde alles für dieses Kapitel. Jetzt machen wir einfach dasselbe noch einmal für die anderen Prozesse.

## 3. Den convertToUpper() Prozess modularisieren

Also lass uns hier eine neue Datei erstellen. Nenne sie Convert to upper nf. Wieder den Shebang kopieren. Und dann den Prozess ausschneiden.

Den Prozessnamen dort kopieren, ein neues include-Statement mit dem neuen Prozessnamen einfügen.

## 4. Den collectGreetings() Prozess modularisieren

Und dann dasselbe für den dritten Prozess machen. Neue Datei, collect greetings,

den Shebang machen. Den Prozess ausschneiden, den Prozess einfügen, und ein neues include-Statement machen.

Jetzt kannst du hier sehen, dass ich eine Fehlerunterstreichung habe, die sagt invalid include source. Und das ist tatsächlich ein echter Fehler, den ich gemacht habe, weil ich ein bisschen zu schnell unterwegs war. Wenn du genau hinsiehst, kannst du sehen, dass ich das T verpasst habe in convert to upper.

Also hat VS Code mir sehr hilfreich mitgeteilt, dass ich dort einen Fehler gemacht habe. Wenn ich diesen Dateinamen korrigiere, verschwindet der Fehler. Es ist ein gutes Beispiel dafür, warum die Fehlerprüfung in VS Code so nützlich ist für das Schreiben von Nextflow-Code. Sonst hätte ich das nicht bemerkt und hätte es erst viel später herausgefunden, wenn ich versucht hätte, den Workflow auszuführen.

Unser Haupt-Pipeline-Skript sieht jetzt viel einfacher aus. Es hat keine Prozesse mehr drin, wir haben nur drei include-Statements und unseren Workflow. Wir haben nichts an der Logik des Workflows geändert. Wir haben nichts am Prozess-Code geändert, also sollte es hoffentlich auf genau die gleiche Weise funktionieren.

## 4.4. Den Workflow ausführen, um zu überprüfen, dass er dasselbe tut wie zuvor

Lass uns das überprüfen. Ich werde ein Terminal öffnen und ich werde genau denselben Befehl wie zuvor ausführen.

Sicher genug, es hat unsere Prozesse ausgeführt, say hello, convert to upper, collect greetings, und uns wieder drei Grüße gegeben.

Also haben wir unseren Code herumverschoben, aber wir haben nichts daran geändert, wie der Workflow ausgeführt wird, und er ist komplett unverändert. Der einzige Unterschied ist, dass wir jetzt saubereren Code haben, einfacher zu warten und einfacher mit anderen zu teilen.

Und das war's. Es war ein kurzes Kapitel. Es ist ein einfaches Konzept, aber es ist sehr mächtig und zentral dafür, wie wir komplexere Nextflow-Workflows schreiben. Also ist es wichtig, dass du es verstehst und dir angewöhnst, es zu verwenden.

Im nächsten Kapitel werden wir das Tempo ein wenig wechseln und aufhören, so viel über die Syntax des Schreibens von Nextflow-Code nachzudenken, und ein wenig darüber nachdenken, wie wir Software in den Prozessen selbst verwenden. Begleite uns in Teil fünf für Hello Containers.

[Nächstes Videotranskript :octicons-arrow-right-24:](05_hello_containers.md)
