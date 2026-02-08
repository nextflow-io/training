# Teil 4: Hello Modules - Videotranskript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../04_hello_modules.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur der Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen

Hallo und herzlich willkommen zurück zu Teil vier von Hello Nextflow. In diesem Abschnitt dreht sich alles um Module, und es ist ein ziemlich kurzer Abschnitt des Kurses. Wir werden eigentlich nicht wirklich viel Code schreiben, es geht mehr darum, wie wir den Code in unserer Pipeline organisieren.

Bis jetzt haben wir alles in eine einzige Datei gepackt, was völlig in Ordnung ist, und so haben wir früher tatsächlich Nextflow-Pipelines gebaut.

Aber wenn die Pipeline wächst, wird das Skript immer länger und länger und länger und immer schwieriger zu navigieren, zu warten, und es bedeutet auch, dass wir den Code nicht wirklich teilen können.

Nextflow-Module ermöglichen es uns, Prozesse aus dem Hauptskript herauszuziehen und sie dann zu importieren. Das bedeutet, dass der Code einfacher zu navigieren ist und es bedeutet auch, dass wir den Modulcode zwischen verschiedenen Pipelines teilen können.

Dieses kleine Diagramm auf der Hauptseite der Dokumentation zeigt das Konzept schön. Anstatt eines riesigen Skripts werden wir diese separaten Moduldateien aus verschiedenen Modulskripten einbinden, und alles wird in den Workflow eingefügt, aber es wird trotzdem auf genau die gleiche Weise ausgeführt.

Also lass uns in GitHub Codespaces springen und uns ein bisschen umsehen. Wie zuvor habe ich meinen Workspace hier etwas aufgeräumt. Die alten Nextflow-Verzeichnisse und das Work-Verzeichnis und so weiter entfernt. Aber es macht nichts, wenn du diese Dateien noch hast.

Ich werde in der hello modules Datei arbeiten, die im Grunde dort ist, wo wir sie am Ende des vorherigen Kapitels gelassen haben. Wir haben unsere drei Prozesse hier drin. Wir haben ein paar Parameter, den Workflow-Block, in dem wir diese drei Prozesse ausführen und sie mit Kanälen zusammenfügen. Dann veröffentlichen wir die Ausgabekanäle und wir haben den Output-Block, der angibt, wie diese Dateien veröffentlicht werden.

## 1. Erstelle ein Verzeichnis zum Speichern von Modulen

Nun, wie ich schon sagte, werden wir nicht wirklich sehr viel Code schreiben oder bearbeiten. Wir werden nur den Code verschieben, den wir bereits haben. Nextflow-Moduldateien haben normalerweise einen einzelnen Prozess darin, und per Konvention bewahren wir sie normalerweise in einem Verzeichnis namens modules auf. Aber du kannst das nennen, wie du willst. Aber ich werde ein modules-Verzeichnis in meinem Repository hier behalten, und dann werde ich eine Datei für jeden Prozess erstellen. Also sage ich neue Datei, sayHello.nf.

## 2. Erstelle ein Modul für sayHello()

Jetzt werde ich meinen Prozess nehmen und ich werde einfach diesen Code auswählen, ihn aus der Haupt-hello-modules-Datei ausschneiden und hier einfügen.

Offensichtlich macht das alleine nichts. Unser Hauptskript braucht immer noch diesen Prozess, also müssen wir ihn irgendwie wieder hineinholen. Und das tun wir mit der include-Anweisung.

Also tippe ich include und einige geschweifte Klammern, und dann nehme ich den Namen des Prozesses. Und ich sage from, und dann gebe ich einen relativen Dateipfad an. Also sagt es, beginnt mit ./ weil es relativ zu dem Ort ist, wo dieses Skript gespeichert ist. Also ist es modules sayHello.nf.

Beachte, dass die VS Code Extension hier ziemlich hilfreich ist. Sie sagt uns, ob sie diese Datei finden kann und ob sie einen Prozess finden kann, den ich benenne. Wenn ich hier absichtlich einen Tippfehler mache, gibt sie mir sofort einen Fehler und sagt mir, dass sie diesen Prozess, den ich zu importieren versuche, nicht finden kann. Also halte Ausschau nach Fehlern, die du findest.

Und das war's im Grunde. Wir haben unseren Prozess immer noch hier. Es sind keine Änderungen hier unten nötig. Der Prozess hat den gleichen Namen und wird auf genau die gleiche Weise ausgeführt. Es ist nur der eigentliche Code des Prozesses, der jetzt in einer separaten Datei ist.

Wir können den Nextflow-Workflow wieder ausführen, er wird genau auf die gleiche Weise funktionieren. Und das ist im Grunde der Rest dieses Kapitels des Kurses - einfach diese drei Prozesse in ihre eigenen Dateien zu verschieben.

Also lass uns das jetzt tun. Ich werde schnell eine neue Moduldatei für den zweiten Prozess erstellen: convertToUpper.nf. Ich werde diesen Code ausschneiden, ihn hier einfügen. Und dann werde ich diesen importieren. Lass uns nur, großartig.

Und dann werde ich eine neue Datei für collectGreetings.nf erstellen. Schneide das aus.

Viel Schneiden, Schneiden und Kopieren und Einfügen.

Und jetzt sieht unser Haupt-Workflow-Skript plötzlich viel, viel kürzer aus, viel zugänglicher und viel einfacher zu lesen.

Und du kannst sehen, wie das Projekt jetzt mit unseren verschiedenen Dateien aufgebaut wird. Wir können ins Detail in den Bereichen eintauchen, die wir wollen. Uns viel einfacher zurechtzufinden, um bestimmte Schritte in der Pipeline zu finden, und schnell einen Überblick darüber zu bekommen, was die Pipeline tut.

## Navigation in Modulen mit VS Code

Natürlich ist der Nachteil dabei, dass du, wenn du eine große Pipeline hast, viele Moduldateien haben wirst und sie in mehreren Unterverzeichnissen oder allen möglichen Dingen organisiert sein könnten. Nun, wieder, ein kleiner Tipp hier. Die VS Code Extension ist ziemlich gut darin, deine Codebasis für dich zu navigieren und dir auch über den Code dort zu erzählen.

Du kannst sehen, dass VS Code versteht, was dieser Prozess ist, und mir einen kleinen Überblick gibt, wenn ich darüber fahre, sodass ich sehen kann, ohne den Quellcode suchen und finden zu müssen, was die Eingaben und die Ausgaben sind, was normalerweise das Wichtigste ist, wenn ich ihn in einem Workflow verwende.

Und auch wenn ich Command halte, ich bin auf einem Mac, und ich auf den Prozessnamen klicke, öffnet es die Datei direkt sofort. Zieht sie rein. Also kann ich einfach direkt dorthin springen, ohne auch nur an die tatsächlichen Dateipfade zu denken. Und das funktioniert überall, ich kann das auch machen, wo Prozesse aufgerufen werden. Das macht es also wirklich schnell.

## 4.4. Führe den Workflow aus

Okay, lass uns einfach überprüfen, dass die Pipeline immer noch so läuft, wie wir es erwarten. Also rufe das Terminal auf. Lass uns "nextflow run hello modules" machen, und sehen, ob es ohne Probleme ausgeführt wird.

Hoffentlich ist der ganze Sinn davon, dass die Pipeline im Grunde unverändert ist, also solltest du eigentlich keine Änderungen im Vergleich zu vorher sehen, als wir sie ausgeführt haben. Die Ausgabe hier sieht genau gleich aus, und du kannst unser results-Verzeichnis mit all den gleichen Dateien sehen, also ist das großartig. Keine Änderung ist gut.

## Ein Hinweis zu nf-core/modules

Bevor wir abschließen, möchte ich kurz auf die Kraft der Zusammenarbeit eingehen, wenn es um Module geht. Diese Dateien befinden sich in meinem Repository, also ist nicht sofort offensichtlich, wie wir bei ihnen zusammenarbeiten könnten. Und es gibt viele verschiedene Möglichkeiten, wie du das tun kannst, aber wahrscheinlich das größte und bekannteste Beispiel dafür ist nf-core.

Wenn ich zur nf-core-Website gehe, gehe ich zu resources, und modules. Du kannst sehen, dass nf-core eine riesige Bibliothek von Modulen hat, fast knapp unter 1700 Module, wenn ich das sehe. Und so kann ich den Namen meiner Lieblingswerkzeuge eingeben, gehen und herausfinden, ob jemand anderes bereits ein Modul dafür geschrieben hat, und diesen vorgeschriebenen Modulprozess hier sehen mit allen Eingaben, den Ausgaben, den Software-Containern, all diesen Informationen, und du kannst hier auf der Seite sehen, wie viele verschiedene nf-core-Pipelines alle diesen einzelnen geteilten Prozess verwenden.

Das ist ein etwas extremes Beispiel, aber du kannst sehen, dass dieser Code wirklich wiederverwendet wird. Und wenn ich zur GitHub-Quelle dafür durchklicke, ist es genau dasselbe wie das, was wir tun. Es ist einfach ein großer Prozess in einer Datei.

Nun, auf der nf-core-Seite machen wir einige Tricks, um diese Dateien teilen zu können und sie in verschiedene Repositories zu bringen. Und wenn du mehr darüber wissen möchtest, schau dir den Kurs an, den wir speziell über die Verwendung und das Erstellen mit nf-core haben. Aber ich wollte dir nur eine Vorstellung davon geben, wie mächtig dieses Konzept der Code-Wiederverwendung sein kann.

## Zusammenfassung

Richtig, das war's für Module. Ich habe dir gesagt, es war ein kurzer Abschnitt des Kurses. Schau dir das Quiz an, stelle sicher, dass du es verstehst und stelle sicher, dass alles noch richtig funktioniert. Und ich sehe dich im nächsten Video wieder, in dem es um Software-Container geht. Vielen Dank.

I.
