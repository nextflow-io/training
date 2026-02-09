# Teil 4: Hello Modules - Videotranskript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../04_hello_modules.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern im Material.

## Willkommen

Hallo und willkommen zurück zu Teil vier von Hello Nextflow. Dieser Abschnitt dreht sich ganz um Module, und es ist ein recht kurzer Abschnitt des Kurses. Wir werden nicht wirklich viel Code schreiben, es geht mehr darum, wie wir den Code in unserer Pipeline organisieren.

Bis jetzt haben wir alles in eine einzige Datei gepackt, was in Ordnung ist, und so haben wir früher tatsächlich Nextflow-Pipelines gebaut.

Aber wenn die Pipeline wächst, wird das Skript länger und länger und länger und schwieriger zu navigieren, zu warten, und es bedeutet auch, dass wir den Code nicht wirklich teilen können.

Nextflow-Module ermöglichen es uns, Prozesse aus dem Hauptskript herauszuziehen und sie dann zu importieren. Das bedeutet, dass der Code einfacher zu navigieren ist und es bedeutet auch, dass wir diesen Modulcode zwischen verschiedenen Pipelines teilen können.

Dieses kleine Diagramm auf der Hauptseite der Dokumentation zeigt das Konzept gut. Anstelle eines riesigen Skripts werden wir diese separaten Moduldateien aus verschiedenen Modulskripten einbinden, und alles wird in den Workflow eingefügt, aber es wird immer noch auf genau die gleiche Weise ausgeführt.

Also lass uns in GitHub Codespaces springen und uns ein bisschen umsehen. Wie zuvor habe ich meinen Arbeitsbereich hier ein wenig aufgeräumt. Die alten Nextflow-Verzeichnisse und das Work-Verzeichnis und so weiter entfernt. Aber es macht nichts, wenn du diese Dateien noch hast.

Ich werde anfangen, in der hello modules-Datei zu arbeiten, die im Grunde dort ist, wo wir sie am Ende eines vorherigen Kapitels verlassen haben. Wir haben unsere drei Prozesse hier drin. Wir haben ein paar Parameter, den Workflow-Block, wo wir diese drei Prozesse ausführen und sie mit Kanälen zusammenfügen. Dann veröffentlichen wir die Ausgabekanäle und wir haben den Output-Block, der angibt, wie diese Dateien veröffentlicht werden sollen.

## 1. Ein Verzeichnis zum Speichern von Modulen erstellen

Nun, wie gesagt, wir werden nicht wirklich viel Code schreiben oder bearbeiten. Wir werden nur den Code verschieben, den wir bereits haben. Nextflow-Moduldateien haben normalerweise einen einzelnen Prozess darin, und per Konvention bewahren wir sie normalerweise in einem Verzeichnis namens modules auf. Aber du kannst das nennen, wie du willst. Aber ich werde ein modules-Verzeichnis in meinem Repository hier behalten, und dann werde ich eine Datei für jeden Prozess erstellen. Also sage ich neue Datei, sayHello.nf.

## 2. Ein Modul für sayHello() erstellen

Jetzt nehme ich meinen Prozess und ich werde einfach diesen Code auswählen, ihn aus der Haupt-hello-modules-Datei ausschneiden und hier einfügen.

Offensichtlich macht das allein noch nichts. Unser Hauptskript braucht immer noch diesen Prozess, also müssen wir ihn irgendwie wieder einbinden. Und das machen wir mit der include-Anweisung.

Also tippe ich include und geschweifte Klammern, und dann nehme ich den Namen des Prozesses. Und ich sage from, und dann gebe ich einen relativen Dateipfad an. Es heißt also, beginnt mit ./, weil es relativ zu dem Ort ist, wo dieses Skript gespeichert ist. Also ist es modules sayHello.nf.

Beachte, dass die VS Code-Erweiterung hier ziemlich hilfreich ist. Sie sagt uns, ob sie diese Datei finden kann und ob sie einen Prozess finden kann, den ich benenne. Wenn ich hier absichtlich einen Tippfehler mache, gibt sie mir sofort einen Fehler und sagt mir, dass sie diesen Prozess, den ich zu importieren versuche, nicht finden kann. Also behalte alle Fehler im Auge, die du findest.

Und das war's wirklich. Wir haben unseren Prozess immer noch hier. Es sind keine Änderungen hier unten nötig. Der Prozess hat den gleichen Namen und wird auf genau die gleiche Weise ausgeführt. Es ist nur so, dass der eigentliche Code des Prozesses jetzt in einer separaten Datei ist.

Wir können den Nextflow-Workflow erneut ausführen, er wird auf genau die gleiche Weise funktionieren. Und das ist im Grunde der Rest dieses Kapitels des Kurses, nur diese drei Prozesse in ihre eigenen Dateien zu verschieben.

Also lass uns das jetzt machen. Ich werde schnell eine neue Moduldatei für den zweiten Prozess erstellen: convertToUpper.nf. Ich werde diesen Code ausschneiden, ihn dort einfügen. Und dann werde ich diesen importieren. Lass uns einfach, großartig.

Und dann werde ich eine neue Datei für collectGreetings.nf erstellen. Das ausschneiden.

Viel Ausschneiden und Kopieren und Einfügen.

Und jetzt sieht unser Haupt-Workflow-Skript plötzlich viel, viel kürzer aus, viel zugänglicher und viel einfacher zu lesen.

Und du kannst sehen, wie sich das Projekt jetzt mit unseren verschiedenen Dateien aufbaut. Wir können an den Stellen, die wir wollen, ins Detail gehen. Uns viel einfacher zurechtfinden, um bestimmte Schritte in der Pipeline zu finden, und schnell einen Überblick darüber bekommen, was die Pipeline macht.

## Navigation in Modulen mit VS Code

Nun, natürlich ist der Nachteil dabei, dass du, wenn du eine große Pipeline hast, viele Moduldateien haben wirst und sie in mehreren Unterverzeichnissen oder allen möglichen Dingen organisiert sein könnten. Nun, wieder ein kleiner Tipp hier. Die VS Code-Erweiterung ist ziemlich gut darin, deine Codebasis für dich zu navigieren und dir auch über den Code dort zu berichten.

Du kannst sehen, dass VS Code versteht, was dieser Prozess ist, und mir einen kleinen Überblick darüber gibt, wenn ich darüber fahre, sodass ich sehen kann, ohne den Quellcode finden zu müssen, was die Eingaben und die Ausgaben sind, was normalerweise das Wichtigste ist, wenn ich ihn in einem Workflow verwende.

Und auch wenn ich Command halte, ich bin auf einem Mac, und ich auf den Prozessnamen klicke, öffnet es die Datei direkt sofort. Zieht sie rein. Also kann ich einfach direkt dorthin springen, ohne auch nur darüber nachzudenken, was die tatsächlichen Dateipfade sind. Und das funktioniert überall, ich kann das auch machen, wo Prozesse aufgerufen werden. Das macht es wirklich schnell.

## 4.4. Den Workflow ausführen

Okay, lass uns einfach überprüfen, dass die Pipeline immer noch so läuft, wie wir es erwarten. Also hole das Terminal hoch. Lass uns "nextflow run hello modules" machen und sehen, ob es ohne Probleme ausgeführt wird.

Hoffentlich ist der ganze Sinn davon, dass die Pipeline im Grunde unverändert ist, also solltest du wirklich keine Änderungen sehen im Vergleich zu vorher, als wir sie ausgeführt haben. Die Ausgabe hier sieht genau gleich aus, und du kannst unser results-Verzeichnis mit all den gleichen Dateien sehen, also ist das großartig. Keine Änderung ist gut.

## Ein Hinweis zu nf-core/modules

Bevor wir abschließen, möchte ich kurz auf die Kraft der Zusammenarbeit eingehen, wenn es um Module geht. Diese Dateien befinden sich in meinem Repository, also ist es nicht sofort offensichtlich, wie wir an ihnen zusammenarbeiten könnten. Und es gibt viele verschiedene Möglichkeiten, wie du das tun kannst, aber wahrscheinlich das größte und bekannteste Beispiel dafür ist nf-core.

Wenn ich auf die nf-core-Website gehe, gehe ich zu Resources und Modules. Du kannst sehen, dass nf-core eine riesige Bibliothek von Modulen hat, fast knapp unter 1700 Module, als ich das ansehe. Und so kann ich den Namen eines meiner Lieblingswerkzeuge eingeben, nachsehen, ob jemand anderes bereits ein Modul dafür geschrieben hat, und diesen vorgeschriebenen Modulprozess hier mit all den Eingaben, den Ausgaben, den Software-Containern, all diesen Informationen sehen, und du kannst auf der Seite hier sehen, wie viele verschiedene nf-core-Pipelines alle diesen einzelnen gemeinsamen Prozess verwenden.

Das ist ein ziemlich extremes Beispiel, aber du kannst sehen, dass dieser Code wirklich wiederverwendet wird. Und wenn ich zur GitHub-Quelle dafür durchklicke, ist es genau das Gleiche wie das, was wir machen. Es ist nur ein großer Prozess in einer Datei.

Nun, auf der nf-core-Seite machen wir einige Tricks, um diese Dateien teilen und in verschiedene Repositories einbringen zu können. Und wenn du mehr darüber wissen willst, schau dir den Kurs an, den wir speziell über die Verwendung und das Bauen mit nf-core haben. Aber ich wollte dir nur eine Vorstellung davon geben, wie mächtig dieses Konzept der Code-Wiederverwendung sein kann.

## Zusammenfassung

Richtig, das war's für Module. Ich habe dir gesagt, es ist ein kurzer Abschnitt des Kurses. Schau dir das Quiz an, stelle sicher, dass du es verstehst und stelle sicher, dass alles noch richtig funktioniert. Und ich sehe dich im nächsten Video wieder, das sich ganz um Software-Container dreht. Vielen Dank.
