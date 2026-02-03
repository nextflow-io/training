# Teil 1: Hello World - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../01_hello_world.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern des Materials.

## Willkommen

Hi, willkommen zu Kapitel Eins von Hello Nextflow.

In diesem ersten Teil eines sechsteiligen Kurses gehen wir auf die grundlegenden Basics von Nextflow ein. Wir beginnen damit, einige Befehle in einem Terminal auszuführen, und dann schauen wir uns an, wie wir diese Bash-Befehle in ein Nextflow-Skript einbauen können.

Wir werden versuchen, diese erste Nextflow-Pipeline auszuführen, sehen, was Nextflow macht, wo es läuft, welche Dateien es erstellt und was der Zweck dieser Dateien ist.

Also gut, legen wir los.

## training.nextflow.io

Als Erstes, geh zu training.nextflow.io. Genau wie vorher ist das gesamte Material hier verfasst, und ich werde es Schritt für Schritt durcharbeiten. Ich zeige meinen Bildschirm, während ich die Schritte des Trainings durchführe, aber alles, was ich sage, steht im Trainingsmaterial, sodass du es in deinem eigenen Tempo verfolgen kannst und dort alles geschrieben findest.

Dieses Video hat auch Videotranskript aktiviert, also fühl dich frei, diese einzuschalten und genau zu verfolgen, was ich sage, während ich es sage.

Okay, gehen wir zu Hello Nextflow. Das ist der Kurs, den wir heute machen werden, und wir haben die Orientierung bereits im ersten Video gemacht, also gehen wir direkt zu Teil eins. Hello World.

Okay, ich verlasse jetzt dieses Trainingsmaterial und wechsle in meine Code Spaces-Umgebung. Das ist das, was wir im ersten Video eingerichtet haben. Hoffentlich hast du etwas, das sehr ähnlich aussieht in deinem eigenen System. Ich verwende VS Code und schaue mir das Trainingsmaterial an, und ich habe das Verzeichnis in das hello Nextflow-Verzeichnis gewechselt.

## 0. Aufwärmen: Hello World direkt ausführen

Okay. Beginnen wir mit ein paar Basics, die hoffentlich jedem vertraut vorkommen. Ich fange einfach an, indem ich einen sehr einfachen Befehl im Terminal schreibe. Hier unten werde ich sagen 'echo Hello World!"' Enter drücken und, keine Überraschungen, das Terminal macht, was ich verlange, und gibt diesen String zurück. Hello World.

Okay, dann drücke ich nach oben, um diesen Befehl zu holen und ihn ein bisschen mehr zu bearbeiten. Lass uns diesmal diese Ausgabe in eine Datei umleiten. Ich werde es stattdessen in output.txt schreiben und Enter drücken, nichts im Terminal dieses Mal, weil die Ausgabe nicht ins Terminal kam. Sie ging in diese Datei.

Ich kann dann diese Datei lesen, indem ich 'cat output.txt' mache, drücke dort Tab, um den Dateinamen automatisch zu erweitern, und da ist es. Die Datei ist da.

Ich kann diese Datei auch in der Seitenleiste im Datei-Explorer in VS Code sehen. Ich kann sie doppelklicken und hier öffnen. Wenn du sie in VS Code öffnen möchtest, ohne etwas anzuklicken, kannst du auch "code" und dann "output.txt" machen, und es macht das Gleiche.

Großartig. Das ist der erste Schritt. Sehr einfach.

## 1. Das Hello World Workflow-Starterskript untersuchen

Okay. Wir werden jetzt genau das Gleiche machen, aber in Nextflow, anstatt direkt im Terminal.

Wir werden das erste Beispielskript verwenden, um zu beginnen, diese Datei heißt Hello World. Ich kann "ls" machen, um es in einem Terminal anzuzeigen, und ich bin auf Mac, also kann ich Befehl-Klick machen, um diese Datei zu öffnen, oder ich hätte einfach in der Seitenleiste hier doppelklicken können.

Es gibt ein paar Dinge, die wir in dieser Datei sehen können. Ganz oben gibt es eine Hash-Anweisung, die sagt, dass dies eine Nextflow-Datei ist und wie sie ausgeführt werden könnte. Es gibt hier einige Kommentare, nur reguläre Code-Kommentare in hellgrau, die die Ausführung nicht beeinflussen und uns nur helfen, das Skript zu lesen.

Und dann gibt es zwei Hauptstrukturen. Es gibt hier einen Process und einen Workflow.

Processes in Nextflow sind die Schritte der Pipeline. Sie sind die Teile, die tatsächlich die Logik ausführen und die Verarbeitung durchführen.

Der Workflow dann am Ende fügt diese Processes zusammen und steuert die Logik des Workflows, wie alles miteinander verbunden ist.

Wir beginnen damit, uns einen Process anzusehen. Wir kommen gleich zum Workflow zurück.

## 1.2 Die Process-Definition

Also beginnt jeder Process mit einem Schlüsselwort process. Hat einen Namen und dann hat es einige geschweifte Klammern, und alles innerhalb dieser geschweiften Klammern ist dieser einzelne Process.

Ein Process muss einen Script-Abschnitt haben, und hier enthalten ist ein Bash-Snippet in einem mehrzeiligen String, der der Teil des Codes ist, der tatsächlich in der Rechenumgebung ausgeführt wird.

Wir haben hier auch eine Output-Anweisung, die Nextflow mitteilt, welche Dateien vom Script erstellt werden sollen. Beachte, dass die Ausgabe hier ein Schlüsselwort path hat, das Nextflow mitteilt, dass dies eine Datei ist, kein Wert oder ein String.

Innerhalb des Script-Blocks ist dies nur eine reguläre Bash-Anweisung, und es ist genau das Gleiche wie das, was wir im Terminal geschrieben haben. Wir echoieren Hello World in eine Datei namens output.txt. Dieses output.txt wird dann von der Output-Definition aufgenommen. Die Output-Definition macht eigentlich nichts. Sie teilt Nextflow nur mit, was zu erwarten ist, und wenn diese Datei nicht erstellt würde, würde Nextflow einen Fehler werfen.

Beachte, dass dieses Beispiel kein großartiges ist, weil wir den Dateinamen hier hartcodiert haben, output.txt und output.txt. Wenn eines davon geändert würde, würde das einen Fehler in unserem Workflow verursachen.

Es gibt einen besseren Weg, dies mit Variablen zu tun, was wir gleich behandeln werden.

## 1.3 Die Workflow-Definition

Okay. Wenn wir nach unten zum Workflow gehen, können wir sehen, dass wir einen Kommentar haben und dann führen wir den Process namens sayHello aus. Dies ist das gleiche Schlüsselwort, das hier oben ist. Dies ist so einfach, wie ein Workflow sein kann. Wir rufen nur einen einzelnen Process ohne variable Eingabe auf, also verbinden wir ihn nicht mit etwas anderem. Im späteren Teil dieses Kurses werden wir darüber sprechen, wie man dies leistungsfähiger macht, indem man variable Eingaben verwendet und Dinge mit Channels verbindet.

## 2. Den Workflow ausführen

Okay, das ist alles, was wir brauchen. Lass uns sehen, ob wir es ausführen können und sehen, was passiert. Ich werde einfach das Terminal löschen und dann werde ich "nextflow run" machen, und ich werde den Dateinamen aufrufen, der hello-world.nf ist. Das ist alles, was wir brauchen, um eine Nextflow-Pipeline auszuführen. Diese Pipeline nimmt keine Eingabe, also brauchen wir keine anderen Argumente.

Lass uns Enter drücken und sehen, was passiert.

Okay. Hoffentlich solltest du eine Ausgabe haben, die so aussieht. Wir haben ein paar Informationen, die uns sagen, dass Nextflow gelaufen ist und welche Version es verwendet hat. Sagt uns, welches Skript gestartet wurde, und es gibt uns einen zufällig generierten Namen für diese bestimmte Workflow-Ausführung. In diesem Fall wurde meiner "gloomy_crick" genannt.

Der wichtigste Teil hiervon ist jedoch, dass es uns sagt, welche Schritte in der Pipeline gelaufen sind. Du kannst sehen, dass unser Process namens sayHello gelaufen ist, und er ist einmal gelaufen und war zu hundert Prozent vollständig.

Dieser Teil hier ist der Hash für diese bestimmte Workflow-Aufgabe. Jeder Process läuft ein oder mehrere Male, und jede dieser Ausführungen wird eine Aufgabe genannt.

## 2.2. Die Ausgabe und Logs im work-Verzeichnis finden

Jede Aufgabe bekommt ihr eigenes isoliertes Verzeichnis, in dem sie läuft, also ist sie vom Rest der Ausführung des Workflows getrennt. Dieser Hash entspricht der Dateistruktur innerhalb des work-Verzeichnisses. Wenn ich "tree work" mache, können wir a0 sehen, und dann eine längere Version eines kurzen Hashs, und dann unsere output.txt-Datei. Du kannst es auch in einer Seitenleiste sehen.

Du kannst in der Seitenleiste sehen, dass es hier einige zusätzliche Dateien gibt. Der Grund, warum diese nicht im Terminal aufgetaucht sind, ist, dass sie versteckte Dateien sind, sie beginnen mit einem Punkt. Und tatsächlich, wenn ich "tree -a" für alle mache, und "work", können wir sie hier sehen.

Diese Dot-Dateien sind in jedem einzelnen work-Verzeichnis vorhanden, das Nextflow erstellt, und jede hat eine etwas andere Aufgabe. Erstens enthält .command.begin nur einige Anweisungen für Nextflow, die die Aufgabe einrichten, bevor sie läuft. .command.run sind die tatsächlichen Anweisungen, die von Nextflow selbst ausgeführt werden. Dann ist .command.sh wahrscheinlich die interessanteste. Dies ist das Skript, das von unserem Process-Block-Script aufgelöst wurde.

Wenn ich es öffne, kannst du sehen, wir haben unser "echo Hello World" in die output.txt-Datei. Dies ist genau das Gleiche wie unser Process in diesem Fall, aber wenn wir irgendwelche Variablen in unserem Nextflow-Code haben, wird jede Aufgabe eine andere .command.sh haben, und du kannst sehen, wie diese Variablen aufgelöst wurden.

Die anderen Dateien haben damit zu tun, wie die Aufgabe ausgeführt wurde. Also sind .command.err, .log und .out die Standardfehler, Standardausgabe und die beiden kombiniert. Und .exitcode teilt Nextflow mit, wie diese Aufgabe mit welchem Exit-Code ausgeführt wurde, ob sie erfolgreich war oder nicht.

Schließlich haben wir unsere output.txt-Datei und tatsächlich, "Hello World", das ist, was wir erwarten, und das ist, was erstellt wurde.

Okay, großartig. Das war dein allererster Nextflow-Lauf. Herzlichen Glückwunsch. Es ist wirklich so einfach.

Als Nächstes werden wir darauf eingehen, wie man dies etwas bequemer macht, sodass wir den Code nicht jedes Mal bearbeiten müssen, wenn wir eine Änderung an der Ausführung der Pipeline vornehmen möchten.

## 3. Workflow-Ausführungen verwalten

Diese Verzeichnisstruktur ist großartig, um alle Aufgaben getrennt und alles organisiert zu halten, aber natürlich ist es nicht sehr bequem, deine Ausgabedateien zu finden. Du möchtest nicht durch Mengen verschachtelter Verzeichnisse graben, um die Ergebnisse deiner Pipeline zu finden.

## 3.1. Ausgaben veröffentlichen

Die gute Nachricht ist, du sollst das nicht. Die work-Verzeichnisse sind wirklich nur für Nextflow, um sie selbst zu verwenden. Also werden wir eine Funktion für Nextflow namens "publishDir" verwenden.

Wir gehen zurück zu unserem Workflow, gehen zum Process. Wir können hier eine neue Anweisung hinzufügen, die eine Direktive genannt wird. Dies ist, was Nextflow diese Dinge am Anfang von Processes nennt, die die Funktionsweise erweitern, und die, die wir verwenden werden, heißt publishDir.

Du kannst sehen, ich habe hier angefangen zu tippen, und die Nextflow-Erweiterung für VS Code hat mir die Direktive vorgeschlagen, also kann ich einfach Enter drücken.

Okay. Ich werde dies mit einem Verzeichnis namens "results" verfolgen, und wir werden ihm sagen, dass es die Ausgabedateien dorthin kopieren soll. Also werde ich sagen mode copy. Großartig. Werde speichern und lass uns den Workflow erneut ausführen.

nextflow run hello-world.nf

Es läuft genau gleich. Beachte aber, wir haben dieses Mal einen etwas anderen Hash. Nextflow wird jedes Mal, wenn du den Workflow ausführst, einen anderen Hash verwenden. Und wir haben infolgedessen einen anderen Satz von work-Verzeichnissen. Bereiche, einer heißt EB stattdessen, aber du kannst sehen, alle Dateien sind gleich. Was jedoch neu ist dieses Mal, ist, dass wir auch ein Verzeichnis namens "results" haben.

Innerhalb von "results" hier haben wir unsere Ausgabedatei. Das ist, was wir Nextflow gesagt haben zu tun. Wir sagten, speichere die Ergebnisdateien in einem Verzeichnis namens "results" und kopiere sie dorthin. Und so ist dies jetzt viel einfacher zu finden. Es ist einfach dort neben dem Ort, wo wir einen Workflow gestartet haben, und all die verschiedenen Dateien können dort organisiert werden, wie auch immer wir wünschen, unabhängig davon, wo oder wie Nextflow die tatsächliche Ausführung durchgeführt hat.

Beachte, dass publishDir Symlinks handhaben kann, was gut ist, wenn du auf einem gemeinsam genutzten Dateisystem arbeitest und Platz sparen möchtest. Und auch musst du nicht alle Dateien definieren, die von einem Process erstellt werden, als Ausgabe.

Nextflow wird nur die Dinge kopieren, die in diesem Output-Block definiert sind. Also wenn du Zwischendateien hast, die vom Schritt erstellt werden, die nicht downstream von diesem Process benötigt werden, definierst du sie einfach nicht in der Ausgabe, und sie werden nicht in publishDir auftauchen. Also ist dies ein Weg, deine Ausgabedateien von einer Pipeline sauber zu halten und Zwischendateien einfach zu löschen, sobald der Arbeitsplatz fertig ist.

Eine kurze Bemerkung hier. Es gibt neue Nextflow-Syntax, die kommt, genannt workflow output definitions, die schließlich publishDir ersetzen wird. Dies gibt uns eine Möglichkeit, alle Ausgaben von einem Workflow auf Pipeline-Ebene unten im Workflow-Block zu definieren. Dies ist in den Nextflow-Docs beschrieben, wenn du es ausprobieren möchtest. Aber für jetzt wird publishDir noch eine Weile da sein, also haben wir das noch in einem Training für 2025.

## 3.2. Einen Workflow mit -resume neu starten

Okay. Ich erwähnte, dass das work-Verzeichnis hier jetzt zwei Sätze von Ergebnissen mit einem anderen Hash von jeder Ausführung des Workflows hat. Das ist gut. Aber manchmal wollen wir Schritte nicht jedes Mal neu berechnen, wenn wir es nicht brauchen.

Vielleicht baust du deinen Workflow iterativ auf und du fügst Schritte hinzu, und du möchtest, dass die ersten Schritte einfach die gecachten Versionen wiederverwenden. Oder vielleicht ist etwas auf deinem Rechensystem in der Mitte deines Workflows schief gegangen, und du möchtest, dass es dort weitermacht, wo es aufgehört hat, aber die Schritte überspringt, die es bereits abgeschlossen hatte.

Nextflow hat eingebaute Funktionalität dafür namens resume. Lass es uns ausprobieren. Also zuerst einmal werde ich mir einfach das work-Verzeichnis ansehen, damit wir uns erinnern können, was dort war.

Und dann werde ich "nextflow run hello-world.nf" machen, und ich werde hier einen einzigen Befehl hinzufügen, "-resume".

Beachte, einzelner Bindestrich, das ist wirklich wichtig. Ich werde es ausführen, und die Ausgabe wird im Grunde genau gleich aussehen, mit ein paar kleinen Unterschieden.

Beachte hier, es steht "cached" in grau. Das bedeutet, dass Nextflow die Aufgabe nicht ausgeführt hat. Dieses Mal fand es etwas, das dem entsprach, was Anforderungen waren, und es verwendete diese Ausgaben direkt wieder, anstatt den Schritt erneut auszuführen.

Und tatsächlich, wenn du dir den Hash hier anschaust, kannst du sehen, dies entspricht dem vorhandenen Hash, den wir von einem vorherigen Lauf hatten.

## 3.3. Ältere work-Verzeichnisse löschen

Okay. Aber wenn du iterativ entwickelst, wirst du eine Menge dieser Workflow-Dateien aufbauen. Das kann ein Problem sein, wenn du möglicherweise wenig Platz hast.

Nextflow kann uns helfen, diese work-Verzeichnisse mit ein paar Hilfsbefehlen aufzuräumen. Wenn ich "nextflow log" mache, das wird mir eine Liste aller verschiedenen Workflow-Läufe geben, die ich in diesem Verzeichnis gemacht habe, und sie haben die Laufnamen hier. Du kannst den gloomy quick sehen, der der erste war, den wir ausgeführt haben, und dann diese beiden neuen.

Wir können jetzt diesen Namen nehmen und diese mit dem "nextflow clean"-Befehl verwenden. Ich kann einen einzelnen Laufnamen angeben. Oder noch besser, ich kann Nextflow sagen, alles vor einem einzelnen Workflow-Namen zu löschen mit "-before", und ich werde "stupefied_shaw" eingeben. Das war mein letzter Lauf, "-n".

Der "-n"-Befehl sagte Nextflow, es als Probelauf zu machen, ohne tatsächlich etwas wirklich zu löschen, und es sagt uns, welche der Hash-Verzeichnisse es entfernt hätte. Tatsächlich ist es nur das eine von der ersten Ausführung. Beide der zweiten Ausführungen verwenden das gleiche Hash-Verzeichnis.

Ich werde es noch einmal ausführen, aber jetzt anstelle von "-n" für Probelauf, werde ich "-f" für force machen, und es hat dieses Hash-Verzeichnis entfernt. Jetzt, wenn ich "tree work" mache, können wir sehen, wir haben nur noch diese Ausgabedatei übrig.

Großartig. Also haben wir es geschafft, dort eine ganze Menge Festplattenspeicher aufzuräumen.

Ein paar Dinge zu beachten beim Löschen von work-Verzeichnissen, wenn du Sachen in dein results-Verzeichnis symlinkt, werden diese Symlink-Quellen jetzt gelöscht, und deine Ergebnisse werden für immer weg sein. Deshalb ist die Verwendung des copy-Modus eine sicherere Sache zu tun, und generell das, was wir empfehlen.

Zweitens hängt die Resume-Funktionalität von Nextflow von diesen work-Verzeichnissen ab. Also wenn du sie löschst und du Nextflow wieder ausführst, wird die Resume-Funktionalität nicht mehr funktionieren. Also liegt es an dir, im Auge zu behalten, welche Dinge du möglicherweise brauchst oder nicht brauchst, und lösche Dinge nur, wenn du sicher bist, dass es sicher ist, dies zu tun.

Die andere Sache, die wir tun können, ist, wir können einfach das gesamte work-Verzeichnis löschen, wenn wir unseren Workflow-Lauf beendet haben und wir sicher sind, dass wir es nicht mehr brauchen.

Also kann ich "rm -r work" machen. Ich weiß, es war nichts Wichtiges darin. Ich habe meine Ergebnisse, die mir wichtig sind, im results-Verzeichnis, wo wir sie kopiert haben. Und so war es sicher, das work-Verzeichnis zu löschen. Es liegt an dir, welchen dieser Ansätze du verwendest.

## 4. Eine variable Eingabe verwenden, die über die Befehlszeile übergeben wird

Okay, was kommt als Nächstes? Ich erwähnte, dass wir einige der Werte in unserem Workflow-Skript hier hartcodiert hatten, die output.txt-Datei, und dass es einen besseren Weg geben könnte, dies zu tun.

Lass uns damit beginnen. Was wir tun werden, sind drei Dinge. Wir werden eine neue Eingabe zum Process hinzufügen. Wir werden dem Process-Script sagen, wie es diese Eingabe verwenden soll, und dann werden wir es im Workflow verdrahten, sodass wir es dynamisch mit einem Befehlszeilen-Flag beim Ausführen von Nextflow verwenden können.

Also zuerst die wichtigsten Dinge. Lass uns hier einen Input-Block hinzufügen. Genau wie bei der Ausgabe. Dies ist ein neuer Abschnitt für den Process, und ich werde sagen, "val greeting".

Beachte hier, ich sage "val", was sagt, dass dies eine Variable ist, kein Pfad.

Ich kann dann nach unten ins Script gehen und dann kann ich diesen hartcodierten Text hier herausnehmen und $greeting machen. Dies funktioniert genau wie jede andere Programmiersprache. Wir definieren hier eine Variable, und wir referenzieren sie innerhalb dieses Script-Blocks. Wenn Nextflow diesen Process ausführt, wird die Variable interpoliert. Und wenn wir uns diese .command.sh-Datei ansehen, werden wir den tatsächlichen hartcodierten String hier stattdessen sehen.

## 4.1.3. Einen CLI-Parameter einrichten und ihn als Eingabe für den Process-Aufruf bereitstellen

Okay, aber wo stellen wir die Variable bereit? Als Nächstes gehen wir zum Workflow-Abschnitt runter, und du kannst sehen, dass die Erweiterung hier sagt, wir erwarten jetzt eine Eingabe, und es hat mir eine Warnung gegeben.

Nun, das Einfachste, was wir tun könnten, ist einfach hartcodieren. Ich könnte "Hello World" schreiben und diesen String-Input für den Process bereitstellen. Aber wieder würde das wirklich keine Probleme lösen. Wir müssten immer noch zurückgehen und den Pipeline-Code jedes Mal bearbeiten, wenn wir etwas ändern wollten, was nicht gut ist.

Die gute Nachricht ist, dass Nextflow ein eingebautes System hat, um Befehlszeilenargumente zu handhaben, genannt Parameter. Also kann ich stattdessen eine dieser speziellen Variablen namens params verwenden, und ich kann sie nennen, wie ich will, aber ich werde greeting sagen, damit es mit der Workflow-Logik übereinstimmt.

Speichern drücken und lass uns sehen, was wir damit machen können.

Also wenn ich zurück zum Terminal gehe. Also machen wir "nextflow run hello-world.nf". Genau wie vorher, aber der entscheidende Unterschied ist, wir machen --greeting

Beachte, es gibt hier zwei Bindestriche, weil dies ein Parameter ist. Als wir den Workflow vorher resumed haben, war das ein einzelner Bindestrich. Das liegt daran, dass resume eine Kern-Nextflow-Option ist, und dies ist ein Parameter, der spezifisch für unsere Pipeline ist.

Verwechsele die beiden nicht. Es ist einfach, das zu tun. Wenn du --resume statt nur einem Bindestrich gemacht hättest, dann wäre das "params.resume", was nichts machen würde. Ebenso, wenn du einen einzelnen Bindestrich hier gemacht hättest, würde Nextflow es nicht als Schlüsselargument erkennen.

Also ist es --greeting, was params.greeting entspricht.

Ich kann diesem jetzt folgen mit welchem Text auch immer ich will. Also bin ich gerade in Schweden, also werde ich sagen, "Hej världen".

Also lass uns es ausführen, sehen, was passiert, Moment der Wahrheit.

Okay, also kannst du sehen, dass der Process erneut gelaufen ist, genau wie vorher, sayHello mit einer einzelnen Ausführung.

Dies wird die Datei überschrieben haben, die im publishDir "results"-Verzeichnis war. Also sei vorsichtig, wenn du die Dateien erneut ausführst, weil Dinge im published air überschrieben werden.

Ich kann jetzt "code results/output.txt" machen, und tatsächlich wurde unsere Ausgabe aktualisiert und sagt jetzt "Hej världen".

## 4.2. Standardwerte für Befehlszeilenparameter verwenden

Okay, das ist großartig. Aber das Problem jetzt ist, dass unser Workflow darauf angewiesen ist, dass wir diesen Parameter immer definieren, und es ist schön, vernünftige Standards zu haben, damit die Dinge auf vernünftige Weise für deinen Workflow laufen, es sei denn, du überschreibst die Standards.

Also der Weg, wie wir das machen, ist, indem wir einen Standardwert für den Parameter in unserem Workflow-Skript setzen.

Also wenn ich zurück zu meiner hello-world.nf-Datei gehe, kann ich ins Skript gehen, direkt über workflow, "params.greeting" tippen und es wie jede andere Variable definieren. Also lass uns hier einen String eingeben und sagen wir "Holà mundo!"

Jetzt hat dieser Parameter einen Standard definiert, der hier verwendet wird, oder wir können ihn immer noch auf der Befehlszeile mit --greeting überschreiben, genau wie wir es vorher gemacht haben.

Also lass uns überprüfen, ob es funktioniert. "nextflow run hello-world.nf"

Dieses Mal keine Befehlszeilenargumente, und überprüfen, ob es das Richtige gemacht hat.

"code results/output.txt". Und da ist es. Wir haben unseren Standard bekommen.

Okay, lass es uns noch einmal versuchen, nur überprüfen, dass ich dir keine Lügen erzähle. Lass uns es noch einmal ausführen, aber --greeting machen, und das Beispiel aus einem Trainingsmaterial verwenden, lass uns sagen "Konnichiwa!"

Führt den Workflow erneut aus, und tatsächlich wurde unsere Ausgabedatei oben gerade mit dem neuen Wert aktualisiert, den wir auf der Befehlszeile bereitgestellt haben.

Großartig. Dies ist ein wirklich zentraler Aspekt beim Schreiben eines jeden Nextflow-Workflows. Vernünftige Standards in deinem Pipeline-Code definieren, aber es sehr einfach machen, für den Endbenutzer zu konfigurieren, indem man Befehlszeilenargumente im Terminal hat.

Beachte, dass der Endbenutzer die Konfiguration an mehreren verschiedenen Stellen überschreiben kann. Du kannst eine Config-Datei in deinem Home-Verzeichnis haben, die auf jeden einzelnen Nextflow-Lauf angewendet wird, den du machst. Du kannst eine Config-Datei in einem Launch-Verzeichnis haben. Du kannst eine Config-Datei in einem Pipeline-Verzeichnis haben. All diese verschiedenen Config-Orte werden in einer bestimmten Reihenfolge geladen, die in den Nextflow-Docs beschrieben ist.

Okay, das ist das Ende von Abschnitt eins. Wir hatten unser allererstes Workflow-Skript in Nextflow mit einem Process und einem Workflow. Wir haben uns Eingaben, Ausgaben, Scripts und Publishing angeschaut und wie man Parameter und einen Input-Channel in unseren Process verdrahtet.

Herzlichen Glückwunsch, dein erster Schritt zum Schreiben von Nextflow-Code ist abgeschlossen.

Mach eine kleine Pause und ich sehe dich in ein paar Minuten zurück für Kapitel zwei.

[Nächstes Videotranskript :octicons-arrow-right-24:](02_hello_channels.md)
