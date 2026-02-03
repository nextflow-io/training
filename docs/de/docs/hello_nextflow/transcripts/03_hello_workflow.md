# Teil 3: Hello Workflow - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../03_hello_workflow.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen

Hallo, willkommen zu Teil drei des "Hello Nextflow" Trainingskurses.

Dieses Kapitel heißt "Hello Workflow".

In Kapitel zwei haben wir einen einfachen Workflow mit einem Prozess erstellt, aber in der Realität sind Pipelines nützlich, weil sie mehrere Analyseschritte miteinander verketten können.

In diesem Kapitel werden wir dieses anfängliche Beispiel nehmen und es etwas realistischer erweitern.

Wir werden einige zusätzliche Schritte hinzufügen und uns ansehen, wie wir Channels verwenden, um diese Schritte zu verbinden.

Wir werden uns mehrere Aufgaben ansehen, die in einen einzigen Prozess zusammengefasst werden können, und wir werden uns Prozesse ansehen, die mehrere Eingaben und mehrere Ausgaben haben können.

Okay, lass uns anfangen.

Also, fangen wir an. Wie zuvor. Lass uns zu training.nextflow.io gehen. Hello Nextflow, Kapitel drei. Hello Workflow. Und lass uns unseren Workspace öffnen. Ich habe alle meine Arbeitsdateien aus den vorherigen Kapiteln aufgeräumt und werde Hello Workflow öffnen.

Das ist jetzt dieselbe Datei, an der wir bisher gearbeitet haben, das sollte also vertraut aussehen. Wir haben unseren say hello Prozess. Wir haben unseren params.greeting mit seiner greetings CSV Datei, und wir haben unseren Workflow unten, der diese CSV Datei lädt, den Channel erstellt und ihn an unseren Prozess übergibt.

## 0. Aufwärmen: hello-workflow.nf ausführen

Wenn du möchtest, können wir das ausprobieren und überprüfen, ob es wie erwartet funktioniert. Öffne ein Terminal für nextflow run hello workflow nf und drücke Enter.

Okay, super. Unsere drei Prozesse laufen. Wir haben unser results Verzeichnis mit unseren drei Ausgaben. Bonjour. Hello. Holà. Also lass uns diese Dateien schließen, das Terminal schließen und zurück zum Skript gehen.

## 1. Einen zweiten Schritt zum Workflow hinzufügen

Okay. Für unser Beispiel bleiben wir einfach und versuchen, domänenunabhängig zu bleiben. Unser zweiter Prozess wird diese Zeichenketten, diese Wörter, einfach auf eine einfache Weise manipulieren. Wir werden den translate Unix-Befehl verwenden, um diese Dateien zu nehmen und sie alle in Großbuchstaben umzuwandeln. Das machen wir mit dem "tr" Befehl.

## 1.1. Den Großschreibungsbefehl definieren und im Terminal testen

Wir können das einfach im Bash-Terminal ausprobieren und sehen, ob es funktioniert. Du machst echo, Hello World, und übergibst das dann mit dem Pipe-Zeichen an tr, und wir geben ihm ein Erkennungsmuster, a bis z und wohin es übersetzt werden soll. A bis Z in Großbuchstaben.

Das ist sehr einfach, weil es buchstäblich die A bis Z Zeichen macht. Es wird also nicht bei allem funktionieren, was Akzente hat oder so. Aber für die Zwecke des Beispiels solltest du das Bild verstehen.

Werde Enter drücken und es druckt ins Terminal, HELLO WORLD in Großbuchstaben. Und genau wie zuvor könnten wir das zu einer Datei umleiten, wenn wir wollten. Outfile.

Okay. Lass uns das aufräumen.

## 1.1. Den Großschreibungsschritt als Nextflow Prozess schreiben

Lass uns zurück zu unserem Skript gehen und einen neuen Prozess schreiben, um diesen Bash-Befehl zu handhaben. Ich werde den vorherigen Prozess kopieren, ihn darunter einfügen und ihn convert to upper nennen. Für Großbuchstaben. Ich werde dasselbe publishDir results verwenden, aber ich werde hier ein paar Änderungen vornehmen. Anstatt ein val zu nehmen, werde ich ein path input file nehmen, und ich werde hier ein Präfix upper haben, damit unsere Ausgabedateien die Ausgabe nicht überschreiben. Und ich werde den Variablennamen aus der Eingabe verwenden. Und dann werde ich das Skript hier unten ändern, und stattdessen werde ich cat für die Eingabedatei verwenden und genau wie wir es in Bash TR gemacht haben, a-z, upper input file .txt. Okay, lass uns speichern.

## 1.2. Einen Aufruf des neuen Prozesses im Workflow-Block hinzufügen

Jetzt, wenn ich nach unten scrolle, müssen wir diesen Prozess tatsächlich aufrufen. Nur den Prozess zum Skript hinzuzufügen reicht nicht aus. Wir müssen Nextflow sagen, dass wir diesen Prozess ausführen müssen und wo das zu tun ist.

Also werde ich hier convert to upper schreiben und

okay, wir bekommen hier einen Fehler, der sagt, es erwartet ein Argument. Sicher, wir müssen etwas an diesen Prozess übergeben, damit er tatsächlich etwas zu tun hat.

## 1.3. Die Ausgabe des ersten Prozesses an den zweiten Prozess übergeben

Was wir tun werden, ist, dass wir die Ausgabe von diesem Prozess nehmen werden. Also nehme ich den Namen, say hello, und wenn ich dot out mache.

Für ein einfaches Beispiel wie dieses, wo wir einen Prozess haben, der nur eine Ausgabe hat, und wir das an einen neuen Prozess übergeben, der also eine Eingabe hat, sollte das alles sein, was wir brauchen. Also werde ich speichern, das Terminal aufrufen und lass uns das nochmal ausführen.

## 1.4. Den Workflow erneut ausführen

Jetzt habe ich mein work Verzeichnis vom letzten Mal, als ich diesen Workflow ausgeführt habe, nicht aufgeräumt. Ich werde es nochmal ausführen und ich werde das als Gelegenheit nutzen, um zu zeigen, wie partielles Caching funktioniert. Also wenn ich einzelnes Minus resume mache. Hoffentlich sollte es die Ausgaben von diesem ersten Prozess wiederverwenden, die genau dieselben waren wie beim letzten Mal, als ich ausgeführt habe. Aber jetzt haben wir hier einen neuen Prozess, der vorher nicht gelaufen ist, der von Grund auf neu läuft. Und tatsächlich kannst du sehen, dass der erste Prozess die Cache-Ausgaben verwendet hat und die zweite Ausgabe drei von drei ausgeführt hat. Du kannst auch sehen, dass wir jetzt beide unsere Prozesse hier haben, unser erster Prozess, say hello, lief dreimal, und unser zweiter Prozess convert to upper lief dreimal.

Wenn ich das nochmal ausführe, zur Erinnerung, mit -ansi-log false, sollten wir sehen, dass sechs verschiedene Prozessaufgaben laufen, drei für jeden von ihnen. Das macht also genau das, was wir erhofft haben. Der erste Prozess läuft dreimal, übergibt diese Ausgaben an einen zweiten Prozess, der dann dreimal läuft.

Also lass uns ins work Verzeichnis schauen und sehen, wie Nextflow diese Dateieingaben handhabt. Wenn ich dieses Hash-Verzeichnis hier vom zweiten Prozess nehme, können wir wieder den tree Befehl mit -a verwenden, nur um diese Dateien anzusehen. Du kannst hier drin sehen, dass wir unsere Eingabedatei haben, die die Bonjour-output.txt Datei ist, und das ist tatsächlich ein Symlink. Das ist es, was uns dieser Pfeil zeigt, und er zeigt auf die Datei im vorherigen work Verzeichnis.

Das macht Sinn. Nextflow handhabt die Ausführung jeder Aufgabe in ihrem eigenen gekapselten Verzeichnis, also ist es vollständig in sich geschlossen. Allerdings muss es die Dateien von einem vorherigen Schritt als Eingabe bereitstellen. Anstatt außerhalb des work Verzeichnisses zu greifen, um diese Dateien zu bekommen, stellt Nextflow sie ins work Verzeichnis bereit.

Wenn wir ein gemeinsames Dateisystem wie hier haben, macht es das mit einem Symlink, sodass es keinen zusätzlichen Dateiplatz verwendet. Wenn wir Cloud-Speicher mit Buckets an verschiedenen Orten verwenden, würde es diese Dateien abrufen und sie tatsächlich ins work Verzeichnis kopieren.

Lass uns die command sh Datei ansehen. Wenn ich code work, command sh mache, kannst du sehen, tatsächlich greift es auf diese Datei aus dem lokalen Verzeichnis zu. Also ist alles sehr in sich geschlossen und sauber.

Wir können auch das results Verzeichnis überprüfen und sicherstellen, dass diese Dateien ordnungsgemäß ausgegeben wurden. Und tatsächlich, in results, können wir alle Ausgabedateien vom ersten Prozess und alle Ausgabedateien vom zweiten sehen. Und sie sind alle in Großbuchstaben, wie wir gehofft hatten.

Hier beginnt die Kraft von Nextflow zu glänzen. Mit sehr minimalem Code hat Nextflow die parallele Ausführung dieser Aufgaben mit sauberer Kapselung in separaten work Verzeichnissen und dem Bereitstellen von Ein- und Ausgabedateien und Dateiveröffentlichung automatisch für uns direkt aus der Box gehandhabt. Du kannst also sehen, wie wertvoll diese Funktionalität wirklich ist, wenn wir die Komplexität unserer Analyse-Workflows skalieren.

## 2. Einen dritten Schritt hinzufügen, um alle Grüße zu sammeln

Okay. Diese Schritte waren eins-zu-eins. Wir hatten eine Ausgabe vom ersten Prozess, die zu einer Eingabe für den zweiten Prozess ging. Als nächstes werden wir darüber sprechen, wie man diese verschiedenen Ausgaben in einer einzigen Prozessaufgabe sammelt, was wiederum eine sehr häufige Sache ist. Also lass uns schnell das Terminal aufrufen und einen Probelauf machen.

## 2.1. Den Sammelbefehl definieren und im Terminal testen

Ich werde schummeln und den Beispiel-Bash-Code aus dem Trainingsmaterial kopieren und einfach Enter drücken.

Was wir hier sehen können, ist, dass wir diesen echo Befehl dreimal zu drei verschiedenen Ausgabedateien ausgeführt haben, die ich hier sehen kann. Und dann den cat Befehl verwendet haben, um die Ausgabe jeder dieser drei verschiedenen Dateien zu drucken und das zu einer einzigen gesammelten Datei umzuleiten.

Und wenn ich "cat COLLECTED-output" mache, kannst du sehen, es hat die Inhalte dieser drei verschiedenen Dateien, jetzt in einer einzigen Datei.

## 2.2. Einen neuen Prozess für den Sammelschritt erstellen

Also lass uns sehen, ob wir dasselbe innerhalb unserer Nextflow Pipeline replizieren können.

Lass uns nach oben scrollen und einen dritten Prozess erstellen. Ich werde diesen vorherigen kopieren, und diesmal werde ich ihn Collect Greetings nennen.

Im Bash-Terminal haben wir es collected output txt genannt. Also werde ich hier denselben path output sagen. Und ich werde die Umleitung hier machen, damit es auf dieselbe Weise gespeichert wird.

Okay. Wir müssen ändern, was am Anfang dieses Befehls passiert, und wir müssen darüber nachdenken, was hier die Eingabedatei ist. Tatsächlich wird dieser Prozess mehrere Eingabedateien nehmen. Ich werde path behalten und ich werde das zu einer neuen Variable namens input files, Plural, ändern.

Ich werde sie dann wieder catten, wie wir es in unserem Bash-Skript gemacht haben. Und ich werde hier die Variable verwenden.

Jetzt könntest du denken, das würde nicht funktionieren. Wir haben zuvor Fehler gesehen, wo ein Array von Zeichenketten oder ein Array von Pfaden an einen Prozess übergeben wurde und das einen Fehler verursacht hat. Aber tatsächlich wird Nextflow das hier automatisch für uns auf die richtige Weise handhaben. Es wird mehrere verschiedene Eingabedateien nehmen, und es wird einfach die verschiedenen Dateipfade hier drucken.

Natürlich hilft es, dass der cat Befehl eine Reihe von Dateinamen wie diesen nehmen kann. Wenn ich einen anderen Befehl verwenden würde, der ein Argument vor jedem Dateipfad oder so erfordert, müssten wir hier etwas mehr Code und Logik haben, um die Iteration dieser Dateipfade handhaben zu können. Aber in diesem Fall sollte es einfach funktionieren.

## 2.3. Den Sammelschritt zum Workflow hinzufügen

Okay, lass uns runter zum Workflow gehen und unseren neuen Prozess hinzufügen. Collect greetings. Und wieder nehmen wir die Ausgabe von convert to upper out. Lass uns das speichern.

Probieren wir es aus. nextflow run hello workflow.

Okay, der Workflow lief, aber etwas ist hier etwas seltsam. Wir haben drei Ausführungen des ersten Schritts, was wir erwarten. Drei Aufgaben für den zweiten, aber wir haben auch drei Aufgaben am Ende, wenn wir erwartet haben, nur eine einzige Aufgabe hier zu haben, die alle Ausgaben zusammenführt.

Wenn wir in unser results Verzeichnis gehen. Wir sehen auch, dass die gesammelte Ausgabe nur einen einzigen Wert hat statt aller drei. Das liegt daran, dass diese Ausgabedatei dreimal mit drei verschiedenen Werten überschrieben wurde.

Das macht Sinn, weil wir hier eine Ausgabe zu einer Eingabe auf dieselbe Weise übergeben haben wie im vorherigen Schritt.

## 2.4. Einen Operator verwenden, um die Grüße in eine einzige Eingabe zu sammeln

Also brauchen wir hier einen Operator, um diesen Channel mit drei Elementen zu nehmen und sie zu einem einzigen Element zusammenzufassen, sodass dieser finale Prozess nur einmal läuft.

Dafür werden wir den collect Operator verwenden. Ich kann das direkt innerhalb des Workflows machen. Ich kann .out machen und hier am Ende zu einem Operator verketten .collect.

Speichern. Und dann für die Zwecke dieses Trainings werde ich auch einige view Operatoren machen, wie wir es zuvor getan haben, damit wir einen Blick auf diesen Channel vor und nach der Verwendung des collect Operators werfen können, damit wir verstehen können, was passiert.

Ich werde diesen Channel nehmen, das collect loswerden und dot view greetings machen, und dann werde ich diese Zeile duplizieren, den collect Operator hinzufügen. Und das zu after ändern.

Das ist getrennt davon, wo wir das aufrufen, aber das ist in Ordnung, weil wir dieselben Operator-Aufrufe auf demselben Ausgabe-Channel verwenden.

Okay, lass uns speichern und im Terminal ausprobieren. Werde nextflow run machen. Hello, workflow. Unser Skript erneut ausführen.

Okay. Das sieht besser aus. Wie zuvor können wir sehen, dass die ersten beiden Prozesse dreimal laufen und jetzt lief unser finaler Prozess nur einmal.

Wenn wir uns ansehen, was vom view Operator gedruckt wurde, hier unten, sagten wir before collect, was diese Ausgabe hier ist, und das wurde dreimal gedruckt. Und du kannst sehen, es gibt einen einzelnen Pfad für jeden von denen. Und dann nach collect kannst du sehen, dass wir dieses Array von drei Pfaden haben. Also ist das wie erwartet.

Okay, lass uns die results Datei überprüfen und sehen, ob es diesmal das ist, was wir erwarten. Tatsächlich, es gibt jetzt drei Zeilen in der Datei - das hat diese drei Ausgaben erfolgreich zu einer einzigen Ausgabedatei zusammengefügt. Fantastisch.

Okay, ich werde aufräumen und lass uns zum nächsten Schritt gehen. Und ich werde diese view Anweisungen löschen, nur um die Dinge sauber zu halten.

## 3. Mehr als eine Eingabe an einen Prozess übergeben, um die finale Ausgabedatei eindeutig zu benennen

Okay. Bisher haben alle unsere Prozesse nur eine einzige Eingabe genommen. Wir werden jetzt eine Übung machen, wo wir mehr als eine Eingabe zu einem Prozess hinzufügen, um zu sehen, wie das funktioniert. Dafür werden wir dieses collect greetings Beispiel verwenden.

Jedes Mal, wenn ich den Workflow ausgeführt habe, hat es diese Datei im results Verzeichnis überschrieben, was vielleicht nicht das ist, was wir wollen.

## 3.1. Den Sammelprozess modifizieren, um einen benutzerdefinierten Namen für die Ausgabedatei zu akzeptieren

Also für dieses Beispiel werden wir einen zusätzlichen Parameter übergeben, damit wir den Ausgabedateinamen anpassen können.

Eine zweite Eingabe zu einem Prozess hinzuzufügen ist sehr einfach. Ich füge einfach eine zweite Zeile im input Block hinzu. Diesmal wird es ein value sein, anstatt ein path, weil wir einen String übergeben wollen, und ich werde es batch underscore name nennen.

Ich kann diese Variable jetzt im script Block verwenden, und ich werde collected dash dollar batch name sagen.

Ich verwende hier geschweifte Klammern um den Variablennamen. Das ist nur, um es vom Rest eines Strings getrennt zu halten, und es wird wahrscheinlich in diesem Fall nicht gebraucht, aber ich denke, es macht es einfacher zu lesen.

Okay. Schließlich, denke daran, den Ausgabe-path zu aktualisieren, weil sich jetzt der Dateiname geändert hat, also werde ich dasselbe tun und den batch name in die Ausgabe von path wie erwartet einfügen.

## 3.2. Einen Batch-Befehlszeilenparameter hinzufügen

Wir müssen jetzt einen batch name von irgendwo übergeben, und ich werde einen zweiten Parameter erstellen, um das zu tun, damit wir es auf der Befehlszeile machen können, wenn wir den Workflow ausführen.

Also werde ich params batch name machen, und standardmäßig, lass uns das test batch nennen. Jetzt kann ich diese spezielle Parameter-Variable unten verwenden, wo wir den Prozess aufrufen.

Und tatsächlich sagt uns VS Code, dass es jetzt nicht genug Argumente für diesen Prozess gibt und dass es eine zweite Eingabe erwartet.

Mache einfach Komma und übergebe unsere neue Variable und der Fehler verschwindet.

Beachte, dass die Reihenfolge der Eingaben hier wirklich wichtig ist. Die erste Prozesseingabe war der path, und die zweite Eingabe ist der Name. Wenn ich die Reihenfolge hier ändere, muss ich auch die Reihenfolge ändern, wenn ich den Prozess aufrufe. Andernfalls. Als nächstes werden wir den falschen Channel zur falschen Eingabe übergeben.

## 3.3. Den Workflow ausführen

Okay, lass uns es ausprobieren und sehen, ob es funktioniert. Lass uns "nextflow run hello- workflow machen. Okay, es lief wie zuvor. Lass uns ins results Verzeichnis schauen.

Tatsächlich, unser Dateiname hier heißt jetzt " collected test batch output txt". Fantastisch.

Und jetzt lass uns sehen, ob wir das überschreiben können, indem wir nochmal ausführen. Diesmal werde ich --batch_name machen, um mit diesem speziellen Parametervariablennamen hier übereinzustimmen. Und ich werde es demo output nennen.

Den Workflow nochmal ausführen und wir werden sehen, ob etwas passiert.

Okay, wir haben jetzt ein collected demo output .txt. Und weil dieser Dateiname anders ist als jener, hat es ihn nicht überschrieben. Beide sind jetzt im results Verzeichnis vorhanden.

## 4. Eine Ausgabe zum Sammelschritt hinzufügen

Okay, dort haben wir also gezeigt, wie man mehrere Eingaben an einen Prozess gibt, aber wie ist es mit mehreren Ausgaben? Für dieses Beispiel werden wir die Anzahl der Grüße berechnen, die verarbeitet werden, und das als sekundäre Ausgabe für diesen collect greeting Schritt ausgeben.

## 4.1. Den Prozess modifizieren, um die Anzahl der Grüße zu zählen und auszugeben

Wir werden hier einen kleinen Trick machen. Nextflow Prozesse haben diesen script Block mit einem mehrzeiligen String, und der wird als Bash-Ausgabe an den dot Befehl dot sh übergeben. Aber wir können tatsächlich beliebigen benutzerdefinierten Code darüber schreiben, und der wird als Teil einer Aufgabe ausgeführt, aber nicht im Bash-Skript enthalten.

Eine der eingebauten Funktionen in der Nextflow Syntax heißt size. Also werde ich die path Eingabe nehmen, und ich werde count underscore greetings sagen, nur um einen Variablennamen zu definieren. Ich werde die input files nehmen und ich werde "size" darauf aufrufen.

Diese Funktion wird die Größe dieses Eingabe-Channels zählen und sie einer Variable zuweisen.

Wir können diese Variable jetzt als Teil des output Blocks zurückgeben. Also sagen wir val, weil es ein Wert ist, keine Datei. Und count greetings.

Jetzt reicht das für sich selbst aus, und wir könnten jetzt auf diese verschiedenen Ausgaben von diesem Prozess zugreifen. Allerdings müssten wir auf sie in einer positionellen Weise zugreifen. Also mit einem Index-Schlüssel wie null und eins.

Um es ein bisschen einfacher zu machen, an die Ausgaben zu kommen, können wir sie benennen und wir machen das, indem wir eine emit Anweisung verwenden.

Also machen wir Komma emit out file oder wie auch immer ich das nennen möchte. Und ich mache hier emit count. Das ist im Grunde nur ein Dekorator, der uns hilft, etwas saubereren Code zu schreiben, damit wir später im workflow Block leicht auf die spezifischen Ausgaben verweisen können.

## 4.2. Die Ausgabe am Ende des Workflows berichten

Okay. Wenn ich runter zum workflow Block scrolle, kann ich jetzt die Ausgaben von collect greetings nehmen, collect greetings machen, dot out, und wir können unsere zwei benannten Ausgaben sehen, die hier von der VS Code Extension vorgeschlagen werden. Sehr praktisch.

Also werde ich dot count machen, um den count Wert zu bekommen, den wir gerade erstellt haben, und ich werde view machen, damit es in der Befehlszeile gedruckt wird. Damit wir es sehen können, wenn wir den Workflow ausführen.

Lass uns hier etwas in die Closure schreiben, nur um es ein bisschen schöner zu machen. num greetings, there were greetings greetings.

Und uns ist die andere Ausgabe eigentlich egal, weil wir sie nicht als Eingabe für irgendwelche anderen Prozesse verwenden. Aber du kannst sehen, wie wir das leicht als Eingabe für einen anderen Prozess übergeben könnten, wenn wir wollten, downstream.

## 4.3. Den Workflow ausführen

Wir werden speichern. Lass uns ins Terminal schauen und es ausprobieren.

Okay, fantastisch. Hier sind wir. There are three greetings. Das ist genau richtig.

Okay, super. Das ist das Ende dieses Kapitels. Wir sind fertig, weil wir es soweit geschafft haben. Du baust jetzt einen ziemlich realistischen Workflow auf, wo wir in der Lage sind, Eingaben und Ausgaben und Logik innerhalb unseres Workflows zu handhaben.

Wenn diese Workflow-Dateien länger werden, beginnen sie ein wenig unhandlich zu werden. Also werden wir im nächsten Kapitel uns ansehen, wie wir Nextflow Code in separate Dateien modularisieren können, sodass es einfacher ist, den Code innerhalb des Workflows zu finden und zu pflegen.

Begleite uns im nächsten Video für Kapitel vier. Hello Modules.

[Nächstes Video-Transkript :octicons-arrow-right-24:](04_hello_modules.md)
