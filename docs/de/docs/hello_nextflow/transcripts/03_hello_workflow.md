# Teil 3: Hello Workflow - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../03_hello_workflow.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern im Material.

## Willkommen und Rückblick

Hallo und willkommen zurück zu Teil drei von Hello Nextflow. Dieser Teil heißt Hello Workflow, und in diesem Teil des Kurses rechtfertigen wir wirklich den Namen Pipeline oder Workflow.

Wir nehmen unser einfaches Pipeline-Skript mit seinem einen Prozess und fügen weitere Prozesse hinzu. Wir schauen uns an, wie Nextflow diese Orchestrierung und den Datenfluss durch die Pipeline handhabt.

Gehen wir zurück zu unseren Code Spaces. Du siehst, ich habe alle .nextflow\*-Verzeichnisse und die work-Verzeichnisse gelöscht, um alles sauber zu halten. Keine Sorge, wenn du diese Dateien aus den vorherigen Kursteilen noch hast.

Wir arbeiten mit einer Datei namens hello-workflow.nf. Wie zuvor repräsentiert diese im Grunde das Skript, das wir bis jetzt aufgebaut haben, und gibt uns einen sauberen Ausgangspunkt. Und wieder unten in der Ausgabe sehen wir, dass der Pfad jetzt hello_workflow ist. Die veröffentlichten Dateien sollten also in ein anderes Unterverzeichnis in deinem results-Ordner gehen.

Um zusammenzufassen, wo wir bisher sind: Wir haben einen einzelnen Prozess hier mit einer Eingabe greeting, einer Ausgabe greeting file. Und dann das einfache Bash-Skript, das nur einen echo-Befehl in eine Datei macht.

Wir haben eine einzelne Workflow-Eingabe, den params-Block hier, wo wir sagen, dass ein Pfad erwartet wird, und der Standard ist data/greetings.csv, was diese Datei hier oben ist.

Dann im Workflow selbst haben wir einen main-Block. Wir erstellen einen Kanal. Wir parsen das CSV in Zeilen und nehmen dann das erste Element jedes Arrays, und wir übergeben diesen Kanal an diesen Prozess, der dann drei Aufgaben generiert, und wir veröffentlichen aus dem Workflow die Ausgaben von diesem Prozess.

Und schließlich im output-Block sagen wir Nextflow, diese Dateien von diesem Kanal in das Verzeichnis namens hello_workflow zu veröffentlichen. Und diese Dateien zu kopieren, anstatt sie per Soft-Link zu verknüpfen.

## 1. Einen zweiten Schritt zum Workflow hinzufügen

Okay, in diesem Teil fügen wir einen zweiten Prozess zu unserem Workflow hinzu. Wir nehmen die Ausgaben des sayHello-Prozesses und verarbeiten sie in einem zweiten Schritt, der alle Buchstaben in diesen Dateien in Großbuchstaben umwandelt - convertToUppercase.

Das ist nur ein einfaches Beispiel, nur etwas einfache String-Verarbeitung, aber es zeigt dir, wie wir die Logik innerhalb des Workflows nutzen können.

Wir werden dafür einen Bash-Befehl namens "tr" verwenden, was für translate steht. Das ist ein Unix-Befehl, der schon ewig existiert. Wenn du damit nicht vertraut bist, kann ich das verstehen. Ich glaube, ich habe ihn vor dem Training nie benutzt, aber du kannst ihn schnell im Terminal ausprobieren. Wenn ich "echo 'hello world'" mache und dann zu 'tr' pipe und dann in Anführungszeichen den Zeichenbereich angebe, also A bis Z, Kleinbuchstaben, und dann A bis Z Großbuchstaben. Und es sagt einfach, übersetze diese Buchstaben in diese Buchstaben.

Und wenn ich Enter drücke, siehst du, dass jetzt alles großgeschrieben ist. Sehr schön, wenn du gerne Leute anschreist.

Das ist also ein sehr einfacher Bash-Befehl, den wir in unserem zweiten Prozess verwenden werden.

## 1.2. Den Großbuchstaben-Schritt als Nextflow-Prozess schreiben

Wenn ich zurück zu meinem Skript gehe, werde ich ein bisschen schummeln und den Code einfach aus der Dokumentation für das Training kopieren. Aber du kannst genau sehen, was passiert.

Wir haben hier einen neuen Prozess. Diesen haben wir convertToUpper genannt, aber wir könnten ihn nennen, wie wir wollen.

Wir haben einen einzelnen Eingabepfad, wie zuvor. Es ist kein Wert-Kanal, es ist ein Pfad-Kanal. Und dann eine einzelne Ausgabe.

Im script-Block machen wir "cat" auf der Eingabedatei. Und wir können das in geschweifte Klammern setzen, wenn wir wollen, was diese Variable nimmt. Und wir führen denselben Bash-Befehl in der Pipe aus und schreiben die Ergebnisse in eine Datei mit diesem Dateinamen, und das wird vom Ausgabepfad aufgenommen.

Wir müssen jetzt etwas mit diesem neuen Prozess machen. Also gehen wir runter zum Workflow, wo wir die verschiedene Logik eines Workflows aufbauen, und nach diesem ersten Prozess werden wir unseren zweiten Prozess ausführen. Also convertToUpper ist der Name des Prozesses hier.

Er nimmt eine Eingabe, also können wir ihn nicht einfach so aufrufen. Wir wollen die Ausgabe des ersten Prozesses verarbeiten. Also genau wie wir das mit diesem sayHello out gemacht haben, wo wir diese Ergebnisse veröffentlichen. Wir wollen dieselben Ergebnisse hier als Eingabe verwenden, also können wir die kopieren und dort einfügen.

Wir wollen den sayHello-Prozess ".out", und Nextflow weiß, dass dies einen einfachen einzelnen Ausgabedatensatz hier bedeutet, was diese Datei ist. Das wird dann als Eingabe an einen zweiten Prozess übergeben.

## 1.5. Die Workflow-Ausgabe-Veröffentlichung einrichten

Okay. Und schließlich, damit wir die Ergebnisse dieses zweiten Prozesses tatsächlich speichern, müssen wir sie auch aus dem Workflow veröffentlichen und dann im output-Block definieren, gleiche Syntax wie zuvor. Also können wir das kopieren und second outputs sagen, oder wie auch immer du es nennen willst.

Nimm den Prozessnamen, an dem wir interessiert sind, convertToUpper out, und dann hier unten im output-Block. Füge das hinzu und wir könnten dieselben Attribute hier machen. Also wollen wir diese Dateien auch im Hello Workflow-Unterverzeichnis, und wir wollen sie auch kopieren.

Großartig. Lass uns versuchen, es auszuführen. Also wenn ich das Terminal hochbringe und "nextflow run hello-workflow.nf" mache, werden wir sehen, was es macht. Schauen wir, ob es anders aussieht als die vorherigen Teile.

Also startet es Nextflow. In der Dokumentation steht, das mit "-resume" zu machen, aber ich habe mein ganzes work-Verzeichnis gelöscht, also hätte es hier keinen Unterschied gemacht. Aber wenn du das getan hast, wird das auch funktionieren.

Und es sieht fast genau gleich aus. Aber du kannst jetzt sehen, dass es eine zweite Ausgabezeile hier gibt, wo du den Namen des zweiten Prozesses sehen kannst, den wir gerade hinzugefügt haben. Und tatsächlich kannst du sehen, dass er dreimal erfolgreich gelaufen ist.

Brillant. Wenn ich meine vorherigen work-Verzeichnisse noch hätte und das mit "-resume" gemacht hätte, wären diese gecacht worden, nur der erste Schritt in der Pipeline. Weil diese Ausgaben genau gleich waren, hätte Nextflow gewusst, diese wieder zu verwenden.

Und so kannst du sehen, wie du -resume verwenden kannst, um deinen Workflow iterativ aufzubauen, Schritt für Schritt, wenn du musst.

Okay, lass uns einen Blick in das results-Verzeichnis hier oben werfen und sehen, ob es funktioniert hat. Wir können sehen, dass wir hier oben einige weitere Dateien haben. Wir haben unsere ursprünglichen Dateien wie zuvor vom ersten Prozess. Und tatsächlich haben wir unsere upper-Dateien und die Buchstaben sind alle großgeschrieben, also hat es funktioniert. Es ist wirklich schön zu sehen.

Es ist auch interessant, einfach in diese work-Verzeichnisse zu schauen. Wie zuvor entspricht der Hash hier den work-Verzeichnissen. Also wenn ich in "ls work" schaue und das dann erweitere, werden wir die verschiedenen Dateien hier sehen.

Wir sehen die Ausgabedatei vom ersten Prozess, die hier als Eingabe hereingezogen wurde. Und wir können die neue Ausgabedatei sehen, die generiert wurde.

Wenn ich das jetzt mit "-la" mache, um alle Dateien aufzulisten und anzuzeigen, werden wir ein paar mehr Dinge sehen. Erstens wirst du sehen, dass diese Datei tatsächlich ein Soft-Link zum ersten Prozess ist. Das ist im Grunde immer ein Soft-Link, wenn es sein kann, um Dateispeicherplatz zu sparen. Wir veröffentlichen die Dateien hier nicht und es referenziert einfach diese Datei von einer ersten Aufgabe in eine zweite Aufgabe, sodass alles in diesem einen Arbeitsverzeichnis gekapselt ist und sicher und isoliert von allem anderen.

Und das muss dort sein, denn wenn wir uns die .command.sh-Datei anschauen, also wenn ich "cat work/b8/56\*" mache, kannst du sehen, dass die Dateipfade hier relativ sind, also cattet es diese Eingabedatei, die in dasselbe Arbeitsverzeichnis per Soft-Link verknüpft wurde.

So wird also jedes work-Verzeichnis aussehen. Wenn du es in Nextflow anschaust, wirst du alle Eingabedateien dort haben, die in dieses Arbeitsverzeichnis gestaged wurden. Und dann hast du auch alle Ausgabedateien, die erstellt wurden. Das ist also großartig. Das sieht so aus, wie wir es erwarten.

## 2.1. Den Sammelbefehl definieren und im Terminal testen

Okay, gehen wir zurück zu unserem Workflow. Was ist der nächste Schritt, den wir machen wollen?

Wir haben jetzt zwei Prozesse und sie nehmen diese eine CSV-Datei, parsen sie und teilen sie auf. Und dann haben wir drei Aufgaben für jeden dieser Prozesse und Nextflow handhabt die Parallelisierung von all dem, sodass alles nebeneinander läuft, wo möglich.

Diese Art, Arbeit aufzuteilen, um Dinge parallel laufen zu lassen, ist sehr üblich. Und das Gegenteil davon ist dann, alles wieder zusammenzuführen. Das werden wir also mit unserem finalen Prozess im Workflow machen - wir haben hier einen dritten, der diese drei verschiedenen Ausgaben nimmt und sie alle in eine einzelne Datei kombiniert.

Wir können das ganz einfach in einem Terminal machen, um ein Gefühl dafür zu bekommen, wie das aussehen wird.

Wenn ich zum results-Ordner gehe. Also "cd results/hello_workflow/", und wir haben hier alle UPPER-Dateien. Ich kann einfach "cat" verwenden, was wir benutzen, um den Inhalt dieser Datei auszugeben, und du kannst "cat" mehrere Dateien geben und es wird eine nach der anderen lesen.

Also kann ich "UPPER-\*" sagen, was mir dieselbe Liste von drei Dateinamen mit Bash-Expansion gibt. Und ich kann combined.txt sagen. Ich glaube, in der Dokumentation sind die genauen Dateinamen aufgelistet, aber es macht dasselbe.

Jetzt, wenn ich "cat combined.txt" verwende, können wir sehen, dass wir den Dateiinhalt von allen drei dieser Dateien haben.

Das ist also im Grunde alles, was dieser Prozess machen wird - wir werden versuchen, ihm alle verschiedenen Ausgabedateien von einem vorherigen Prozess in einer einzelnen Prozessaufgabe zu geben, und dann werden wir sie zusammen "catten" und die Ausgabedatei speichern.

## 2.2. Einen neuen Prozess für den Sammelschritt erstellen

Okay, also fügen wir unseren neuen Prozess hinzu. Ich werde das aus den Trainingsmaterialien einfügen, und du kannst sehen, es hat uns hier ein bisschen eine Übung für den Leser mit diesen Fragezeichen gelassen. Aber du kannst den allgemeinen Umriss des Prozesses sehen - es ist im Grunde das, was wir gerade im Terminal gemacht haben, wo wir "cat" von einem Haufen Eingabedateien machen und es in eine Ausgabedatei hier namens collected schreiben, und dann erwartet die Ausgabe wieder diesen einzelnen Pfad.

Also brauchen wir hier irgendeine Art von Eingabe und es wird ein Set von Pfaden sein. Also wieder definieren wir einen Eingabepfad-Kanal und nennen wir ihn input_files. Jetzt hat uns das vorher einen einzelnen Pfad hier gegeben, aber ein Pfad kann auch mehrere Dateien hier haben, auch wenn es immer noch eine einzelne Deklaration ist.

Ich werde das hier runter kopieren, weil wir diese Dateien "catten" wollen. Und du denkst vielleicht, dass wir hier einige Probleme mit dem Ausgeben eines Arrays oder so haben, aber Nextflow ist generell ziemlich vernünftig, wenn es darum geht. Und wenn es einen Kanal mit mehreren Dateien darin wie diesem bekommt, wird es sie alle mit Leerzeichen-Trennzeichen zusammenfügen. Das wird uns also die korrekte Syntax geben.

Das ist großartig. Also jetzt verdrahten wir unseren neuen Prozess. Ich gehe runter zum Workflow. Ich sage combine the outputs, der neue Prozessname, und genau wie zuvor. Ich nehme diesen vorherigen Prozess, convertToUpper und mache ".out".

Großartig. Probieren wir es aus und schauen, ob es im Terminal funktioniert. Wenn ich einfach ein paar Verzeichnisse zurückgehe und dann den Nextflow-Befehl erneut ausführe, werden wir sehen, was passiert.

Also der Workflow wurde gestartet und jetzt kannst du sehen, dass wir drei verschiedene Prozessnamen haben, was großartig ist. Die ersten beiden sehen beide gleich aus wie zuvor, und der dritte neue läuft, was gut ist.

Allerdings gibt es hier etwas Seltsames. Wir wollten diese Ausgabedateien in eine einzelne Datei kombinieren, und doch können wir sehen, dass dieser Prozess dreimal gelaufen ist, nicht einmal.

Tatsächlich, wenn wir in eines dieser work-Verzeichnisse gehen. Und "cat work/" "collected" machen, dann werden wir sehen. Es gibt hier nur ein einzelnes Wort, nicht drei.

Und was passiert ist, ist, dass Nextflow diese Parallelisierung genau wie in den vorherigen Schritten fortgesetzt hat. Und dieser Prozess gab uns einen Kanal mit drei Elementen, und diese drei Kanalelemente wurden an unseren nachgelagerten Prozess übergeben, der drei Prozessaufgaben generierte.

Es hat im Grunde versucht, dreimal separat zu sammeln, und jedes Mal hatte es nur eine einzelne Datei, also hat es einfach cat einzelne Datei zu einer Ausgabe gemacht, und tatsächlich können wir das auch in der .command.sh-Datei sehen.

Wenn ich .command.sh mache, können wir sehen, dass es hier nur einen einzelnen Dateinamen hat und nur eine einzelne Datei in dieses Arbeitsverzeichnis gestaged wurde.

## 2.3. Den Sammelschritt zum Workflow hinzufügen

Also müssen wir Nextflow irgendwie sagen, dass es alle diese Ausgaben von einem vorherigen Prozess zusammenbringen und sie diesem nachgelagerten Prozess als ein einzelnes Kanalelement geben soll, anstatt drei.

Wir machen das mit einem Kanal-Operator namens _collect_.

Das ist ein super nützlicher Operator, den du in Nextflow-Pipelines die ganze Zeit sehen wirst. Das ist hier ein Kanal, dieser Ausgabekanal, genau wie der, den wir oben erstellt haben. Und so können wir Kanal-Operatoren daran anhängen, genau wie wir es zuvor gemacht haben. Wir können einfach Punkt machen, und dann in diesem Fall collect, Klammern.

Und das ist alles, was wir brauchen. Das wird dann diesen Kanal manipulieren, bevor er an diesen Prozess übergeben wird.

Wenn du sehen willst, was damit passiert, können wir es auch hier anzeigen. Also hier, das ist überhaupt nicht mit dem Ausführen dieses Prozesses verbunden, also könnte ich es an jedem Punkt nach dem Ausführen dieses Prozesses setzen. Aber wir nehmen denselben Ausgabekanal, und wir schauen ihn mit .view an, und dann schauen wir ihn wieder mit .collect.view an.

Und wenn wir das ausführen, wird es uns die zwei verschiedenen Strukturen dieses Kanals zeigen, vor und nach collect. Also probieren wir das jetzt aus. Okay, ich habe gerade ein bisschen rausgezoomt, weil einige der Ausgaben ziemlich lang sind, aber wenn ich die Pipeline ausführe, werden wir sehen, ob es funktioniert.

Ich hoffe, ein dritter Prozess wird nur einmal laufen, weil er die Ausgaben sammelt, und tatsächlich kannst du sehen, collectGreetings als eins von eins. Das hat also nur eine Aufgabe ausgeführt.

Und dann, wenn wir uns die view-Anweisungen anschauen, haben wir drei view-Anweisungen für die drei Elemente von vorher, mit einem Dateipfad in jedem.

Und dann nach dieser collect-Anweisung wurde das nur einmal ausgelöst, weil es ein einzelnes Element in diesem Kanal gibt. Und jetzt haben wir diese Liste von drei verschiedenen Dateipfaden.

Das ist genau das, was wir erhofft haben. Und du kannst hoffentlich sehen, das ist im Grunde das Gegenteil von diesem "map"-Operator, den wir gemacht haben, um von den CSV-Arrays in separate Kanalelemente zu gehen. Jetzt nehmen wir separate Kanalelemente und fügen sie zurück in ein einzelnes Array.

Großartig, wir können diese view-Anweisungen aufräumen. Wir brauchen diese nicht mehr. Wir können zum nächsten Schritt übergehen.

Bevor ich weitermache, und bevor ich es vergesse, werde ich hier eine neue publish-Anweisung hinzufügen. Third output. Du kannst das in deinem Workflow semantischer und beschreibender nennen. Und dann werde ich das wieder zum output-Block hinzufügen und path 'hello_workflow' mode 'copy' sagen. Nur damit die Ausgabedatei, die von diesem Prozess generiert wird, in unserem results-Ordner hier oben gespeichert wird.

Nur um schnell zu überprüfen, dass das funktioniert. Sollte jetzt ein bisschen sauberer sein, weil wir diese view-Anweisungen nicht haben. Und wir werden sehen, ob wir unsere neue Ausgabedatei hier oben bekommen. Eins von, eine Aufgabe lief, bekam eine neue Datei namens collected, und jetzt haben wir alle drei dieser Wörter. Fantastisch. Was kommt als Nächstes?

## 3. Zusätzliche Parameter an einen Prozess übergeben

Okay. Als Nächstes schauen wir uns an, wie man mehrere Eingaben in einen einzelnen Prozess handhabt. Bisher kannst du sehen, dass alle unsere Prozesse nur eine Sache als Eingabe nehmen. Sie haben alle eine einzelne Zeile unter ihrer Eingabe.

Wir werden das demonstrieren, indem wir Nextflow erlauben, eine andere Batch-Kennung anzugeben, sodass du vielleicht diesen Workflow mehrmals ausführst und ihm jedes Mal eine andere Batch-ID geben kannst.

Ich werde einfach eine zweite Zeile in der Eingabe hier für collectGreetings hinzufügen. Und ich werde es "val" nennen, weil das ein String ist. Jetzt ist es ein Wert, kein Pfad, und ich werde es "batch_name" nennen.

Dann werde ich das Skript hier unten bearbeiten, um diese Variable zu verwenden, und ich werde versuchen, es an derselben Stelle wie im Trainingsmaterial zu platzieren. Also setze ich es in die Mitte dieses Dateipfads COLLECTED-$\{batch_name\}-output.

Noch nicht ganz fertig. Denk daran, dass wir Nextflow sagen müssen, wie die Ausgabedateinamen sein werden. Also müssen wir auch dasselbe hier oben machen: COLLECTED-$\{batch_name\}-output.txt".

Fantastisch. Nextflow bekommt jetzt eine zweite Variable-Eingabe und interpoliert diese in das Skript und die Ausgabe.

Eine letzte Sache, wir müssen jetzt finden, wo das aufgerufen wird, und wir müssen die zweite Eingabe an den Prozess übergeben. Das ist genau wie jede andere Eingabe in eine Funktion in jeder anderen Sprache.

Genau wie wir es früher im Training gemacht haben, werde ich hier das spezielle "params" verwenden, und wir werden es "params.batch" nennen, sodass wir eine --batch CLI-Option haben können. Und jetzt kannst du sehen, dass unser Prozess hier zwei separate Eingaben hat, nur durch Komma getrennt, die übergeben werden.

Es ist wirklich wichtig, die Reihenfolge richtig zu bekommen, also die Reihenfolge der Argumente hier für channel und dann den param muss übereinstimmen. Der channel und der batch name dort. Das ist nur positionelles Matching.

Okay. Ich kann diese Pipeline jetzt direkt mit --batch ausführen, aber lass uns zuerst das Richtige tun und es in der Eingabe hier in Params definieren. Also werde ich es zu batch hinzufügen und dann werden wir sagen, es ist ein String und geben wir ihm einen Standard. Also nennen wir es einfach batch. Okay? Jetzt versuchen wir, den Workflow auszuführen.

--batch Trio. Ich glaube, das steht im Trainingsmaterial, aber wir könnten jeden String verwenden, den wir dort wollen. Und hoffentlich werden wir sehen, dass diese Ergebnisausgabedatei hier oben erscheint.

Und tatsächlich, COLLECTED-trio-output - das hat richtig funktioniert. Es hat unsere Datei umbenannt. Und du kannst dir jetzt vorstellen, das ist nützlich, denn wenn ich das nochmal mit einem anderen Batch-Namen ausführe, wie replicate_two, dann wird es uns hier oben einen anderen Batch-Namen geben.

Und und es wird dann die Ausgabedateien in diesem Fall nicht überschreiben. Das ist also schön.

## 4. Eine Ausgabe zum Sammelschritt hinzufügen

Okay, wir haben jetzt also mehrere Eingaben zu unserem Prozess hier. Aber was passiert, wenn wir mehrere Ausgaben erstellen wollen? Unser Beispiel hier ist dann, dass wir einen Bericht für diesen Prozess erstellen werden, der einfach sagt, so viele Dateien wurden gesammelt.

Und wir machen das mit einem echo-Befehl hier. Also können wir echo sagen. There were, ich werde das aus einem Trainingsmaterial kopieren, damit du mir nicht beim Tippen zusehen musst.

There were $\{count_greetings\} greetings in this batch, und das in eine neue Datei jetzt namens $\{batch_name\} speichern, also dieselbe Variable, wir können die so oft wiederverwenden, wie wir wollen, report.txt.

## 4.1.1. Die Anzahl der gesammelten Begrüßungen zählen

Wir müssen das irgendwie berechnen. Wir könnten diese Logik im Bash-Skript machen, wenn wir wollten, mit Bash-Logik. Allerdings können wir auch direkt innerhalb des Nextflow-Codes skripten, solange es innerhalb des script-Blocks im Prozess und über dem zitierten Abschnitt ist.

Alles hier wird nicht im finalen gerenderten Skript enthalten sein, und es wird einfach von Nextflow ausgeführt, wenn es eine Aufgabe rendert.

Also machen wir hier einfach etwas Logik. Wir erstellen eine neue Variable namens count_greetings. Wir nehmen den input files-Kanal hier, und wir rufen .size() darauf auf.

Okay, diese Funktion wird mir hier eine Nummer in diese Variable geben, und jetzt ist unsere Warnung weg, weil diese Variable definiert wird.

Okay, also erstellen wir diese zweite Datei im work-Verzeichnis, aber wir müssen Nextflow sagen, dass es sie als veröffentlichte Ausgabe dieses Prozesses erwarten soll. Also machen wir das mit genau derselben Syntax wie für die erste Datei.

Wir sagen path, weil es, wieder, wir könnten hier eine Variable veröffentlichen, wenn wir wollten mit "val", aber wir werden "path" sagen. Und dann den erwarteten Dateinamen. Beachte, es ist hier nicht hervorgehoben. Das liegt daran, dass ich einfache Anführungszeichen verwendet habe. Ich muss doppelte Anführungszeichen verwenden.

## 4.1.2. Die Berichtsdatei ausgeben und Ausgaben benennen

Okay, das ist großartig. Und wir könnten jetzt anfangen, auf diese Ausgaben hier unten zuzugreifen, genau wie ich es hier gemacht habe. Aber jetzt ist es ein Array von verschiedenen Objekten, also könnte ich collectGreetings.out[0] machen, um das erste zu bekommen, oder eins, um das zweite zu bekommen, was unser neuer Bericht ist.

Aber ich mag das nicht wirklich sehr, weil es ziemlich einfach ist, das Index-Zählen durcheinander zu bringen. Und du sitzt da und zählst viel Zeilen und du fügst eine neue Ausgabe hinzu und plötzlich bricht alles zusammen. Also

ist es viel schöner, stattdessen alles nach Namen zu referenzieren. Und wir können das mit einem speziellen Schlüssel hier namens "emit" machen.

Also können wir das nennen, wie wir wollen. Nennen wir es emit outfile, und emit reports. Wenn du diese definierst und du kannst es bei einem oder vielen machen, liegt es an dir. Jetzt kann ich hier runtergehen und stattdessen kann ich dot out dot reports gehen und es einfach nach Namen aufrufen, was viel einfacher ist, deinen Code zu verstehen, wenn du ihn liest, und es ist sicherer bei Änderungen im Code.

Ich habe das .out.report hier hinzugefügt, aber eigentlich muss ich zwei verschiedene Ausgaben veröffentlichen lassen. Also werde ich es als etwas Interessanteres umbenennen wie collected und report und ist das, wie ich es genannt habe? Ich habe es out file genannt, sorry. Also dieser emit-Name hier outfile und report. Weil wir zwei verschiedene Ausgabekanäle veröffentlichen und wir müssen also beide im publish-Block referenzieren.

Dann müssen wir diese auch im output-Block definieren. Also habe ich das collected umbenannt, und wieder für reports, ein bisschen verbose hier, aber es ist wirklich nützlich, wenn du reinkommst, um einen neuen Workflow zu lesen, alle verschiedenen Ausgaben hier zu sehen, alle verschiedenen Kanäle nebeneinander aufgelistet, und es gibt Wege, das weniger verbose zu machen, was wir später berühren werden.

Okay, probieren wir es aus und führen unseren Workflow aus und schauen, was passiert.

Hoffentlich sollte es jetzt im Grunde gleich laufen wie zuvor. Und wir werden eine neue Ausgabedatei hier oben namens replicate_two, report bekommen. Und da haben wir es. Es hat sich geöffnet und es sagt, es gibt drei Begrüßungen im Batch, was wir erwartet haben, also ist es perfekt.

Wenn ich hier ins work-Verzeichnis gehe, nur um dir zu beweisen, dass es im Nextflow-Code ausgeführt wurde und nicht im Bash-Skript, kann ich zu cat work/ command.sh gehen, und du wirst hier sehen, dass es einfach diesen String direkt echot. There were three greetings in this batch, und so wurde diese Variable von Nextflow interpoliert. Sie wurde im script-Block berechnet, bevor es die .command.sh-Datei geschrieben hat. Also ist die resultierende Variablenberechnung im Grunde hart in diese codiert, bevor sie auf deiner Compute-Umgebung in diesem Fall ausgeführt wird.

Und so kannst du diese Trennung zwischen dem script-Block hier und allem darüber sehen. Ich hoffe, das macht Sinn.

## Fazit und Quiz

Okay, das ist das Ende dieses Teils von Hello Nextflow. Also wie zuvor, geh und schau dir das Quiz an. Mach es auf der Webseite oder in der CLI, geh durch einige der Fragen und überprüf einfach, ob du etwas vom Material verstanden hast, das wir behandelt haben. Schau, ob es dort etwas gibt, das etwas hervorhebt, das du vielleicht nicht verstanden hast. Nicht zu viele Fragen. Schön und einfach zu machen. Oder du kannst es auch hier unten auf der Webseite machen.

Und mach eine kleine Pause, einen kleinen Spaziergang und komm zurück und schließ dich uns in Teil vier von Hello Nextflow an, wo wir über Module sprechen werden. Vielen Dank.
