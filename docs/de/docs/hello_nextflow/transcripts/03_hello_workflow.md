# Teil 3: Hello Workflow - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../03_hello_workflow.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern aus dem Material.

## Willkommen und Zusammenfassung

Hallo und willkommen zurück zu Teil drei von Hello Nextflow. Dieser Teil heißt Hello Workflow, und in diesem Teil des Kurses werden wir wirklich anfangen, den Namen Pipeline oder Workflow zu rechtfertigen.

Wir werden unser einfaches Pipeline-Skript, das bisher nur einen Prozess hat, nehmen und weitere Prozesse hinzufügen. Wir werden sehen, wie Nextflow diese Orchestrierung und den Datenfluss durch die Pipeline handhabt.

Lass uns zurück zu unseren Codespaces gehen. Du siehst, ich habe alle meine .nextflow\*-Verzeichnisse und die work-Verzeichnisse und alles gelöscht, um es sauber zu halten. Mach dir keine Sorgen, wenn du diese Dateien noch von früheren Kursteilen herumliegen hast.

Wir werden mit einer Datei namens hello-workflow.nf arbeiten. Wie zuvor repräsentiert diese im Grunde das Skript, das wir bis zu diesem Punkt aufgebaut haben, und gibt uns einen sauberen Ausgangspunkt. Und wieder unten in der Ausgabe können wir sehen, dass der Pfad jetzt hello_workflow ist. Die veröffentlichten Dateien sollten also in ein anderes Unterverzeichnis in deinem results-Ordner gehen.

Um zusammenzufassen, wo wir bisher sind: Wir haben einen einzelnen Prozess hier, mit einer Eingabe greeting, einer Ausgabe greeting file. Und dann das einfache Bash-Skript, das nur einen echo-Befehl in eine Datei macht.

Wir haben eine einzelne Workflow-Eingabe, den params-Block hier, wo wir sagen, dass ein Pfad erwartet wird, und der Standard ist data/greetings.csv, was diese Datei hier oben ist.

Dann im Workflow selbst haben wir einen main-Block. Wir erstellen einen Kanal. Wir parsen das CSV in Zeilen und nehmen dann das erste Element jedes Arrays, und wir übergeben diesen Kanal an diesen Prozess, der dann drei Aufgaben generiert, und wir veröffentlichen aus dem Workflow die Ausgaben von diesem Prozess.

Und schließlich im output-Block sagen wir Nextflow, dass diese Dateien aus diesem Kanal in das Verzeichnis namens hello_workflow veröffentlicht werden sollen. Und diese Dateien zu kopieren, anstatt sie zu soft-linken.

## 1. Einen zweiten Schritt zum Workflow hinzufügen

Okay, in diesem Teil werden wir einen zweiten Prozess zu unserem Workflow hinzufügen. Wir nehmen die Ausgaben des sayHello-Prozesses und verarbeiten sie in einem zweiten Schritt, der alle Buchstaben innerhalb dieser Dateien in Großbuchstaben umwandelt - convertToUppercase.

Das ist nur ein albernes Beispiel, es ist wieder nur einfache String-Verarbeitung, aber es zeigt dir, wie wir die Logik innerhalb des Workflows nehmen können.

Wir werden dafür einen Bash-Befehl namens "tr" verwenden, was kurz für translate ist. Es ist ein Unix-Befehl, der schon ewig existiert. Wenn du nicht damit vertraut bist, kann ich das verstehen. Ich glaube nicht, dass ich ihn jemals vor dem Training verwendet habe, aber du kannst ihn sehr schnell im Terminal ausprobieren. Wenn ich "echo 'hello world'" mache und dann zu 'tr' pipe und dann in Anführungszeichen den Zeichenbereich angebe, also A bis Z, klein, und dann willst du A bis Z groß machen. Und es sagt einfach, übersetze diese Buchstaben in diese Buchstaben.

Und wenn ich Enter drücke, siehst du, dass jetzt alles großgeschrieben ist. Sehr schön, wenn du gerne Leute anschreist.

Das ist also ein sehr einfacher Stil von Bash-Befehl, den wir in unserem zweiten Prozess verwenden werden.

## 1.2. Den Großschreib-Schritt als Nextflow-Prozess schreiben

Wenn ich also zurück zu meinem Skript gehe, werde ich ein bisschen schummeln und einfach den Code aus den, aus den Docs für das Training kopieren. Aber du kannst genau sehen, was vor sich geht.

Wir haben hier einen neuen Prozess. Diesen haben wir convertToUpper genannt, aber wir könnten ihn nennen, wie wir wollen.

Wir haben einen einzelnen Eingabepfad, wie zuvor. Es ist kein value channel, es ist ein path channel. Und dann eine einzelne Ausgabe.

Im script-Block machen wir "cat" auf die Eingabedatei. Und wir können das in geschweifte Klammern setzen, wenn wir wollen. Und das nimmt diese Variable. Und wir führen denselben Bash-Befehl in der Pipe aus und schreiben die Ergebnisse in eine Datei mit diesem Dateinamen, und das wird vom output path aufgenommen.

Wir müssen jetzt etwas mit diesem neuen Prozess machen. Also gehen wir runter zum Workflow, wo wir die verschiedene Logik eines Workflows aufbauen, und nach diesem ersten Prozess werden wir unseren zweiten Prozess ausführen. Also convertToUpper ist hier der Name des Prozesses.

Er nimmt eine Eingabe, also können wir ihn nicht einfach so aufrufen. Wir wollen die Ausgabe des ersten Prozesses verarbeiten. Also genau wie wir das mit diesem sayHello out gemacht haben, wo wir diese Ergebnisse veröffentlichen. Wir wollen dieselben Ergebnisse hier als Eingabe verwenden, also können wir die kopieren und dort einfügen.

Wir wollen den sayHello-Prozess ".out", und Nextflow weiß, dass das einen einfachen einzelnen Ausgabedatensatz hier bedeutet, was diese Datei ist. Das wird dann als Eingabe an einen zweiten Prozess übergeben.

## 1.5. Die Workflow-Ausgabeveröffentlichung einrichten

Okay. Und schließlich, damit wir die Ergebnisse dieses zweiten Prozesses tatsächlich speichern, müssen wir sie auch aus dem Workflow veröffentlichen und dann im output-Block definieren, gleiche Syntax wie zuvor. Also können wir das kopieren und second outputs sagen, oder wie auch immer du es nennen willst.

Nimm den Prozessnamen, der uns interessiert, convertToUpper out, und dann hier unten im output-Block. Füge das hinzu und wir könnten dieselben Attribute hier machen. Also wollen wir auch diese Dateien im Hello Workflow-Unterverzeichnis, und wir wollen sie auch kopieren.

Großartig. Lass uns versuchen, es auszuführen. Wenn ich also das Terminal aufrufe und "nextflow run hello-workflow.nf" mache, werden wir sehen, was es tut. Schauen wir, ob es anders aussieht als in den vorherigen Teilen.

Also startet es Nextflow. In den Docs heißt es, das mit "-resume" zu machen, aber ich habe alle meine work-Verzeichnisse gelöscht, also hätte es hier keinen Unterschied gemacht. Aber wenn du es getan hast, dann wird das auch funktionieren.

Und es sieht fast genauso aus. Aber du kannst jetzt sehen, dass es eine zweite Ausgabezeile hier gibt, wo du den Namen des zweiten Prozesses sehen kannst, den wir gerade hinzugefügt haben. Und tatsächlich kannst du sehen, dass er dreimal erfolgreich gelaufen ist.

Brilliant. Wenn ich meine vorherigen work-Verzeichnisse noch hätte und das mit "-resume" gemacht hätte, wären diese zwischengespeichert worden, nur der erste Schritt in der Pipeline. Denn diese Ausgaben waren genau dieselben, also hätte Nextflow gewusst, diese wieder zu verwenden.

Und so kannst du sehen, wie du -resume verwenden kannst, um deinen Workflow iterativ aufzubauen, Schritt für Schritt, wenn du musst.

Okay, lass uns einen Blick in das results-Verzeichnis hier oben werfen und sehen, ob es funktioniert hat. Wir können sehen, dass wir hier oben einige weitere Dateien haben. Wir haben unsere ursprünglichen Dateien wie zuvor vom ersten Prozess. Und tatsächlich haben wir unsere upper-Dateien und die Buchstaben sind alle groß, also hat es funktioniert. Es ist wirklich schön zu sehen.

Es ist auch interessant, einfach in diese work-Verzeichnisse zu schauen. Wie zuvor entspricht der Hash hier den work-Verzeichnissen. Wenn ich also in "ls work" schaue und das dann erweitere, sehen wir die verschiedenen Dateien hier.

Wir sehen die Ausgabedatei vom ersten Prozess, die hier als Eingabe gezogen wurde. Und wir können die neue Ausgabedatei sehen, die generiert wurde.

Wenn ich das nun mit "-la" mache, um alle Dateien aufzulisten und zu zeigen, werden wir ein paar mehr Dinge sehen. Erstens siehst du, dass diese Datei tatsächlich ein Soft-Link zum ersten Prozess ist. Das ist grundsätzlich immer ein Soft-Link, wenn es sein kann, um Dateispeicher zu sparen. Wir veröffentlichen die Dateien hier nicht und es referenziert nur diese Datei von einer ersten Aufgabe in eine zweite Aufgabe, sodass alles innerhalb dieses einen Arbeitsverzeichnisses gekapselt ist und sicher und isoliert von allem anderen.

Und das muss dort sein, denn wenn wir uns die .command.sh-Datei ansehen, wenn ich also "cat work/b8/56\*" mache, kannst du sehen, dass die Dateiteile hier relativ sind, also catted es diese Eingabedatei, die ins selbe Arbeitsverzeichnis soft-gelinkt wurde.

So wird also jedes work-Verzeichnis aussehen. Wenn du es in Nextflow anschaust, hast du alle Eingabedateien dort in dieses Arbeitsverzeichnis staged. Und dann hast du auch alle Ausgabedateien, die erstellt wurden. Das ist also großartig. Das sieht so aus, wie wir es erwarten.

## 2.1. Den Sammelbefehl definieren und im Terminal testen

Okay, lass uns zurück zu unserem Workflow gehen. Was ist der nächste Schritt, den wir machen wollen?

Wir haben jetzt zwei Prozesse und sie nehmen diese eine CSV-Datei, parsen sie und teilen sie auf. Und dann haben wir drei Aufgaben für jeden dieser Prozesse und Nextflow handhabt die Parallelisierung von all dem, sodass alles nebeneinander läuft, wo möglich.

Diese Art, Arbeit aufzuteilen, um Dinge parallel laufen zu lassen, ist sehr üblich. Und das Gegenteil davon ist dann, alles wieder zusammenzuführen. Das ist also, was wir mit unserem finalen Prozess im Workflow machen werden - wir haben einen dritten hier, der diese drei verschiedenen Ausgaben nimmt und sie alle in eine einzelne Datei kombiniert.

Wir können das ziemlich einfach in einem Terminal machen, nur um ein Gefühl dafür zu bekommen, wie das aussehen wird.

Wenn ich zum results-Ordner gehe. Also, "cd results/hello_workflow/", und wir haben hier alle UPPER-Dateien. Ich kann einfach "cat" verwenden, was wir benutzen, um den Inhalt dieser Datei zu drucken, und du kannst mehrere Dateien an "cat" geben und es liest eine nach der anderen.

Also kann ich "UPPER-\*" sagen, was mir dieselbe Liste von drei Dateinamen mit Bash-Expansion gibt. Und ich kann combined.txt sagen. Ich glaube in den Docs werden die genauen Dateinamen aufgelistet, aber es macht dasselbe.

Jetzt, wenn ich "cat combined.txt" verwende, können wir sehen, dass wir den Dateiinhalt von allen drei dieser Dateien haben.

Das ist also im Grunde alles, was dieser Prozess tun wird - wir werden versuchen, ihm alle verschiedenen Ausgabedateien von einem vorherigen Prozess in einer einzelnen Prozessaufgabe zu geben, und dann werden wir sie zusammen "catten" und die Ausgabedatei speichern.

## 2.2. Einen neuen Prozess für den Sammelschritt erstellen

Okay, also lass uns unseren neuen Prozess hinzufügen. Ich werde das aus dem Trainingsmaterial einfügen, und du kannst sehen, dass es uns hier ein bisschen als Übung für den Leser mit diesen Fragezeichen überlassen hat. Aber du kannst den allgemeinen Umriss des Prozesses sehen - im Grunde das, was wir gerade im Terminal gemacht haben, wo wir "cat" von einem Haufen Eingabedateien machen und es in eine Ausgabedatei hier namens collected schreiben, und dann erwartet die Ausgabe wieder diesen einzelnen Pfad.

Wir brauchen also hier irgendeine Eingabe und es wird ein Set von Pfaden sein. Also wieder definieren wir einen input path channel und nennen wir ihn input_files. Jetzt hat uns das vorher einen einzelnen Pfad hier gegeben, aber ein Pfad kann auch mehrere Dateien hier haben, obwohl es immer noch eine einzelne Deklaration ist.

Ich werde das hier runterkopieren, weil wir diese Dateien "catten" wollen. Und du denkst vielleicht, dass wir hier einige Probleme haben mit dem Drucken eines Arrays oder solchen Sachen, aber Nextflow ist generell ziemlich vernünftig, wenn es darum geht. Und wenn es einem Kanal mit mehreren Dateien darin wie diesem gegeben wird, wird es sie alle mit Leerzeichen-Trennzeichen zusammenfügen. Das gibt uns also die korrekte Syntax.

Das ist großartig. Also lass uns jetzt unseren neuen Prozess verdrahten. Ich gehe runter zum Workflow. Ich werde sagen, kombiniere die Ausgaben, der neue Prozessname, und genau wie zuvor. Ich nehme diesen vorherigen Prozess, convertToUpper und mache ".out".

Großartig. Lass es uns ausprobieren und sehen, ob es im Terminal funktioniert. Wenn ich einfach ein paar Verzeichnisse hochgehe und dann den Nextflow-Befehl erneut ausführe, werden wir sehen, was passiert.

Also der Workflow wurde gestartet und jetzt kannst du sehen, dass wir drei verschiedene Prozessnamen haben, was großartig ist. Die ersten beiden sehen beide gleich aus wie zuvor, und der dritte neue läuft, was gut ist.

Allerdings gibt es hier etwas Seltsames. Wir wollten diese Ausgabedateien in eine einzelne Datei kombinieren, und doch können wir sehen, dass dieser Prozess dreimal gelaufen ist, nicht einmal.

Tatsächlich, wenn wir in eines dieser work-Verzeichnisse gehen. Und machen "cat work/" "collected", dann werden wir sehen. Es gibt nur ein einzelnes Wort hier drin, nicht drei.

Und was also passiert ist, ist dass Nextflow diese Parallelisierung fortgesetzt hat, genau wie es in den vorherigen Schritten getan hat. Und dieser Prozess gab uns einen Kanal mit drei Elementen, und diese drei Kanalelemente wurden an unseren nachgelagerten Prozess übergeben, der drei Prozessaufgaben generierte.

Er hat im Grunde versucht, dreimal separat zu sammeln, und jedes Mal hatte er nur eine einzelne Datei, also hat er einfach cat einzelne Datei zu einer Ausgabe gemacht, und tatsächlich können wir das auch in der .command.sh-Datei sehen.

Wenn ich .command.sh mache, können wir sehen, dass es nur einen einzelnen Dateinamen hier hat und nur eine einzelne Datei in dieses Arbeitsverzeichnis staged wurde.

## 2.3. Den Sammelschritt zum Workflow hinzufügen

Also müssen wir Nextflow irgendwie sagen, dass es alle diese Ausgaben von einem vorherigen Prozess zusammenbringt und sie diesem nachgelagerten Prozess als ein einzelnes Kanalelement gibt, anstatt drei.

Wir machen das mit einem Kanal-Operator namens _collect_.

Das ist ein super nützlicher Operator, den du in Nextflow-Pipelines die ganze Zeit sehen wirst. Das ist hier ein Kanal, dieser Ausgabekanal, genau wie der, den wir oben erstellt haben. Und so können wir Kanal-Operatoren daran anhängen, genau wie wir es zuvor getan haben. Wir können einfach dot machen, und dann in diesem Fall collect, Klammern.

Und das ist alles, was wir brauchen. Das wird dann diesen Kanal manipulieren, bevor er an diesen Prozess übergeben wird.

Wenn du sehen willst, was damit passiert, können wir es auch hier viewen. Also hier, das ist überhaupt nicht damit verbunden, diesen Prozess auszuführen, also könnte ich es an jedem Punkt nach Ausführung dieses Prozesses setzen. Aber wir nehmen denselben Ausgabekanal, und wir schauen ihn mit .view an, und dann schauen wir ihn wieder mit .collect.view an.

Und wenn wir das ausführen, wird es uns die zwei verschiedenen Strukturen dieses Kanals zeigen, vor und nach collect. Also lass uns das jetzt versuchen. Okay, ich habe nur ein bisschen rausgezoomt, weil einige der Ausgaben ziemlich lang sind, aber wenn ich die Pipeline ausführe, werden wir sehen, ob es funktioniert.

Ich hoffe, ein dritter Prozess wird nur einmal laufen, weil er die Ausgaben sammelt und tatsächlich kannst du sehen, collectGreetings als eins von eins. Also ist das nur eine Aufgabe gelaufen.

Und dann, wenn wir uns die view-Statements ansehen, haben wir drei view-Statements für die drei Elemente von vorher, mit einem Dateipfad in jedem.

Und dann nach diesem collect-Statement wurde das nur einmal ausgelöst, weil es ein einzelnes Element in diesem Kanal gibt. Und jetzt haben wir diese Liste von drei verschiedenen Dateipfaden.

Das ist genau das, was wir erhofft haben. Und du kannst sehen, hoffentlich, das ist im Grunde das Gegenteil von diesem "map"-Operator, den wir gemacht haben, um von den CSV-Arrays in separate Kanalelemente zu gehen. Jetzt nehmen wir separate Kanalelemente und setzen sie zurück in ein einzelnes Array.

Großartig, wir können diese view-Statements aufräumen. Wir brauchen diese nicht mehr. Wir können zum nächsten Schritt übergehen.

Bevor ich weitermache, und bevor ich es vergesse, werde ich hier ein neues publish-Statement hinzufügen. Third output. Du kannst das in deinem Workflow etwas semantischer und beschreibender nennen. Und dann werde ich das wieder zum output-Block hinzufügen und path 'hello_workflow' mode 'copy' sagen. Nur damit die von diesem Prozess generierte Ausgabedatei in unserem results-Ordner hier oben gespeichert wird.

Nur um schnell zu überprüfen, dass das funktioniert. Sollte jetzt ein bisschen sauberer sein, weil wir diese view-Statements nicht haben. Und wir werden sehen, ob wir unsere neue Ausgabedatei hier oben bekommen. Eins von einer Aufgabe lief, bekam eine neue Datei namens collected, und jetzt haben wir alle drei dieser Wörter. Fantastisch. Was kommt als nächstes?

## 3. Zusätzliche Parameter an einen Prozess übergeben

Okay. Als nächstes werden wir uns ansehen, wie man mehrere Eingaben in einen einzelnen Prozess handhabt. Bisher kannst du sehen, dass alle unsere Prozesse nur eine Sache als Eingabe nehmen. Sie alle haben eine einzelne Zeile unter ihrer Eingabe.

Wir werden das demonstrieren, indem wir Nextflow erlauben, eine andere Batch-Kennung anzugeben, sodass du vielleicht diesen Workflow mehrmals ausführst und ihm jedes Mal eine andere Batch-ID geben kannst.

Ich werde einfach eine zweite Zeile in der Eingabe hier für collectGreetings hinzufügen. Und ich werde es "val" nennen, weil das ein String ist. Jetzt ist es ein Wert, kein Pfad, und ich werde es "batch_name" nennen.

Dann werde ich das Skript hier unten bearbeiten, um diese Variable zu verwenden, und ich werde versuchen, es an derselben Stelle wie im Trainingsmaterial zu setzen. Also setze ich es in die Mitte dieses Dateipfads COLLECTED-$\{batch_name\}-output.

Noch nicht ganz fertig. Denk daran, dass wir Nextflow sagen müssen, wie die Ausgabedateinamen sein werden. Also müssen wir auch dasselbe hier oben machen: COLLECTED-$\{batch_name\}-output.txt".

Fantastisch. Nextflow bekommt jetzt eine zweite Variableneingabe und interpoliert diese in das Skript und die Ausgabe.

Eine letzte Sache, wir müssen jetzt finden, wo das aufgerufen wird, und wir müssen die zweite Eingabe an den Prozess übergeben. Das ist genau wie jede andere Eingabe in eine Funktion in jeder anderen Sprache.

Genau wie wir es früher im Training gemacht haben, werde ich hier das spezielle "params" verwenden, und wir werden es "params.batch" nennen, sodass wir eine --batch CLI-Option haben können. Und jetzt kannst du sehen, dass unser Prozess hier zwei separate Eingaben hat, nur durch Komma getrennt, die übergeben werden.

Es ist wirklich wichtig, die Reihenfolge richtig zu bekommen, also die Reihenfolge der Argumente hier für channel und dann den param muss übereinstimmen. Der Kanal und der batch name dort. Das ist nur positionelles Matching.

Okay. Ich kann diese Pipeline jetzt direkt mit --batch ausführen, aber lass uns zuerst das Richtige tun und es in der Eingabe hier in Params definieren. Also werde ich es zu batch hinzufügen und dann werden wir sagen, es ist ein String und lass uns einen Standard geben. Also nennen wir es einfach batch. Okay? Jetzt lass uns versuchen, den Workflow auszuführen.

--batch Trio. Ich glaube, das steht im Trainingsmaterial, aber wir könnten dort jeden String verwenden, den wir wollen. Und hoffentlich werden wir diese Ergebnisausgabedatei hier oben sehen.

Und tatsächlich, COLLECTED-trio-output - das hat richtig funktioniert. Es hat unsere Datei umbenannt. Und du kannst dir jetzt vorstellen, das ist nützlich, denn wenn ich das nochmal mit einem anderen Batch-Namen ausführe, wie replicate_two, dann wird es uns hier oben einen anderen Batch-Namen geben.

Und und es wird dann die Ausgabedateien in diesem Fall nicht überschreiben. Das ist also schön.

## 4. Eine Ausgabe zum Sammelschritt hinzufügen

Okay, wir haben jetzt also mehrere Eingaben für unseren Prozess hier. Aber was passiert, wenn wir mehrere Ausgaben erstellen wollen? Unser Beispiel hier ist dann, dass wir einen Report für diesen Prozess erstellen werden, der nur sagt, so viele Dateien wurden gesammelt.

Und wir werden das mit einem echo-Befehl hier machen. Also können wir echo sagen. Es gab, ich werde das vom Trainingsmaterial kopieren, damit du mir nicht beim Tippen zusehen musst.

Es gab $\{count_greetings\} Begrüßungen in diesem Batch, und speichere das in eine neue Datei jetzt namens $\{batch_name\}, also dieselbe Variable, wir können die so oft wiederverwenden, wie wir wollen, report.txt.

## 4.1.1. Die Anzahl der gesammelten Begrüßungen zählen

Wir müssen das irgendwie berechnen. Wir könnten diese Logik im Bash-Skript machen, wenn wir wollten, mit Bash-Logik. Allerdings können wir auch direkt innerhalb des Nextflow-Codes skripten, solange es innerhalb des script-Blocks im Prozess ist und über dem zitierten Abschnitt.

Alles hier wird nicht im finalen gerenderten Skript enthalten sein, und es wird einfach von Nextflow ausgeführt, wenn es eine Aufgabe rendert.

Also machen wir hier einfach etwas Logik. Wir erstellen eine neue Variable namens count_greetings. Wir nehmen den input files channel hier, und wir rufen .size() darauf auf.

Okay, diese Funktion wird mir hier eine Zahl in diese Variable geben, und jetzt ist unsere Warnung weg, weil diese Variable definiert wird.

Okay, wir erstellen also diese zweite Datei im work-Verzeichnis, aber wir müssen Nextflow sagen, dass es sie als veröffentlichte Ausgabe dieses Prozesses erwarten soll. Wir machen das mit genau derselben Syntax, wie wir es für die erste Datei gemacht haben.

Wir sagen path, denn es ist, wieder, wir könnten hier eine Variable veröffentlichen, wenn wir wollten mit "val", aber wir werden "path" sagen. Und dann der erwartete Dateiname. Beachte, es ist hier nicht hervorgehoben. Das liegt daran, dass ich einfache Anführungszeichen verwendet habe. Ich muss doppelte Anführungszeichen verwenden.

## 4.1.2. Die Report-Datei ausgeben und Ausgaben benennen

Okay, das ist großartig. Und wir könnten jetzt anfangen, auf diese Ausgaben hier unten zuzugreifen, genau wie ich es hier gemacht habe. Aber jetzt ist es ein Array von verschiedenen Objekten, also könnte ich collectGreetings.out[0] machen, um die erste zu bekommen, oder eins, um die zweite zu bekommen, was unser neuer Report ist.

Aber ich mag das nicht wirklich sehr, weil es ziemlich einfach ist, die Index-Zählung durcheinander zu bringen. Und du sitzt da und zählst viel Zeilen und du fügst eine neue Ausgabe hinzu und plötzlich bricht alles zusammen. Also

ist es viel schöner, stattdessen alles per Name zu referenzieren. Und wir können das mit einem speziellen Schlüssel hier namens "emit" machen.

Wir können das also nennen, wie wir wollen. Nennen wir es emit outfile, und emit reports. Wenn du diese definierst und du kannst es bei einem oder vielen machen, das liegt an dir. Jetzt kann ich hier runtergehen und stattdessen kann ich dot out dot reports gehen und es einfach per Name aufrufen, was viel einfacher ist, deinen Code zu verstehen, wenn du ihn liest, und es ist sicherer bei Änderungen im Code.

Ich habe das .out.report hier hinzugefügt, aber eigentlich muss ich zwei verschiedene Ausgaben veröffentlichen lassen. Also werde ich es als etwas Interessanteres umbenennen wie collected und report und ist das, was ich es genannt habe? Ich nannte es out file, sorry. Also dieser emit-Name hier outfile und report. Denn wir veröffentlichen zwei verschiedene Ausgabekanäle und so müssen wir beide im publish-Block referenzieren.

Dann müssen wir diese auch im output-Block definieren. Also habe ich das collected umbenannt, und wieder für reports, ein bisschen verbose hier, aber es ist wirklich nützlich, wenn du reinkommst, um einen neuen Workflow zu lesen, alle verschiedenen Ausgaben hier zu sehen, alle verschiedenen Kanäle nebeneinander aufgelistet, und es gibt Wege, das weniger verbose zu machen, was wir später berühren werden.

Okay, lass es uns ausprobieren und unseren Workflow ausführen und sehen, was passiert.

Hoffentlich sollte es jetzt im Grunde genauso laufen wie zuvor. Und wir werden eine neue Ausgabedatei hier oben bekommen namens replicate_two, report. Und da haben wir es. Es ist geöffnet und es sagt, es gibt drei Begrüßungen im Batch, was wir erwartet haben, also ist es perfekt.

Wenn ich in das work-Verzeichnis hier gehe, nur um dir zu beweisen, dass es im Nextflow-Code ausgeführt wurde, anstatt im Bash-Skript, kann ich zu cat work/ command.sh gehen, und du wirst hier sehen, dass es nur diesen String direkt echoed. Es gab drei Begrüßungen in diesem Batch, und so wurde diese Variable von Nextflow interpoliert. Sie wurde im script-Block berechnet, bevor es die .command.sh-Datei geschrieben hat. Also ist die resultierende Variablenberechnung im Grunde in diese hartcodiert, bevor sie auf deiner Compute-Umgebung in diesem Fall ausgeführt wird.

Und so kannst du diese Trennung zwischen dem script-Block hier und allem darüber sehen. Ich hoffe, das macht Sinn.

## Fazit und Quiz

Okay, das ist das Ende dieses Teils von Hello Nextflow. Also wie zuvor, geh und schau dir das Quiz an. Mach es auf der Webseite oder in der CLI, geh durch einige der Fragen und überprüfe einfach, ob du einiges vom Material verstanden hast, das wir behandelt haben. Schau, ob es dort etwas gibt, das etwas hervorhebt, das du vielleicht nicht verstanden hast. Nicht zu viele Fragen. Schön und einfach zu machen. Oder du kannst es auch hier unten auf der Webseite machen.

Und mach eine kleine Pause, einen kleinen Spaziergang und komm zurück und schließ dich uns in Teil vier von Hello Nextflow an, wo wir über Module sprechen werden. Vielen Dank.
