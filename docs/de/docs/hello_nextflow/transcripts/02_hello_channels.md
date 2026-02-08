# Teil 2: Hello Channels - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../02_hello_channels.md) zurück.

    Die im Transkript gezeigten Abschnittsnummern dienen nur zur Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus dem Material.

## Willkommen

Hallo und willkommen zurück zu Teil 2 von Hello Nextflow. Dieses Kapitel heißt Hello Channels.

Kanäle sind wie der Klebstoff in deiner Nextflow-Pipeline. Sie sind die Teile, die all die verschiedenen Prozesse zusammenhalten, die Nextflow nutzt, um alle Informationen weiterzugeben und deinen Workflow zu orchestrieren.

Es gibt noch einen weiteren Teil von Kanälen, das sind Operatoren. Das sind im Grunde Funktionen, die wir auf Kanälen verwenden können, um die Inhalte zu ändern. Lass uns in VS Code eintauchen und sehen, wo wir stehen.

Ich bin sehr nah rangezoomt an diesem VS Code. Um die Dinge sauber und aufgeräumt zu halten, habe ich alle _.nextflow\*_-Dateien und das _work/_-Verzeichnis und die results/ und alles aus Kapitel Eins entfernt. Und ich fange hier einfach frisch an. Aber mach dir keine zu großen Sorgen. Wenn du nicht möchtest, kannst du diese Dateien herumliegen lassen. Sie werden keine Probleme verursachen.

Wir werden für dieses Kapitel mit _hello-channels.nf_ arbeiten, und wenn ich das öffne, sollte es der Datei sehr ähnlich aussehen, an der wir zuvor gearbeitet haben. Es kann sein, dass verschiedene Teile an verschiedenen Stellen des Skripts sind, aber alles sollte im Grunde gleich sein.

Eine Sache, die anders ist, ist dass der Pfad im Ausgabeblock hier jetzt _hello_channels_ für diesen Teil ist, was bedeutet, dass die Ergebnisdateien in einem anderen Unterverzeichnis in deinen results gespeichert werden, falls du das dort noch hast. Es sollte also ein schöner und sauberer Ort zum Starten sein, ohne durch Ausgaben verwirrt zu werden.

Okay, lass uns kurz daran erinnern, was dieses Skript tut, wenn wir diesen Workflow ausführen. Wir machen _"nextflow run hello-channels.nf"_. Wir können _"--input myinput"_ machen, und wenn wir das ausführen, wird es diesen Parameter params.input verwenden, der als Variable für den sayHello-Prozess hier oben übergeben wurde, der in greeting geht und in output.txt gespeichert wird. Und das können wir in der Ergebnisdatei sehen. Großartig.

## 1. Variable Eingaben explizit über einen Kanal bereitstellen

Das ist schön. Aber es ist, es ist ziemlich simpel. Wir haben eine Variable in diesem Parameter, die in einen Prozess geht, der einmal läuft und nicht wirklich skaliert. Und wir können ihm nicht viele verschiedene Dateien zum Erstellen geben. Wir können ihm nicht viele verschiedene Begrüßungen geben. Wir haben nur eine.

In Wirklichkeit geht es bei Nextflow darum, deine Analyse zu skalieren. Du möchtest also wahrscheinlich, dass es mehr als eine Sache tut. Und das tun wir mit _Kanälen_.

Kanäle sind für viele Leute, die Nextflow aufgreifen, ein ziemlich einzigartiges Konzept. Es kommt aus diesen Konzepten der funktionalen Programmierung, und es kann ein bisschen Zeit brauchen, um den Kopf darum zu bekommen, aber sobald es klick macht, entfesseln sie wirklich die Kraft von Nextflow und es ist der Schlüssel dazu, wie du deine Workflows schreibst.

## 1.1. Einen Eingabekanal erstellen

Lass uns damit anfangen, dieses Skript zu nehmen und es einen _Kanal_ verwenden zu lassen, anstatt nur einen _Parameter_.

Wir gehen runter zum Workflow, wo all unsere Workflow-Logik ist, um Dinge zusammenzufügen. Und ich werde hier reingehen und ich werde einen neuen Kanal erstellen.

Erstelle einen neuen Kanal.

Und ich werde ihn "_greeting_ch_" nennen. Das ist eine Konvention, "_\_ch_" so zu machen, einfach damit du dich erinnern kannst, dass diese Variable ein Kanal ist. Aber du kannst ihn nennen, wie du möchtest.

Und dann werde ich sagen equals, und ich werde _"Channel.of"_ machen.

Channel ist wie der Namensraum für alles, was mit Kanälen zu tun hat. Kleinbuchstabe "c", falls du Nextflow schon vorher benutzt hast. Und das _".of"_ ist etwas, das Channel Factory genannt wird, was im Grunde eine Möglichkeit ist, einen Kanal zu erstellen.

Es gibt viele verschiedene Channel Factories. Wenn ich hier einfach "." mache, kannst du sehen, dass VS Code viele davon vorschlägt, aber _".of"_ ist die einfachste und nimmt einfach eine Eingabe hier.

Also kann ich Klammern machen und ich werde sagen _"Hello Channels!"_.

Großartig. Ich habe einen Kanal. Fantastisch. Ich kann speichern drücken, ich könnte es wieder ausführen, aber nichts Interessantes wird passieren. VS Code hat mir eine orangefarbene Warnlinie hier unten gegeben und mir gesagt, dass das eingerichtet ist: Du hast das erstellt, aber du hast es nie wirklich für irgendetwas benutzt. Dieser Kanal wird nicht konsumiert.

Okay, wie benutzen wir ihn? Ganz einfach. Ich werde das kopieren, und ich werde _params.input_ löschen und ich werde stattdessen _"greeting_ch"_ hier einfügen. Also werden wir diesen Kanal als Eingabe an sayHello übergeben.

Beachte, dass ich diesen String vorerst fest codiert habe. Das ist ein bisschen ein Rückschritt nach unserem schönen Parameter, den wir am Ende des letzten Kapitels verwendet haben, aber es hält die Dinge einfach, um zu starten, damit du die Logik sehen kannst.

Okay, ich werde in mein Terminal gehen und ich werde den Workflow wieder ausführen. Dieses Mal ohne _"--input"_, und es wird laufen und es wird diesen Kanal verwenden, den wir erstellt haben, und hoffentlich sollten wir eine Datei hier oben in _results/hello_channels/_ haben und es steht jetzt "Hello Channels!" drin. Fantastisch. Das ist also das, was wir von unserem Kanal hier erwarten. Großartig.

## 1.4. view() verwenden, um die Kanalinhalte zu inspizieren

Eine weitere Sache, die ich hier hinzufügen möchte, ist nur eine kurze Einführung in eine weitere Funktion, die wir auf Kanälen verwenden können, genannt "_.view_".

Das ist analog zum _print_-Befehl in Python oder anderen Sprachen, die du vielleicht gewohnt bist, und es gibt einfach den Inhalt dieses Kanals ins Terminal aus, wenn wir es ausführen.

Also mache "_.view_", und dann, wenn ich den Workflow wieder ausführe, sollte es ins Terminal drucken, was der Inhalt dieses Kanals ist, zu dem Zeitpunkt, zu dem wir ihn erstellt haben.

Sicher genug, du kannst sehen, dass es hier ins Terminal gedruckt wurde. _"Hello Channels!"_.

Beachte, dass du diese Dinge über Zeilen brechen kannst, wenn du möchtest, und tatsächlich wird der automatische Nextflow-Formatierer versuchen, das für dich zu tun. Leerzeichen sind hier nicht wirklich wichtig, also kannst du diese Dinge eines nach dem anderen verketten.

## 2. Den Workflow so modifizieren, dass er mit mehreren Eingabewerten läuft

Okay, unser Kanal hat also eine Sache drin, was schön ist, aber es ist im Grunde das Gleiche wie vorher. Also lass es uns ein bisschen komplizierter machen. Lass uns ein paar mehr Sachen in unseren Kanal einfügen.

Die "_.of()_"-Channel Factory kann mehrere Elemente annehmen, also lass uns ein paar mehr schreiben. Wir machen _Hello, Bonjour, Hej_. Und dann können wir diesen Workflow wieder ausführen und sehen, was passiert.

Sollte wieder laufen. Und wir haben jetzt gedruckt. _"Hello", "Bonjour"_ und _"Hej"_ ins Terminal mit unserer view-Anweisung. Fantastisch.

## 2.1.2. Den Befehl ausführen und die Log-Ausgabe ansehen

Du könntest denken, dass wir an dieser Stelle fertig sind. Aber tatsächlich gibt es hier eine kleine Falle, die uns stolpern lassen wird. Wenn wir uns unsere Ausgabedatei hier ansehen. Du kannst sehen, dass _"Hello"_ drin steht, aber es hat keine der anderen Ausgaben. Tatsächlich ist es nur diese eine.

Wenn wir diesen Workflow mehrmals ausführen, könnten wir sogar sehen, dass es manchmal _"Bonjour"_ hat, manchmal _"Hej"_. Es ist ein bisschen zufällig.

Wenn wir uns das Terminal ansehen, können wir sehen, dass es dreimal gelaufen ist und wir können die verschiedenen view-Ausgaben sehen. Aber wenn ich zum work-Verzeichnis gehe, kann ich _"cat work"_ machen. Diesen Hash einfügen und das erweitern und _output.txt_. Du kannst sehen, dass diese Datei im work-Verzeichnis anders ist als das results-Verzeichnis, und diese hier ist _"Hej"_. Also funktioniert hier etwas nicht ganz richtig.

Und der Schlüssel ist, dass wir drei Tasks hatten, die gelaufen sind. Die Nextflow-Ausgabe versucht das zusammenzufassen, während die Verarbeitung weitergeht, damit sie nicht komplett dein gesamtes Terminal übernimmt, und dass ANSI Logging ANSI-Escape-Codes verwendet, hat im Grunde die anderen Tasks überschrieben. Also zeigt es dir nur den letzten, der zufällig aktualisiert wurde.

## 2.1.3. Den Befehl erneut mit der Option -ansi-log false ausführen

Es gibt ein paar Dinge, die wir tun können, um das tatsächlich ein bisschen besser zu verstehen. Wir können in das work-Verzeichnis selbst schauen und du kannst all die verschiedenen Work Dirs dort sehen, aber das ist ein bisschen verwirrend, weil es mit verschiedenen Nextflow-Ausführungen vermischt wird.

Oder wir können Nextflow sagen, keine ANSI-Escape-Codes zu verwenden.

Also wenn ich den Befehl wieder ausführe, aber dieses Mal sage ich _"-ansi-log false"_, um es auszuschalten, könnte ich auch die Umgebungsvariablen _$NO_COLOR_ oder _"$NXF_ANSI_LOG=false"_ verwenden. Dann verwendet es die Art altmodischeren Stil des Nextflow-Loggings ohne diese Escape-Codes. Es druckt einfach direkt ins Terminal ohne clevere Aktualisierungen.

Und jetzt können wir alle drei dieser Prozesse sehen, die gelaufen sind. Und jeder von ihnen hat seinen eigenen Task-Hash. Und wenn wir in diese work-Verzeichnisse gehen, werden wir die drei verschiedenen Begrüßungen sehen, die wir angegeben haben.

Das macht jetzt ein bisschen mehr Sinn. Hoffentlich verstehst du, dass Nextflow das getan hat, es war nur ein bisschen clever mit dem, was es dir im Terminal mit diesen work-Verzeichnissen gezeigt hat.

Allerdings ist das für ein Problem mit den work-Verzeichnissen behoben, aber es hat ein Problem mit der Ausgabedatei nicht behoben. Wir haben immer noch nur eine Ausgabedatei, die _"Hello"_ sagt.

## 2.2. Sicherstellen, dass die Ausgabedateinamen eindeutig sein werden

Um das zu verstehen, müssen wir zu unserem Workflow-Skript zurückgehen. Wir generieren unseren Kanal hier, wir übergeben ihn an unseren Prozess, und wenn wir uns den Prozess ansehen, schreiben wir die Begrüßung in eine Datei namens _"output.txt"_ und geben diese Ausgabedatei zurück an den Ausgabeblock hier unten, veröffentlichen sie.

Allerdings, jedes Mal, wenn dieser Prozess diese drei verschiedenen Tasks läuft. Sie alle generieren eine Datei namens _"output.txt"_, alle diese Ausgabedateien werden ins results-Verzeichnis veröffentlicht, und sie überschreiben sich alle gegenseitig. Also welche Ergebnisdatei du dort auch bekommst, ist nur die letzte, die generiert wurde, aber alle anderen wurden überschrieben. Das ist nicht wirklich das, was wir wollen.

## 2.2.1. Einen dynamischen Ausgabedateinamen konstruieren

Es gibt verschiedene Möglichkeiten, damit umzugehen, aber die einfachste ist im Moment, einfach verschiedene eindeutige Dateinamen zu erstellen. Also jedes Mal, wenn der Task mit einer anderen Begrüßung läuft, wird er eine andere Ausgabedatei generieren, die nicht mehr kollidiert, wenn sie veröffentlicht wird. Und dann bekommen wir drei eindeutige Ausgabedateien.

Wir tun das auf genau die gleiche Weise. Wir können diese Variable überall innerhalb des Skriptblocks verwenden und wir können sie mehrmals verwenden.

Also kann ich das hier einfügen, _"$\{greeting\}\_output.txt"_, und dann muss ich es auch hier oben einfügen, weil wir nicht mehr eine Datei namens output.txt erstellen. Wenn ich das also nicht aktualisiere, wird Nextflow mit einem Fehler abstürzen und sagen, dass es eine Datei erwartet hat, die nie generiert wurde.

Also muss ich das Gleiche dort machen und ich muss doppelte Anführungszeichen verwenden, nicht einfache Anführungszeichen, damit diese Variable verstanden wird.

Okay, lass es uns ausprobieren und sehen, ob es funktioniert hat. Wir werden den Workflow wieder ausführen. Hoffentlich wird es uns die drei verschiedenen Tasks innerhalb der drei verschiedenen work-Verzeichnisse zeigen. Und sicher genug, du kannst hier oben im results-Ordner auf der linken Seite sehen. Wir haben jetzt drei verschiedene Dateien mit drei verschiedenen Dateinamen und jede mit den verschiedenen Inhalten, die wir erwarten. Also überschreiben die Dateien sich nicht mehr gegenseitig, und alles ist da, wie wir es erwarten.

Das ist ein bisschen ein triviales Setup, das wir hier durchgemacht haben, aber es unterstreicht einige der Schlüsselkonzepte, die du verstehen musst, wie Datei-Publishing funktioniert, und einige der Dinge, in die du als Fallen geraten könntest. Also hoffentlich kannst du das in deinen eigenen Workflows vermeiden.

Es ist auch erwähnenswert, dass das, was wir hier getan haben, im echten Leben ein bisschen unpraktisch ist. Wir haben einige Eingabedaten genommen und wir verwenden diese Daten, aber wir benennen auch die Datei nach diesen Daten, was du normalerweise nicht tun kannst.

Also in echten reiferen Nextflow-Pipelines wirst du oft ein Meta-Objekt mit allen Metadaten, die mit einer bestimmten Probe verbunden sind, herumreichen. Du kannst dann dynamische Dateinamen basierend darauf erstellen, was viel praktischer ist.

Wenn du daran interessiert bist, wie man das mit Best Practices macht, gibt es eine Side Quest auf _training.nextflow.io_, die speziell über Metadaten und Meta Maps ist, also kannst du dort für mehr Details eintauchen.

## 3. Mehrere Eingaben über ein Array bereitstellen

Okay. Als Nächstes werden wir ein bisschen darüber erkunden, wie Kanäle strukturiert sind und wie sie sich von anderen Arten von Datenstrukturen in der Programmiersprache unterscheiden. Und ich werde ein bisschen darüber nachdenken, wie ich potenziell ein Array verwenden könnte, was ein vertrautes Konzept sein könnte, wenn du von anderen Sprachen kommst.

Kann ich ein Array in einem Kanal verwenden? Lass es uns versuchen. Ich werde ein Array erstellen, und ich habe das aus den Docs kopiert, _"greetings_array"_ und _"Hello", "Bonjour"_ und _"Holà"_. Und dann werde ich das hier anstelle meiner fest codierten Strings einfügen. Also werde ich "Channel.of" _"greetings_array"_ sagen, dieses Array in einen Kanal übergeben. Lass es uns versuchen.

Terminal hochholen und die Pipeline ausführen.

Okay. Du kannst sehen, dass die view-Anweisung hier unser Array wie erwartet gedruckt hat, aber dann all dieser rote Text, oder es wird nicht rot sein, wenn du immer noch _"-ansi-log"_ aus hast, aber all dieser rote Text sagt uns, dass etwas schief gelaufen ist.

Wir haben hier kein schönes grünes Häkchen mehr. Wir haben ein rotes Kreuz, und wenn ich das nur ein bisschen breiter mache, damit es leichter zu lesen ist, sagt uns Nextflow, was schief gelaufen ist.

Also lass uns das Abschnitt für Abschnitt aufschlüsseln. Es sagt, der Fehler wurde verursacht durch, und dann der Grund für den Fehler, welcher ist fehlende Ausgabedateien. Also im Grunde sagte dieser Ausgabeblock, dass diese Datei erstellt werden sollte und sie wurde es nicht. Als Nächstes sagt es, das ist der Befehl, der ausgeführt wurde. Also ist das im Grunde der Inhalt dieser _.command.sh_-Datei. So sah es aus, nachdem all diese Variablen eingefügt wurden.

Und du kannst hier sehen, unser echo-Befehl wurde tatsächlich nur einmal ausgeführt und er hat das gesamte Array verwendet, aber in einer String-Darstellung, was nicht wirklich das ist, was wir wollten.

Und dann ist der Befehl so beendet worden, und das war das work-Verzeichnis, wo wir hingehen und die Dateien sehen können, um ein bisschen mehr zu verstehen.

Okay. Was also passiert ist, war. Nextflow hat einfach dieses gesamte Array als einzelnes Kanalelement an den Prozess übergeben, was bedeutete, dass der Prozess nur einmal lief. Er hatte einen Task und er hat die Daten nicht in einer Struktur verwendet, die wir erwartet haben.

## 3.2. Einen Operator verwenden, um Kanalinhalte zu transformieren

Also müssen wir diesem Kanal zuerst etwas tun, bevor er verwendet werden kann. Und das bereitet die Bühne für die Verwendung von Operatoren, was spezielle Funktionen sind, die wir auf Kanälen verwenden können, um Kanalinhalte zu manipulieren.

In diesem Fall werden wir etwas verwenden, das _flatten_ genannt wird. Das wir am Ende des Kanals hier übergeben. Also erstellen wir den Kanal und dann führen wir _flatten_ aus. Und wieder, wenn wir darüber hovern, zeigt es uns die Dokumentation für diesen Befehl direkt in VS Code, was sehr hilfreich ist. Du kannst all diese Docs auch auf der Nextflow-Website finden, die Dokumentation.

Ich könnte diesen Code jetzt einfach ausführen und sehen, ob er funktioniert, aber es ist auch eine schöne Gelegenheit einzuführen, wie man dynamischen Code innerhalb von Operatoren und innerhalb von Nextflow-Code macht, die Closures genannt werden.

Also werde ich hier einen view-Befehl wieder einfügen, bevor wir _flatten_ ausführen. Und hier hat dieser diese geschweiften Klammern, was die dynamische Closure ist. Und es gibt nur etwas beliebigen Code hier drin, der ausgeführt wird, innerhalb des Kontexts eines view-Operators.

Hier sagt das, nimm die Begrüßung, was die Eingaben des view-Operators sind, und das ist hier. Ich könnte das nennen, wie ich wollte, ich könnte das _"foo"_ nennen und ich muss nur später darauf als _"foo"_ verweisen. Und dann sage ich damit, gib das zurück.

Und dann setze ich eine String zurück, der sagt vor dem flatten für eine Variable. Sehr einfach.

Ich werde jetzt einen weiteren davon hinzufügen, genau das Gleiche, aber ich werde sagen nach _flatten_.

Also was das macht, weil das in Sequenz läuft, wirst du sehen, wie der Kanal aussieht, bevor wir _flatten_ ausführen, und dann wieder nach wir _flatten_ ausführen.

Und dann wird dieser greeting-Kanal immer noch erstellt, also wird er immer noch an den Prozess übergeben. Und hoffentlich wird der Workflow jetzt laufen. Lass es uns ausprobieren.

Großartig. Also das Erste zuerst ist, dass die Pipeline dieses Mal nicht abgestürzt ist. Wir hatten drei Prozesse, die richtig liefen und wir haben ein kleines Häkchen. Und dann können wir sehen, dass unsere view-Anweisungen funktioniert haben.

Wir haben vor _flatten_, was dieses Array ist, das wir vorher vom Fehler gesehen haben, und dann haben wir dreimal das nach _flatten_ aufgerufen wurde, wo wir _"Hello", "Bonjour"_ haben, und all die anderen drei separaten Elemente im Array, die jetzt wie gehofft drei separate Elemente im Kanal sind.

Und du kannst sehen, dass der _view_-Operator dreimal ausgeführt wurde. Und das ist, weil dieser Kanal nach _flatten_ jetzt drei Elemente hat. Und also wird der Operator dreimal aufgerufen.

Sehr schnell würde ich nur erwähnen, dass als ich vorher Channel Factories erstellt habe, ich _"."_ gemacht habe, und dann sahen wir, dass es viele verschiedene Möglichkeiten gab, Kanäle zu erstellen, und eine davon heißt "_fromList_". Und das ist tatsächlich speziell dafür entwickelt, diese gleiche Operation zu tun. Also hätten wir einfach fromList greetings away machen können, und das wird funktionieren. Es ist etwas sauberer und schönere Syntax. Aber für die Zwecke dieser Demonstration wollten wir es ein bisschen schrittweiser machen, damit du sehen kannst, wie der Kanal manipuliert wird und wie verschiedene Operatoren ändern können, was im Inhalt eines Kanals ist.

## 4. Eingabewerte aus einer CSV-Datei lesen

Okay, wie können wir das ein bisschen realistischer machen? Du wirst wahrscheinlich nicht viel Code in deiner Nextflow-Pipeline mit fest codierten Arrays erstellen wollen. Du wirst wahrscheinlich die Daten von außen nehmen wollen, wenn du startest, und diese Daten werden mit ziemlicher Sicherheit in Dateien sein.

Also das Nächste, was wir tun werden, ist, dass wir das replizieren werden, aber anstatt die Daten von einem einzelnen CLI-Parameter oder von einem fest codierten String oder Array zu nehmen, werden wir es von einer Datei nehmen.

Also lass uns unser greetings away loswerden. Und jetzt werden wir diese Channel Factory wieder ändern. Ich habe gerade gesagt, es gibt eine Menge zur Auswahl und es gibt eine namens _".fromPath"_. Und ich werde ihr sagen, dass sie in diesem Fall _params.input_ nimmt, was zu unserer Eingabe zurückgeht, die wir früher verwendet haben.

Jetzt ist dieser Parameter noch nicht wirklich bereit, verwendet zu werden. Wir sagen immer noch, dass es ein String ist und er ist hier mit einem Standard fest codiert, aber wir könnten diesen String überschreiben. Wir wollen jetzt, dass das eine Datei ist. Also ist der Typ anders. Es ist nicht mehr ein _String_. Es ist ein _Path_.

Und dann können wir den Standard setzen, wenn wir wollen, wieder zu einem Path. Und wenn ich in explore auf der linken Seite schaue, kannst du sehen, in diesem Repository, in diesem Arbeitsverzeichnis, habe ich ein Verzeichnis namens data. Ich habe eine Datei dort namens _"greetings.csv"_.

Also kann ich einfach den Standard hier setzen zu _"data/greetings.csv"_. Jetzt, wenn ich diese Pipeline wieder ohne irgendwelche Befehlszeilenoptionen ausführe, wird sie diesen Standardwert verwenden. Sie weiß, dass es ein Path ist, also weiß sie, dass sie das als Path behandeln sollte und nicht als String.

Und dann wird sie das in eine Channel Factory von diesem _params.input_ übergeben und unseren Kanal erstellen, der dann in diesem Prozess namens _sayHello_ verwendet wird. Lass es uns ausprobieren.

Okay. Fehlgeschlagen. Keine Sorge. Das war erwartet. Und wenn du dem Trainingsmaterial folgst, wirst du sehen, dass es dort auch erwartet wurde. Lass uns sehen, was hier passiert.

Es hat versucht, die Pipeline auszuführen. Es hat versucht, den Prozess auszuführen, und es hat einen ziemlich ähnlichen Fehler bekommen wie den, den wir vorher gesehen haben.

Hier sagt es: Wir haben versucht, _echo_ auszuführen, aber anstatt den Inhalt dieser CSV-Datei zu echoen, hat es nur den Pfad geechoed. Und du kannst sehen, es ist der volle absolute Pfad hier zu dieser CSV-Datei.

Und dann sicher genug, weil es versucht hat, das in diesen wirklich komplizierten Pfad zu schreiben, wusste es nicht wirklich, was es tun sollte. Und es war außerhalb des Geltungsbereichs des Prozess-work-Verzeichnisses.

Ich habe am Anfang erwähnt, dass Nextflow jeden ausgeführten Task in einem speziellen work-Verzeichnis kapselt. Und wenn du versuchst, in Daten zu schreiben, die außerhalb dieses work-Verzeichnisses sind, wird Nextflow dich als Sicherheitsvorkehrung stoppen. Und das ist hier passiert. Wir haben versucht, in einen absoluten Pfad zu schreiben und Nextflow ist fehlgeschlagen und hat uns verhindert.

## 4.2. Den splitCsv()-Operator verwenden, um die Datei zu parsen

Okay, lass uns diesen Kanal ansehen und sehen, wie er aussieht. Wir können _".view"_ machen, und ich habe das von der Website kopiert. Also _.view_, und wir haben eine dynamische Closure hier und wir sagen einen Variablennamen "_csv_" als Eingabe. Also das sind die Kanalinhalte, und wir sagen vor splitCsv, und so sieht es aus.

Wenn ich es wieder ausführe, wird es immer noch fehlschlagen, aber es wird uns zeigen, was in diesem Kanal ist. Es ist nicht besonders aufregend. Es ist diese _Path_-Variable. Du kannst also sehen, es ist nur ein String hier, weil er ins Terminal gedruckt wird, aber es ist ein _Path_-Objekt, das die Informationen und Metadaten über diese Datei enthält.

Wir wollen nicht die Metadaten der Datei an die Eingabe übergeben. Wir wollen den Inhalt dieser Datei übergeben. Wenn wir uns die _greetings.csv_-Datei ansehen, kannst du hier sehen, dass sie diese verschiedenen Variablen hier hat. _Hello, Bonjour, Holà_ wieder. Und das sind die Dinge, die wir wirklich an unseren Prozess übergeben wollen, nicht nur die Datei selbst als einzelnes Objekt.

Also müssen wir diese CSV-Datei parsen. Wir müssen sie auspacken, an den Inhalt der CSV-Datei gelangen und dann den Inhalt innerhalb des Kanals an den Prozess übergeben.

Wie du wahrscheinlich von der Log-Nachricht erkennen kannst, wollen wir _splitCsv_ verwenden, was ein weiterer Operator ist, ein weiterer Kanaloperator. Also wenn ich "_dot_" "_s_" mache, und dann kannst du sehen, es wurde automatisch vorgeschlagen. Ups, _splitCsv_ und ein paar Klammern.

Und dann nach _splitCsv_ werde ich eine weitere _view_-Anweisung einfügen, nur damit wir sehen können, wie es danach aussieht. Lass uns die Pipeline ausführen und sehen, was wir bekommen haben.

Okay. Es ist immer noch fehlgeschlagen, aber auf eine neue und aufregende Weise, was Fortschritt ist.

Dieses Mal haben wir wieder ein Problem mit unserem Skript, das gerendert wurde. Jetzt haben wir nicht mehr den finalen Pfad, aber wir haben ein Array von Variablen, was sehr ähnlich dem Fehler aussieht, den wir früher hatten, als wir ein Array als feste Eingabe übergeben haben.

Mit unserem Logging vom view-Operator können wir sehen, vor _splitCsv_ war der Pfad. Und sicher genug, nach _splitCsv_ haben wir drei verschiedene Ausgaben und jede dieser Ausgaben sieht sehr ähnlich aus wie jede der Zeilen aus der _greetings.csv_-Datei, was Sinn macht.

Also was hier passiert ist, ist dass Nextflow diese CSV-Datei geparst hat und uns drei Objekte gegeben hat, ein Array für jede Zeile der CSV-Datei. Also dann haben wir dreimal ein Array von Variablen an den Kanal übergeben anstelle eines einzelnen String-Werts.

Okay, also letztes Mal hatten wir dieses Problem, haben wir _flatten_ verwendet. Lass uns das sehr schnell versuchen. Versuche flatten und sieh, was passiert.

Ich kann diese Variablen nennen, wie auch immer. Also werde ich es _myarray_ nennen, weil es nicht mehr wirklich ein CSV ist. Lass uns versuchen, es wieder auszuführen und sehen, was mit _flatten_ passiert.

Also dieses Mal werden wir ausführen, wir haben das CSV in drei Array-Objekte geparst, und dann haben wir es geflattened. Und dieses Mal wurde es übergeben. Und die Nextflow-Pipeline lief. Allerdings kannst du sehen, dass _flatten_ wirklich zur Sache geht und alles flattened. Und also bekommen wir drei unabhängige Array-Einträge für jede Zeile. Und also lief der Prozess dreimal für jede Zeile eines CSV. Und jetzt haben wir eine ganze Menge Ergebnisdateien, und 123, 456, und alle möglichen Dinge, nicht nur diese erste Spalte des CSV, was wir wirklich wollten.

## 4.3. Den map()-Operator verwenden, um die Begrüßungen zu extrahieren

Wie also kommen wir an nur die erste Spalte? Wenn flatten hier zu simpel ist, brauchen wir einen komplexeren Operator, wo wir tatsächlich anpassen und ihm sagen können, was wir vom CSV wollen.

Um das zu tun, werden wir _map_ verwenden. Im Grunde sagt _map_ einfach, führe etwas Code, eine Funktion aus über jedes Element, das ich bekomme, und mache eine Art Transformation darauf. Und weil es so flexibel ist, wirst du sehen, dass es ständig in Nextflow-Code auftaucht.

Von sich aus tut es nichts. Also wollen wir keine regulären Klammern, wir wollen hier eine Closure und wir müssen ihr sagen, was sie tun soll. Also werde ich sagen _"row"_, weil das Zeilen vom CSV bekommt, also ist es ein logischer Variablenname. Ist die Eingabe. Und ich will nur das erste Element dieses Arrays zurückgeben.

Arrays in Nextflow sind nullbasiert, also werden wir sagen, nur das erste Element, was Zeile Null ist. Wenn wir die zweite Spalte wollten, könnte ich eins sein oder die dritte Spalte zwei sein, und so weiter. Wir können hier zurückgeben, was wir wollen, aber ich werde nur den ersten Wert zurückgeben.

Und jetzt können wir die Pipeline wieder ausführen und sehen, ob sie das tut, was wir erwarten.

Sicher genug, nach _splitCsv_ haben wir unsere Arrays, und dann nach dem _map_ haben wir unsere schönen sauberen Strings, nur _"Hello", "Bonjour"_ und _"Holà"_. Und die Pipeline tut jetzt das, was wir wollen. Fantastisch.

Also können wir alle diese view-Befehle jetzt loswerden. Wir brauchen sie nicht mehr.

## Zusammenfassung

Wir haben unser Debugging beendet und das ist der Code, mit dem wir enden. Wir nehmen unseren CLI-Parameter namens _input_, der als _Path_ klassifiziert ist. Nextflow findet den Pfad, lädt ihn und versteht die CSV-Datei. Gibt alle verschiedenen Zeilen zurück. Und dann mappen wir nur das erste Element dieser Zeile in den Kanal, der uns die Kanalinhalte gibt, die an den Prozess übergeben werden.

Und der Prozess läuft über jedes Element im Kanal, was drei ist. Und er führt den Prozess dreimal aus und gibt ihm drei Tasks. Und diese Ergebnisse werden dann vom Workflow veröffentlicht, vom Prozess-Output aufgenommen. Vom Workflow veröffentlicht und im Ausgabeblock in einem Unterverzeichnis namens _"hello_channels"_ gespeichert.

Ziemlich cool. Wir kommen jetzt zu etwas, das näher an einer echten Nextflow-Pipeline ist, die du für eine echte Analyse ausführen könntest.

## Fazit

Okay. Hoffentlich bekommst du jetzt ein Gefühl dafür, was Nextflow-Kanäle und Operatoren sind und wie Operatoren auf Kanälen arbeiten und wie du sie erstellen kannst.

Kanäle, wie ich am Anfang dieses Videos sagte, sind der Klebstoff von Nextflow. Und du kannst hier sehen, dass wir verschiedene Eingaben nehmen und sie manipulieren und diese Daten nehmen und sie dann in nachgelagerte Workflow-Logik übergeben können.

Und dieser Workflow-Block hier ist wirklich, wo du all diese Parallelisierung und all die clevere Logik aufbaust und Nextflow erklärst, wie du deinen Workflow-DAG aufbaust und wie du deine Pipeline orchestrierst.

Kanäle sind nicht das einfachste Konzept, um den Kopf darum zu bekommen. Also mache eine Pause, denke ein bisschen darüber nach, vielleicht lies das Material nochmal durch und stelle wirklich sicher, dass du diese Konzepte verstanden hast, weil das der Schlüssel zu deinem Verständnis von Nextflow ist und je besser du Kanäle und die verschiedenen Kanaloperatoren und die verschiedenen Channel Factories verstehst. Desto mehr Spaß wirst du beim Schreiben von Nextflow haben und desto leistungsfähiger werden deine Pipelines sein.

Das ist nicht dasselbe wie reguläres Programmieren in Python oder anderen Sprachen. Wir verwenden hier keine _if_-Anweisungen, das ist funktionale Flow-Programmierung mit Kanälen und Operatoren. Also ist es ein bisschen anders, aber es ist auch super mächtig.

Das ist das Ende dieses Kapitels. Geh und mach eine kurze Pause und ich sehe dich im nächsten Video für Teil Drei, wo wir durch Hello Workflow gehen werden und ein bisschen mehr über die Workflows sprechen werden.

Genau wie im vorherigen Kapitel gibt es ein paar Quiz-Fragen am Ende der Webseite hier, also kannst du diese schnell durchgehen und sicherstellen, dass du alle verschiedenen Teile des Materials verstehst, das wir gerade gemacht haben. Und abgesehen davon werde ich dich im nächsten Video sehen. Vielen Dank.

Okay.
