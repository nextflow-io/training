# Teil 2: Hello Channels - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../02_hello_channels.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern im Material.

## Willkommen

Hallo und willkommen zurück zu Teil 2 von Hello Nextflow. Dieses Kapitel heißt Hello Channels.

Kanäle sind wie der Klebstoff in deiner Nextflow-Pipeline. Sie sind die Teile, die alle verschiedenen Prozesse zusammenhalten, die Nextflow verwendet, um alle Informationen weiterzugeben und deinen Workflow zu orchestrieren.

Es gibt noch einen weiteren Teil von Kanälen, nämlich Operatoren. Das sind im Grunde Funktionen, die wir auf Kanäle anwenden können, um deren Inhalte zu verändern. Lass uns in VS Code eintauchen und sehen, wo wir stehen.

Ich bin sehr stark in VS Code hineingezoomt, also habe ich, um die Dinge sauber und ordentlich zu halten, alle _.nextflow\*_-Dateien und das _work/_-Verzeichnis sowie _results/_ und alles aus Kapitel Eins entfernt. Und ich fange hier einfach frisch an. Aber mach dir keine Sorgen. Wenn du nicht willst, kannst du diese Dateien auch liegen lassen. Sie werden keine Probleme verursachen.

Wir werden für dieses Kapitel mit _hello-channels.nf_ arbeiten, und wenn ich das öffne, sollte es der Datei, an der wir zuvor gearbeitet haben, sehr ähnlich sehen. Es kann sein, dass verschiedene Teile an verschiedenen Stellen im Skript sind, aber alles sollte im Grunde gleich sein.

Ein Unterschied ist, dass der Pfad im Output-Block hier jetzt _hello_channels_ für diesen Teil ist, was bedeutet, dass die Ergebnisdateien in einem anderen Unterverzeichnis in deinen Ergebnissen gespeichert werden, falls du das noch dort hast. Es sollte also ein schöner und sauberer Ort sein, um zu starten, ohne durch Ausgaben verwirrt zu werden.

Okay, lass uns kurz daran erinnern, was dieses Skript macht, wenn wir diesen Workflow ausführen. Wir machen _"nextflow run hello-channels.nf"_. Wir können _"--input myinput"_ machen, und wenn wir das ausführen, wird es diesen Parameter, params.input, verwenden, der als Variable für den sayHello-Prozess hier oben übergeben wurde, der in greeting geht und in output.txt gespeichert wird. Und das können wir in der Ergebnisdatei sehen. Großartig.

## 1. Variable Eingaben explizit über einen Kanal bereitstellen

Das ist schön. Aber es ist ziemlich simpel. Wir haben eine Variable in diesem Parameter, die in einen Prozess geht, der einmal läuft, und das skaliert nicht wirklich. Und wir können ihm nicht viele verschiedene Dateien zum Erstellen geben. Wir können ihm nicht viele verschiedene Begrüßungen geben. Wir haben nur eine.

In Wirklichkeit geht es bei Nextflow darum, deine Analyse zu skalieren. Du möchtest also wahrscheinlich, dass es mehr als eine Sache macht. Und das machen wir mit _Kanälen_.

Kanäle sind für viele Leute, die Nextflow lernen, ein etwas einzigartiges Konzept. Es kommt aus diesen Konzepten der funktionalen Programmierung, und es kann ein bisschen Zeit brauchen, bis man es versteht, aber sobald es klick macht, erschließen sie wirklich die Macht von Nextflow und sind der Schlüssel dazu, wie du deine Workflows schreibst.

## 1.1. Einen Eingabekanal erstellen

Lass uns damit beginnen, dieses Skript zu nehmen und es einen _Kanal_ verwenden zu lassen, anstatt nur einen _Parameter_.

Wir gehen runter zum Workflow, wo unsere gesamte Workflow-Logik zum Zusammenfügen von Dingen ist. Und ich werde hier reingehen und einen neuen Kanal erstellen.

Erstelle einen neuen Kanal.

Und ich werde ihn "_greeting_ch"_ nennen. Es ist Konvention, "_\_ch"_ so zu machen, nur damit du dich daran erinnern kannst, dass diese Variable ein Kanal ist. Aber du kannst ihn nennen, wie du willst.

Und dann sage ich gleich, und ich mache _"channel.of"._

Channel ist wie der Namensraum für alles, was mit Kanälen zu tun hat. Kleinbuchstabe "c", falls du Nextflow schon mal benutzt hast. Und das _".of"_ ist etwas, das Channel Factory genannt wird, was im Grunde eine Möglichkeit ist, einen Kanal zu erstellen.

Es gibt viele verschiedene Channel Factories. Wenn ich hier nur "." mache, kannst du sehen, dass VS Code viele davon vorschlägt, aber _".of"_ ist die einfachste und nimmt einfach eine Eingabe hier.

Also kann ich ein paar Klammern machen und ich sage _"Hello Channels!"_.

Großartig. Ich habe einen Kanal. Fantastisch. Ich kann speichern, ich könnte es nochmal ausführen, aber nichts Interessantes wird passieren. VS Code hat mir eine orange Warnlinie darunter gegeben und mir gesagt, dass dies eingerichtet ist: du hast das erstellt, aber du hast es nie wirklich für irgendetwas verwendet. Dieser Kanal wird nicht konsumiert.

Okay, wie verwenden wir ihn also? Ganz einfach. Ich nehme das, kopiere es, und ich lösche _params.input_ und ich setze stattdessen _"greeting_ch"_ hier ein. Also werden wir diesen Kanal als Eingabe an sayHello übergeben.

Beachte, dass ich diesen String vorerst fest codiert habe. Das ist ein bisschen ein Rückschritt nach unserem schönen Parameter, den wir am Ende des letzten Kapitels verwendet haben, aber es hält die Dinge einfach, um zu beginnen, damit du die Logik sehen kannst.

Okay, ich gehe in mein Terminal und ich führe den Workflow nochmal aus. Ohne _"--input"_ dieses Mal, und es wird laufen und es wird diesen Kanal verwenden, den wir erstellt haben, und hoffentlich sollten wir eine Datei hier oben in _results/hello_channels/_ haben und es steht jetzt "Hello Channels!" drin. Fantastisch. Das ist also das, was wir von unserem Kanal hier erhofft haben. Großartig.

## 1.4. view() verwenden, um die Kanalinhalte zu inspizieren

Noch eine Sache, die ich hier hinzufügen möchte, nur eine kurze Einführung in eine weitere Funktion, die wir auf Kanälen verwenden können, genannt "_.view"_.

Das ist analog zum _print_-Befehl in Python oder anderen Sprachen, die du vielleicht kennst, und es gibt einfach den Inhalt dieses Kanals ins Terminal aus, wenn wir es ausführen.

Also mache "_.view"_, und dann, wenn ich den Workflow nochmal ausführe, sollte es ins Terminal ausgeben, was der Inhalt dieses Kanals ist, zu dem Zeitpunkt, als wir ihn erstellt haben.

Tatsächlich kannst du sehen, dass es hier ins Terminal ausgegeben wurde. _"Hello Channels!"_.

Beachte, dass du diese Dinge über Zeilen hinweg aufbrechen kannst, wenn du willst, und tatsächlich wird der automatische Nextflow-Formatierer versuchen, das für dich zu tun. Leerzeichen sind hier nicht wirklich wichtig, also kannst du diese Dinge nacheinander verketten.

## 2. Den Workflow modifizieren, um mit mehreren Eingabewerten zu laufen

Okay, unser Kanal hat also eine Sache drin, was schön ist, aber es ist im Grunde das Gleiche wie vorher. Also lass es uns ein bisschen komplizierter machen. Lass uns ein paar mehr Dinge in unseren Kanal hinzufügen.

Die "_.of()"_-Channel Factory kann mehrere Elemente aufnehmen, also lass uns ein paar mehr schreiben. Wir machen _Hello, Bonjour, Hej_. Und dann können wir diesen Workflow nochmal ausführen und sehen, was passiert.

Sollte nochmal laufen. Und wir haben jetzt ausgegeben. _"Hello", "Bonjour"_ und _"Hej"_ ins Terminal mit unserer view-Anweisung. Fantastisch.

## 2.1.2. Den Befehl ausführen und die Log-Ausgabe ansehen

Du denkst vielleicht, dass wir an diesem Punkt fertig sind. Aber tatsächlich gibt es hier einen kleinen Fallstrick, der uns stolpern lassen wird. Wenn wir uns unsere Ausgabedatei hier ansehen. Du kannst sehen, dass _"Hello"_ drin steht, aber keine der anderen Ausgaben. Tatsächlich ist es nur diese eine.

Wenn wir diesen Workflow mehrmals ausführen, sehen wir vielleicht sogar, dass manchmal _"Bonjour"_ drin steht, manchmal _"Hej"_. Es ist ein bisschen zufällig.

Wenn wir uns das Terminal ansehen, können wir sehen, dass es dreimal gelaufen ist und wir können die verschiedenen view-Ausgaben sehen. Aber wenn ich zum work-Verzeichnis gehe, kann ich _"cat work"_ machen. Diesen Hash eingeben und das erweitern und _output.txt_. Du kannst sehen, dass diese Datei im work-Verzeichnis anders ist als im results-Verzeichnis, und diese hier ist _"Hej"._ Also funktioniert hier etwas nicht ganz richtig.

Und der Schlüssel ist, dass drei Aufgaben gelaufen sind. Die Nextflow-Ausgabe versucht das zusammenzufassen, während die Verarbeitung läuft, damit es nicht dein gesamtes Terminal übernimmt, und dieses ANSI-Logging verwendet ANSI-Escape-Codes, hat im Grunde die anderen Aufgaben überschrieben. Also zeigt es dir nur die letzte, die zufällig aktualisiert wurde.

## 2.1.3. Den Befehl nochmal mit der Option -ansi-log false ausführen

Es gibt ein paar Dinge, die wir tun können, um das tatsächlich besser zu verstehen. Wir können uns das work-Verzeichnis selbst ansehen und du kannst all die verschiedenen work-Verzeichnisse dort sehen, aber das ist ein bisschen verwirrend, weil es mit verschiedenen Nextflow-Ausführungsläufen vermischt wird.

Oder wir können Nextflow sagen, die ANSI-Escape-Codes nicht zu verwenden.

Also wenn ich den Befehl nochmal ausführe, aber dieses Mal sage ich _"-ansi-log false"_, um es auszuschalten, könnte ich auch die Umgebungsvariablen _$NO_COLOR_ oder _"$NXF_ANSI_LOG=false"_ verwenden. Dann verwendet es die eher altmodische Art des Nextflow-Loggings ohne diese Escape-Codes. Es gibt einfach direkt ins Terminal aus ohne clevere Updates.

Und jetzt können wir alle drei dieser Prozesse sehen, die gelaufen sind. Und jeder von ihnen hat seinen eigenen Aufgaben-Hash. Und wenn wir in diese work-Verzeichnisse gehen, sehen wir die drei verschiedenen Begrüßungen, die wir angegeben haben.

Das macht jetzt ein bisschen mehr Sinn. Hoffentlich verstehst du, dass Nextflow das gemacht hat, es war nur ein bisschen clever mit dem, was es dir im Terminal mit diesen work-Verzeichnissen gezeigt hat.

Allerdings ist dies für ein Problem mit den work-Verzeichnissen behoben, aber es hat ein Problem mit der Ausgabedatei nicht behoben. Wir haben immer noch nur eine Ausgabedatei, die _"Hello"_ sagt.

## 2.2. Sicherstellen, dass die Ausgabedateinamen eindeutig sein werden

Um das zu verstehen, müssen wir zurück zu unserem Workflow-Skript gehen. Wir generieren unseren Kanal hier, wir übergeben ihn an unseren Prozess, und wenn wir uns den Prozess ansehen, schreiben wir die Begrüßung in eine Datei namens _"output.txt"_ und übergeben diese Ausgabedatei zurück an den Output-Block hier unten, der sie veröffentlicht.

Allerdings generieren alle drei Male, wenn dieser Prozess läuft, diese drei verschiedenen Aufgaben, alle eine Datei namens _"output.txt"_, und alle diese Ausgabedateien werden ins results-Verzeichnis veröffentlicht, und sie überschreiben sich alle gegenseitig. Also ist die Ergebnisdatei, die du dort bekommst, einfach die letzte, die generiert wurde, aber alle anderen wurden überschrieben. Das ist nicht wirklich das, was wir wollen.

## 2.2.1. Einen dynamischen Ausgabedateinamen konstruieren

Es gibt verschiedene Möglichkeiten, damit umzugehen, aber die einfachste für jetzt ist einfach, verschiedene eindeutige Dateinamen zu erstellen. Jedes Mal, wenn die Aufgabe mit einer anderen Begrüßung läuft, wird sie eine andere Ausgabedatei generieren, die nicht mehr kollidiert, wenn sie veröffentlicht wird. Und dann bekommen wir drei eindeutige Ausgabedateien.

Wir machen das auf genau die gleiche Weise. Wir können diese Variable überall innerhalb des script-Blocks verwenden und wir können sie mehrmals verwenden.

Also kann ich sie hier einfügen, _"$\{greeting\}\_output.txt"_, und dann muss ich sie auch hier oben einfügen, weil wir nicht mehr eine Datei namens _output.txt_ erstellen. Wenn ich das also nicht aktualisiere, wird Nextflow mit einem Fehler abstürzen, der besagt, dass es eine Datei erwartet hat, die nie generiert wurde.

Also muss ich das Gleiche dort machen und ich muss doppelte Anführungszeichen verwenden, nicht einfache Anführungszeichen, damit diese Variable verstanden wird.

Okay, lass es uns ausprobieren und sehen, ob es funktioniert hat. Wir werden den Workflow nochmal ausführen. Hoffentlich zeigt es uns die drei verschiedenen Aufgaben innerhalb der drei verschiedenen work-Verzeichnisse. Und tatsächlich kannst du hier oben im results-Ordner auf der linken Seite sehen. Wir haben jetzt drei verschiedene Dateien mit drei verschiedenen Dateinamen und jede mit den verschiedenen Inhalten, die wir erwarten. Also überschreiben sich die Dateien nicht mehr gegenseitig, und alles ist da, wie wir es erwarten.

Das ist ein bisschen ein triviales Setup, das wir hier durchgegangen sind, aber es unterstreicht einige der Schlüsselkonzepte, die du über die Funktionsweise der Datei-Veröffentlichung verstehen musst, und einige der Dinge, in die du als Fallen tappen könntest. Hoffentlich kannst du das in deinen eigenen Workflows vermeiden.

Es ist auch erwähnenswert, dass das, was wir hier gemacht haben, in realen Situationen ein bisschen unpraktisch ist. Wir haben einige Eingabedaten genommen und wir verwenden diese Daten, aber wir benennen die Datei auch nach diesen Daten, was du normalerweise nicht tun kannst.

In realen, ausgereifteren Nextflow-Pipelines wirst du also oft ein Meta-Objekt mit allen Metadaten, die mit einer bestimmten Probe verbunden sind, herumreichen. Du kannst dann dynamische Dateinamen basierend darauf erstellen, was viel praktischer ist.

Wenn du daran interessiert bist, wie man das mit Best Practices macht, gibt es eine Nebenquest auf _training.nextflow.io_, die sich speziell mit Metadaten und Meta-Maps befasst, also kannst du dort für mehr Details eintauchen.

## 3. Mehrere Eingaben über ein Array bereitstellen

Okay. Als Nächstes werden wir ein bisschen darüber erkunden, wie Kanäle strukturiert sind und wie sie sich von anderen Arten von Datenstrukturen in der Programmiersprache unterscheiden. Und ich werde ein bisschen darüber nachdenken, wie ich potenziell ein Array verwenden könnte, was ein vertrautes Konzept sein könnte, wenn du von anderen Sprachen kommst.

Kann ich ein Array in einem Kanal verwenden? Lass es uns versuchen. Ich werde ein Array erstellen, und ich habe das aus den Docs kopiert, _"greetings_array"_ und _"Hello", "Bonjour"_ und _"Holà"_. Und dann werde ich das hier anstelle meiner fest codierten Strings einfügen. Also sage ich "channel.of" _"greetings_array"_, dieses Array in einen Kanal übergeben. Lass es uns versuchen.

Terminal hochbringen und die Pipeline ausführen.

Okay. Du kannst sehen, dass die view-Anweisung hier unser Array wie erwartet ausgegeben hat, aber dann all dieser rote Text, oder er wird nicht rot sein, wenn du immer noch _"-ansi-log"_ aus hast, aber all dieser rote Text sagt uns, dass etwas schief gelaufen ist.

Wir haben hier kein schönes grünes Häkchen mehr. Wir haben ein rotes Kreuz, und wenn ich das nur ein bisschen breiter mache, damit es leichter zu lesen ist, sagt uns Nextflow, was schief gelaufen ist.

Also lass uns das Abschnitt für Abschnitt aufschlüsseln. Es sagt, der Fehler wurde verursacht durch, und dann der Grund für den Fehler, nämlich fehlende Ausgabedateien. Also im Grunde sagte dieser Output-Block, dass diese Datei erstellt werden sollte und sie wurde nicht erstellt. Als Nächstes sagt es, dies ist der Befehl, der ausgeführt wurde. Das ist also im Grunde der Inhalt dieser _.command.sh_-Datei. So sah es aus, nachdem all diese Variablen eingefügt wurden.

Und du kannst hier sehen, dass unser echo-Befehl tatsächlich nur einmal ausgeführt wurde und er das gesamte Array verwendet hat, aber in einer String-Darstellung, was nicht wirklich das ist, was wir wollten.

Und dann ist der Befehl so beendet worden, und das war das work-Verzeichnis, wo wir hingehen und die Dateien ansehen können, um ein bisschen mehr zu verstehen.

Okay. Was also passiert ist, war. Nextflow hat einfach dieses gesamte Array als einzelnes Kanalelement an den Prozess übergeben, was bedeutete, dass der Prozess nur einmal lief. Er hatte eine Aufgabe und er hat die Daten nicht in einer Struktur verwendet, die wir erwartet haben.

## 3.2. Einen Operator verwenden, um Kanalinhalte zu transformieren

Wir müssen also etwas mit diesem Kanal machen, bevor er verwendet werden kann. Und das bereitet die Bühne für die Verwendung von Operatoren, die spezielle Funktionen sind, die wir auf Kanälen verwenden können, um Kanalinhalte zu manipulieren.

In diesem Fall werden wir etwas verwenden, das _flatten_ heißt. Das wir am Ende des Kanals hier anhängen. Also erstellen wir den Kanal und dann führen wir _flatten_ aus. Und wieder, wenn wir darüber schweben, zeigt es uns die Dokumentation für diesen Befehl direkt in VS Code, was sehr hilfreich ist. Du kannst all diese Docs auch auf der Nextflow-Website finden, in der Dokumentation.

Ich könnte diesen Code jetzt einfach ausführen und sehen, ob er funktioniert, aber es ist auch eine schöne Gelegenheit, einzuführen, wie man dynamischen Code innerhalb von Operatoren und innerhalb von Nextflow-Code macht, die Closures genannt werden.

Also werde ich hier einen view-Befehl wieder hinzufügen, bevor wir _flatten_ ausführen. Und hier hat dieser diese geschweiften Klammern, was die dynamische Closure ist. Und es gibt nur etwas beliebigen Code hier drin, der ausgeführt wird, im Kontext eines view-Operators.

Hier sagt das, nimm die Begrüßung, was die Eingaben des view-Operators sind, und das ist hier. Ich könnte das nennen, wie ich wollte, ich könnte das _"foo"_ nennen und ich muss es nur später als _"foo"_ bezeichnen. Und dann sage ich damit, gib das zurück.

Und dann setze ich einen String zurück, der sagt vor dem flatten für eine Variable. Sehr einfach.

Ich werde jetzt noch einen davon genau gleich hinzufügen, aber ich werde sagen nach _flatten_.

Was das also macht, weil das in Sequenz läuft, wirst du sehen, wie der Kanal aussieht, bevor wir _flatten_ ausführen, und dann wieder nachdem wir _flatten_ ausgeführt haben.

Und dann wird dieser greeting-Kanal immer noch erstellt, also wird er immer noch an den Prozess übergeben. Und hoffentlich wird der Workflow jetzt laufen. Lass es uns ausprobieren.

Großartig. Also zuerst ist die Pipeline dieses Mal nicht abgestürzt. Wir hatten drei Prozesse, die ordnungsgemäß gelaufen sind, und wir haben ein kleines Häkchen. Und dann können wir sehen, dass unsere view-Anweisungen funktioniert haben.

Wir haben vor _flatten_, was dieses Array ist, das wir vorher vom Fehler gesehen haben, und dann haben wir dreimal das nach _flatten_ aufgerufen wurde, wo wir _"Hello", "Bonjour"_ haben, und all diese anderen drei separaten Elemente im Array, die jetzt, wie wir gehofft haben, drei separate Elemente im Kanal sind.

Und du kannst sehen, dass der _view_-Operator dreimal ausgeführt wurde. Und das liegt daran, dass dieser Kanal nach _flatten_ jetzt drei Elemente hat. Und so wird der Operator dreimal aufgerufen.

Ganz kurz würde ich nur erwähnen, dass als ich vorher Channel Factories erstellt habe, ich _"."_ gemacht habe, und dann haben wir gesehen, dass es viele verschiedene Möglichkeiten gibt, Kanäle zu erstellen, und eine davon heißt "_fromList"_. Und das ist tatsächlich speziell dafür konzipiert, diese gleiche Operation zu machen. Wir hätten also einfach fromList greetings_array machen können, und das wird funktionieren. Es ist etwas sauberer und schönere Syntax. Aber für die Zwecke dieser Demonstration wollten wir es ein bisschen schrittweiser machen, damit du sehen kannst, wie der Kanal manipuliert wird und wie verschiedene Operatoren ändern können, was im Inhalt eines Kanals ist.

## 4. Eingabewerte aus einer CSV-Datei lesen

Okay, wie können wir das ein bisschen realistischer machen? Du wirst wahrscheinlich nicht viel Code in deiner Nextflow-Pipeline mit fest codierten Arrays erstellen wollen. Du wirst wahrscheinlich die Daten von außen nehmen wollen, wenn du startest, und diese Daten werden mit ziemlicher Sicherheit in Dateien sein.

Das Nächste, was wir also tun werden, ist, dass wir das replizieren werden, aber anstatt die Daten von einem einzelnen CLI-Parameter oder von einem fest codierten String oder Array zu nehmen, werden wir sie aus einer Datei nehmen.

Also lass uns unser greetings*array loswerden. Und jetzt werden wir diese Channel Factory wieder ändern. Ich habe gerade gesagt, es gibt eine Menge zur Auswahl und es gibt eine namens *".fromPath"_. Und ich werde ihr sagen, in diesem Fall \_params.input_ zu nehmen, was zurück zu unserer Eingabe geht, die wir früher verwendet haben.

Jetzt ist dieser Parameter noch nicht wirklich bereit, verwendet zu werden. Wir sagen immer noch, dass es ein String ist und er ist hier mit einem Standard fest codiert, aber wir könnten diesen String überschreiben. Wir wollen jetzt, dass das stattdessen eine Datei ist. Also ist der Typ anders. Es ist nicht mehr ein _String_. Es ist ein _Path_.

Und dann können wir den Standard setzen, wenn wir wollen, wieder auf einen Path. Und wenn ich im Explorer auf der linken Seite schaue, kannst du sehen, dass ich in diesem Repository, in diesem Arbeitsverzeichnis, ein Verzeichnis namens data habe. Ich habe dort eine Datei namens _"greetings.csv"._

Also kann ich den Standard hier einfach auf _"data/greetings.csv"_ setzen. Wenn ich diese Pipeline jetzt wieder ohne Befehlszeilenoptionen ausführe, wird sie diesen Standardwert verwenden. Sie weiß, dass es ein Pfad ist, also weiß sie, dass sie das als Pfad behandeln sollte und nicht als String.

Und dann wird sie das in eine Channel Factory von diesem _params.input_ übergeben und unseren Kanal erstellen, der dann in diesem Prozess namens _sayHello_ verwendet wird. Lass es uns ausprobieren.

Okay. Fehlgeschlagen. Keine Sorge. Das war erwartet. Und wenn du dem Trainingsmaterial folgst, wirst du sehen, dass es dort auch erwartet wurde. Lass uns sehen, was hier passiert.

Es hat versucht, die Pipeline auszuführen. Es hat versucht, den Prozess auszuführen, und es hat einen ziemlich ähnlichen Fehler wie den bekommen, den wir vorher gesehen haben.

Hier steht: wir haben versucht, _echo_ auszuführen, aber anstatt den Inhalt dieser CSV-Datei zu echoen, hat es nur den Pfad geechot. Und du kannst sehen, es ist der vollständige absolute Pfad hier zu dieser CSV-Datei.

Und dann, weil es versucht hat, das in diesen wirklich komplizierten Pfad zu schreiben, wusste es nicht wirklich, was es tun sollte. Und es war außerhalb des Geltungsbereichs des Prozess-work-Verzeichnisses.

Ich habe am Anfang erwähnt, dass Nextflow jede ausgeführte Aufgabe innerhalb eines speziellen work-Verzeichnisses kapselt. Und wenn du versuchst, Daten zu schreiben, die außerhalb dieses work-Verzeichnisses sind, wird Nextflow dich als Sicherheitsvorkehrung stoppen. Und das ist hier passiert. Wir haben versucht, in einen absoluten Pfad zu schreiben und Nextflow ist fehlgeschlagen und hat uns daran gehindert.

## 4.2. Den splitCsv()-Operator verwenden, um die Datei zu parsen

Okay, lass uns einen Blick auf diesen Kanal werfen und sehen, wie er aussieht. Wir können _".view"_ machen, und ich habe das von der Website kopiert. Also _.view_, und wir haben eine dynamische Closure hier und wir sagen einen Variablennamen "_csv"_ als Eingabe. Das sind also die Kanalinhalte, und wir sagen vor splitCsv, und so sieht es aus.

Wenn ich es nochmal ausführe, wird es immer noch fehlschlagen, aber es wird uns zeigen, was in diesem Kanal ist. Es ist nicht besonders aufregend. Es ist diese _path_-Variable. Du kannst also sehen, es ist nur ein String hier, weil es ins Terminal ausgegeben wird, aber es ist ein _path_-Objekt, das die Informationen und Metadaten über diese Datei enthält.

Wir wollen nicht die Metadaten der Datei an die Eingabe übergeben. Wir wollen den Inhalt dieser Datei übergeben. Wenn wir uns die _greetings.csv_-Datei ansehen, kannst du hier sehen, dass sie diese verschiedenen Variablen hier hat. _Hello, Bonjour, Holà_ wieder. Und das sind die Dinge, die wir wirklich an unseren Prozess übergeben wollen, nicht nur die Datei selbst als einzelnes Objekt.

Wir müssen also diese CSV-Datei parsen. Wir müssen sie auspacken, an den Inhalt der CSV-Datei kommen und dann den Inhalt innerhalb des Kanals an den Prozess übergeben.

Wie du wahrscheinlich aus der Log-Nachricht erkennen kannst, wollen wir _splitCsv_ verwenden, was ein weiterer Operator ist, ein weiterer Kanaloperator. Also wenn ich "_dot" "s"_ mache, und dann kannst du sehen, es wird automatisch vorgeschlagen. Ups, _splitCsv_ und ein paar Klammern.

Und dann nach _splitCsv_ werde ich eine weitere _view_-Anweisung einfügen, nur damit wir sehen können, wie es danach aussieht. Lass uns die Pipeline ausführen und sehen, was wir haben.

Okay. Es ist immer noch fehlgeschlagen, aber auf eine neue und aufregende Weise, was Fortschritt ist.

Dieses Mal haben wir wieder ein Problem mit unserem Skript, das gerendert wurde. Jetzt. Wir haben nicht mehr den finalen Pfad, aber wir haben ein Array von Variablen, was sehr nach dem Fehler aussieht, den wir früher hatten, als wir ein Array als feste Eingabe übergeben haben.

Mit unserem Logging vom view-Operator können wir sehen, vor _splitCsv_ war der Pfad. Und tatsächlich, nach _splitCsv_, haben wir drei verschiedene Ausgaben und jede dieser Ausgaben sieht sehr nach jeder der Zeilen aus der _greetings.csv_-Datei aus, was Sinn macht.

Was also hier passiert ist, ist, dass Nextflow diese CSV-Datei geparst hat, uns drei Objekte gegeben hat, ein Array für jede Zeile der CSV-Datei. Also haben wir dann dreimal ein Array von Variablen an den Kanal übergeben anstelle eines einzelnen String-Werts.

Okay, also das letzte Mal, als wir dieses Problem hatten, haben wir _flatten_ verwendet. Lass uns ganz schnell. Flatten ausprobieren und sehen, was passiert.

Ich kann diese Variablen nennen, wie auch immer. Also werde ich es _myarray_ nennen, weil es nicht mehr wirklich ein CSV ist. Lass uns versuchen, es nochmal auszuführen und sehen, was mit _flatten_ passiert.

Also dieses Mal werden wir laufen, wir haben das CSV in drei Array-Objekte geparst, und dann haben wir es geflattened. Und dieses Mal ist es durchgegangen. Und die Nextflow-Pipeline ist gelaufen. Allerdings kannst du sehen, dass _flatten_ wirklich zur Sache geht und alles flattened. Und so bekommen wir drei unabhängige Array-Einträge für jede Zeile. Und so lief der Prozess dreimal für jede Zeile eines CSV. Und jetzt haben wir eine ganze Menge Ergebnisdateien, und 123, 456, und alle möglichen Dinge, nicht nur diese erste Spalte des CSV, was wir wirklich wollten.

## 4.3. Den map()-Operator verwenden, um die Begrüßungen zu extrahieren

Wie kommen wir also nur an die erste Spalte? Wenn flatten hier zu simpel ist, brauchen wir einen komplexeren Operator, wo wir tatsächlich anpassen und ihm sagen können, was wir vom CSV wollen.

Um das zu tun, werden wir _map_ verwenden. Im Grunde sagt _map_ einfach, führe etwas Code aus, eine Funktion über jedes Element, das ich bekomme, und mache irgendeine Art von Transformation darauf. Und weil es so flexibel ist, wirst du es ständig in Nextflow-Code auftauchen sehen.

Für sich allein macht es nichts. Also wollen wir keine regulären Klammern, wir wollen hier eine Closure und wir müssen ihm sagen, was zu tun ist. Also sage ich _"row"_, weil das Zeilen vom CSV bekommt, also ist es ein logischer Variablenname. Ist die Eingabe. Und ich möchte nur das erste Element dieses Arrays zurückgeben.

Arrays in Nextflow sind nullbasiert, also werden wir nur das erste Element sagen, was Zeile null ist. Wenn wir die zweite Spalte wollten, könnte ich eins sein oder die dritte Spalte zwei sein, und so weiter. Wir können hier zurückgeben, was wir wollen, aber ich werde nur den ersten Wert zurückgeben.

Und jetzt können wir die Pipeline nochmal ausführen und sehen, ob sie das macht, was wir erwarten.

Tatsächlich, nach _splitCsv_ haben wir unsere Arrays, und dann nach dem _map_ haben wir unsere schönen sauberen Strings, nur _"Hello", "Bonjour"_ und _"Holà"_. Und die Pipeline macht jetzt das, was wir wollen. Fantastisch.

Wir können jetzt all diese view-Befehle loswerden. Wir brauchen sie nicht mehr.

## Zusammenfassung

Wir haben unser Debugging beendet und das ist der Code, mit dem wir enden. Wir nehmen unseren CLI-Parameter namens _input_, der als _Path_ klassifiziert ist. Nextflow findet den Pfad, lädt ihn und versteht die CSV-Datei. Gibt alle verschiedenen Zeilen zurück. Und dann mappen wir nur das erste Element dieser Zeile in den Kanal, der uns sozusagen die Kanalinhalte gibt, die an den Prozess übergeben werden.

Und der Prozess läuft über jedes Element im Kanal, was drei sind. Und er führt den Prozess dreimal aus, was ihm drei Aufgaben gibt. Und diese Ergebnisse werden dann vom Workflow veröffentlicht, vom Prozess-Output aufgenommen. Vom Workflow veröffentlicht und im Output-Block in einem Unterverzeichnis namens _"hello_channels"_ gespeichert.

Ziemlich cool. Wir kommen jetzt zu etwas, das einer realen Nextflow-Pipeline, die du für eine echte Analyse ausführen könntest, näher kommt.

## Fazit

Okay. Hoffentlich bekommst du jetzt ein Gefühl dafür, was Nextflow-Kanäle und -Operatoren sind und wie Operatoren auf Kanälen arbeiten und wie du sie erstellen kannst.

Kanäle, wie ich am Anfang dieses Videos gesagt habe, sind der Klebstoff von Nextflow. Und du kannst hier sehen, dass wir verschiedene Eingaben nehmen und sie manipulieren und diese Daten nehmen und sie dann in nachgelagerte Workflow-Logik übergeben können.

Und dieser Workflow-Block hier ist wirklich, wo du all diese Parallelisierung und all die clevere Logik aufbaust und Nextflow erklärst, wie es deinen Workflow-DAG bauen soll und wie es deine Pipeline orchestrieren soll.

Kanäle sind nicht das einfachste Konzept, um den Kopf herumzubekommen. Also mach eine Pause, denk ein bisschen darüber nach, lies vielleicht das Material nochmal durch und stelle wirklich sicher, dass du diese Konzepte verstanden hast, denn das ist der Schlüssel zu deinem Verständnis von Nextflow, und je besser du Kanäle und die verschiedenen Kanaloperatoren und die verschiedenen Channel Factories verstehst. Desto mehr Spaß wirst du beim Schreiben von Nextflow haben und desto mächtiger werden deine Pipelines sein.

Das ist nicht dasselbe wie reguläres Programmieren in Python oder anderen Sprachen. Wir verwenden hier keine _if_-Anweisungen, das ist funktionale Flow-Programmierung mit Kanälen und Operatoren. Es ist also ein bisschen anders, aber es ist auch super mächtig.

Das ist das Ende dieses Kapitels. Geh und mach eine kurze Pause und ich sehe dich im nächsten Video für Teil drei, wo wir Hello Workflow durchgehen werden und ein bisschen mehr über die Workflows sprechen.

Genau wie im vorherigen Kapitel gibt es ein paar Quiz-Fragen am Ende der Webseite hier, also kannst du diese schnell durchgehen und sicherstellen, dass du alle verschiedenen Teile des Materials verstehst, das wir gerade gemacht haben. Und abgesehen davon sehe ich dich im nächsten Video. Vielen Dank.

Okay.

​
