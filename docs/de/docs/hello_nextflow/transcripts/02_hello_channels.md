# Teil 2: Hello Channels - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../02_hello_channels.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen

Hi, willkommen zu Teil zwei von Hello Nextflow.

Dieses Kapitel heißt Hello Channels. Wir werden über diesen grundlegenden Teil von Nextflow sprechen.

Channels sind die Dinge, die die verschiedenen Schritte in deiner Pipeline verbinden, die Art und Weise, wie deine Daten und Logik durch deinen Workflow fließen.

Okay, lass uns eintauchen.

Beginnen wir damit, zu training.nextflow.io zu gehen

Hello Nextflow in der Seitenleiste und klicke auf Teil zwei. Hello Channels.

Das gesamte Material ist hier aufgeschrieben, sodass du in deinem eigenen Tempo folgen und alles nachholen kannst, was du vielleicht verpasst hast.

Sobald du die Website geöffnet hast, kannst du Codespaces laden und wir machen dort weiter, wo wir am Ende des letzten Kapitels aufgehört haben.

## 0. Aufwärmen: hello-channels.nf ausführen

Für dieses Kapitel werden wir eine andere Datei bearbeiten. Diese heißt Hello Channels, du findest sie also in der Seitenleiste, doppelklicke darauf, um sie zu öffnen.

Wenn du gerade von Kapitel eins kommst, wird dir diese Datei sehr bekannt vorkommen. Der Ausgangspunkt hier ist im Grunde dort, wo wir Kapitel eins beendet haben, mit unserem process namens sayHello, unserer Eingabe, Ausgabe, unserem publishDir und unserem params.greeting und unserem einfachen workflow.

Wir beginnen mit einer neuen Datei, also ist es für alle eine ebene Spielfläche, aber du kannst mit deiner vorherigen Datei fortfahren, wenn du möchtest.

Beachte, ich habe auch alle .nextflow\*-Dateien und die work-Verzeichnisse hier gelöscht, nur damit es ein sauberer Ausgangspunkt ist. Es ist nicht wichtig, ob du das tust oder nicht, das liegt bei dir.

Okay. Lass uns zunächst überprüfen, dass diese Pipeline noch wie erwartet funktioniert. Ich werde das Terminal hier aufrufen.

Gib "nextflow run hello-channels.nf" ein und drücke Enter.

Es wird diesen kleinen Workflow ausführen, unseren sayHello-Schritt ausführen, ein work-Verzeichnis mit diesem Hash generieren, und hier ist unser results-Ordner und da ist unsere Ausgabedatei, genau wie wir es von unserem Standard-params.greeting erwartet haben.

Das ist großartig. Genau wie in Kapitel eins, funktioniert wie erwartet.

## 1. Variable Eingaben explizit über einen channel bereitstellen

In Kapitel eins hast du bereits channels verwendet, du hast es nur nicht bemerkt. Als wir hier einen String angegeben haben, hat Nextflow automatisch einen channel um diesen String für uns erstellt, einfach weil es wusste, dass wir einen process aufrufen, also brauchten wir einen Eingabe-channel.

Das erste, was wir tun werden, ist, es explizit zu machen, indem wir den channel tatsächlich selbst aufschreiben.

## 1.1. Einen Eingabe-channel erstellen

Ich gehe also zum workflow hier am Ende des Skripts, und ich werde greeting_ch sagen. Das ist eine Konvention, die wir oft in Nextflow-Code verwenden, ein Unterstrich ch am Ende eines Variablennamens zu haben, wenn es ein channel ist, nur damit es leicht zu erkennen ist, dass es ein channel ist, aber du musst das nicht tun. Gleich channel of Hello Channels.

Was wir gerade verwendet haben, wird in der Nextflow-Sprache als "Channel Factory" bezeichnet. Das ist dieses Ding hier, wir setzen diese Variable auf einen neuen channel, und diese channel factory hier erstellt einen channel für uns auf eine bestimmte Weise.

Es gibt eine Handvoll verschiedener channel factories, die Nextflow hat, um channels aus verschiedenen Arten von Eingaben zu erstellen. Dot of ist die einfachste und nimmt einfach alle Strings, die wir ihr geben.

Beachte, wenn ich in VS Code über diese Wörter fahre, gibt mir die Nextflow-Erweiterung ein Popup, das erklärt, was diese Syntax macht, und es gibt auch einen "Mehr lesen"-Text am unteren Rand dieses Popup-Fensters.

Wenn ich darauf klicke, öffnet es die Nextflow-Dokumentation in einem neuen Tab und führt mich direkt zur Dokumentation für diese spezifische Sache. In diesem Fall für channel.of.

## 1.2. Den channel als Eingabe zum process-Aufruf hinzufügen

Beachte, dass die Erweiterung uns auch eine Warnung gibt, die besagt, dass wir hier einen neuen channel erstellt haben, aber er wird von nichts verwendet.

Also, lass uns das beheben. Ich nehme den neuen channel-Namen und werde dieses params.greeting durch unseren neuen channel ersetzen.

Beachte, dass wir jetzt das Kommandozeilen-Flag --greeting nicht mehr verwenden, params.greeting wird nicht verwendet, wir gehen zurück zum Hardcodieren dieses Strings. Das ist okay. Ich versuche nur, die Dinge einfach zu halten. Wir kommen später zurück und verwenden die params wieder.

## 1.3. Den workflow-Befehl erneut ausführen

Okay, lass uns nur noch einmal überprüfen, dass dies funktioniert. Rufe das Terminal auf und beachte nochmal. Nextflow run hello channels. Überprüfe output.txt, und da ist es.

Tolles Beispiel, das genau dasselbe tut wie zuvor, aber jetzt ist die Logik zumindest ein bisschen klarer. Wir sind explizit beim Schreiben eines neuen channels.

Wir haben effektiv einfach mehr Code geschrieben, um dasselbe zu tun. Aber das wird mehr Sinn ergeben, wenn wir mit der Art und Weise, wie wir unsere channels erstellen, etwas komplizierter werden.

## 2. Den Workflow modifizieren, um mit mehreren Eingabewerten zu arbeiten

Okay, lass uns das ein bisschen interessanter machen. Es ist sehr selten, dass du eine Nextflow-Pipeline auf einer einzigen Eingabe ausführen möchtest, also lass uns ihr mehrere Eingaben geben.

## 2.1. Mehrere Grüße in den Eingabe-channel laden

Aus der Dokumentation hier. Ich werde diese verschiedenen Strings kopieren, drei davon. Hello, Bonjour, Olà. Oh, hoffe. Copilot schlägt ein paar andere vor. Also lass uns diese mit Tab Enter eingeben.

Die Nextflow-Dokumentation hier sagt uns, dass wir diesem Operator mehrere Werte geben können, also sollte es funktionieren, aber lass es uns ausprobieren und sehen, was passiert.

## 2.1.2. Den Befehl ausführen und die Log-Ausgabe betrachten

Nun. Ja und nein. Lass uns sehen. Es sagt, dass fünf von fünf Aufgaben hier ausgeführt wurden, aber es zeigt uns nur einen Hash, was ein bisschen seltsam ist. Das ist okay. Alles ist wie erwartet hier. Standardmäßig verwendet Nextflow eine spezielle Art der Ausgabe an ein Terminal namens ANSI-Kontrollcodes, was bedeutet, dass es bestimmte Zeilen überschreibt, um eine schöne komprimierte Ansicht aller verschiedenen Prozesse zu geben, die ausgeführt werden.

Das macht viel mehr Sinn, wenn du größere Workflows hast und Hunderte oder Tausende verschiedener Proben ausführst. Du kannst einfach so viel Ausgabe auf dem Terminal generieren, dass es unmöglich ist, sie anzuschauen, während diese aktualisierende Ansicht dir einen Echtzeit-Fortschritt gibt.

## 2.1.3. Den Befehl erneut mit der Option -ansi-log false ausführen

Wenn du möchtest, kannst du es erneut ausführen, und diesmal werde ich ein zusätzliches Nextflow-Kernargument mit einem einzelnen Bindestrich verwenden und sagen: "-ansi-log false". Dies verwendet die vorherige Version der Nextflow-Log-Ausgabe. Und hier kannst du alle einzelnen Prozesse sehen, die gestartet wurden.

Es liegt bei dir, ob du dies tust oder nicht. Die Ausgabe von Nextflow ist in beiden Fällen genau gleich.

## 2.2. Sicherstellen, dass die Ausgabedateinamen eindeutig sind

Okay, schauen wir uns die Ausgabedateien an, dann gehen wir zu results. Aber es gibt nur eine einzige Ausgabedatei. Was ist passiert? Wir haben gesehen, dass der process viele Male ausgeführt wurde. Wir können ins work-Verzeichnis gehen und alle verschiedenen Hashes sehen, alle Aufgaben wurden ordnungsgemäß ausgeführt. Aber wenn du dich an unseren process hier erinnerst, speichern wir alles in einer output.txt-Datei und veröffentlichen diese dann in diesem Verzeichnis.

Die gleiche Datei wurde also fünfmal erstellt und dann fünfmal überschrieben. Und wir haben nur, welche Aufgabe auch immer zuletzt ausgeführt wurde.

## 2.2.1. Einen dynamischen Ausgabedateinamen konstruieren

Die Art, wie wir das beheben, ist die Verwendung eines dynamischen Ausgabedateinamens. Hier haben wir bereits eine Variable namens greeting innerhalb des process, also können wir diese im Ausgabedateinamen verwenden. Ich kopiere das und mache $greeting-output.txt.

Ich werde das in Anführungszeichen setzen, nur damit bash nicht durch irgendwelche Leerzeichen verwirrt wird, die hier hineinschleichen könnten. Und dann nehme ich denselben Dateinamen und aktualisiere die Ausgabe hier.

Es ist wirklich wichtig, dass die Ausgabe damit übereinstimmt, denn sonst wird diese Datei nicht gefunden und Nextflow wird abstürzen.

Ich werde noch eine wirklich wichtige Änderung vornehmen, nämlich diese einfachen Anführungszeichen durch doppelte Anführungszeichen ersetzen. Beachte, dass sich die Farbe des Codes geändert hat, als ich das tat. Diese Variable wird nur erweitert, wenn wir doppelte Anführungszeichen verwenden. Wenn ich hier einfache Anführungszeichen verwende, wird es als literaler Wert verwendet, und ich würde eine einzelne Datei namens $greeting-output bekommen, was nicht das ist, was ich will.

## 2.2.2. Den Workflow ausführen

Also lass uns die doppelten Anführungszeichen zurücksetzen und es versuchen.

Ich werde nur mein Verzeichnis aufräumen, bevor ich beginne, damit es einfach ist, die neuen Dateien zu sehen. Ich werde alles löschen, was .nextflow, work und results heißt.

Und ich werde diesen Nextflow-Befehl erneut ausführen und sehen, welche Dateien erstellt werden. Es führt also die fünf Prozesse dort aus. Wenn du sehr genau zugeschaut hast, hast du vielleicht gesehen, dass sich diese Zeile aktualisiert hat, während sie lief.

Und jetzt können wir ins results-Verzeichnis gehen, und tatsächlich haben wir fünf verschiedene Ausgaben, und sie sind alle mit dem jeweiligen Gruß vorangestellt.

Wenn ich jede davon öffne, werden wir sehen, dass sie jeweils den entsprechenden Gruß enthalten. Fantastisch. Das ist es, was wir wollen.

## 3. Einen Operator verwenden, um den Inhalt eines channels zu transformieren

Okay, jetzt wissen wir, was channels sind und was channel factories sind. Was ist mit Operatoren? Dies ist ein weiterer Begriff für einen Teil der Nextflow-Sprache, der eine Reihe von Funktionen ist, die es uns ermöglichen, auf channels zu operieren, um bestimmte Dinge mit ihnen zu tun. Nextflow kommt mit einer Reihe von Operatoren, die es uns ermöglichen, channels auf verschiedene Arten zu manipulieren.

## 3.1. Ein Array von Werten als Eingabe für den channel bereitstellen

Lass uns das mit einem Beispiel durcharbeiten. Nehmen wir an, wir möchten diese Eingabe-Strings nehmen, aber anstatt sie direkt in eine channel factory zu stecken, möchten wir sie als Array definieren.

## 3.1.1. Die Eingabevariable einrichten

Ich werde diese nehmen und das als neue Zeile darüber machen und sagen, greetings, array.

Da haben wir's. Ich nehme diese Array-Variable und stecke sie in das channel.of und speichere.

## 3.1.3. Den Workflow ausführen

Jetzt, lass uns sehen, was passiert. Zurück zu meinem Terminal. Ich werde nur all diese temporären Dateien wieder aufräumen. Und lass uns den Workflow ausführen.

Nicht gut. Okay. Es ist kaputt gegangen. Das ist okay. Ich habe erwartet, dass es diesmal kaputt geht. Das Debuggen dessen, was schief läuft, wenn ein Nextflow-Workflow fehlschlägt, ist ein wichtiger Teil davon, Nextflow-Entwickler\*in zu sein. Das wird oft passieren, und es ist wichtig zu verstehen, was die Fehlermeldung sagt und wie man damit umgeht.

Die Nextflow-Fehlermeldungen sind tatsächlich ziemlich strukturiert. Sie sagen uns, welcher process schief gelaufen ist. Sie geben uns eine Fehlermeldung für einen Grund. Sie sagen, was der Befehl war, den es zu führen versuchte innerhalb dieser speziellen Aufgabe, was der Exit-Status war, was die Ausgabe war und wo dieses Aufgaben-work-Verzeichnis war.

Beachte, dass ich dies in VS Code mit Option anklicken kann und es in einer Seitenleiste öffnet, sodass ich direkt dorthin gehen und alle diese versteckten Dateien ansehen kann, über die wir im vorherigen Kapitel gesprochen haben, einschließlich der .command.sh-Datei. Diese kannst du sehen, ist dieselbe wie die Befehle, die hier ausgeführt wurden.

Indem wir uns diese Datei ansehen, können wir ein Gefühl dafür bekommen, was hier schief gelaufen sein könnte, anstatt eine einzelne Aufgabe für jedes Element im Array auszuführen, wie es beim letzten Mal tat, hat es einfach das gesamte Array auf einmal als String bereitgestellt. Also müssen wir dieses Array in einzelne Werte entpacken, bevor wir es in den channel übergeben. Lass uns zurückgehen und sehen, ob wir das mit einem Operator tun können.

## 3.2. Einen Operator verwenden, um channel-Inhalte zu transformieren

In diesem Fall werden wir das Array nicht ändern, bevor wir es in den channel übergeben. Wir werden den channel so anpassen, dass er sich so verhält, wie wir es erwarten. Wir werden das tun, indem wir den flatten-Operator verwenden, kannst du dot anfangen zu tippen, und wir können sehen, dass die VS Code-Erweiterung beginnt, alle verschiedenen verfügbaren Operatoren vorzuschlagen.

## 3.2.1. Den flatten()-Operator hinzufügen

Und ich wähle flatten aus. Beachte, dass Leerzeichen in diesem Kontext für Nextflow keine Rolle spielen. Du kannst diese Operatoren also auf eine neue Zeile setzen, wenn du möchtest. Ich kann das also hier ablegen und einrücken, sodass es unter ".of" sitzt, und du wirst sehen, dass die Leute oft viele Operatoren wie diesen auf einen channel verketten und ihn auf diese Weise einrücken, damit er leichter zu lesen ist.

Du kannst auch sehen, wie zuvor kann ich darüber fahren und lesen, was der flatten-Operator macht, und auch einem Link zur Dokumentation folgen, wenn ich möchte.

Dieser Operator nimmt also diesen channel, der ein einzelnes Array darin hat, und trennt die Array-Werte.

## 3.2.2. view() hinzufügen, um channel-Inhalte zu inspizieren

Wir können in die channels hineinschauen, indem wir den speziellen view-Operator verwenden, und ich werde hier ein paar davon hinzufügen. Das ist ein bisschen wie die Verwendung von Print-Anweisungen in anderen Sprachen. Ich werde also dot view machen und dann diese geschweiften Klammern verwenden.

Das wird closure genannt. Das gibt im Grunde zusätzlichen Code an den view-Operator, den er auf jedem Element innerhalb des channels ausführt. In diesem Fall werde ich sagen greeting before flatten. Greeting.

Ich definiere hier eine Variable, die nur innerhalb des Geltungsbereichs dieser closure liegt. Diese Variable wird also nur hier verwendet, und ich könnte sie nennen, wie ich wollte. Es ist nicht wirklich wichtig. Ich verwende nur greeting, um es leicht lesbar zu machen.

In einigen Nextflow-Pipelines siehst du vielleicht, dass Leute eine spezielle implizite Variable namens "$it" verwenden. So. Das ist eine spezielle Variable innerhalb von Nextflow-Code, die eine Abkürzung ist, sodass du die kleine Definition einer Variablen nicht machen musst. Im Laufe der Zeit denken wir jedoch, dass dies für Leute, die neu bei Nextflow sind, nicht sehr klar ist, und wir raten von der Verwendung von "$it" jetzt ab.

Also werde ich beim vorherigen Verhalten von greeting bleiben und es so verwenden, weil das expliziter ist und klarer ist, was vor sich geht.

Ich werde dann diese Zeile kopieren und genau dasselbe noch einmal nach den flatten-Argumenten machen. Der view-Operator ist ein bisschen speziell, weil er etwas mit den Elementen macht, aber er gibt sie auch einfach an den nächsten Operator weiter, sodass wir ihn in der Mitte einer Kette von Operationen wie dieser verketten können, und er wird den Status dort drucken und weitermachen. Hoffentlich wird uns das zeigen, wie der channel vor und nach dem flatten-Operator aussieht.

## 3.2.3. Den Workflow ausführen

Lass es uns ausprobieren. Räume alles im Workspace auf. Führe die Pipeline erneut aus.

Okay, wir können sehen, dass es unsere fünf Prozesse erneut ausgeführt hat. Es ist nicht mit einem Fehler abgestürzt, das ist definitiv gut. Und jetzt haben wir das before flatten und es hat sicher genug unser Array und wir haben after flatten, fünfmal gedruckt, einmal für jedes Element des Arrays. Das ist genau das, worauf wir gehofft hatten. Das sind wirklich gute Nachrichten. Und das passt genau zu dem, was wir vom Code erwarten würden.

Wir brauchen diese Debug-Anweisungen nicht mehr, also kann ich sie entweder auskommentieren oder löschen. Ich werde sie löschen, nur um meinen Code schön und sauber zu halten. Okay, toll. Dieses Beispiel funktioniert jetzt schön und wir können beginnen zu sehen, wie channels etwas kompliziertere Logik machen können.

## 4. Einen Operator verwenden, um Eingabewerte aus einer CSV-Datei zu parsen

Jetzt werden wir versuchen, dies mit einer Datei mit einer Reihe von Eingaben zu tun. Das ist eine sehr übliche Art, Nextflow-Pipelines zu schreiben, indem man ein Sample Sheet oder ein CSV mit Metadaten verwendet.

## 4.1. Das Skript modifizieren, um eine CSV-Datei als Quelle der Grüße zu erwarten

Wenn ich zur Seitenleiste gehe, kannst du greetings.csv im Beispiel-Repository sehen, und das ist eine sehr, sehr einfache CSV-Datei, die nur drei Zeilen mit drei verschiedenen Grüßen enthält. Lass uns sehen, ob wir diese CSV-Datei innerhalb unseres Workflows verwenden können.

Ich werde jetzt wieder params verwenden, wie wir es in Kapitel eins getan haben, damit wir eine Kommandozeilen-Eingabe haben können.

Ich werde dieses greetings-Array löschen.

## 4.1.1. Den Eingabeparameter auf die CSV-Datei umstellen

Ich setze params greeting auf den Dateinamen, der greetings.csv ist, und ich werde diese spezielle Variable verwenden, um den channel zu generieren. Ich werde das da reinstecken, und die Fehler verschwinden. Denk daran, dass dies diese Variable jetzt standardmäßig setzt. Wenn ich also die Pipeline ohne Argumente ausführe, wird sie greetings.csv verwenden, aber ich könnte --greeting machen, um diese Variable zu überschreiben, wenn ich wollte.

## 4.1.2. Zu einer channel factory wechseln, die für eine Datei ausgelegt ist

Okay, wir übergeben jetzt eine Datei statt eines Strings oder eines Arrays von Strings, also brauchen wir wahrscheinlich eine andere channel factory.

Wir werden "of" loswerden, das wir bisher verwendet haben, und stattdessen .fromPath verwenden. Dies tut genau das, wonach es klingt. Es erstellt einen channel mit Pfaden anstelle von Werten, unter Verwendung eines String-Dateinamens oder Globs. Ich werde auch den flatten-Operator entfernen, da wir diesen nicht mehr benötigen, jetzt, da wir eine Datei übergeben.

## 4.1.3. Den Workflow ausführen

Ich werde speichern, das Terminal öffnen, den Workflow ausführen und dann sehen, was passiert.

Okay. Es ist wieder abgestürzt. Keine Sorge. Ich habe diesen auch erwartet. Schauen wir uns die Fehlermeldung an und sehen, ob wir herausfinden können, was schief läuft. Hier können wir den ausgeführten Befehl sehen, und ein bisschen wie zuvor, wo wir das ganze Array gedruckt hatten. Jetzt wird der Dateipfad in den Befehl geechot, anstatt durch den Inhalt der Datei zu gehen.

## 4.2. Den splitCsv()-Operator verwenden, um die Datei zu parsen

Um also den Inhalt der Datei stattdessen zu verwenden, brauchen wir einen anderen Operator. Der Operator, den wir dafür verwenden werden, heißt splitCsv. Macht Sinn, weil es eine CSV-Datei ist, die wir laden.

## 4.2.1. splitCsv() auf den channel anwenden

Ok, also splitCsv. Klammer schließen. Wir brauchen hier keine Argumente. Und wieder werde ich einige view-Operatoren verwenden, um einen Einblick zu geben, was hier vor sich geht.

.view csv after splitCsv. Before split Csv.

## 4.2.2. Den Workflow erneut ausführen

Okay, lass uns das ausführen und sehen, was passiert.

Okay, wir haben diesmal etwas mehr Ausgabe, aber es ist immer noch fehlgeschlagen. Wir können uns die view-Anweisungen ansehen, und hier kannst du before split CSV sehen, und wir haben einen Dateipfad, wie wir in der vorherigen Fehlermeldung gesehen haben. After split CSV haben wir jetzt drei Werte, die den drei Zeilen in der CSV-Datei entsprechen.

Du kannst jedoch sehen, dass jeder dieser Werte von eckigen Klammern umgeben ist. Jeder von ihnen war also ein eigenes Array, und das hat uns dasselbe Gebiet gegeben, das wir vorher hatten, wo es versucht, ein Array zu echoen, anstatt nur einen einzelnen String.

Wenn wir an eine CSV-Datei denken, macht das irgendwie Sinn. Typischerweise hat eine CSV-Datei Zeilen und Spalten, also macht split CSV ein zweidimensionales Array. Die erste Dimension des Arrays ist jede Zeile, und dann gibt es eine zweite Dimension, die jede Spalte für jede Zeile ist.

Hier haben wir also nur einen einzigen Wert auf jeder Zeile, also haben wir eine einzige Spalte, also haben wir ein Ein-Element-Array für jede Zeile der Datei.

Das ist in Ordnung. Wir brauchen nur einen weiteren Operator, um dieses Array für jede Zeile der geparsten CSV-Datei zu kollabieren. Lass uns das aufräumen. Werde das Terminal los und sehen, was wir tun können.

## 4.3. Den map()-Operator verwenden, um die Grüße zu extrahieren

Jetzt könnten wir den flatten-Operator wieder verwenden, den wir vorher verwendet haben. Wir haben gesehen, wie dieser ein Array in eine Reihe von Werten kollabieren kann, was hier sehr gut funktionieren würde. Aber ich werde die Gelegenheit nutzen, um einen anderen Operator zu demonstrieren, der innerhalb von Workflows sehr häufig ist, genannt der map-Operator.

## 4.3.1. map() auf den channel anwenden

Ich werde dot map machen und ich werde item item[0] machen.

Wenn du viel Code in anderen Sprachen schreibst, bist du vielleicht bereits mit dem map-Operator vertraut. Er nimmt ein iterierbares Objekt, wie ein Array oder einen channel, und er führt eine Operation auf jedem Wert davon aus.

Hier sagen wir, dass wir eine Variable namens item innerhalb des Geltungsbereichs dieser closure definieren sollten, und dann wollen wir nur den ersten Wert in diesem Array zurückgeben. Also item Index null.

Dies flacht das Array effektiv ab. Du kannst sehen, wie wir dies erweitern könnten, um komplexer zu sein: Wenn unsere CSV-Datei sechs Spalten hätte, wir aber nur an der vierten Spalte interessiert sind, könnten wir hier auf einen bestimmten Index zugreifen. Oder jede andere Art von Operation auf dem Wert durchführen, bevor wir ihn an die nachgelagerte Verarbeitung weitergeben.

Der map-Operator ist also extrem flexibel und sehr leistungsfähig zum Modifizieren von channels im Flug. Lass uns eine weitere view-Anweisung einbauen, nur damit wir sehen können, was sie in unserer Ausführung macht. Kann diese Zeile nehmen und nach unten verschieben. Und after map.

## 4.3.2. Den Workflow ein letztes Mal ausführen

Lass uns das Terminal aufrufen und versuchen, den Workflow auszuführen.

Okay, diesmal keine Fehler. Das ist ein gutes Zeichen. Wir können jetzt alle diese verschiedenen Ausgaben von den view-Anweisungen durchgehen. Before split CSV hatten wir einen einzelnen Pfad. After split CSV hatten wir die Ein-Wert-Arrays, und dann after map haben wir nur die Werte ohne Array-Syntax. Lass uns zum results-Verzeichnis gehen, und hier sind unsere Ausgabedateien, die sich genau so verhalten, wie wir es wollten.

Es gibt hier einen kleinen Bonus. Du kannst tatsächlich sehen, dass die view-Operatoren in der Reihenfolge, in der sie die Ausgabe gemacht haben, leicht durcheinander sind. Dies liegt daran, dass Nextflow die Parallelisierung dieser verschiedenen Aufgaben durchführt. Nachdem es also das CSV gesplittet hat, gibt es drei Elemente in diesem channel, und es handhabt die Verarbeitung dieser drei Elemente automatisch parallel. Das bedeutet, dass die Reihenfolge der Ausgaben stochastisch ist und variieren kann. In diesem Fall geschah es nur, dass einige der view-Operatoren zurückkehrten, nachdem der nachfolgende Schritt abgeschlossen war, und so kam es in dieser Reihenfolge.

Wenn ich denselben Workflow erneut ausführe. Dann ist es tatsächlich in einer anderen Reihenfolge gekommen, und diesmal haben wir die split CSVs und die maps in der Reihenfolge, die wir erwarten würden.

Bedenke also, du kannst dich nicht auf die Reihenfolge der Ausgaben einer process-Aufgabe verlassen, weil Nextflow diese Parallelisierung automatisch für dich handhabt. Nextflow macht das für dich mit seiner Datenfluss-Logik, und das ist die wahre Stärke von Nextflow.

Okay, das ist wahrscheinlich eines der wichtigsten Kapitel des gesamten Trainings. Sobald du channels, channel factories und Operatoren verstehst, beginnst du, die Stärke von Nextflow zu verstehen und was es als Programmiersprache einzigartig macht. Diese Funktionalität ermöglicht es Nextflow, alle deine Workflows für dich zu parallelisieren und extrem komplexe Workflow-Logik mit einer sehr sauberen Syntax und einem Push-Datenfluss-Modell zu generieren. Es kann zunächst ein etwas seltsames Konzept sein, aber sobald du dich daran gewöhnt hast, Code wie diesen zu schreiben, wird es sich schnell natürlich anfühlen, und bevor du es weißt, wirst du fantastische Workflows schreiben.

Mach eine Pause, eine Tasse Tee, geh ein bisschen herum und lass uns zu Kapitel drei übergehen, wo wir beginnen, diese Konzepte in komplexere Workflows zu erweitern. Wir sehen uns im nächsten Video.

[Nächstes Video-Transkript :octicons-arrow-right-24:](03_hello_workflow.md)
