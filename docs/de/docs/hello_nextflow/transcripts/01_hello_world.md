# Teil 1: Hello World - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../01_hello_world.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern des Materials.

## Willkommen

Hallo und willkommen zurück.

Du bist jetzt in Teil 1 des "Hello Nextflow"-Kurses mit dem Titel "Hello World". In diesem Kapitel werden wir beginnen, die absoluten Grundlagen von Nextflow zu verstehen.

Hoffentlich bist du jetzt in Codespaces oder etwas Ähnlichem mit VS Code eingerichtet, und du hast deinen Hello Nextflow-Ordner im Workspace im Explorer mit all diesen verschiedenen Dateien hier.

Wir fangen damit an, einfach ein paar grundlegende Dinge im Terminal mit Bash zu machen, und dann sehen wir, ob wir die gleichen Dinge in Nextflow machen können, damit du ein Gefühl für die Syntax bekommst.

## 0. Aufwärmen

Also, fangen wir wirklich einfach an. Beginnen wir einfach mit "echo", um etwas ins Terminal zu drucken. "Hello World". Ich drücke Enter und das geht ins Terminal. Hello World. Hoffentlich ist das für niemanden, der diesen Kurs anschaut, eine Überraschung.

Okay, lass uns etwas damit machen. Anstatt es nur ins Terminal zu drucken, schreiben wir es in eine Datei. Ich drücke die Aufwärtspfeiltaste meiner Tastatur, die durch den Bash-Verlauf blättert, also gibt sie mir meinen letzten Befehl, und ich hänge am Ende ein kleines Größer-als-Symbol dran, das die Ausgabe von diesem Befehl in eine Datei umleitet, und ich nenne sie output.txt.

Enter nochmal, um den Befehl auszuführen, diesmal nichts im Terminal, aber wir können auf der linken Seite sehen, dass die neue Datei hier aufgetaucht ist, genannt output.txt.

Wir können das in einem Terminal mit etwas wie cat ansehen. Also cat output.txt und tatsächlich steht dort "Hello World". Wir können auch doppelt drauf klicken und es öffnet sich im Code-Editor in VS Code.

## 1.1. Den Code untersuchen

Gut. Ich sagte dir, es war einfach. Was kommt als Nächstes? Lass uns versuchen, diesen Prozess zu wiederholen, aber diesmal machen wir es in Nextflow.

Wie ich sagte, beginnen alle verschiedenen Kapitel in diesem Kurs mit einem Skript und dieses heißt Hello World. Also werde ich Hello World finden. Es zeigt eine Vorschau, wenn ich es einmal anklicke, ich klicke doppelt, um es im Editor hier zu öffnen. Und ich werde kurz das Terminal loswerden.

Nun, das ist ein sehr einfaches Skript, also etwa so einfach wie es nur geht. Es ist nur 22 Zeilen lang und macht im Grunde das Gleiche. Tatsächlich sollte einiges davon vertraut aussehen. Es ist das, was wir gerade getippt haben. Wir können unseren Bash-Befehl sehen, der in eine Datei umleitet.

Okay. Was noch? Außerdem können wir in dieser Datei einige der Kernkonzepte von Nextflow sehen. Wir haben hier einen Prozess in Rot und einen Workflow. Das sind die speziellen Schlüsselwörter und spezielle Terminologie in Nextflow.

## 1.1.1. Die Prozess-Definition

Verschiedene Prozesse innerhalb eines Workflows umschließen verschiedene logische Einheiten deines Workflows. Jeder Prozess macht eine Sache.

Wenn wir es ausführen, generiert es eine Aufgabe oder mehrere Aufgaben, die tatsächliche Schritte einer Pipeline ausführen. Alle Prozesse werden dann innerhalb eines Workflow-Blocks orchestriert, den wir unten sehen, und in diesem Fall wird einfach dieser eine Prozess ausgeführt.

Der Prozessname folgt auf dieses Schlüsselwort hier, und das kann im Grunde alles sein. Und dann sind die Inhalte des Prozesses innerhalb dieser geschweiften Klammern.

Es gibt wirklich nur eine Anforderung für einen Prozess, nämlich dass er irgendeine Art von Script- oder Exec-Block enthält. Das ist in den dreifachen Anführungszeichen hier, und das ist das Bash-Skript, das in das Arbeitsverzeichnis geschrieben wird, wenn wir die Pipeline ausführen, und das ist das, was tatsächlich auf deinem Computer oder Server läuft.

Das ist typischerweise Bash, aber du kannst auch hier oben ein anderes Hashbang einfügen, und es könnte ein Python-Skript oder ein R-Skript sein. Es spielt keine Rolle. Was auch immer in diesem Skript steht, wird ausgeführt.

Es gibt noch eine andere Sache, die wir diesem Prozess hier hinzugefügt haben, nämlich die Output-Deklaration. Das teilt Nextflow mit, dass dieser Prozess eine Ausgabedatei namens output.txt erwartet. Es heißt, es ist ein Pfad, also sollte es wie eine Datei behandelt werden, nicht sagen wir, wenn das val wäre, würde es sagen, es ist wie eine Variable oder ein Wert.

Beachte, dass dies die Datei nicht erstellt. Es generiert sie nicht tatsächlich. Das wird vom Skript hier unten gemacht. Es teilt Nextflow nur mit, dass es eine Ausgabedatei mit diesem Dateinamen erwarten soll.

## 1.1.2. Die Workflow-Definition

Okay. Und dann haben wir unten einen Workflow, und wieder haben wir eine Deklaration. Diese heißt Main. Das ist das Workflow-Äquivalent zu einem Script-Block, sozusagen. Es ist der Teil des Workflows, der etwas tut. Und in diesem Fall sagen wir, rufe den Prozess namens sayHello auf.

Normalerweise wird deine Pipeline natürlich viel komplexer aussehen als das. Du wirst wahrscheinlich mehr als einen Prozess haben, und du wirst Kanäle verwenden, um den Datenfluss zwischen ihnen zu orchestrieren. Darauf kommen wir in den nächsten Teilen dieses Kurses, aber fürs Erste ist das genug. Das ist eine gültige Pipeline, die funktionieren sollte.

Ich kann sogar hier in VS Code auf "DAG-Vorschau" klicken. Der DAG ist eine Darstellung einer Datenflussstruktur in der Pipeline, und wir können ihn seitlich als Mermaid-Diagramm gerendert sehen. In diesem Fall ist es sehr einfach. Es gibt eine Box, die den Workflow darstellt, und einen Prozess namens sayHello, aber das könnte interessanter aussehen, wenn wir weitermachen.

## 1.2. Den Workflow ausführen

Okay, lass uns versuchen, diesen Workflow auszuführen und sehen, was passiert.

Ich hole das Terminal wieder unten hervor, lösche die Ausgabe und tippe Nextflow Run. Und dann tippe ich einfach den Skriptnamen ein, der hello-world.nf ist. Und ich drücke Enter.

Okay, es gibt einige Standardsachen oben, die uns sagen, dass Nextflow lief und welche Version lief und wie der Skriptname war und alles.

Und wirklich das Wichtige, wonach wir hier suchen, ist _hier_, eine Zusammenfassung der verschiedenen Aufgaben, die ausgeführt wurden.

Wenn deines so aussieht mit einem kleinen grünen Häkchen, dann gut gemacht. Du hast gerade deine erste Pipeline ausgeführt. Fantastisch.

Es sagt uns hier den Namen des Prozesses, der ausgeführt wurde, der Say Hello hieß, und es sagt uns, dass er einmal lief und dass er erfolgreich war. Das aktualisiert sich während du weitergehst, also wenn du eine größere Pipeline ausführst, siehst du den Fortschritt hier dargestellt. Aber weil das so winzig ist, läuft es im Grunde sofort.

## 1.2.2. Die Ausgabe und Logs im Work-Verzeichnis finden

Wenn du jetzt eine Nextflow-Pipeline ausführst, wird jeder dieser Prozesse zusammengefügt, und jeder Prozess, wie ich vorhin sagte, kann Aufgaben generieren, eine oder mehrere. In diesem Fall hatten wir also eine einzelne Aufgabe aus diesem Prozess. Sie wurde nur einmal ausgeführt und das geschah unter diesem Aufgaben-_Hash_.

Nextflow arbeitet nicht direkt mit den Dateien in deinem Arbeitsverzeichnis, es erstellt einen speziellen Ordner namens work. Und wenn ich "ls" mache, werden wir sehen, dass er hier erschienen ist: _work_, und darin gibt es Unterverzeichnisse für jede einzelne Aufgabe, die läuft. Und das stimmt mit diesem Hash überein. Du kannst sehen, wenn ich zu "ls work/c4" gehe, und dann ist es abgeschnitten, aber es beginnt mit 203, und das ist das Arbeitsverzeichnis, das von diesem Prozess erstellt wurde, als wir die Pipeline ausführten. Und du kannst es auch an der Seite sehen.

Wenn ich diese Dateien aufliste, kannst du sehen, dass die Datei output.txt generiert wurde. Du kannst sie auch hier sehen. Und es gibt eine Reihe versteckter Dateien, die bei meinem normalen "ls" nicht angezeigt werden.

Wenn ich auf output.txt klicke, haben wir tatsächlich unsere Ausgabe. Fantastisch. Also hat die Pipeline funktioniert.

Es mag wie ziemlich viel Boilerplate erscheinen, um im Wesentlichen ein einzeiliges Bash-Skript auszuführen, aber es wird mehr Sinn ergeben, wenn unsere Prozesse komplizierter werden. Und dieses Work-Verzeichnis mit Nextflow und diese Dateien, die erstellt werden, sind wirklich das Rückgrat dessen, was Nextflow so leistungsfähig macht.

Jede Aufgabe, jedes Element einer Pipeline ist von jeder anderen Aufgabe isoliert. Es ist reproduzierbar. Sie kollidieren nicht miteinander, und alles kann parallel laufen. Es ist tatsächlich eine wirklich schöne Art, wenn man sich daran gewöhnt hat, wegen dieser Isolation, dass du hineingehen und genau sehen kannst, was für eine einzelne Aufgabe passiert ist, und debuggen kannst.

Lass uns einen kurzen Blick auf diese anderen Dateien im Work-Verzeichnis werfen. Von oben nach unten haben wir eine Datei namens _.command.begin_. Diese ist leer. Es ist nur eine sogenannte Sentinel-Datei, die von Nextflow erstellt wurde und sagt, okay, ich starte die Aufgabe. Nichts Interessantes da.

Dann gibt es _.command.error_, _.command.log_ und _.command.out_. Das sind alles Ausgaben vom Bash-Befehl oder diesem Skript, das lief. Das ist Standard-Error. Das ist Standard-Out, und das sind die beiden zusammen, wie sie herauskamen. Du bekommst also die logische Reihenfolge.

Okay, die waren auch alle leer, also nicht sehr interessant, aber es wird interessanter, wenn du zu _.command.run_ kommst.

Das ist typischerweise ein sehr langes Skript. Und das ist, was Nextflow tatsächlich ausführt. Wenn du hier hineingehst, wirst du die ganze innere Logik von Nextflow sehen und sehen, was es tut und wie es deinen Prozess ausführt. Das hängt davon ab, wo du läufst, ob wir lokal laufen oder es als Job an SLURM senden, in welchem Fall wir SLURM-Header oben haben werden. All diese verschiedenen Setups.

Generell musst du aber nie wirklich in diese Datei schauen. Sie wird von Nextflow automatisch generiert und es gibt nichts wirklich besonders Einzigartiges für deine Pipeline, was darin ist. Aber das ist wirklich der Kern dessen, was läuft.

Die nächste ist viel interessanter. _.command.sh_ ist das generierte Skript, das aus deinem Prozess kam, und hier kannst du sehen, dass Nextflow den Bash-Header hinzugefügt hat, und dann hat es unseren Befehl ausgeführt, der in unserem Script-Block war.

Und das ist alles, was die _.command.run_-Datei macht, sie führt einfach diese _.command.sh_-Datei aus.

Das ist eine wirklich nützliche, die normalerweise am meisten angeschaut wird, wenn du versuchst, etwas zu debuggen und zu überprüfen, dass die Logik deiner Nextflow-Pipeline das tut, was du erwartest.

Schließlich haben wir eine Datei namens _.exitcode_, und die erfasst einfach den Exit-Code von einer Aufgabe, der in diesem Fall erfolgreich war. Also war der Exit-Code null.

Wenn etwas schief geht, dir der Speicher ausgeht oder etwas anderes und es fehlschlägt, dann ist das sehr nützlich zu verstehen, was schief gelaufen ist.

## 1.3. Den Workflow erneut ausführen

Noch eine Sache, die man über Work-Verzeichnisse verstehen sollte, ist, dass wenn ich diese Pipeline wiederholt ausführe, also wenn ich _"nextflow run hello-world.nf"_ mache, wird sie genau das Gleiche tun, aber diesmal wird sie eine neue Task-ID haben. Du kannst sehen, dass dieser Hash hier anders ist, und wenn ich jetzt in work schaue, gibt es zwei Hash-Verzeichnisse. Und diese sind wieder getrennt voneinander.

Also jedes Mal, wenn du einen Nextflow-Workflow ausführst, es sei denn, du verwendest resume, das den Cache nutzt, darauf kommen wir später, wird es diese Prozesse in neuen Work-Verzeichnissen erneut ausführen, die voneinander getrennt sind. Du bekommst keine Dateinamenskollisionen, du wirst keine Probleme wie das haben. Alles ist isoliert und sauber.

Und wenn wir in dieses Verzeichnis gehen, kannst du all die gleichen Dateien und die gleiche _output.txt_ sehen, die von Grund auf neu erstellt wurde.

## 2. Ausgaben veröffentlichen

Okay, das ist großartig für Nextflow selbst, während es deine Pipeline ausführt, damit alle Dinge voneinander getrennt und sauber sind und verwaltet werden können.

Aber es ist nicht besonders nützlich, wenn du eine Person bist, die versucht, deine Ergebnisse zu erkunden. Du möchtest nicht wirklich durch Tausende und Abertausende verschiedener Work-Verzeichnisse wühlen, um deine Ergebnisdateien zu finden. Und das sollst du auch nicht. Die Work-Verzeichnisse sind nicht als endgültiger Zustand gedacht, wo deine Dateien erstellt werden.

Wir tun dies, indem wir unsere Dateien veröffentlichen.

## 2.1.1. Die Ausgabe des sayHello-Prozesses deklarieren

Wenn ich also zurück zu unserem Skript gehe, werden wir in unserem Workflow-Block hier arbeiten. Wir werden ihm sagen, welche Dateien zu erwarten sind, welche Dateien uns wichtig sind, und dann werden wir einen neuen Block darunter erstellen, der Output-Block genannt wird.

Das ist die neue Syntax, die mit dem Syntax-Parser kam und standardmäßig in Version 26.04 von Nextflow ist. Wenn du also Nextflow schon ein bisschen benutzt hast, ist das eine der Sachen, die neu sind.

Also haben wir den Main-Block, und als Nächstes sage ich publish und ich sage Nextflow, was vom Veröffentlichen zu erwarten ist. Wir nennen es _first_output_, und wir nennen es _sayHello.out_.

Ich habe dort versehentlich einen Tippfehler gemacht, aber das ist eine gute Gelegenheit, auch auf einige der Funktionen der Nextflow VS Code-Extension hinzuweisen. Du kannst sehen, dass es mir sofort eine kleine gewellte rote Linie darunter gegeben hat, die sagt, dass etwas nicht stimmt. Und wenn ich darüber fahre, wird es mir sagen, dass diese Variable nicht definiert ist. Ich weiß nicht, was es ist.

Es ist in diesem Fall ziemlich offensichtlich, ich habe einen Tippfehler gemacht. Ich wollte sayHello tippen, und dann verschwindet die gewellte Linie.

Jetzt ist es lila. Der Nextflow-Syntax-Parser weiß, dass das ein Prozess ist, und wenn ich darüber fahre, gibt er mir eine reduzierte Darstellung, wie dieser Prozess aussieht. Ich kann also sehr schnell auf einen Blick sehen, dass er keine Eingaben nimmt und uns diese Ausgabe gibt. Also gibt dir die Arbeit in VS Code mit dieser Extension viele kontextuelle Informationen, während du Code schreibst.

Beachte, dass wir auf die Ausgabe von diesem Prozess mit der _.out_-Syntax verweisen können. Und im Moment können wir das nennen, wie wir wollen, es ist nur ein beliebiger Variablenname.

## 2.1.2. Einen output:-Block zum Skript hinzufügen

Wo es wichtig wird, ist, wenn wir unseren neuen Block hier machen, und das ist jetzt unter dem Workflow-Block, wir sind nicht mehr innerhalb des Workflows. Wieder geschweifte Klammern. Und hier sagen wir Nextflow einfach, wohin alle Dateien gelegt werden sollen, die vom Workflow erstellt werden.

Jetzt nehme ich diesen Variablennamen, den ich hier erstellt habe, und ich werde ihn dort hinsetzen und einige geschweifte Klammern dafür setzen. Und ich sage Nextflow, einen Pfad zu verwenden. Ups. Path, in Anführungszeichen. Und ich werde einen Punkt verwenden. Das sagt Nextflow einfach, die Datei in die Wurzel des Ergebnisverzeichnisses zu legen. Also keine Unterverzeichnisse oder so.

Lass uns unseren Workflow noch einmal ausführen. Wenn ich _"nextflow run hello-world.nf"_ mache, sollte es hoffentlich im Grunde genau gleich aussehen. Nichts hat sich hier wirklich bei Nextflow geändert. Es führt die gleichen Dinge aus. Es macht sie nur wieder in Work-Verzeichnissen.

Aber jetzt, wenn ich _"ls results/"_ mache, wirst du sehen, dass es hier ein neues Verzeichnis gibt, das erstellt wurde und results heißt, das das Standard-Basisverzeichnis für Workflow-Veröffentlichungen ist. Und darin ist eine Datei namens _output.txt_.

Wenn ich _"ls -l results"_ mache, wirst du sehen, dass das tatsächlich eine Softlink zum Work-Verzeichnis ist. Das ist also keine echte Datei, sie ist mit dem Work-Verzeichnis verlinkt und hat alle Dateien dort für uns gesammelt.

## 2.2. Einen benutzerdefinierten Speicherort festlegen

"Results" ist der Standardname für diesen Pfad. Wenn ich den Workflow erneut ausführe, und diesmal mache ich _dash_ einzelner Bindestrich, das ist, weil es eine Kern-Nextflow-Option ist. _"-output-dir **my**results"_. Könnte auch einfach _"-o"_ kurz machen. Dann wird es ein anderes Basisverzeichnis festlegen, wo die Dateien gespeichert werden, und noch einmal, hier oben in _myresults/_, haben wir jetzt eine _output.txt_.

Das ist großartig, aber wir wollen wahrscheinlich nicht alle Dateien nur in der Wurzel. Wir wollen etwas Organisation, also können wir hier auch ein Unterverzeichnis erstellen, das wir nennen, wie wir wollen. Sagen wir _"path 'hello_world'"_, und ich führe das einfach nochmal aus. _"nextflow run hello-world.nf"_. Es sollte ins Ergebnisverzeichnis in ein Unterverzeichnis gehen und tatsächlich, jetzt unter results hier oben haben wir _hello_world/_ und wir haben _output.txt_.

Wichtig zu bemerken, die alte _output.txt_-Datei ist immer noch da. Das Ergebnisverzeichnis wird nicht gelöscht, wenn du das machst. Nur neue Dateien werden dort hineinkopiert. Sie überschreiben Dateien, die bereits da sind, wenn sie den gleichen Dateinamen haben, aber sie löschen keine alten. Du musst also ein bisschen vorsichtig sein, wenn du Pipelines erneut ausführst. Wenn du nicht willst, dass sie über den Dateien sind, die bereits da sind, stelle sicher, dass du ein leeres Verzeichnis verwendest.

## 2.3. Den Veröffentlichungsmodus auf copy setzen

Okay, ich habe erwähnt, dass diese Dateien Softlinks sind, also wenn ich _"ls -l results/hello_world/"_ mache, kannst du sehen, dass es auf das Work-Verzeichnis soft-linkt. Das ist generell eine gute Sache, wenn du an etwas wie HPC arbeitest, und das sind wirklich riesige Dateien und du willst sie nicht duplizieren, weil es bedeutet, dass die Dateien nur einmal im Dateisystem gespeichert werden.

Allerdings bedeutet es, dass wenn du das Work-Verzeichnis löschst: wenn ich _"rm -r work"_ mache und alle diese Zwischendateien lösche, die erstellt wurden. Jetzt, wenn ich versuche, diese Datei _"results/hello_world/"_ zu lesen. Es wird als Softlink auf eine Datei zeigen, die nicht mehr existiert, und die Daten sind für immer verloren und unwiederbringlich, was vielleicht nicht großartig ist.

Also generell sage ich, es ist gute Praxis, die Dateien zu kopieren anstatt soft-zu-linken, wenn du kannst, weil es sicherer ist. Sei dir nur bewusst, dass es doppelt so viel Festplattenspeicher verwenden wird, es sei denn, du löschst diese Work-Verzeichnisse.

Um das mit dem Output-Block zu machen, gehe ich zum ersten Output hier. Ich habe vorher den Pfad gesetzt und jetzt werde ich den Modus setzen und du kannst sehen, während ich tippe, schlägt die VS Code-Extension Sachen vor, sie weiß, dass es eine Output-Direktive hier ist. Und ich sage copy. Ich speichere.

Lass uns den Workflow neu ausführen. Es wird die Dateien wieder erstellen, neues Work-Verzeichnis.

Jetzt, wenn ich zu _"ls -l results/hello_world/"_ gehe, kannst du sehen, das ist eine echte Datei und es ist kein Softlink mehr, und Nextflow hat das kopiert. Gut zu wissen. Also path und mode sind Dinge, die du ziemlich oft schreiben wirst.

Nun, natürlich ist das sehr einfach. Wir werden das komplexer und mächtiger machen, wenn wir weitergehen, und du wirst sehen, wie man diese Dinge dynamisch macht und nicht zu ausführlich.

## 2.4. Hinweis zu publishDir-Direktiven auf Prozessebene

Nun, ich sagte, als wir damit anfingen, dass dies eine ziemlich neue Form der Syntax ist. Sie ist nur in den neuesten Versionen von Nextflow verfügbar, während ich das aufnehme, und sie heißt Workflow Outputs.

Wenn du das verwendest, ist es großartig. Es schaltet viele andere coole Funktionen innerhalb von Nextflow frei, wie zum Beispiel Nextflow Lineage, um die Herkunft dieser Dateien zu verfolgen, während sie erstellt werden, und bald wird das der Standard in 26.04 sein. Und zu einem späteren Zeitpunkt in der Zukunft wird dies die einzige Möglichkeit sein, deine Workflows zu schreiben.

Da wir uns jedoch gerade in dieser Übergangsphase befinden, könntest du durchaus Pipelines in freier Wildbahn sehen, die etwas namens publishDir verwenden, was der alte Weg ist, es zu tun, und das ist nicht auf Workflow- und Output-Ebene definiert, sondern auf Prozessebene.

Und diese Deklaration sagt im Grunde das Gleiche. Sie sagt, veröffentliche die Ergebnisdateien in ein Verzeichnis namens results, und verwende einen Copy-Modus. Du kannst also sehen, dass die Syntax sehr ähnlich ist. Aber wenn du jetzt neue Pipelines schreibst, versuche nicht, diese publishDir-Direktive zu verwenden, selbst wenn du sie in KI-Ergebnissen oder in Dokumentation oder anderen Pipelines siehst, denn das ist der alte Weg, es zu tun.

Im Jahr 2026 sollten wir alle Workflow Outputs verwenden.

Das ist alles dokumentiert, wenn du das machst und Nextflow schon mal benutzt hast, kannst du zu den Nextflow-Docs hier gehen, nextflow.io/docs/. Und wenn ich zu Tutorials scrolle, gibt es ein Tutorial namens _Migrating to Workflow Outputs_.

Es ist wirklich gut. Es geht durch die ganze Syntax, wie es der alten Syntax entspricht, warum wir es geändert haben, und hat eine Zeitleiste und alles. Und es geht durch alle verschiedenen Szenarien mit Unmengen von Beispielen. Du kannst also bestehenden Nextflow-Code leicht zur neuen Syntax konvertieren.

## 3.1. Den sayHello-Prozess ändern, um eine variable Eingabe zu erwarten

Okay, wir haben also unser einfaches Skript, das einen Prozess ausführt, eine Datei erstellt, Nextflow mitteilt, dass es eine Ausgabe ist, und dann sagen wir Nextflow, wo diese Datei gespeichert werden soll. Das ist ein guter Start.

Aber es wäre interessanter, wenn es nicht alles hartcodiert wäre. Als Nächstes überlegen wir also, wie wir Nextflow mitteilen können, dass dieser Prozess eine variable Eingabe annehmen kann, die wir zur Laufzeit kontrollieren können, wenn wir einen Workflow starten.

Wir müssen ein paar verschiedene Dinge tun, um das möglich zu machen.

Zuerst müssen wir diesem Prozess mitteilen, dass er eine Eingabevariable akzeptieren kann, und wir tippen hier _input_ als neuen Deklarationsblock. Und wir nennen das _"val greeting"_.

Das val-Bit ist das Äquivalent zu einem path hier unten. Es sagt Nextflow, dass dies eine Variable ist, in diesem Fall wie ein String. Und wenn du darüber fährst, sagt dir die Extension wieder, was das bedeutet.

Als Nächstes sagen wir Nextflow, was damit zu tun ist. Es reicht nicht, nur zu sagen, dass es eine Variable gibt. Du musst im Skript sagen, wie diese Variable zu verwenden ist. Und ich werde also diesen hartcodierten String hier loswerden und eine Variable einfügen.

Ich mache es schnell ohne geschweifte Klammern, nur um dir zu zeigen, dass das erlaubt ist, und das ist die alte Art, es zu tun. Aber jetzt mit der neuen Syntax empfehlen wir wirklich, es in geschweifte Klammern wie diese zu setzen, und es macht es wirklich klar, dass das hier von Nextflow interpoliert wird.

Großartig. Also _"input greeting"_ geht in _$\{greeting\}._ Die letzte Sache ist, dass wir Nextflow auf Workflow-Ebene mitteilen müssen, dass dieser Prozess jetzt eine Eingabe nimmt. Und um das zu tun, geben wir ihm im Grunde eine Variable.

## 3.2. Einen Befehlszeilenparameter einrichten, um Benutzereingaben zu erfassen

Wir könnten es wieder hartcodieren, wie Hello World, und das würde gut funktionieren, aber offensichtlich gibt es uns keinen wirklichen Vorteil. Wir wollten es zur Laufzeit konfigurieren können, also wollen wir es auf der CLI machen können, wenn du Nextflow startest.

Und die Art, wie wir das machen, ist ein spezielles Nextflow-Konzept namens _params_. Wir nennen das _params.input_.

Was das macht, ist, dass es diese Eingabevariable auf der CLI verfügbar macht, und da verwenden wir einen doppelten Bindestrich, wenn wir Nextflow starten.

Ich kann das nennen, wie ich will, ich kann es _hello, greeting_ nennen. Spielt keine Rolle. Was auch immer ich dort mache, wird als CLI-Option verfügbar gemacht, wenn wir eine Pipeline starten. Und das ist ein echter Zaubertrick von Nextflow, weil es bedeutet, dass du dein Workflow-Skript sehr schnell mit diesen Parametern aufbauen kannst, und du baust im Wesentlichen eine benutzerdefinierte CLI für deine Pipeline auf, was es wirklich einfach macht, verschiedene Optionen spontan anzupassen, wenn du startest.

Also. Probieren wir es aus. Zurück zu unserem Terminal. Wir haben unseren _"nextflow run"_-Befehl hier. Und jetzt mache ich _"--input"_, was mit dem _"params.input"_ übereinstimmt, das wir vorher gesehen haben. Ich denke, in den Docs ist es auf Französisch. Geraldine spricht gerne Französisch. Ich mache es auf Schwedisch, weil ich in Schweden lebe. Also sage ich "_Hej Världen_" und drücke Enter.

Kann einzelne oder doppelte Anführungszeichen verwenden, es beeinflusst nur, wie Bash es interpretiert.

Es führt die Nextflow-Pipeline genau auf die gleiche Weise aus. Du kannst sehen, dass das Arbeitsverzeichnis und alles gleich ist. Aber jetzt, wenn ich zu _"results/hello_world/output"_ gehe. Können wir unser schönes Schwedisch hier stattdessen sehen.

Wir haben also dynamisch eine Eingabe von einer CLI an einen Parameter übergeben. Wir haben das als Eingabe an den Prozess übergeben, und der Prozess hat das interpretiert und in einen Script-Block eingefügt, der dann dynamisch die Ausgabe dieses Skriptergebnisses geändert hat. Ziemlich cool.

Ziemlich komplexe Logik mit sehr wenig Syntax hier. Und du kannst hoffentlich sehen, wie das jetzt skaliert. Und so bauen wir wirklich die Logik und die Anpassbarkeit unserer Pipelines in das Nextflow-Skript ein.

## 3.4. Standardwerte für Befehlszeilenparameter verwenden

Okay, das ist großartig. Das Problem ist jetzt jedoch, dass ich jedes Mal, wenn ich diese Pipeline ausführe, dash input machen muss, damit sie läuft.

Wenn ich versuche, ohne diesen Parameter auszuführen, wird Nextflow jetzt einen Fehler werfen und sagen, dass es diesen Parameter brauchte und er nicht gesetzt wurde. Und es wusste also nicht, was zu tun ist.

Das ist übrigens eine coole neue Sache. In der Vergangenheit hätte Nextflow einfach mit einem leeren String gelaufen, und du hättest alle Arten von seltsamen Fehlern gehabt, die schwer zu verstehen gewesen wären. Aber im neuen Nextflow-Syntax-Parser ist es etwas vorsichtiger und sagt es dir sofort.

Wir wollen also nicht immer jede einzelne Option angeben. Es ist gute Praxis, vernünftige Standards festzulegen. Wie machen wir das also in unserem Skript?

Du wirst bemerken, dass wir, als wir das schrieben, einfach _params.input_ direkt dort eingefügt haben, wo wir es verwenden. Die offensichtliche Lösung ist also, dass wir einen Standard definieren, und das machen wir hier oben im Skript in einem speziellen params-Block im Workflow. Das ist im Workflow-Skript hier.

Wieder etwas neue Syntax hier, also pass auf. Das ist wirklich cooles Zeug. Wir haben den Namen des Parameters, der hier erwartet wird.

Und dann nach diesem Doppelpunkt-Zeichen definieren wir einen Typ der Variable. Du musst das nicht tun, du kannst es einfach leer lassen, aber es ist wirklich schön. Es sagt Nextflow, dass wir einen String erwarten und ihn entsprechend behandeln.

Wenn wir stattdessen zum Beispiel eine Zahl wollen, könnten wir float schreiben, und das würde sagen, wir wollen eine Fließkommazahl. Und wenn wir versuchen, damit auszuführen, wird es einen Fehler werfen. Wenn wir ihm einen String geben, der kein Float ist. Und es wird ihn auch entsprechend übergeben. Bei string weiß es, dass es ein String ist. Und selbst wenn es führende Nullen hat und vollständig numerisch ist, wird es es trotzdem als tatsächlichen String übergeben.

Diese Typsicherheit ist also eine sehr neue Funktion von Nextflow, aber wirklich leistungsstark, um deinen Code sicherer zu schreiben und auszuführen.

Dann danach haben wir ein Gleichheitszeichen und dann den Standardwert hier. Nextflow wurde ursprünglich in Barcelona geschrieben, also scheint es angemessen, dass wir etwas Spanisches hier haben, _"Holà mundo!"_ als Standard.

Richtig, ich speichere das Skript, gehe zurück, führe das Skript wieder ohne _--input_ aus. Und diesmal sollte es laufen und es wird unsere neue Datei in _results_ erstellen. Und in dieser Datei steht jetzt _"Holà mundo!"_.

Das ist aber nur ein Standard, also bedeutet das nicht, dass wir nicht immer noch das Gleiche wie vorher tun können. Wenn ich zurückgehe und mein altes Skript hier finde, _"Hej Världen"_, weil ich _--input_ auf der Befehlszeile mache, wird das diesen Standard überschreiben und das wieder in der output.txt-Datei verwenden.

Das im Skript ist also nur der Standardwert, den ich setze.

Während wir unseren Workflow aufbauen, komplexer machen und mehr Parameter einfügen, wird sich dieser params-Block oben im Skript anfangen, sie alle an einem Ort zu sammeln.

Und du endest mit dieser ziemlich schönen Symmetrie in deinem Skript, wo du effektiv alle deine Workflow-Eingaben hier hast und deine Workflow-Ausgaben unten. Und es ist sehr klar, was die Schnittstelle deines Workflows zur Außenwelt ist. Du kannst also sehr schnell eine neue Pipeline mit der neuen Syntax aufgreifen und verstehen, wie man sie benutzt.

Eine letzte coole Sache. Wir müssen keinen Standardwert damit setzen. Wenn wir params input machen, aber keinen Standardwert setzen, dann sagt das Nextflow, dass dieser Parameter erforderlich ist, und wieder wird die Pipeline ohne ihn nicht laufen, aber sie gibt dir eine nützlichere Fehlermeldung als etwas über null sein.

Also sagt es, wir erwarten, seine Eingabe ist erforderlich, aber sie wurde nicht auf der Befehlszeile angegeben. Sehr schön.

Okay, hoffentlich ist es jetzt also klar, wie man seine Nextflow-Pipeline mit variablen Eingaben und Parametern einrichtet, wie man den Standard setzt, die Typen setzt, es könnte ein Boolean true false Flag sein oder eine Ganzzahl oder verschiedene Typen hier. Wie man sie in deinen Workflow übergibt, wo es durchgeht, und dann in deinen Prozess interpoliert. Und du weißt auch, wie man diese auf der Befehlszeile anpasst, wenn du Nextflow startest. Das fängt an, interessanter auszusehen als unser einfacher Bash-Befehl.

## 4. Workflow-Ausführungen verwalten

Okay. Was kommt als Nächstes? Für den letzten Teil dieses Kapitels werden wir ein wenig darüber sprechen, wie man all die verschiedenen Workflow-Ausführungen verwaltet. Wenn du in meiner Seitenleiste hier und im Explorer unter work schaust, wirst du sehen, dass ich eine Reihe verschiedener Pipelines ausgeführt habe und diese Work-Verzeichnisse ziemlich lang werden, es gibt viele davon.

Und die andere Sache ist, wie ich vorher gesagt habe, jedes Mal, wenn ich diese Pipeline erneut ausführe, erstellt sie einen neuen Satz von Work-Verzeichnissen, und sie führt alle Prozesse von Grund auf neu aus, was eine gute Sache ist. Das ist beabsichtigtes Verhalten. Es ist reproduzierbar und regeneriert alles frisch. Aber es ist offensichtlich, wenn du sehr lange laufende Prozesse ausführst, ärgerlich, immer deine Pipeline von Anfang an starten zu müssen, wenn sie auf halbem Weg abgestürzt ist oder wenn du etwas am Ende der Pipeline änderst.

## 4.1. Einen Workflow mit -resume neu starten

Zum Glück ist Nextflow wirklich gut darin zu wissen, was vorher gelaufen ist und was verfügbar ist, und diese alten Ergebnisse wiederzuverwenden ist sehr einfach. Wir fügen einfach ein neues Flag am Ende des Befehls hinzu _"-resume"_.

Beachte jetzt, dass es zwei Bindestriche bei input gibt, weil das der Parameter ist. Es gibt nur einen Bindestrich bei resume, weil das eine Kern-Nextflow-Option ist.

Das verwirrt die Leute die ganze Zeit, selbst wenn du Nextflow schon lange benutzt hast. Also denk immer daran, ein oder zwei Bindestriche. Hängt davon ab, ob es eine Kern-Nextflow-Option ist.

Okay, also jetzt mache ich _-resume_ und ich führe genau den gleichen Workflow wieder aus. Und diesmal sollte es ziemlich genau gleich aussehen mit einem wichtigen Unterschied.

In der Ausgabe hier kannst du sehen, dass die Ergebnisse gecacht wurden. Und tatsächlich ist dieser Task-Hash hier genau der gleiche wie beim vorherigen Lauf, und es hat einfach dieses Work-Verzeichnis in seiner Gesamtheit wiederverwendet. Die Eingaben und die Ausgaben und das Skript waren alle unverändert. Und so nimmt es einfach diese Datei davon und wenn es nachfolgende Schritte im Prozess gibt, würde es sie an den nächsten Schritt in der Pipeline übergeben.

Es führt also immer noch die gesamte Pipeline von Anfang bis Ende aus, aber es verwendet gecachte Ergebnisse für jeden dieser Tasks, wo es kann.

Wenn du jetzt _-resume_ machst, nimmt es einfach den letzten Pipeline-Lauf in deinem Arbeitsverzeichnis wieder auf, was auch immer das war. Aber du kannst tatsächlich von jedem vorherigen Lauf fortsetzen, den du dort gemacht hast. Und wir haben jetzt ziemlich viele gemacht.

## 4.2. Das Log vergangener Ausführungen inspizieren

Um alle anzusehen, können wir _"nextflow log"_ statt _"nextflow run"_ machen, und das gibt uns eine schöne Ausgabe, die all diese verschiedenen.. Ich muss meinen Bildschirm etwas kleiner machen, damit wir es sehen können, all diese verschiedenen Läufe zeigt, wann wir sie gemacht haben, die Session-ID, den Befehl und alles.

Und wir können hier reinschauen und wir können den Laufnamen von jedem davon nehmen und dann einen von diesen spezifischen fortsetzen. Ich kann also zurückgehen und den mit dem Namen _hungry_ekeblad_ fortsetzen. Und ich setze das einfach nach dem _resume_.

Wenn du übrigens neugierig bist, all diese Adjektive und Wissenschaftlernamen sind im Nextflow-Quellcode. Es ist ein wirklich guter Weg, deinen allerersten Pull-Request an Nextflow zu bekommen, indem du hingehst und es findest und deinen Lieblingswissenschaftler hinzufügst.

Und jedenfalls, ich habe das gemacht und es ging zurück und schaute sich die gecachten Ergebnisse von diesem Workflow-Lauf an, realisierte, dass es sie immer noch wiederverwenden konnte, und tat es. Ich bekam also wieder die gecachten Ergebnisse.

## 4.3. Ältere Work-Verzeichnisse löschen

Das ist großartig. Was ist, wenn ich diese Work-Verzeichnisse aufräumen möchte? Es gibt hier Unmengen davon. Es gibt Unmengen von Dateien. Vielleicht weiß ich sicher, dass ich von den letzten paar Pipeline-Läufen fortsetzen möchte, aber mir sind alle davor egal.

Dann kann ich einen hier auswählen und ich kann einen anderen Nextflow-Befehl verwenden, der _"nextflow clean"_ ist, und ich kann _"nextflow clean"_ machen, ich mache _"-before"_, und den bestimmten Laufnamen, der in diesem Fall _reverent_pike_ war, und ich mache _"-n"_, was Nextflow sagt, nur einen Trockenlauf zu machen. Es sagt mir also nur, was es löschen würde. Ohne tatsächlich etwas zu tun, also würde es diese Work-Verzeichnisse entfernen.

Das sieht vernünftig aus. Also mache ich den gleichen Befehl nochmal, aber statt _"-n"_ mache ich _"-f"_, um das Aufräumen tatsächlich zu machen. Und diesmal hat es tatsächlich all diese Verzeichnisse entfernt. Und wenn ich reingehe und mir die Work-Verzeichnisse anschaue, sieht es jetzt viel leichter aus. Fantastisch.

So räumst du also all deine lokalen Work-Verzeichnisse auf ziemlich sichere Weise auf, ohne den Cache komplett zu zerstören. Du kannst also immer noch fortsetzen, wenn du willst.

Wenn du jemals vergisst, was diese Flags für jeden Nextflow-Befehl sind, kannst du _"nextflow help"_ machen, und dann den Namen des Befehls. Also wenn ich _"nextflow help clean"_ mache, kannst du alle verschiedenen Optionen sehen: _-after, -before, -but_, alle verschiedenen Arten, dieses Aufräumverhalten zu konfigurieren. Ziemlich cool.

## Fazit

Okay, das ist das Ende von Teil 1 von Hello Nextflow. Es ist ein ziemlich intensiver Start in den Kurs, aber hoffentlich hast du jetzt ein ziemlich gutes Verständnis davon, wie ein Nextflow-Skript aussieht; mit verschiedenen Schlüsselteilen, den Prozessen, den Workflows, den Ausgaben und den Parametern. Du weißt, wie man sie mit grundlegenden Überschreibungen von der Befehlszeile konfiguriert, wie man einen dynamischen Input-Block mit einem dynamischen Skript macht, und du weißt, wie man all deine Workload-Ausführungen verwaltet: zu sehen, was du bereits ausgeführt hast, fortzusetzen, aufzuräumen. Es gibt viel Zeug. Du bist weit gekommen. Wenn du also eine Pause machen und kurz spazieren gehen und eine Tasse Tee trinken willst, ist jetzt wahrscheinlich ein guter Zeitpunkt. Du hast es dir verdient.

Von hier an bauen wir im Grunde auf dieser Grundlage auf. Wie können wir das komplexer, leistungsfähiger machen? Wie können wir es flexibler machen? Die Dinge tun, die wir für unsere Analyse im großen Maßstab tun wollen.

## Quiz

Wenn du jetzt zu Teil 1, Hello World, auf der Webseite scrollst, wirst du ein kleines Quiz sehen, und das ist etwas Neues, das wir für diese Version des Nextflow-Trainings gemacht haben. Und du kannst durchgehen und dich selbst quizzen, um zu überprüfen, dass du das ganze Material verstanden hast, das wir in diesem Kapitel gemacht haben.

Das wird nicht an uns gesendet oder so, es ist nur in deinem Browser gespeichert. Wir wissen also nicht, was deine Antworten sind, aber es ist nur eine kleine Selbstüberprüfung, um sicherzustellen, dass du nichts verpasst oder missverstanden hast. Und du kannst es so oft versuchen, wie du willst.

Wenn du wie ich bist, möchtest du vielleicht im Terminal in deiner VS Code-Instanz bleiben, in diesem Fall kannst du den _quiz_-Befehl eingeben und dann einfach sagen, in welchem Kapitel du bist. Also machen wir _"Hello World"_, und dann kannst du genau die gleichen Quizfragen machen, die im Webbrowser sind, aber einfach in deinem Terminal.

Cool. Okay. Hoffe, das gefällt dir. Hab ein bisschen Spaß und wir sehen uns im nächsten Kapitel in nur einer Minute, um alles über Nextflow-Kanäle zu sprechen.
