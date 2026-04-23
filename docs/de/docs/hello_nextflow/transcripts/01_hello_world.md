# Teil 1: Hello World - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../01_hello_world.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen

Hallo und willkommen zurück.

Du bist jetzt in Teil Eins des "Hello Nextflow"-Kurses mit dem Titel "Hello World". In diesem Kapitel werden wir beginnen, ein Verständnis für die absoluten Grundlagen von Nextflow aufzubauen.

Hoffentlich bist du jetzt in Codespaces oder einer vergleichbaren Umgebung mit VS Code eingerichtet, und du hast deinen Hello Nextflow-Ordner im Workspace im Explorer mit all diesen verschiedenen Dateien hier.

Wir werden damit beginnen, einige sehr grundlegende Dinge im Terminal mit Bash zu machen, und dann schauen wir, ob wir die gleichen Dinge in Nextflow machen können, damit du ein Gefühl für die Syntax bekommst.

## 0. Aufwärmen

Also fangen wir wirklich einfach an. Beginnen wir einfach mit "echo", um etwas ins Terminal auszugeben. "Hello World". Ich drücke Enter und das geht ins Terminal. Hello World. Hoffentlich ist das für niemanden, der diesen Kurs schaut, eine Überraschung.

Okay, lass uns etwas damit machen. Anstatt es nur ins Terminal auszugeben, schreiben wir es in eine Datei. Ich drücke die Pfeiltaste nach oben auf meiner Tastatur, die durch die Bash-Historie scrollt, sodass ich meinen letzten Befehl bekomme, und ich hänge am Ende ein kleines Größer-als-Symbol an, das die Ausgabe dieses Befehls in eine Datei umleitet, und ich nenne sie output.txt.

Nochmal Enter, um den Befehl auszuführen, diesmal nichts im Terminal, aber wir können auf der linken Seite sehen, dass die neue Datei hier aufgetaucht ist, genannt output.txt.

Wir können das in einem Terminal mit etwas wie cat anzeigen. Also cat output.txt und tatsächlich steht da "Hello World". Wir können auch darauf doppelklicken und es öffnet sich im Code-Editor in VS Code.

## 1.1. Den Code untersuchen

Gut. Ich sagte dir, es war einfach. Was kommt als Nächstes? Lass uns versuchen, diesen Prozess zu nehmen und ihn nochmal zu machen, aber diesmal machen wir es innerhalb von Nextflow.

Wie ich sagte, beginnen alle verschiedenen Kapitel in diesem Kurs mit einem Skript und dieses heißt Hello World. Also werde ich Hello World finden. Es zeigt eine Vorschau, wenn ich einmal darauf klicke, ich werde doppelklicken, um es im Editor hier zu öffnen. Und ich werde schnell das Terminal loswerden.

Nun, das ist ein sehr einfaches Skript, also so einfach wie es nur geht. Es ist nur 22 Zeilen lang und macht im Grunde das Gleiche. Tatsächlich sollte dir einiges davon bekannt vorkommen. Wir können unseren Bash-Befehl sehen, der in eine Datei umleitet.

Okay. Was noch? Auch in dieser Datei können wir einige der Kernkonzepte von Nextflow sehen. Wir haben hier einen Prozess in Rot und einen Workflow. Das sind die speziellen Schlüsselwörter und spezielle Terminologie in Nextflow.

## 1.1.1. Die Prozessdefinition

Verschiedene Prozesse innerhalb eines Workflows kapseln verschiedene logische Einheiten deines Workflows. Jeder Prozess macht eine Sache.

Wenn wir ihn ausführen, erzeugt er eine Aufgabe oder mehrere Aufgaben, die tatsächliche Ausführungsschritte einer Pipeline sind. Alle Prozesse werden dann innerhalb eines Workflow-Blocks orchestriert, den wir unten sehen, und in diesem Fall wird nur dieser eine Prozess ausgeführt.

Der Prozessname folgt diesem Schlüsselwort hier, und das kann im Grunde alles sein. Und dann befinden sich die Inhalte des Prozesses innerhalb dieser geschweiften Klammern.

Es gibt wirklich nur eine Anforderung für einen Prozess, nämlich dass er irgendeine Art von Script- oder Exec-Block enthält. Das ist hier in den dreifachen Anführungszeichen, und das ist das Bash-Skript, das in das Arbeitsverzeichnis geschrieben wird, wenn wir die Pipeline ausführen, und das ist das, was tatsächlich auf deinem Computer oder Server läuft.

Das ist typischerweise Bash, aber du kannst auch einen anderen Shebang hier oben einfügen, und es könnte ein Python-Skript oder ein R-Skript sein. Es spielt keine Rolle. Was auch immer in diesem Skript steht, wird ausgeführt.

Es gibt noch eine andere Sache, die wir in diesen Prozess hier eingefügt haben, nämlich die Output-Deklaration. Das sagt Nextflow, dass dieser Prozess eine Ausgabedatei namens output.txt erwartet. Es sagt, es ist ein path, also sollte es wie eine Datei behandelt werden, nicht etwa, wenn das val wäre, würde es sagen, es ist wie eine Variable oder ein Wert.

Beachte, dass dies diese Datei nicht erstellt. Es erzeugt sie nicht wirklich. Das wird vom Skript hier unten gemacht. Es sagt Nextflow nur, dass es eine Ausgabedatei mit diesem Dateinamen erwarten soll.

## 1.1.2. Die Workflow-Definition

Okay. Und dann haben wir unten einen Workflow hier, und wieder haben wir eine Deklaration. Diese heißt Main. Das ist das Workflow-Äquivalent eines Script-Blocks, wenn du so willst. Es ist der Teil des Workflows, der etwas tut. Und in diesem Fall sagen wir, rufe den Prozess namens sayHello auf.

Normalerweise wird deine Pipeline natürlich viel komplexer aussehen als das. Du wirst wahrscheinlich mehr als einen Prozess haben, und du wirst Kanäle verwenden, um den Datenfluss zwischen ihnen zu orchestrieren. Darauf kommen wir in den nächsten Teilen dieses Kurses, aber für jetzt reicht das. Das ist eine gültige Pipeline, die funktionieren sollte.

Ich kann sogar hier in VS Code auf "Preview DAG" klicken. Der DAG ist eine Darstellung einer Datenflussstruktur in der Pipeline, und wir können ihn auf der Seite als Mermaid-Diagramm gerendert sehen. In diesem Fall ist es sehr einfach. Es gibt eine Box, die der Workflow ist, und einen Prozess, der sayHello heißt, aber das könnte interessanter aussehen, wenn wir weitermachen.

## 1.2. Den Workflow ausführen

Okay, lass uns versuchen, diesen Workflow auszuführen und sehen, was passiert.

Ich werde das Terminal unten wieder hochholen, die Ausgabe löschen, und ich werde Nextflow Run eingeben. Und dann werde ich einfach den Skriptnamen eingeben, der hello-world.nf ist. Und ich werde Enter drücken.

Okay, es hat einige Standardsachen oben, die uns sagen, dass Nextflow gelaufen ist und welche Version lief und wie der Skriptname war und alles.

Und wirklich das Wichtige, wonach wir hier suchen, ist _hier_, was eine Zusammenfassung der verschiedenen Aufgaben ist, die ausgeführt wurden.

Wenn deins so aussieht mit einem kleinen grünen Häkchen, dann gut gemacht. Du hast gerade deine erste Pipeline ausgeführt. Fantastisch.

Es sagt uns hier den Namen des Prozesses, der lief, der sayHello hieß, und es sagte uns, dass er einmal lief und dass er erfolgreich war. Das wird aktualisiert, während du weitergehst, also wenn du eine größere Pipeline ausführst, wirst du den Fortschritt hier dargestellt sehen. Aber weil das so klein ist, läuft es im Grunde sofort.

## 1.2.2. Die Ausgabe und Logs im Work-Verzeichnis finden

Wenn du eine Nextflow-Pipeline ausführst, wird jeder dieser Prozesse zusammengefügt, und jeder Prozess kann, wie ich vorher sagte, Aufgaben erzeugen, eine oder mehrere. In diesem Fall hatten wir also eine einzelne Aufgabe von diesem Prozess. Sie lief nur einmal und das wurde unter diesem Aufgaben-_Hash_ gemacht.

Nextflow arbeitet nicht direkt mit den Dateien in deinem Arbeitsverzeichnis, es erstellt einen speziellen Ordner namens work. Und wenn ich "ls" mache, werden wir sehen, dass er hier erschienen ist: _work_, und darin sind Unterverzeichnisse für jede einzelne Aufgabe, die läuft. Und das passt zu diesem Hash. Du kannst also sehen, wenn ich zu "ls work/c4" gehe, und dann ist es abgeschnitten, aber es beginnt mit 203, und das ist das Arbeitsverzeichnis, das von diesem Prozess erstellt wurde, als wir die Pipeline ausgeführt haben. Und du kannst es auch auf der Seite sehen.

Wenn ich diese Dateien aufliste, kannst du sehen, dass die output.txt-Datei erzeugt wurde. Du kannst sie auch hier sehen. Und es gibt eine Reihe versteckter Dateien, die nicht mit meinem regulären "ls" angezeigt werden.

Wenn ich auf output.txt klicke, haben wir tatsächlich unsere Ausgabe. Fantastisch. Die Pipeline hat also funktioniert.

Es mag wie ziemlich viel Boilerplate für das Ausführen dessen erscheinen, was im Wesentlichen ein einzeiliges Bash-Skript war, aber es wird mehr Sinn ergeben, wenn unsere Prozesse komplizierter werden. Und dieses Work-Verzeichnis mit Nextflow und diese Dateien, die erstellt werden, sind wirklich das Rückgrat dessen, was Nextflow so mächtig macht.

Jede Aufgabe, jedes Element einer Pipeline ist von jeder anderen Aufgabe isoliert. Es ist reproduzierbar. Sie kollidieren nicht miteinander, und alles kann parallel laufen. Es ist tatsächlich eine wirklich schöne Art, wenn du dich daran gewöhnt hast, wegen dieser Isolation, dass du hineingehen und genau sehen kannst, was für eine einzelne Aufgabe passiert ist und debuggen kannst.

Lass uns einen kurzen Blick auf diese anderen Dateien im Work-Verzeichnis werfen. Von oben nach unten haben wir eine Datei namens _.command.begin_. Diese ist leer. Es ist nur eine sogenannte Sentinel-Datei, die von Nextflow erstellt wurde und sagt, okay, ich starte die Aufgabe. Nichts Interessantes dort.

Dann gibt es _.command.error_, _.command.log_ und _.command.out_. Das sind alles Ausgaben vom Bash-Befehl oder diesem Skript, das lief. Das ist Standard-Error. Das ist Standard-Out, und das sind die beiden zusammen, wie sie herauskamen. Du bekommst also die logische Reihenfolge.

Okay, die waren alle auch leer, also nicht sehr interessant, aber die Dinge werden interessanter, wenn du zu _.command.run_ kommst.

Das ist typischerweise ein sehr langes Skript. Und das ist es, was Nextflow tatsächlich ausführt. Wenn du hier hineingehst, wirst du anfangen, die gesamte innere Logik von Nextflow zu sehen und zu sehen, was es tut und wie es deinen Prozess ausführt. Das hängt davon ab, wo du läufst, ob wir lokal laufen oder es als Job an SLURM übermitteln, in welchem Fall wir SLURM-Header oben haben werden. All diese verschiedenen Setups.

Generell musst du aber nie wirklich in diese Datei schauen. Sie wird von Nextflow automatisch generiert und es gibt nichts wirklich Besonderes für deine Pipeline, was darin ist. Aber das ist wirklich der Kern dessen, was läuft.

Die nächste ist viel interessanter. _.command.sh_ ist das generierte Skript, das von deinem Prozess kam, und hier kannst du sehen, dass Nextflow den Bash-Header hinzugefügt hat, und dann hat es unseren Befehl ausgeführt, der in unserem Script-Block war.

Und das ist alles, was die _.command.run_-Datei macht, sie führt einfach diese _.command.sh_-Datei aus.

Das ist eine wirklich nützliche, die du normalerweise am meisten anschaust, wenn du versuchst, etwas zu debuggen und zu überprüfen, dass die Logik deiner Nextflow-Pipeline das tut, was du erwartest.

Schließlich haben wir eine Datei namens _.exitcode_, und diese erfasst nur den Exit-Code von einer Aufgabe, der in diesem Fall erfolgreich war. Der Exit-Code war also null.

Wenn etwas schiefgeht, dir der Speicher ausgeht oder etwas anderes und es fehlschlägt, dann ist das sehr nützlich, um zu verstehen, was schiefgelaufen ist.

## 1.3. Den Workflow erneut ausführen

Noch eine Sache, die man über Work-Verzeichnisse verstehen sollte: Wenn ich diese Pipeline wiederholt ausführe, also wenn ich _"nextflow run hello-world.nf"_ mache, wird es genau das Gleiche tun, aber diesmal wird es eine neue Task-ID haben. Du kannst sehen, dass dieser Hash hier anders ist, und wenn ich jetzt in work schaue, gibt es zwei Hash-Verzeichnisse. Und diese sind wieder voneinander getrennt.

Jedes Mal, wenn du also einen Nextflow-Workflow ausführst, es sei denn, du verwendest resume, das den Cache verwendet, darauf kommen wir später, wird es diese Prozesse in neuen Work-Verzeichnissen erneut ausführen, die voneinander getrennt sind. Du wirst keine Dateinamen-Kollisionen bekommen, du wirst keine Probleme wie das haben. Alles ist isoliert und sauber.

Und wenn wir in dieses Verzeichnis gehen, kannst du alle gleichen Dateien und die gleiche _output.txt_ sehen, die von Grund auf neu erstellt wurde.

## 2. Ausgaben veröffentlichen

Okay, das ist großartig für Nextflow selbst, während es deine Pipeline ausführt, sodass alle Dinge voneinander getrennt und sauber sind und verwaltet werden können.

Aber es ist nicht super nützlich, wenn du eine Person bist, die versucht, deine Ergebnisse zu erkunden. Du willst nicht wirklich durch Tausende und Abertausende verschiedener Work-Verzeichnisse graben, um deine Ergebnisdateien zu finden. Und das sollst du auch nicht. Die Work-Verzeichnisse sind nicht als finaler Zustand gedacht, wo deine Dateien erstellt werden.

Wir tun dies, indem wir unsere Dateien veröffentlichen.

## 2.1.1. Die Ausgabe des sayHello-Prozesses deklarieren

Wenn ich also zurück zu unserem Skript gehe, werden wir in unserem Workflow-Block hier arbeiten. Wir werden ihm sagen, welche Dateien zu erwarten sind, welche Dateien uns wichtig sind, und dann werden wir einen neuen Block darunter erstellen, genannt Output-Block.

Das ist die neue Syntax, die mit dem Syntax-Parser kam und standardmäßig in Version 26.04 von Nextflow ist. Wenn du also Nextflow schon ein bisschen benutzt hast, ist das eine der Sachen, die neu sind.

Also haben wir den Main-Block, und als Nächstes werde ich publish sagen und ich werde Nextflow sagen, was vom Publishing zu erwarten ist. Wir werden es _first_output_ nennen, und wir werden es _sayHello.out_ nennen.

Ich habe dort versehentlich einen Tippfehler gemacht, aber das ist eine gute Gelegenheit, auch auf einige der Funktionen der Nextflow VS Code-Erweiterung hinzuweisen. Du kannst sehen, dass es mir sofort eine kleine rote Wellenlinie darunter gegeben hat, die sagt, etwas ist falsch. Und wenn ich darüber fahre, wird es mir sagen, diese Variable ist nicht definiert. Ich weiß nicht, was es ist.

Es ist in diesem Fall ziemlich offensichtlich, ich habe einen Tippfehler gemacht. Ich wollte sayHello eingeben, und dann verschwindet die Wellenlinie.

Jetzt ist es lila. Der Nextflow-Syntax-Parser weiß, dass dies ein Prozess ist, und wenn ich darüber fahre, gibt er mir eine reduzierte Darstellung davon, wie dieser Prozess aussieht. Ich kann also sehr schnell auf einen Blick sehen, dass er keine Eingaben nimmt und uns diese Ausgabe gibt. Das Arbeiten in VS Code mit dieser Erweiterung gibt dir also viele kontextuelle Informationen, während du Code schreibst.

Beachte, dass wir auf die Ausgabe von diesem Prozess mit der _.out_-Syntax verweisen können. Und im Moment können wir das nennen, wie wir wollen, es ist nur ein beliebiger Variablenname.

## 2.1.2. Einen output:-Block zum Skript hinzufügen

Wo es wichtig wird, ist, wenn wir unseren neuen Block hier machen, und das ist jetzt unterhalb des Workflow-Blocks, wir sind nicht mehr innerhalb von workflow. Geschweifte Klammern wieder. Und hier sagen wir Nextflow einfach, wo alle Dateien hingelegt werden sollen, die vom Workflow erstellt werden.

Jetzt werde ich diesen Variablennamen nehmen, den ich hier erstellt habe, und ich werde ihn dort einfügen und einige geschweifte Klammern dafür setzen. Und ich werde Nextflow sagen, einen path zu verwenden. Hoppla. Path, in Anführungszeichen. Und ich werde einen Punkt verwenden. Das sagt Nextflow einfach, die Datei in das Root des Results-Verzeichnisses zu legen. Also keine Unterverzeichnisse oder so.

Lass uns versuchen, unseren Workflow erneut auszuführen. Wenn ich _"nextflow run hello-world.nf"_ mache, dann sollte es hoffentlich im Grunde genau gleich aussehen. Nichts hat sich wirklich mit Nextflow hier geändert. Es führt die gleichen Dinge aus. Es macht sie nur wieder in Work-Verzeichnissen.

Aber jetzt, wenn ich _"ls results/"_ mache, wirst du sehen, dass hier ein neues Verzeichnis erstellt wurde, genannt results, was das Standard-Basisverzeichnis für Workflow-Publishing ist. Und darin ist eine Datei namens _output.txt_.

Wenn ich _"ls -l results"_ mache, wirst du sehen, dass dies tatsächlich ein Soft-Link zum Work-Verzeichnis ist. Das ist also keine echte Datei, sie ist mit dem Work-Verzeichnis verlinkt und hat alle Dateien dort für uns gesammelt.

## 2.2. Einen benutzerdefinierten Speicherort festlegen

"Results" ist der Standardname für diesen Pfad. Wenn ich den Workflow erneut ausführe, und diesmal mache ich _dash_ einzelner Bindestrich, das ist, weil es eine Kern-Nextflow-Option ist. _"-output-dir **my**results"_. Könnte auch einfach _"-o"_ kurz machen. Dann wird es ein anderes Basisverzeichnis festlegen, wo die Dateien gespeichert werden, und wieder, hier oben in _myresults/_, haben wir jetzt eine _output.txt_.

Das ist großartig, aber wir wollen wahrscheinlich nicht alle Dateien nur im Root. Wir wollen etwas Organisation, also können wir auch ein Unterverzeichnis hier erstellen, genannt wie wir wollen. Sagen wir _"path 'hello_world'"_, und ich führe das einfach nochmal aus. _"nextflow run hello-world.nf"_. Es sollte ins Results-Verzeichnis in ein Unterverzeichnis gehen und tatsächlich, jetzt unter results hier oben haben wir _hello_world/_ und wir haben _output.txt_.

Wichtig zu bemerken, die alte _output.txt_-Datei ist immer noch da. Das Results-Verzeichnis wird nicht gelöscht, wenn du das machst. Nur neue Dateien werden dort hineinkopiert. Sie überschreiben Dateien, die bereits dort sind, wenn sie den gleichen Dateinamen haben, aber sie werden alte nicht löschen. Du musst also ein bisschen vorsichtig sein, wenn du Pipelines erneut ausführst. Wenn du nicht willst, dass sie auf den Dateien sind, die bereits dort sind, stelle sicher, dass du ein leeres Verzeichnis verwendest.

## 2.3. Den Publish-Modus auf copy setzen

Okay, ich erwähnte, dass diese Dateien Soft-Links sind, also wenn ich _"ls -l results/hello_world/"_ mache, kannst du sehen, dass es zum Work-Verzeichnis soft-linkt. Das ist generell eine gute Sache, wenn du an etwas wie HPC arbeitest, und das sind wirklich riesige Dateien und du willst sie nicht duplizieren, weil es bedeutet, dass die Dateien nur einmal im Dateisystem gespeichert werden.

Es bedeutet jedoch, dass wenn du das Work-Verzeichnis löschst: wenn ich _"rm -r work"_ mache und all diese Zwischendateien lösche, die erstellt wurden. Jetzt, wenn ich versuche, diese Datei zu lesen _"results/hello_world/"_. Sie wird als Soft-Link auf eine Datei zeigen, die nicht mehr existiert, und die Daten sind für immer weg und nicht wiederherstellbar, was vielleicht nicht großartig ist.

Also generell sage ich, es ist gute Praxis, die Dateien zu kopieren, anstatt sie zu soft-linken, wenn du kannst, weil es sicherer ist. Sei dir nur bewusst, dass es doppelt so viel Festplattenspeicher verwenden wird, es sei denn, du löschst diese Work-Verzeichnisse.

Um das mit dem Output-Block zu machen, werde ich zum first_output hier gehen. Ich habe vorher den Pfad gesetzt und jetzt werde ich den Modus setzen und du kannst sehen, während ich tippe, schlägt die VS Code-Erweiterung Sachen vor, sie weiß, dass es hier eine Output-Direktive ist. Und ich werde copy sagen. Ich drücke Speichern.

Lass uns den Workflow erneut ausführen. Er wird die Dateien wieder erstellen, neues Work-Verzeichnis.

Jetzt, wenn ich zu _"ls -l results/hello_world/"_ gehe, kannst du sehen, dass dies eine echte Datei ist und kein Soft-Link mehr, und Nextflow hat das kopiert. Gut zu wissen. Also path und mode sind Dinge, die du ziemlich oft schreiben wirst.

Nun, natürlich ist das sehr einfach. Wir werden das komplexer und mächtiger machen, während wir weitermachen, und du wirst sehen, wie man diese Dinge dynamisch macht und nicht zu ausführlich.

## 2.4. Hinweis zu publishDir-Direktiven auf Prozessebene

Nun, ich sagte, als wir damit anfingen, dass dies eine ziemlich neue Form der Syntax ist. Sie ist nur in den neuesten Versionen von Nextflow verfügbar, während ich das aufnehme, und sie heißt Workflow Outputs.

Wenn du das verwendest, ist es großartig. Es schaltet viele andere coole Funktionen innerhalb von Nextflow frei, wie zum Beispiel Nextflow Lineage, um die Herkunft dieser Dateien zu verfolgen, während sie erstellt werden, und wird bald der Standard in 26.04 sein. Und zu einem späteren Zeitpunkt in der Zukunft wird dies die einzige Art sein, deine Workflows zu schreiben.

Da wir uns jedoch gerade in dieser Übergangsphase befinden, könntest du durchaus Pipelines in freier Wildbahn sehen, die etwas namens publishDir verwenden, was der alte Weg ist, es zu tun, und das wird nicht auf Workflow- und Output-Ebene definiert, sondern auf Prozessebene.

Und diese Deklaration sagt im Grunde das Gleiche. Sie sagt, veröffentliche die Ergebnisdateien in ein Verzeichnis namens results und verwende einen Copy-Modus. Du kannst also sehen, dass die Syntax sehr ähnlich ist. Aber wenn du jetzt neue Pipelines schreibst, versuche diese publishDir-Direktive nicht zu verwenden, selbst wenn du sie in KI-Ergebnissen oder in Dokumentation oder anderen Pipelines siehst, weil das der alte Weg ist, es zu tun.

Im Jahr 2026 sollten wir alle Workflow Outputs verwenden.

Das ist alles dokumentiert, wenn du das machst und Nextflow schon vorher benutzt hast, kannst du zu den Nextflow-Docs hier gehen, nextflow.io/docs/. Und wenn ich zu Tutorials runterscrolle, gibt es ein Tutorial namens _Migrating to Workflow Outputs_.

Es ist wirklich gut. Es geht durch die gesamte Syntax, wie sie zur alten Syntax äquivalent ist, warum wir sie geändert haben, und hat eine Zeitleiste und alles. Und es geht durch all die verschiedenen Szenarien mit vielen, vielen Beispielen. Du kannst also bestehenden Nextflow-Code leicht zur neuen Syntax konvertieren.

## 3.1. Den sayHello-Prozess ändern, um eine variable Eingabe zu erwarten

Okay, wir haben also unser einfaches Skript, das einen Prozess ausführt, eine Datei erstellt, Nextflow sagt, dass es eine Ausgabe ist, und dann sagen wir Nextflow, wo diese Datei gespeichert werden soll. Das ist ein guter Anfang.

Aber es wäre interessanter, wenn nicht alles fest codiert wäre. Als Nächstes überlegen wir also, wie wir Nextflow sagen können, dass dieser Prozess eine variable Eingabe annehmen kann, die etwas ist, das wir zur Laufzeit steuern können, wenn wir einen Workflow starten.

Wir müssen ein paar verschiedene Dinge tun, damit das funktioniert.

Erstens müssen wir diesem Prozess sagen, dass er eine Eingabevariable akzeptieren kann, und wir geben hier _input_ als neuen Deklarationsblock ein. Und wir werden das _"val greeting"_ nennen.

Das val-Bit ist das Äquivalent zu einem path hier unten. Es sagt Nextflow, dass dies eine Variable ist, wie ein String in diesem Fall. Und wenn du wieder darüber fährst, sagt es dir von der Erweiterung, was das bedeutet.

Als Nächstes werden wir Nextflow sagen, was damit zu tun ist. Es reicht nicht, nur zu sagen, dass es eine Variable gibt. Du musst im Skript sagen, wie diese Variable verwendet werden soll. Und ich werde diesen fest codierten String hier loswerden, und ich werde eine Variable einfügen.

Ich werde es schnell ohne geschweifte Klammern machen, nur um dir zu zeigen, dass das erlaubt ist, und das ist die alte Art, es zu tun. Aber jetzt mit der neuen Syntax empfehlen wir wirklich, es in geschweifte Klammern wie diese zu setzen, und es macht wirklich klar, dass dies von Nextflow hier interpoliert wird.

Großartig. Also _"input greeting"_ geht in _$\{greeting\}._ Das Letzte ist, wir müssen Nextflow auf Workflow-Ebene sagen, dass dieser Prozess jetzt eine Eingabe nimmt. Und um das zu tun, werden wir ihm im Grunde eine Variable geben.

## 3.2. Einen Kommandozeilen-Parameter einrichten, um Benutzereingaben zu erfassen

Wir könnten es wieder fest codieren, wie Hello World, und das würde gut funktionieren, aber offensichtlich gibt es uns keinen wirklichen Vorteil. Wir wollten in der Lage sein, das zur Laufzeit zu konfigurieren, also wollen wir in der Lage sein, es auf der CLI zu tun, wenn du Nextflow startest.

Und die Art, wie wir das tun, ist ein spezielles Nextflow-Konzept namens _params_. Wir werden das _params.input_ nennen.

Was das macht, ist, dass es diese Input-Variable auf der CLI freilegt, und da verwenden wir einen doppelten Bindestrich, wenn wir Nextflow starten.

Ich kann das nennen, wie ich will, ich kann es _hello, greeting_ nennen. Spielt keine Rolle. Was auch immer ich dort mache, wird als CLI-Option freigelegt, wenn wir eine Pipeline starten. Und das ist ein echter Zaubertrick von Nextflow, denn es bedeutet, dass du dein Workflow-Skript sehr schnell mit diesen Parametern aufbauen kannst, und du baust im Wesentlichen eine benutzerdefinierte CLI für deine Pipeline auf, was es wirklich einfach macht, verschiedene Optionen spontan anzupassen, wenn du startest.

Also. Lass es uns ausprobieren. Zurück zu unserem Terminal. Wir haben unseren _"nextflow run"_-Befehl hier. Und jetzt werde ich _"--input"_ machen, was dem _"params.input"_ entspricht, das wir vorher gesehen haben. Ich glaube, in den Docs ist es auf Französisch. Geraldine spricht gerne Französisch. Ich werde es auf Schwedisch machen, weil ich in Schweden lebe. Also werde ich sagen, "_Hej Världen_" und Enter drücken.

Kann einfache oder doppelte Anführungszeichen verwenden, es beeinflusst nur, wie Bash es interpretiert.

Es führt die Nextflow-Pipeline genau auf die gleiche Weise aus. Du kannst sehen, dass das Arbeitsverzeichnis und alles gleich ist. Aber jetzt, wenn ich zu _"results/hello_world/output"_ gehe. Können wir unser schönes Schwedisch hier stattdessen sehen.

Wir haben also dynamisch eine Eingabe von einer CLI an einen Parameter übergeben. Wir haben das als Eingabe an den Prozess übergeben und der Prozess hat das interpretiert und in einen Script-Block eingefügt, der dann dynamisch die Ausgabe dieses Skriptergebnisses geändert hat. Ziemlich cool.

Ziemlich komplexe Logik mit sehr wenig Syntax hier. Und du kannst hoffentlich sehen, wie das jetzt skaliert. Und so bauen wir wirklich die Logik und die Anpassbarkeit unserer Pipelines in das Nextflow-Skript ein.

## 3.4. Standardwerte für Kommandozeilen-Parameter verwenden

Okay, das ist großartig. Das Problem ist jetzt aber, dass ich jedes Mal, wenn ich diese Pipeline ausführe, dash input machen muss, damit sie läuft.

Wenn ich versuche, ohne diesen Parameter zu laufen, wird Nextflow jetzt einen Fehler werfen und sagen, dass es diesen Parameter brauchte und er nicht gesetzt wurde. Und es wusste also nicht, was zu tun ist.

Das ist übrigens eine coole neue Sache. In der Vergangenheit wäre Nextflow einfach mit einem leeren String gelaufen, und du hättest alle möglichen seltsamen Fehler gehabt, die schwer zu verstehen gewesen wären. Aber im neuen Nextflow-Syntax-Parser ist es etwas vorsichtiger und sagt es dir sofort.

Wir wollen also nicht immer jede einzelne Option angeben. Es ist gute Praxis, vernünftige Standards anzugeben. Wie machen wir das also in unserem Skript?

Du wirst bemerken, dass wir, als wir das geschrieben haben, einfach _params.input_ direkt dort eingefügt haben, wo wir es verwenden. Die offensichtliche Lösung ist also, dass wir einen Standard definieren, und wir tun das oben im Skript hier in einem speziellen Params-Block im Workflow. Das ist im Workflow-Skript hier.

Wieder etwas neue Syntax hier, also pass auf. Das ist wirklich cooles Zeug. Wir haben den Namen des Parameters, der hier erwartet wird.

Und dann nach diesem Doppelpunkt-Zeichen definieren wir einen Typ der Variable. Du musst das nicht tun, du kannst es einfach leer lassen, aber es ist wirklich schön. Es sagt Nextflow, dass wir einen String erwarten und ihn als solchen behandeln.

Wenn wir stattdessen zum Beispiel eine Zahl wollen, könnten wir float schreiben, und das würde sagen, wir wollen eine Gleitkommazahl. Und wenn wir versuchen, damit zu laufen, wird es einen Fehler werfen. Wenn wir ihm einen String geben, der kein Float ist. Und es wird ihn auch als solchen übergeben. Wenn wir string machen, dann weiß es, dass es ein String ist. Und selbst wenn es führende Nullen hat und ganz numerisch ist, wird es ihn trotzdem als tatsächlichen String übergeben.

Diese Typsicherheit ist also eine sehr neue Funktion von Nextflow, aber wirklich mächtig, um deinen Code sicherer zu schreiben und auszuführen.

Dann danach haben wir ein Gleichheitszeichen und dann den Standardwert hier. Nextflow wurde ursprünglich in Barcelona geschrieben, also scheint es angemessen, dass wir etwas Spanisches hier haben, _"Holà mundo!"_ als Standard.

Richtig, ich werde das Skript speichern, zurückgehen, das Skript wieder ohne _--input_ ausführen. Und diesmal sollte es laufen und es wird unsere neue Datei oben in _results_ erstellen. Und in dieser Datei steht jetzt _"Holà mundo!"_.

Das ist aber nur ein Standard, also bedeutet es nicht, dass wir nicht immer noch das Gleiche wie vorher tun können. Wenn ich zurückgehe und mein altes Skript hier finde, _"Hej Världen"_, weil ich _--input_ auf der Kommandozeile mache, wird das diesen Standard überschreiben und das wieder in der output.txt-Datei verwenden.

Das im Skript ist also nur der Standardwert, den ich setze.

Während wir unseren Workflow aufbauen, um komplexer zu sein und mehr Parameter einzuschließen, wird dieser Params-Block oben im Skript anfangen, sie alle an einem Ort zu sammeln.

Und du endest mit dieser ziemlich schönen Symmetrie in deinem Skript, wo du effektiv alle deine Workflow-Eingaben hier und deine Workflow-Ausgaben unten hast. Und es ist sehr klar, was die Schnittstelle deines Workflows zur Außenwelt ist. Du kannst also eine neue Pipeline sehr schnell mit der neuen Syntax aufnehmen und verstehen, wie man sie verwendet.

Eine letzte coole Sache. Wir müssen keinen Standardwert damit setzen. Wenn wir params input machen, aber keinen Standardwert setzen, dann sagt es Nextflow, dass dieser Parameter erforderlich ist, und wieder wird die Pipeline nicht ohne ihn laufen, aber sie wird dir eine nützlichere Fehlermeldung geben, anstatt etwas darüber, dass es null ist.

Es sagt also, wir erwarten, dass seine Eingabe erforderlich ist, aber sie wurde nicht auf der Kommandozeile angegeben. Sehr schön.

Okay, hoffentlich ist jetzt klar, wie man seine Nextflow-Pipeline mit variablen Eingaben und Parametern einrichtet, wie man den Standard setzt, die Typen setzt, es könnte ein Boolean true-false-Flag oder eine Ganzzahl oder verschiedene Typen hier sein. Wie man sie in deinen Workflow übergibt, wo es durchgeht, und dann in deinen Prozess interpoliert. Und du weißt auch, wie man diese auf der Kommandozeile anpasst, wenn du Nextflow startest. Das fängt an, interessanter auszusehen als unser einfacher Bash-Befehl.

## 4. Workflow-Ausführungen verwalten

Okay. Was kommt als Nächstes? Für den letzten Teil dieses Kapitels werden wir ein bisschen darüber sprechen, wie man all die verschiedenen Workflow-Ausführungen verwaltet. Wenn du in meiner Seitenleiste hier und im Explorer unter work schaust, wirst du sehen, dass ich eine Reihe verschiedener Pipelines ausgeführt habe und diese Work-Verzeichnisse ziemlich lang werden, es gibt viele davon.

Und die andere Sache ist, wie ich vorher sagte, jedes Mal, wenn ich diese Pipeline erneut ausführe, erstellt sie einen neuen Satz von Work-Verzeichnissen, und sie führt alle Prozesse von Grund auf neu aus, was eine gute Sache ist. Das ist beabsichtigtes Verhalten. Es ist reproduzierbar und regeneriert alles frisch. Aber es ist offensichtlich, wenn du sehr lang laufende Prozesse ausführst, ärgerlich, deine Pipeline immer von Anfang an starten zu müssen, wenn sie auf halbem Weg abgestürzt ist, oder wenn du etwas am Ende der Pipeline änderst.

## 4.1. Einen Workflow mit -resume neu starten

Glücklicherweise ist Nextflow wirklich gut darin zu wissen, was zuvor ausgeführt wurde und was verfügbar ist, und diese alten Ergebnisse wiederzuverwenden ist sehr einfach. Wir fügen einfach ein neues Flag am Ende des Befehls hinzu _"-resume"_.

Beachte nun, dass es zwei Bindestriche bei input gibt, weil das der Parameter ist. Es gibt nur einen Bindestrich bei resume, weil das eine Kern-Nextflow-Option ist.

Es bringt Leute die ganze Zeit durcheinander, selbst wenn du Nextflow schon lange benutzt hast. Also denk immer daran, ein oder zwei Bindestriche. Hängt davon ab, ob es eine Kern-Nextflow-Option ist.

Okay, also jetzt mache ich _-resume_ und ich führe genau den gleichen Workflow wieder aus. Und diesmal sollte es ziemlich genau gleich aussehen mit einem wichtigen Unterschied.

In der Ausgabe hier kannst du sehen, dass die Ergebnisse gecacht wurden. Und tatsächlich ist dieser Task-Hash hier genau der gleiche wie beim vorherigen Lauf, und er hat einfach dieses Work-Verzeichnis in seiner Gesamtheit wiederverwendet. Die Eingaben und die Ausgaben und das Skript waren alle unverändert. Und so nimmt es einfach diese Datei davon und wenn es nachgelagerte Schritte im Prozess gibt, würde es sie an den nächsten Schritt in der Pipeline weitergeben.

Es führt also immer noch die gesamte Pipeline von Anfang bis Ende aus, aber es verwendet gecachte Ergebnisse für jeden dieser Tasks, wo es kann.

Wenn du jetzt _-resume_ machst, nimmt es einfach den letzten Pipeline-Lauf in deinem Arbeitsverzeichnis wieder auf, was auch immer das war. Aber du kannst tatsächlich von jedem vorherigen Lauf wieder aufnehmen, den du dort gemacht hast. Und wir haben jetzt ziemlich viele gemacht.

## 4.2. Das Log vergangener Ausführungen inspizieren

Um sie alle anzuschauen, können wir _"nextflow log"_ statt _"nextflow run"_ machen, und das wird uns eine schöne Ausgabe geben, die all diese verschiedenen... Ich muss meinen Bildschirm ein bisschen kleiner machen, damit wir es sehen können, all diese verschiedenen Läufe zeigt, wann wir sie gemacht haben, die Session-ID, den Befehl und alles.

Und wir können hier reinschauen und wir können den Run-Namen von jedem davon nehmen und dann einen dieser spezifischen wieder aufnehmen. Ich kann also zurückgehen und ich kann den namens _hungry_ekeblad_ wieder aufnehmen. Und ich setze das einfach nach dem resume.

Wenn du übrigens neugierig bist, all diese Adjektive und Wissenschaftlernamen sind im Nextflow-Quellcode. Es ist ein wirklich guter Weg, deinen allerersten Pull Request zu Nextflow zu bekommen, indem du hingehst und ihn findest und deinen Lieblingswissenschaftler hinzufügst.

Und jedenfalls, ich habe das gemacht und es ging zurück und schaute sich die gecachten Ergebnisse von diesem Workflow-Lauf an, erkannte, dass es sie immer noch wiederverwenden konnte, und tat es. Ich bekam also wieder die gecachten Ergebnisse.

## 4.3. Ältere Work-Verzeichnisse löschen

Das ist großartig. Was ist, wenn ich diese Work-Verzeichnisse aufräumen will? Es gibt viele davon hier. Es gibt viele Dateien. Vielleicht weiß ich mit Sicherheit, dass ich von den letzten paar Pipeline-Läufen wieder aufnehmen will, aber mir sind alle davor egal.

Dann kann ich hier einen auswählen und ich kann einen anderen Nextflow-Befehl verwenden, der _"nextflow clean"_ ist, und ich kann _"nextflow clean"_ machen, ich werde _"-before"_ machen, und den bestimmten Run-Namen, der in diesem Fall _reverent_pike_ war, und ich werde _"-n"_ machen, was Nextflow sagt, nur einen Dry-Run zu machen. Es sagt mir also nur, was es löschen würde. Ohne tatsächlich etwas zu tun, würde es also diese Work-Verzeichnisse entfernen.

Das sieht vernünftig aus. Ich werde also den gleichen Befehl nochmal machen, aber statt _"-n"_ werde ich _"-f"_ machen, um das Cleanup tatsächlich durchzuführen. Und diesmal hat es tatsächlich all diese Verzeichnisse entfernt. Und wenn ich reingehe und mir die Work-Verzeichnisse anschaue, sieht es jetzt viel leichter aus. Fantastisch.

So räumst du also all deine lokalen Work-Verzeichnisse auf ziemlich sichere Weise auf, ohne den Cache komplett zu zerstören. Du kannst also immer noch wieder aufnehmen, wenn du willst.

Falls du jemals vergisst, was diese Flags für jeden Nextflow-Befehl sind, kannst du _"nextflow help"_ machen, und dann den Namen des Befehls. Wenn ich also _"nextflow help clean"_ mache, kannst du all die verschiedenen Optionen sehen: _-after, -before, -but_, all die verschiedenen Wege, dieses Cleanup-Verhalten zu konfigurieren. Ziemlich cool.

## Fazit

Okay, das ist das Ende von Teil Eins von Hello Nextflow. Es ist ein ziemlich intensiver Start in den Kurs, aber hoffentlich hast du jetzt ein ziemlich gutes Verständnis davon, wie ein Nextflow-Skript aussieht; mit verschiedenen Schlüsselteilen, den Prozessen, den Workflows, den Ausgaben und den Parametern. Du weißt, wie man sie mit grundlegenden Überschreibungen von der Kommandozeile konfiguriert, wie man einen dynamischen Input-Block mit einem dynamischen Skript macht, und du weißt, wie man all deine Workload-Ausführungen verwaltet: sehen, was du bereits ausgeführt hast, wieder aufnehmen, aufräumen. Das ist eine Menge Zeug. Du bist weit gekommen. Wenn du also eine Pause machen und eine kleine Runde gehen und eine Tasse Tee trinken willst, ist jetzt wahrscheinlich ein guter Zeitpunkt. Du hast es dir verdient.

Von hier an bauen wir im Grunde auf dieser Grundlage auf. Wie können wir das komplexer, mächtiger machen? Wie können wir es flexibler machen? Die Dinge tun, die wir für unsere Analyse im großen Maßstab tun wollen.

## Quiz

Wenn du jetzt zu Teil Eins, Hello World, auf der Webseite runterscrollst, wirst du ein kleines Quiz sehen, und das ist etwas Neues, das wir für diese Version des Nextflow-Trainings gemacht haben. Und du kannst durchgehen und dich selbst abfragen, um zu überprüfen, dass du das gesamte Material verstanden hast, das wir in diesem Kapitel gemacht haben.

Das wird nicht an uns gesendet oder so, es wird nur in deinem Browser gespeichert. Wir wissen also nicht, was deine Antworten sind, aber es ist nur eine kleine Selbstüberprüfung, um sicherzustellen, dass du nichts verpasst oder missverstanden hast. Und du kannst es so oft versuchen, wie du willst.

Wenn du wie ich bist, möchtest du vielleicht im Terminal in deiner VS Code-Instanz bleiben, in diesem Fall kannst du den _quiz_-Befehl eingeben und dann einfach sagen, in welchem Kapitel du bist. Wir machen also _"Hello World"_, und dann kannst du genau die gleichen Quizfragen machen, die im Webbrowser sind, aber einfach in deinem Terminal.

Cool. Okay. Hoffe, du genießt das. Hab ein bisschen Spaß und wir sehen uns im nächsten Kapitel in nur einer Minute, um alles über Nextflow-Kanäle zu sprechen.
