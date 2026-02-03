# Teil 6: Hello Config - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../06_hello_config.md) zurück.

    Die im Transkript angezeigten Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern im Material.

## Willkommen

Hallo, willkommen zu Teil sechs des Hello Nextflow Trainingskurses.

Dieses Kapitel heißt Hello Config und ist der letzte Teil unseres Trainingskurses.

In diesem Kapitel werden wir über Nextflow-Konfiguration sprechen. Die Nextflow-Konfiguration ist wirklich leistungsstark. Sie ermöglicht es uns, dieselbe Pipeline auf mehreren verschiedenen Computing-Infrastrukturen mit unterschiedlicher Software-Bereitstellung und verschiedenen Optionen in der Pipeline selbst auszuführen.

Das bedeutet, dass du Nextflow-Pipelines, die von anderen Leuten erstellt wurden, auf deinem System ausführen kannst, obwohl sie möglicherweise für eine völlig andere Infrastruktur erstellt wurden. Diese Fähigkeit, Nextflow zu konfigurieren, macht Workflows wirklich portabel und gemeinsam nutzbar.

In diesem Kapitel werden wir den Workflow verwenden, den wir in früheren Teilen erstellt haben, aber wir werden den Workflow-Code überhaupt nicht bearbeiten. Wir werden uns nur unsere Nextflow-Konfigurationsdatei ansehen und sehen, wie das Ändern der Konfiguration die Art und Weise verändert, wie Nextflow läuft.

Okay, lass uns anfangen.

Genau wie zuvor, lass uns damit beginnen, zu training.nextflow.io zu gehen. Gehe links zu Hello Nextflow und Kapitel sechs, Hello config. Ich werde jetzt in meine GitHub Codespaces-Umgebung gehen und das Skript überprüfen, das wir verwenden werden.

## 0. Aufwärmen: Überprüfe, dass Docker aktiviert ist und führe den Hello Config Workflow aus

Dieses hier heißt Hello Config und beginnt dort, wo wir vorher waren. Also sieht es genau gleich aus mit unseren drei Parametern: Greetings für die CSV-Datei, batch für den Ausgabe-Batch-Namen und character für den Cowpy-Namen. Wir haben unsere vier Imports der verschiedenen Prozesse, und dann haben wir einen Workflow, in dem wir sie miteinander verknüpfen.

Ich werde diese Datei jetzt tatsächlich schließen, weil wir die Nextflow-Datei in diesem Kapitel überhaupt nicht anfassen werden. Wir werden rein innerhalb der Konfigurationsdatei arbeiten. Wenn ich in die Nextflow.config-Datei schaue, die wir im vorherigen Kapitel fünf kurz betrachtet haben, können wir sehen, dass wir hier eine einzelne Anweisung haben: Docker enabled equals true, die Nextflow mitteilt, Docker zu verwenden, wenn dieser Workflow ausgeführt wird.

Ich verwende Nextflow.config im Pipeline-Stammverzeichnis hier, die automatisch geladen wird, wenn ich Nextflow ausführe. Aber denk daran, Nextflow kann Konfigurationsdateien von mehreren Orten laden.

Wenn ich mit Nextflow docs zu Configuration gehe, kannst du eine Liste dieser Orte und eine Priorität sehen, in der sie geladen werden.

Okay. Lass uns überprüfen, dass unser Workflow wie erwartet ausgeführt wird. Öffne ein Terminal. Mache Nextflow run hello-config und drücke Enter. Wir sollten diese vier Prozesse laufen haben, die mit einem Cowpy-Befehl enden. Tatsächlich hat das richtig funktioniert. Ich hatte Docker aktiviert, es hat Docker heruntergeladen und Cowpy für mich ausgeführt, genau wie am Ende von Kapitel fünf.

## 1. Bestimme, welche Software-Packaging-Technologie verwendet werden soll

Okay. Sagen wir, ich arbeite auf einem HPC und ich habe Docker nicht installiert. Das Beste in diesem Szenario wäre, Singularity oder Apptainer zu verwenden. Wenn ich das tun würde, würde ich in das Modul cowpy gehen und diesen Container ändern, um das Singularity-Image zu verwenden, wie ich im vorherigen Kapitel gezeigt habe, mit einem oras://, das du auch von Seqera Containers bekommen kannst.

Ich würde dann zu Nextflow.config gehen, Docker enabled auf false setzen und singularity enabled equals true machen. Oder, wenn Apptainer verwendet wird, apptainer enabled equals true und das würde funktionieren.

Nextflow unterstützt auch andere Technologien neben Containern, etwas mit dem du vielleicht vertraut bist, ist conda. Hier können wir conda enabled equals true machen und Docker auf false setzen. conda verwendet nicht dieselbe container-Direktive. Stattdessen können wir hier eine neue hinzufügen, die conda heißt. Wir spezifizieren dann das conda-Paket, das wir verwenden möchten. Es ist gute Praxis, so spezifisch wie möglich zu sein, um zu versuchen, die Pipeline so reproduzierbar wie möglich zu machen. Also werde ich den conda-Channel spezifizieren, conda-forge, und dann cowpy, und die exakte Version, die 1.1.5 war.

Ich könnte auch einfach cowpy schreiben, wenn ich wollte, aber das könnte bei verschiedenen Ausführungen der Pipeline zu einer anderen Version von cowpy führen.

Das Schöne daran ist, dass ich die Docker-Direktive überhaupt nicht angefasst habe. Dieses Docker-Image ist immer noch da. Ich stelle jetzt nur zwei Alternativen bereit, und diese können ein- oder ausgeschaltet werden, indem man allein eine Konfigurationsdatei verwendet.

## 1.3. Führe den Workflow aus, um zu überprüfen, dass er Conda verwenden kann

Conda ist jetzt aktiviert, also lass es uns ausprobieren.

Großartig. Es läuft und du kannst sehen, dass hier eine Nachricht von Nextflow ist, die besagt, dass Nextflow eine Conda-Umgebung für mich erstellt, und es verwendet diesen Cache-Speicherort.

Im Hintergrund führt Nextflow "conda create"-Befehle für mich aus, um eine neue isolierte Conda-Umgebung mit nur den Paketen zu erstellen, die ich möchte, und installiert und holt dann diese Conda-Pakete, damit es den Prozess ausführen kann.

Du kannst sehen, dass es dort ein bisschen Zeit genommen hat, weil es die Umgebung erstellt und die Software zum ersten Mal installiert hat. Allerdings ist diese Umgebung gecacht, also wenn ich denselben Nextflow-Befehl erneut ausführe, sollte es viel schneller sein, weil es dieselbe Conda-Umgebung wiederverwenden wird.

Eines der coolen Dinge daran ist, dass diese Direktiven auf Prozessebene spezifiziert werden können, nicht nur für den gesamten Workflow. Also wenn du möchtest, kannst du mischen und kombinieren, welche Technologie für verschiedene Prozesse verwendet wird.

## 2. Weise Rechenressourcen mit Prozessdirektiven zu

Die Nextflow-Konfigurationsdatei kann viel mehr als nur Software-Packaging. Wir können Nextflow auch sagen, wie die Schritte in der Pipeline tatsächlich ausgeführt werden sollen. Ein Beispiel ist, einem Host-System mitzuteilen, welche Ressourcen für jede ausgeführte Aufgabe verfügbar gemacht werden sollten.

Standardmäßig gibt Nextflow nicht sehr viel. Es gibt eine einzelne CPU und nur zwei Gigabyte Speicher für jeden Prozess.

Das ist wahrscheinlich etwas, das wir ändern möchten, damit Prozesse, die lange laufen, mehr Ressourcen haben können und schneller laufen, aber es kann schwierig sein zu wissen, was man einem Prozess zuweisen soll. Nextflow hat einige nette Tricks auf Lager, um dir dabei zu helfen.

## 2.1. Führe den Workflow aus, um einen Ressourcennutzungsbericht zu generieren

Lass uns den Workflow erneut ausführen. Dieses Mal werde ich ein zusätzliches Argument hinzufügen, nämlich -with-report. Es ist eine Nextflow-Kernoption, also ist es ein einzelner Bindestrich. Und dann ein beliebiger Dateiname. In diesem Fall werde ich es report-config-one.html nennen.

Ich werde den Workflow erneut ausführen. Er wird genau wie zuvor laufen, aber er wird mir einen zusätzlichen Hilfsbericht geben, der, wie du sehen kannst, jetzt hier in der Seitenleiste aufgetaucht ist.

Ich werde auf diese Datei rechtsklicken, auf Download klicken, was sie von GitHub Codespaces auf mein lokales System herunterlädt, damit ich sie dann leicht im Webbrowser hier oben ansehen kann.

Dieser Bericht kann für jeden Nextflow-Durchlauf generiert werden und enthält viele Informationen. Er beginnt oben mit einigen Metadaten darüber, welcher Befehl verwendet wurde, wann der Workflow lief, wie lange es dauerte, aber wenn du nach unten scrollst, bekommen wir detailliertere Informationen über die Ressourcen, die von jedem Schritt in der Pipeline verwendet wurden.

Da jeder Prozess mehrmals für verschiedene Aufgaben läuft, haben wir ein Boxplot, das die Variation der Ressourcen zeigt, die wir für jeden Prozess verwendet haben.

Wenn ich ein bisschen weiter nach unten scrolle, sehe ich ähnliche Informationen über verwendeten Speicher und Aufgabendauer. Auch Lese-/Schreibvorgänge auf der Festplatte.

Du kannst dir vorstellen, dass dies bei einer großen Pipeline mit lang laufenden Aufgaben sehr informativ sein kann, wie man die Konfiguration der Ressourcen, die du anforderst, fein abstimmt, damit du nicht zu viel anforderst, aber auch genug bereitstellst, damit es schnell läuft.

Wenn ich den Bericht weiter nach unten scrolle, sehen wir auch eine Aufgabentabelle, die uns detaillierte Informationen über jede einzelne Aufgabe zeigt, die im Workflow ausgeführt wurde. Dies beinhaltet Informationen wie das aufgelöste Skript, das ausgeführt wurde.

Okay, lass uns zurück zu unserer Konfigurationsdatei gehen. Wir haben gesehen, dass wir für unseren Workflow wirklich nicht viel benötigt haben, also lass uns Nextflow mitteilen, dass wir nur ein Gigabyte Speicher für jeden Prozess im Workflow benötigen.

Wenn wir es jetzt so auf Prozessebene definieren, wird dies auf jeden einzelnen Prozess in der Pipeline angewendet.

## 2.3. Setze Ressourcenzuweisungen für einen einzelnen Prozess

Lass uns der Argumentation halber so tun, als würde cowpy wirklich viel schwere Arbeit leisten und mehr Ressourcen als die anderen Aufgaben benötigen. Wir können hier einen zusätzlichen Konfigurationsblock definieren, der nur auf diesen Prozess angewendet wird, indem wir withName cowpy verwenden.

Das nennt man einen Konfigurationsselektor, und wir können hier verschiedene Muster definieren, um verschiedene Prozesse zu matchen. Zum Beispiel könnte ich cow\* machen. Ich folge dem dann mit einigen geschweiften Klammern und lass uns ihm zwei Gigabyte Speicher statt einem geben und sagen wir zwei CPUs.

Jetzt wird Nextflow jedem Prozess im Workflow ein Gigabyte geben, außer dieser Anforderung, die spezifischer ist. Also überschreibt sie es. Und nur für alle Prozesse, die cowpy heißen, werden zwei Gigs Speicher und zwei CPUs bekommen.

Beachte, dass Nextflow clever in Bezug auf Ressourcennutzung ist. Also wenn du anfängst, diese Zahlen auf höhere Werte zu setzen, wirst du sehen, dass Nextflow beginnt, Job-Einreichungen nacheinander in die Warteschlange zu stellen, anstatt sie alle parallel auszuführen, damit es nicht mehr als die verfügbaren Ressourcen anfordert.

## 2.4. Führe den Workflow mit der geänderten Konfiguration aus

Lass uns versuchen, den Workflow erneut auszuführen und lass uns dieses Mal einen neuen Bericht speichern.

Okay, wir können diese Datei herunterladen und einen Blick darauf werfen.

Ja, unüberraschend sieht es im Grunde genau gleich aus, weil dies ein Dummy-Workflow ist, der nichts Echtes macht. Aber du kannst dir vorstellen, wie dieser iterative Ansatz, Limits zu definieren und echte Workflows mit dieser Art von Berichterstattung durchzuführen, es dir ermöglicht, einen evidenzbasierten Ansatz zum Setzen angemessener Konfiguration zu verfolgen und wirklich das Beste aus den dir zur Verfügung stehenden Rechenressourcen zu machen.

Du kannst anfangen, wirklich clever zu sein. Nextflow hat eine eingebaute Fähigkeit, Fehler erneut zu versuchen, und du kannst das in deiner Konfigurationsdatei nutzen, indem du eine Closure wie diese verwendest und dynamisch die Ressourcen setzt, die verfügbar gemacht werden. Also habe ich hier Nextflow gesagt, diese zwei Gigabyte mit dem retry-Versuch zu multiplizieren. Also wird der zweite Versuch vier Gigs bekommen, der dritte Versuch sechs Gigs und so weiter. Das geht ein bisschen über den Rahmen dieses Trainingskurses hinaus, aber wenn du interessiert bist, schau dir die Nextflow-Dokumentation an, die einen schönen Abschnitt über dynamische Retry-Logik hat.

## 2.5. Füge Ressourcenlimits hinzu

Nun, eine Sache, die dir daran auffallen könnte, ist, dass diese Art von Sache es ziemlich einfach machen kann, versehentlich über die auf deinem System verfügbaren Ressourcen hinauszugehen. Wenn du mehr Ressourcen anforderst, als verfügbar sind, wird Nextflow einen Fehler über deine Konfiguration werfen und den Durchlauf anhalten. Um das zu vermeiden, kannst du etwas verwenden, das Ressourcenlimits genannt wird.

Unter dem process-Scope in unserem Workflow können wir Ressourcenlimits wie dieses definieren, das ein Array nimmt, und wir können den maximalen Speicher, CPUs und Zeit spezifizieren, die auf diesem System verfügbar sind.

Das Setzen hoher Werte hier erhöht nicht die Menge an Ressourcen, die angefordert werden. Wir werden immer noch ein Gigabyte in unseren Anforderungen verwenden, aber es bedeutet, dass wenn eine dieser Anforderungen 750 erreicht, sie diese Obergrenze erreichen wird und nicht mehr als das angefordert wird, was bedeutet, dass Nextflow weiter laufen wird und nicht wegen nicht verfügbarer Ressourcen abstürzen wird.

Also ist das eine nette Schutzmaßnahme, die man verwenden kann, besonders wenn du dynamische Logik mit deiner Ressourcenzuweisung verwendest.

Die andere Situation, wo dies wirklich nützlich ist, ist, wenn du Pipelines verwendest, die öffentlich sind und nicht von dir kontrolliert werden. Sie könnten mit Konfigurationsstandards kommen, und Nextflow wird automatisch den richtigen Ansatz verfolgen, alle Ressourcenanforderungen zu begrenzen, um auf deinem System zu laufen.

Okay, großartig. Wir haben über Software gesprochen. Wir haben über Ressourcenzuweisung gesprochen und wir haben verschiedene Scopes von Konfiguration beschrieben, sowohl für alle Prozesse als auch für spezifische Prozesse.

## 3. Verwende eine Parameterdatei, um Workflow-Parameter zu speichern

Okay, als Nächstes werden wir unsere Aufmerksamkeit auf Parameter richten. Wir können Parameter in der Konfigurationsdatei definieren, genau wie wir es zuvor im Nextflow-Skript getan haben. Also params.greeting equals hello oder verwende params-Scope und setze foo equals bar.

Und das ist großartig, um Standardwerte für deinen Workflow zu setzen. Wenn du jedoch Pipelines ausführst, kann es schön sein, Parameter in einer JSON- oder einer YAML-Datei zu spezifizieren.

Die Verwendung einer Datei wie dieser ist viel besser als das Spezifizieren von Kommandozeilenoptionen mit --. Denn wenn du einen Workflow ausführst, musst du möglicherweise viele Parameter spezifizieren und es kann mühsam sein, sie alle in einer einzigen CLI zu schreiben und fehleranfällig. Außerdem ist es unwahrscheinlich, dass du dich an alle Parameter erinnern wirst, die du verwendet hast, also wenn du das in eine Datei kodierst, ist es einfacher, den Workflow in Zukunft erneut mit denselben Parametern zu starten.

Wir haben hier eine Beispieldatei namens test-params, und du kannst sehen, dass diese die drei Parameter spezifiziert, die wir in unserem Workflow haben, mit drei verschiedenen Werten. Persönlich finde ich YAML einfacher zu schreiben als JSON. Also nur um zu demonstrieren, dass es funktioniert, werde ich eine neue Datei namens test.yaml erstellen und diese hineinkopieren, die Anführungszeichen loswerden und speichern.

Diese JSON- und YAML-Dateien können einfacher zu schreiben sein, da sie eine vertrautere Syntax haben. Aber beachte, dass diese nur für Parameter sind und sie nur Schlüssel-Wert-Syntax wie diese akzeptieren.

## 3.1. Führe den Workflow mit einer Parameterdatei aus

Lass es uns ausprobieren. Mache denselben Befehl wie zuvor. Entferne den Bericht und ich werde -params-file test-params.yaml machen.

Nein, das ist eine Nextflow-Kernoption, also ist es ein einzelner Bindestrich.

Okay. Es hat den Workflow ausgeführt und es hat die Parameter in dieser YAML-Datei verwendet, anstatt dass ich sie alle in der Kommandozeile spezifiziere. Mag wie Overkill nur für dieses einfache Beispiel erscheinen, aber du kannst dir vorstellen, wenn du 10 oder 20 verschiedene Parameter hast, kann es mühsam sein, sie manuell einzutippen, und das ist einfach viel einfacher in einem Code-Editor zu bearbeiten und zum Zweck der Reproduzierbarkeit festzuhalten.

## 3. Bestimme, welche Executors verwendet werden sollen, um die Arbeit zu erledigen

Okay. Wir haben über Software-Packaging mit Docker und conda gesprochen. Wir haben über Prozessressourcenanforderungen mit CPUs und Speicher gesprochen. Und wir haben ein bisschen darüber gesprochen, wie man Parameter beim Ausführen von Workflows spezifiziert.

Die letzten Teile der Konfiguration sind wirklich die Ausführung, die zugrunde liegende Computing-Infrastruktur selbst, und das ist das wahre Juwel in der Krone von Nextflow: dass wir denselben Workflow über mehrere verschiedene Computing-Infrastrukturen hinweg ausführen können.

Ich werde tatsächlich für eine Sekunde zum geschriebenen Trainingsmaterial wechseln. Unter diesem Teil des Trainings können wir einige verschiedene Beispiele sehen, wie verschiedene Executors, in diesem Fall HPC-Scheduler, die Ressourcenanforderungen definieren, die benötigt werden, um einen Job einzureichen.

Also für Slurm hast du diese SBATCH-Header, die --mem und die CPU-Nummer definieren. Wenn du PBS verwendest, hast du andere Header, und wenn du Grid Engine verwendest, hast du wieder andere Header.

Du kannst dir vorstellen, dass es noch unterschiedlicher ist, wenn du in der Cloud ausführen möchtest, sei es AWS Batch, Google Cloud, Azure oder mehr.

Jede dieser zugrunde liegenden Computing-Infrastrukturen wird Executor genannt und Nextflow weiß, wie man mit all diesen verschiedenen Executors kommuniziert, um Jobs mit der korrekten Syntax einzureichen.

Die gute Nachricht ist, dass du darüber nicht Bescheid wissen musst. Alles, was du tun musst, ist Nextflow zu sagen, welchen Executor es verwenden soll.

## 3.1. Auf ein anderes Backend abzielen

Wir gehen zurück zu unserer Konfigurationsdatei und dem Prozess machen wir executor, und ich werde local tippen.

Local ist tatsächlich der Standard, wenn du keinen anderen Executor spezifizierst, wird local verwendet, und das bedeutet einfach dein Host-System, wo auch immer du Nextflow gestartet hast.

Ich könnte stattdessen Slurm spezifizieren. Und das würde Slurm-Jobs einreichen, oder ich könnte AWS Batch sagen, und das würde Jobs an AWS Batch einreichen.

Du benötigst in einigen Fällen zusätzliche Konfiguration, zum Beispiel benötigt das Ausführen in der Cloud bestimmte Anmeldeinformationen, aber wirklich ist das der Kern davon, und es kann so einfach sein wie eine oder zwei Zeilen Konfiguration, um deinen Workflow in einer völlig anderen Computing-Umgebung auszuführen.

Obwohl wir auf einem einfachen System innerhalb von Codespaces laufen, kann ich trotzdem ein bisschen damit herumspielen und so tun, als würden wir auf Slurm laufen. Wenn ich dann den Workflow erneut starte, Nextflow run hello-config. Es wird fehlschlagen, weil es nicht in der Lage sein wird, Jobs an Slurm einzureichen. Aber wir können immer noch in die Arbeitsverzeichnisse gehen und sehen, was Nextflow getan hat. Also wenn wir zu diesem Arbeitsverzeichnis gehen und uns .command.run ansehen. Du kannst oben in dieser Datei sehen, dass wir jetzt diese sbatch-Header-Zeilen haben, die versucht haben, die für den Slurm-Job benötigten Ressourcen zu spezifizieren.

## 4. Verwende Profile, um voreingestellte Konfigurationen auszuwählen

Okay, wir sind fast da. Der letzte Teil dieses Kapitels ist das Sprechen über Konfigurationsprofile. Wenn du deine Pipeline auf mehreren verschiedenen Systemen ausführst, könnte es ärgerlich sein, all diese verschiedenen Nextflow-Konfigurationsdateien zu haben, die du jedes Mal spezifizieren musst.

Stattdessen kannst du Gruppierungen von Konfiguration innerhalb deiner Nextflow-Konfigurationsdatei kodieren und diese Gruppen ein- und ausschalten, indem du ein Profil-Flag verwendest. Lass uns sehen, wie das aussieht.

## 4.1. Erstelle Profile zum Wechseln zwischen lokaler Entwicklung und Ausführung auf HPC

Wir werden in unserem Beispiel hier zwei Profile erstellen, eines für meinen Laptop und eines für ein schwereres HPC-System. Ich werde ein bisschen schummeln und einfach den Code aus dem Trainingsmaterial kopieren und hier einfügen.

Wir haben einen neuen Scope namens profiles, und dann haben wir einen Namen für jedes Profil, der beliebig sein kann. Und darin haben wir Konfiguration, die genau gleich aussieht wie die Top-Level-Konfiguration, die wir bereits geschrieben haben. Also haben wir wieder process-Scope, Docker-Scope.

Beim Profil namens my_laptop sage ich, dass es mit dem local-Executor ausgeführt werden soll, also auf meinem Host-System, und Docker verwendet werden soll.

Beim university_hpc-Profil hier sage ich, Slurm zu verwenden, um Jobs einzureichen, conda statt Docker zu verwenden, und ich spezifiziere verschiedene Ressourcenlimits, die zur Systemgröße von Nodes auf dem von mir verwendeten HPC passen könnten.

Standardmäßig wird keine dieser Konfigurationen verwendet, wenn ich Nextflow ausführe, ich muss spezifizieren, dass ich eines dieser Profile verwenden möchte.

## 4.2. Führe den Workflow mit einem Profil aus

Lass uns nextflow run hello-config machen. Und ich werde -profile machen, einzelner Bindestrich, weil es eine Nextflow-Kernoption ist. Und dann der Name, den ich ihm gegeben habe, nämlich my_laptop. Nextflow sollte jetzt den Konfigurationsblock verwenden, der innerhalb dieses Konfigurationsprofils spezifiziert wurde, und ihn anwenden, wenn es Nextflow ausführt. Wenn ich den anderen Konfigurationsblock verwenden wollte, müsste ich nur diesen Profilnamen wechseln. Viel einfacher zu merken. Viel einfacher zu verwenden.

## 4.3. Erstelle ein Test-Profil

Beachte, dass die Profile jede Art von Konfiguration haben können, also muss es nicht mit deiner Ausführungsumgebung zusammenhängen. Zum Beispiel, lass uns hier ein neues Profil erstellen, das einen Satz von Parametern hat. Wir können das zu tux ändern und zu my profile ändern, und jetzt, wenn wir profile test machen, wird es diese Parameter spezifizieren, die die Parameter überschreiben werden, die auf der obersten Ebene des Workflows spezifiziert sind.

Wenn du Nextflow ausführst, kannst du mehrere Profile verketten und sie werden nacheinander angewendet.

## 4.4. Führe den Workflow lokal mit dem Test-Profil aus

Also kann ich den vorherigen Befehl nehmen und Komma test machen. Das wird zuerst die my_laptop-Konfiguration anwenden und dann die test-Konfiguration anwenden. Wenn es eine Überschneidung gibt, wird das Profil auf der rechten Seite jede Konfiguration in vorherigen Profilen überschreiben. Wenn ich Enter drücke, lass uns sehen, was passiert.

Okay, wir haben hier eine neue Ergebnisdatei. Du kannst das My Profile sehen, das ich als eine der Optionen spezifiziert habe. Und wir können auch cowpy, my profile sehen, und tatsächlich, da ist tux. Also hat das funktioniert.

## Zusammenfassung

Okay! Fantastisch. Das war's. Du hast es bis zum Ende des Kurses geschafft. Du bekommst ein bisschen Feier-Konfetti. Gut gemacht für das Beenden dieses Kapitels.

[Nächstes Video-Transkript :octicons-arrow-right-24:](07_next_steps.md)
