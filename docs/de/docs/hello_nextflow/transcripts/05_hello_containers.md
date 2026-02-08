# Teil 5: Hello Containers - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../05_hello_containers.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur der Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen und Hintergrund

Hallo und willkommen zurück zu Hello Nextflow. Das ist Teil 5, genannt Hello Containers. In diesem Teil des Kurses sprechen wir darüber, wie man die Softwareanforderungen einer Pipeline so kapselt, dass die Nutzer\*innen der Pipeline nicht über die Installation der Software nachdenken müssen.

Falls du schon genauso lange in der Bioinformatik arbeitest wie ich, erinnerst du dich vielleicht an das, was ich oft die schlechten alten Zeiten nenne. Wenn du damals die Pipeline von jemand anderem ausführen oder deren Arbeit replizieren wolltest, hast du Stunden oder Tage damit verbracht, alle verschiedenen Software-Tools zu installieren, die verwendet wurden – in den gleichen Versionen, sie auf deinem Rechner zu kompilieren. Das war ein Albtraum. Es war wirklich schwierig.

Wenn du auf einem HPC gearbeitet hast, hast du vielleicht Environment-Module verwendet, bei denen die Systemadministrator\*innen versucht haben, Software für dich zu installieren. Das war okay, aber immer noch nicht perfekt.

Aber jetzt haben wir bessere Möglichkeiten. Nextflow hat integrierte Unterstützung für verschiedene Software-Container-Technologien. Docker ist die gängigste. Die werden wir heute verwenden. Sie funktioniert gut in Codespaces, gut auf deinem lokalen Computer und gut in der Cloud.

Aber auch Singularity oder Apptainer, die auf HPC-Systemen sehr verbreitet sind und im Grunde genommen genauso funktionieren. Oder Podman, Shifter – es gibt eine Reihe anderer, die alle sehr ähnlich sind.

Die eine zusätzliche Technologie, die ähnlich, aber nicht ganz gleich ist und die Nextflow unterstützt, ist Conda. Nextflow kann Conda-Umgebungen für dich prozessspezifisch verwalten, was viel besser ist, als eigene Conda-Umgebungen zu erstellen. Und auch das kann mit einer Pipeline ausgeliefert werden.

Wir starten dieses Kapitel mit einem Überblick über Container-Technologien, Docker und wie sie funktionieren. Die erste Hälfte machen wir manuell in Docker, damit du verstehst, was unter der Haube passiert und wie das funktioniert. Denn das ist wirklich wichtig, um zu verstehen, was Nextflow tut und wie dein Workflow bei der Ausführung funktioniert.

Springen wir also rüber zu unseren Codespaces. Ich habe wieder alles aufgeräumt, aber wenn wir zu Hello Containers gehen, solltest du sehen, dass alle unsere Skripte und alles da sind – genauso wie am Ende des Modules-Kapitels. Wir haben also unsere verschiedenen Module hier, die ich im modules-Verzeichnis erstellt habe.

Sie sind immer noch da. Sie müssen da sein, damit es laufen kann. Der Workflow und die Ausgabe sind alle gleich, außer dass wir den Ausgabe-Publishing-Pfad auf Hello Containers geändert haben, damit deine Dateien in diesem Verzeichnis landen.

Wir können das jetzt ausführen, um zu prüfen, ob es funktioniert, wenn du möchtest, oder wir können direkt im Terminal weitermachen.

## 1. Container 'manuell' verwenden

Wir werden Docker verwenden, um unsere Container zu verwalten. Ich kann prüfen, ob es in meinem Codespaces installiert ist, indem ich "docker -v" eingebe, was mir die installierte Version zeigt und dass alles ordnungsgemäß funktioniert.

Container und Docker haben zwei Konzepte, die wirklich wichtig sind. Eines nennt sich Image, das andere Container. Das Image ist sozusagen der Snapshot des gesamten Dateisystems, das du verwenden wirst, und der Container ist die laufende Umgebung. Du erstellst einen Container aus einem Image.

Sobald du in diesem Container bist, funktioniert er normalerweise wie ein komplettes Betriebssystem. Er ist von der Außenwelt abgeschnitten. Er ist von allem anderen getrennt, und das ist gut so. So erreichen wir eine so gute Reproduzierbarkeit mit Nextflow.

Denn für Aufgaben, die innerhalb eines Containers laufen, werden sie nicht durch Konfigurationsdateien auf deinem lokalen System beeinflusst. Keine anderen externen Einflüsse – sie laufen in ihrer eigenen kleinen Sandbox. Die Dateien werden dann auf sehr, sehr reproduzierbare Weise erzeugt, weil du die gleichen zugrunde liegenden Bibliotheken verwendest, alle gleichen Abhängigkeiten, genau die gleiche Software für jede Person, die auf jeder verschiedenen Rechenumgebung läuft. Was ehrlich gesagt fantastisch ist, und ich finde es immer noch erstaunlich, dass das funktioniert. Und selbst bis heute fasziniert es mich, dass das möglich ist.

## 1.1. Container-Image pullen

Wir werden also einige Docker-Images und Docker ausprobieren. Wenn du es auf deinem System ausführst, hat Docker eine Docker-Registry auf deinem Computer oder in diesem Fall in dem Codespace, die alle verschiedenen Images verfolgt, die in der Vergangenheit heruntergeladen und verwendet wurden, und die verschiedenen Layer, aus denen sie aufgebaut sind.

Wir können sehen, welche Images wir lokal mit Docker haben, indem wir "docker image ls" eingeben. In diesem Fall siehst du hier eine Reihe von Docker-Images, die alle mit der Einrichtung dieses Codespaces zu tun haben. Alles mit Dev-Containern und so weiter. Du brauchst dir also keine großen Sorgen um sie zu machen, aber während wir mehr Images hinzufügen und herunterladen, während dieser Kurs weitergeht, kannst du diese Liste überprüfen und du wirst sehen, dass die lokale Registry all diese Dinge verfolgt, die wir gepullt haben.

Aber wir holen uns ein neues, indem wir "docker pull" machen. Das sagt Docker, dass es ein neues Image aus dem Internet holen soll.

Dann geben wir die URI für diesen Container ein. Das könnte ein Docker-Image sein, das du lokal erstellt und dann ins Internet gepusht hast. Es könnte ein Image sein, das jemand anders erstellt hat. Es gibt viele, viele verschiedene Möglichkeiten, Docker-Images zu erstellen, aber eine der einfachsten Möglichkeiten ist, das outzusourcen und jemand anderen das für dich machen zu lassen.

Und was wir in diesem Tutorial verwenden werden, ist ein Service von Seqera namens Seqera Containers.

Seqera Containers ist völlig kostenlos und verwendet ein Open-Source-Tool, das wir entwickelt haben, genannt Wave, das gebaut wurde, um Container ergänzend zu Nextflow zu verwalten. Es behandelt viele der häufigen Anwendungsfälle, mit denen wir uns bei Nextflow beschäftigen.

Es ist sehr üblich, dass die Software, die wir brauchen, in Conda, in Bioconda oder Conda-Forge-Channels oder anderen domänenspezifischeren Channels verpackt ist. Und Wave und Seqera Containers sind wirklich gut darin, Images daraus zu erstellen.

Ich kann also zu dieser Web-Oberfläche gehen und wir werden mit dem Paket "cowpy" herumspielen. Ich gebe den Namen des Pakets ein, das ich möchte. Es sucht, es hat es im Python Package Index gefunden, also kann ich das verwenden. Oder wenn ich ein bisschen länger warte, durchsucht es Bioconda und Conda-Forge. Und du siehst, ich kann hier jeden Conda-Channel angeben. Wenn du also einen Nvidia-Channel oder etwas anderes finden willst, sollte das auch funktionieren.

Und dann kann ich angeben, ob ich möchte, dass es ein Docker-Image oder ein Singularity-Image für mich baut, und auch welche CPU-Architektur ich verwenden möchte – amd64 oder arm64.

Sobald die Bioconda-Ergebnisse aufgelistet sind, kann ich jetzt auch alle verschiedenen verfügbaren Versionen sehen. Ich werde das einfügen. Und jetzt könnte ich weitersuchen und mehr Pakete von Conda bekommen, wenn ich wollte, und diesen Container nach Belieben zusammenstellen, aber ich will nur dieses eine. Also klicke ich auf Get Container.

Jemand anders hat bereits den gleichen Container angefordert und er wird aus einer Registry zurückgegeben, also bekommen wir ihn sofort. Aber wenn noch nie jemand nach diesem Softwarepaket oder dieser Kombination von Softwarepaketen gefragt hätte, würden Wave und Seqera Containers es spontan für uns bauen.

Wir können diese URL kopieren und wir können auch die Build-Details ansehen. Und das zeigt uns, was der Service im Backend gemacht hat. Er hat eine Conda-Environment-Datei erstellt, ein Dockerfile, und dann ist das hier der Docker-Build-Prozess. Er hat auch einen Scan durchgeführt, einen Sicherheitsscan, sodass du alle CVEs sehen kannst. Und es sagt dir, wann das erstellt wurde.

Wave und Seqera Containers können viel mehr als das, aber das ist eine Art einfacher Anwendungsfall, der am häufigsten vorkommt. Und ich sollte sagen, dass diese Images für mindestens fünf Jahre gehostet werden. Du kannst also diese URLs in deine Pipelines einbauen und weißt, dass sie nicht so bald verschwinden werden.

Ich habe also meine URL für mein Docker-Image für cowpy.

Ich kann jetzt "docker pull" mit dieser URL machen, und es wird alle verschiedenen Layer holen und dieses Image herunterladen, sodass es für mich lokal verfügbar ist.

## 1.2. Container verwenden, um cowpy als einmaligen Befehl auszuführen

Okay, lass uns jetzt versuchen, es tatsächlich zu benutzen. Jetzt werde ich einen "docker run"-Befehl verwenden anstelle von "docker pull", und ich werde das Flag "--rm" verwenden, das Docker einfach sagt, diesen Container herunterzufahren, sobald er fertig ist mit dem, was ich von ihm verlangt habe. Dann gebe ich den Identifier für den Container ein, der einfach eine URI ist.

Und dann am Ende spezifiziere ich den Befehl, den Docker innerhalb des aus diesem Image erzeugten Containers ausführen soll. Ich sage einfach cowpy, was der Name des Tools ist, das von Conda-Forge installiert wurde und innerhalb des Images verfügbar ist.

Ich drücke Enter und da haben wir es. Wir haben cowpy auf einem System ausgeführt. Wir haben eine kleine Kuh, die uns einige Informationen gibt.

Beachte, dass cowpy nicht auf meinem lokalen System installiert ist. Wenn ich es also ohne all das Docker-Zeug ausführe, sagt es "command not found". Dies hat also ein Image gepullt. Es hat mit Docker einen Container erstellt und ist dann in diesen Container gegangen und hat diesen Befehl für uns ausgeführt und die Ausgabe an unser Terminal zurückgegeben. Sehr, sehr cool.

## 1.3. Container verwenden, um cowpy interaktiv auszuführen

Okay, wir gehen jetzt noch einen Schritt weiter und führen diesen Container interaktiv aus und schauen uns ein bisschen um, damit wir sehen können, was innerhalb des Containers passiert.

Wenn ich also zurückgehe und meinen Run-Befehl nehme, werde ich cowpy am Ende entfernen, denn ich will eigentlich nicht cowpy ausführen. Ich will ein Bash-Terminal ausführen.

Und dann gehe ich hierher zurück und mache "-it", was für Interactive und Terminal oder TTY steht, und ich drücke Enter.

Und jetzt siehst du, der Prompt, der Teil bevor ich tippe, hat sich geändert. Das war der Codespaces-Prompt, wo das Verzeichnis stand, und jetzt steht da base und root und tmp. Ich bin jetzt also innerhalb des Containers, und wenn ich "ls" mache, wirst du sehen, dass die Dateien, die ich in diesem Verzeichnis sehe, anders sind als die Dateien, die ich in meinem Workspace habe.

Und tatsächlich kann ich keine der Dateien aus meinem lokalen Codespaces-Workspace oder meiner lokalen Festplatte innerhalb des Docker-Containers sehen. Die Docker-Container-Laufzeit ist vollständig isoliert und kann keine Dateien aus einem Host-Dateisystem außerhalb lesen oder schreiben.

Ich kann jedoch die Software sehen, die innerhalb des Containers installiert ist, und sie ausführen. Ich kann also cowpy ausführen und wir können ein bisschen mehr darüber sehen, wie man cowpy benutzt. Hier kann ich "cowpy 'Hello World'" machen und das sagt ihm, dass es mein Zitat in eine kleine Sprechblase setzen soll. Und man kann auch verschiedene Arten von Kühen ausführen, es muss also keine Kuh sein. Man kann ein "-c" machen. Und ich bin in Schweden, also wähle ich einen Elch. Sehr schön. Hat ihm ein Geweih gegeben.

Und es gibt eine ganze Reihe verschiedener, mit denen du herumspielen kannst, die in den Trainingsdokumenten beschrieben sind.

## 1.3.4. Daten in den Container mounten

Okay. Es wäre schön, wenn wir cowpy auf den Dateien in unserem Dateisystem ausführen könnten.

Natürlich ist es nicht besonders nützlich, nur den Container zu haben und gar keinen Zugriff auf irgendetwas. Es mag sicher und reproduzierbar sein, aber es ist nicht sehr nützlich.

Wie machen wir das also? Ich werde aus diesem Docker-Container aussteigen, indem ich exit eingebe, und du siehst, der Prompt sagt uns, dass wir jetzt wieder in unserem normalen Codespaces sind.

Und ich werde den gleichen Befehl nochmal ausführen. Aber diesmal werde ich hier hinten ein paar zusätzliche Flags hinzufügen. Und das wichtige ist "-v", was für das Mounten eines Volumes steht, was im Grunde wie ein Teil eines Speicherplatzes ist.

Das "-v" nimmt zwei Teile: es gibt einen String, dann einen Doppelpunkt und einen String. Und der erste Teil ist das lokale Dateisystem, das in den Container gemountet werden soll. Und dann ist der zweite Teil, wo das innerhalb des Containers landen soll.

Ich will hier einfach mein ganzes lokales Dateisystem laden. Also "." ist das aktuelle Arbeitsverzeichnis. Ich mache also einfach "." und dann ":", und dann werden wir das in ein neues Verzeichnis innerhalb des Containers namens "my_project" setzen. Das könnte wirklich alles heißen.

Und dann führe ich wieder aus.

Im Arbeitsverzeichnis, in dem ich abgelegt werde, was /tmp ist, sind die Dateien nicht da. Aber wenn ich "ls my_project" mache, da haben wir es: All die gleichen Dateien, die wir lokal auf unseren Codespaces hatten, sind jetzt innerhalb des Containers unter diesem Pfad verfügbar.

Das ist Lese- und Schreibzugriff, sodass ich neue Dateien in diesem Verzeichnis erstellen kann und sie werden auf meinem Host-Dateisystem erscheinen. Dieses spezielle Verzeichnis verhält sich dann genau so, als wäre ich außerhalb des Containers, sodass ich jetzt lesen und schreiben und Dinge tun kann.

## 1.3.5. Die gemounteten Daten verwenden

Okay, lass uns einfach beweisen, dass wir das tun können. Ich mache "cat /my_project/data/greetings.csv". Wenn du dich erinnerst, sieht der Inhalt dieser Datei so aus. Ich kann das jetzt zu cowpy pipen und die Kuh wird die verschiedenen Ausgaben dieser Datei in ihrer kleinen Sprechblase ausdrucken, was irgendwie lustig ist.

Du siehst also, wir können jetzt die Software im Container verwenden, um mit den Dateien auf unserem Host-System zu interagieren.

Okay, lass uns zurückspringen und mit dem Rest des Trainingsmaterials weitermachen.

## 2. Container in Nextflow verwenden

Das ist also wirklich cool, Container zu verwenden. Hoffentlich macht das Sinn. Und du kannst den Wert dieser Container sehen und warum das nützlich ist für die Ausführung von Analysesoftware.

Aber wie machen wir diesen ganzen gleichen Prozess innerhalb von Nextflow? Wir wollen nicht selbst eine Menge Docker-Befehle ausführen. Wir wollen einfach, dass Nextflow das alles für uns handhabt.

Also lass uns das durcharbeiten. Wir werden einen neuen Prozess zu unserer Pipeline hinzufügen, um cowpy auszuführen. Okay, lass uns also ein neues Modul für unseren neuen Prozess erstellen. Also gehen wir in modules, nennen wir es cowPy.nf, und dann kopiere ich den Code aus dem Trainingsmaterial hier.

Aber du siehst, der Prozess ist sehr einfach. Er sieht aus wie die, die wir bisher gemacht haben, wir haben einen Input-Block mit einem Path, der unsere Eingabedatei ist, und auch einen Value hier, sodass das ein Zeichen sein wird, damit wir wieder einen Elch verwenden können, wenn wir wollen.

Und dann eine Ausgabe, die eine einzelne Datei hier ist, ein Path, und dann ein Script. Und wir machen das Gleiche, was wir interaktiv innerhalb des Containers gemacht haben: wir machen "cat", um die Eingabedatei zu lesen. Wir pipen diesen Inhalt zu cowpy. Wir wählen ein bestimmtes Zeichen basierend auf dieser Eingabe, wir schreiben in eine Ausgabedatei namens cowpy input file, die dann zur Ausgabe geechot wird.

Großartig. Lass uns das einbinden. Also include \{ cowpy \} from "./modules/cowpy.nf", habe ich es cowpy genannt? Ja.

Und dann lass uns unseren neuen Prozess hier unten im main-Block des Workflows aufrufen. Also lass uns cowpy ausführen. Und wir nehmen unseren neuen cowpy-Prozess und wir sagen collectGreetings.out.

Und wenn du dich erinnerst, gab es zwei Ausgaben für dieses Modul. Eine namens outfile und eine namens report. Die VS-Code-Extension schlägt diese automatisch für uns vor und wir wollen .outfile.

Du kannst immer in diesen Prozess hier springen. Du hoverst entweder darüber und es sollte dir schnell zeigen, was die Ausgaben waren. Und wir können auch mit Befehl-Klick hineinklicken und es wird die Moduldatei öffnen, wenn du mehr Details sehen willst.

Also hier. Das ist das outfile da, und das ist der Path. Das wird jetzt also die Eingabedatei für unseren cowpy-Prozess sein. Fantastisch.

Wenn du dich erinnerst, hat ein cowpy-Prozess zwei Eingaben. Wir hatten auch den Value-Kanal für das Zeichen. Also können wir "params.character" hier hinzufügen. Ich hätte das hart codieren können, wenn ich wollte, aber lass es uns zu einer CLI-Option machen, sodass wir --character machen können.

Richtig. Ich muss jetzt den Eingabeparameter definieren, den wir gerade aufgerufen haben, und ihm einen Standardwert geben. Also character, String. Und ich mag den Elch, also setze ich ihn standardmäßig auf moose.

Richtig, lass uns versuchen, es auszuführen. Also wenn ich Nextflow run hello containers mache, werden wir sehen, was passiert.

Ich hätte -resume verwenden können, wenn ich die alten Work-Verzeichnisse herumliegen hätte. Und wieder wären diese ersten Prozesse gecacht worden und es wäre ein bisschen schneller gewesen, aber es sollte im Grunde das Gleiche sein.

Jetzt sehen wir sofort, dass es einen Fehler geworfen hat, als es zu unserem neuen Prozess kam, es sagt uns hier, dass es einen Fehler beim Ausführen des cowpy-Prozesses gab und es mit einem Exit-Status 127 beendet wurde. Das ist der Befehl, den es versucht hat auszuführen. Es sieht richtig aus, es sieht aus, wie wir es erwartet haben. Es nimmt diesen Ausgabedateinamen, der etwa richtig aussieht, es führt ihn mit einem moose-Zeichen aus und versucht zu speichern.

Aber du siehst, der Command-Error hier sagt, cowpy command's not found. Und das macht Sinn, weil wir Nextflow noch nicht gesagt haben, dass es einen Container verwenden soll. Wir haben ihm nur den cowpy-Befehl gegeben. Und wie ich vorher sagte, ist cowpy nicht auf unserem lokalen System installiert. Als es also versuchte, es auszuführen, ist es fehlgeschlagen.

## 2.3.1. Einen Container für cowpy angeben

Wir müssen Nextflow sagen, dass ein Container verfügbar ist und dass es ihn verwenden kann. Wie machen wir das?

Wenn wir in unser Modul hier hineingehen, werden wir oben eine neue Deklaration namens "container" hinzufügen. Und wir werden das dann auf einen String setzen.

Wenn du dich erinnerst, kann ich in Seqera Containers diese URL kopieren und ich lasse sie einfach hier in Anführungszeichen fallen.

Jetzt zurückgehen und versuchen, es nochmal auszuführen.

Mal sehen, ob es diesmal funktioniert.

Leider schlägt es auf genau die gleiche Weise fehl, obwohl wir jetzt einen Container für den Prozess definiert haben. Um unser Docker-Image zu verwenden, müssen wir Nextflow sagen, dass es die Docker-Nutzung aktivieren soll, wenn wir den Workflow ausführen.

Und wir werden das tun, indem wir eine neue Config-Datei erstellen. Ich werde also sagen touch nextflow.config.

Das ist ein spezieller Dateiname, wo, wenn er im Arbeitsverzeichnis ist, während ich die Pipeline starte, er automatisch geladen wird. Wenn ich also in diese nextflow.config-Datei gehe, siehst du, sie existiert tatsächlich schon, was ich vergessen hatte. Und wir haben docker.enabled hier bereits drin, aber es ist auf false gesetzt, was der Standard ist.

Wenn ich das also auf equals True ändere, docker.enabled. Und es gibt Referenzdokumente für all diese Config-Scopes in den Nextflow-Docs. Und auch, wenn ich mit der VS-Code-Extension darüberhover, zieht sie die Docs spezifisch dafür herein und sagt mir, was es bedeutet und wie man es setzt.

Also haben wir es jetzt auf true gesetzt, und wenn ich Nextflow nochmal ausführe, wird Nextflow jetzt wissen, dass es Docker-Image für uns holen soll, wenn wir es noch nicht lokal haben, und dann diesen Prozess mit dieser Container-Umgebung ausführen.

Und wir können sehen, dass es erfolgreich ausgeführt wurde und wir haben ein kleines Häkchen neben cowpy. Fantastisch. Wenn ich hochgehe und im results-Verzeichnis nachschaue, ist die Datei noch nicht da. Und das liegt daran, dass wir diese Ausgabedatei noch publishen müssen, genau wie alle anderen.

Also gehen wir zum published-Block innerhalb des Workflows, sagen mycowpy equals cowpy.out.

Und dann hier unten im output-Block, mycowpy, geschweifte Klammern path. Hoppla. Hello containers. Mode, copy.

Wenn ich jetzt nochmal ausführe, sollte es auf genau die gleiche Weise laufen. Ich hätte -resume verwenden können und ich vergesse es jedes Mal. Und dann gehe ich hoch und jetzt haben wir eine neue Datei erstellt namens cowpy-COLLECTED, und da ist mein Elch, der BONJOUR, HELLO, HOLà sagt. Fantastisch.

Jetzt könnte ich natürlich auch "--character" übergeben. Was sind die verschiedenen Optionen? Ich glaube, es gibt einen Turkey? Also kann ich character Turkey verwenden. Es wird auf genau die gleiche Weise laufen. Ich habe eine weitere Gelegenheit verpasst, -resume zu verwenden, und jetzt, wenn wir unsere Datei laden, haben wir jetzt einen Truthahn. Fantastisch.

## 2.3.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat

Okay. Eine letzte kleine Sache. Lass uns diesen Befehl nochmal schnell ausführen, diesmal resume, und einen kurzen Blick in das work-Verzeichnis werfen, um zu sehen, was Nextflow unter der Haube macht, damit das alles für uns funktioniert.

Diesmal ist es super schnell, lass uns in dieses work-Verzeichnis gehen, cd work/. Wenn du dich erinnerst, haben wir hier eine Reihe von Punktdateien und die, die uns in diesem Fall interessiert, ist die, von der ich sagte, dass wir sie fast nie anschauen müssen, genannt .command.run.

Wenn ich code .command.run mache, wird es sie im Editor öffnen. Und ich kann in dieser Datei suchen und wenn ich runterscrolle, sollte ich Docker run sehen. Und du siehst, Nextflow macht den docker run-Befehl für uns, wenn Docker im Config aktiviert ist. Es hat eine ganze Reihe verschiedener Flags und Sachen hier, aber du siehst das "-v"-Flag, das wir selbst verwendet haben, als wir gelaufen sind. Und du siehst, es mountet das lokale Workspace-Verzeichnis in den Container, sodass der Container auf unsere Eingabedateien zugreifen und die Ausgaben speichern kann. Und dann am Ende führt es auch .command.sh aus, was das generierte Script ist, das den cowpy-Befehl drin hat.

Und du siehst, dass Nextflow die Workflow-Logik nimmt, die das Zeug ist, das uns tatsächlich interessiert, das spezifisch für unsere Analyse ist, und es macht all das clevere Zeug hinter den Kulissen, damit Docker auf unserem System funktioniert.

Und es macht das auf eine wirklich portable Weise, sodass ein Endnutzer der Pipeline die Technologie, die er verwendet, austauschen kann: Docker, Singularity, Apptainer, Conda. Das spielt für die Pipeline-Logik keine Rolle, aber Nextflow wird alle zugrunde liegenden Infrastrukturbedürfnisse handhaben, sodass es überall läuft.

Und das ist wirklich die Superkraft von Nextflow: Reproduzierbarkeit und Portabilität. Und mit Nextflow kannst du deinen Workflow tatsächlich teilen und andere Leute können ihn auf ihren Systemen ausführen und es wird einfach funktionieren.

Das ist eine wirklich, wirklich schwierige Sache, und jetzt weißt du auch, wie man das mit deinen Workflows macht.

Okay, das war's für dieses Kapitel. Wenn du zum Ende des Kurses gehst, findest du wieder ein Quiz über Container. Hoffentlich hat das alles Sinn gemacht. Es ist ein wirklich cooler Weg, mit Analysen zu arbeiten. Und wenn du neu bei Containern bist, hoffe ich, dass ich dich überzeugt habe, dass das der richtige Weg ist, und du nie zurückschauen wirst.

Aber damit, mach vielleicht eine kleine Pause, und du triffst mich in ein paar Minuten, um den finalen Teil 6 von Hello Nextflow durchzugehen, der ganz über Konfiguration geht.

Vielen Dank.
