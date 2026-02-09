# Teil 5: Hello Containers - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../05_hello_containers.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern im Material.

## Willkommen und Hintergrund

Hallo und willkommen zurück zu Hello Nextflow. Dies ist Teil fünf mit dem Titel Hello Containers. In diesem Teil des Kurses sprechen wir darüber, wie du die Software-Anforderungen für eine Pipeline kapseln kannst, sodass Personen, die die Pipeline ausführen, sich keine Gedanken über die Installation der Software machen müssen.

Wenn du schon so lange in der Bioinformatik arbeitest wie ich, erinnerst du dich vielleicht an das, was ich oft die schlechten alten Zeiten nenne. Wenn du die Pipeline von jemand anderem ausführen oder deren Arbeit replizieren wolltest, hast du Stunden oder Tage damit verbracht, all die verschiedenen Software-Tools zu installieren, die sie verwendet haben, in denselben Versionen, und versucht, sie auf deinem Rechner zu kompilieren – es war ein Albtraum. Es war wirklich schwierig.

Wenn du auf einem HPC gearbeitet hast, hast du vielleicht Umgebungsmodule verwendet, bei denen die Systemadministrator\*innen versucht haben, Software für dich zu installieren, was okay war, aber immer noch unvollkommen.

Aber jetzt haben wir bessere Möglichkeiten, dies zu tun. Nextflow hat integrierte Unterstützung für verschiedene Software-Container-Technologien. Docker ist die gängigste. Das werden wir heute verwenden. Es funktioniert gut in Codespaces, auf deinem lokalen Computer und in der Cloud.

Aber auch Singularity oder Apptainer, die auf HPC-Systemen sehr verbreitet sind und im Grunde genau gleich funktionieren. Oder Podman, Shifter – es gibt eine Reihe anderer, die alle sehr ähnlich sind.

Die eine zusätzliche Technologie, die ähnlich, aber nicht ganz gleich ist und die Nextflow unterstützt, ist Conda. Nextflow kann Conda-Umgebungen für dich auf Prozessbasis verwalten, was viel besser ist, als deine eigenen Conda-Umgebungen zu erstellen. Und auch das kann mit einer Pipeline ausgeliefert werden.

Wir beginnen dieses Kapitel damit, ein wenig über Container-Technologien und Docker zu sprechen und wie sie funktionieren. Die erste Hälfte machen wir manuell in Docker, damit du verstehst, was unter der Haube passiert und wie das funktioniert. Denn das ist wirklich wichtig, um zu verstehen, was Nextflow macht und wie dein Workflow bei der Ausführung funktioniert.

Also, springen wir zu unseren Codespaces. Ich habe wieder alles aufgeräumt, aber wenn du zu Hello Containers gehst, solltest du sehen, dass alle unsere Skripte und alles andere genauso da sind wie am Ende des Modul-Kapitels. Wir haben also unsere verschiedenen Module hier, die ich im modules-Verzeichnis erstellt habe.

Sie sind immer noch da. Sie müssen da sein, damit es laufen kann. Der Workflow und die Ausgabe sind alle gleich, außer dass wir den Ausgabe-Veröffentlichungspfad auf Hello Containers geändert haben, sodass deine Dateien in diesem Verzeichnis landen.

Wir können das jetzt ausführen, um zu prüfen, ob es funktioniert, wenn du möchtest, oder wir können mit dem Terminal weitermachen.

## 1. Einen Container 'manuell' verwenden

Wir werden Docker verwenden, um unsere Container zu verwalten, und ich kann überprüfen, ob es in meinen Codespaces installiert ist, indem ich "docker -v" eingebe, was mir die installierte Version zeigt und dass alles ordnungsgemäß funktioniert.

Container und Docker haben zwei Konzepte, die wirklich wichtig sind. Eines heißt Image und eines heißt Container. Das Image ist sozusagen der Snapshot des gesamten Dateisystems, das du verwenden wirst, und der Container ist die laufende Umgebung. Du erstellst einen Container mit einem Image.

Sobald du in diesem Container bist, funktioniert er normalerweise wie ein ganzes Betriebssystem. Er ist von der Außenwelt abgeschnitten. Er ist von allem anderen getrennt, und das ist gut so. So erreichen wir eine so gute Reproduzierbarkeit mit Nextflow.

Denn für Aufgaben, die innerhalb eines Containers ausgeführt werden, werden sie nicht durch Konfigurationsdateien auf deinem lokalen System beeinträchtigt. Keine anderen externen Einflüsse – sie laufen in ihrer eigenen kleinen Sandbox. Die Dateien werden dann auf sehr, sehr reproduzierbare Weise erstellt, weil du dieselben zugrunde liegenden Bibliotheken, alle dieselben Abhängigkeiten, genau dieselbe Software für jede Person verwendest, die in jeder unterschiedlichen Rechenumgebung läuft. Was ehrlich gesagt fantastisch und erstaunlich ist, dass es funktioniert. Und selbst heute noch irgendwie meinen Verstand sprengt, dass das möglich ist.

## 1.1. Das Container-Image pullen

Wir werden also versuchen, einige Docker-Images und Docker zu verwenden. Wenn du es auf deinem System ausführst, gibt es eine Docker-Registry auf deinem Computer oder in diesem Fall im Codespace, die alle verschiedenen Images verfolgt, die in der Vergangenheit heruntergeladen und verwendet wurden, und die verschiedenen Schichten, aus denen sie aufgebaut sind.

Wir können sehen, welche Images wir lokal mit Docker haben, indem wir "docker image ls" eingeben. In diesem Fall siehst du eine Reihe von Docker-Images hier, die alle mit der Einrichtung dieser Codespaces zu tun haben. Alles zu tun mit Dev-Containern und so. Du musst dir darüber nicht zu viele Gedanken machen, aber wenn wir im Laufe dieses Kurses weitere Images hinzufügen und herunterladen, kannst du diese Liste überprüfen und sehen, dass die lokale Registry all diese Dinge verfolgt, die wir gepullt haben.

Aber wir werden ein neues holen, indem wir "docker pull" eingeben. Das sagt Docker, ein neues Image aus dem Web zu holen.

Wir geben dann die URI für diesen Container ein. Das könnte ein Docker-Image sein, das du lokal erstellt und dann ins Internet gepusht hast. Es könnte ein Image sein, das jemand anderes erstellt hat. Es gibt viele, viele verschiedene Möglichkeiten, Docker-Images zu erstellen, aber wohl eine der einfachsten Möglichkeiten ist, das auszulagern und jemand anderen das für dich machen zu lassen.

Und was wir in diesem Tutorial verwenden werden, ist ein Service von Seqera namens Seqera Containers.

Seqera Containers ist völlig kostenlos und verwendet ein Open-Source-Software-Stück, das wir entwickeln, genannt Wave, das gebaut wurde, um Container auf komplementäre Weise zu Nextflow zu verwalten. Es behandelt viele der häufigen Anwendungsfälle, mit denen wir uns bei Nextflow befassen.

Es ist sehr üblich, dass die Software, die wir benötigen, in Conda, in Bioconda- oder conda-forge-Kanälen oder anderen domänenspezifischeren Kanälen verpackt ist. Und Wave und Seqera Containers sind wirklich gut darin, Images daraus zu erstellen.

Ich kann also zu dieser Web-UI gehen und wir werden mit dem Paket namens "cowpy" herumspielen. Ich gebe den Namen des Pakets ein, das ich möchte. Es sucht, es hat es im Python Package Index gefunden, also kann ich das verwenden. Oder wenn ich etwas länger warte, durchsucht es bioconda und conda-forge. Und du siehst, ich kann hier jeden Conda-Kanal angeben. Wenn du also einen Nvidia-Kanal oder etwas anderes finden möchtest, sollte das auch funktionieren.

Und dann kann ich angeben, ob ich möchte, dass es ein Docker-Image oder ein Singularity-Image für mich erstellt, und auch welche CPU-Architektur ich verwenden möchte. Also amd64 oder arm64.

Und sobald die bioconda-Ergebnisse aufgelistet sind, kann ich jetzt auch alle verschiedenen verfügbaren Versionen sehen. Ich werde das einfügen. Und jetzt könnte ich weiter suchen und mehr Pakete aus Conda holen, wenn ich möchte, und diesen Container zusammenstellen, wie ich möchte, aber ich möchte nur dieses eine. Also klicke ich auf Get Container.

Jemand anderes hat bereits denselben Container angefordert und er wird aus einer Registry zurückgegeben, also bekommen wir ihn sofort. Aber wenn niemand sonst jemals dieses Softwarepaket oder diese Kombination von Softwarepaketen angefordert hätte, würden Wave und Seqera Containers es spontan für uns erstellen.

Wir können diese URL kopieren und wir können auch die Build-Details anzeigen. Und das zeigt uns, was der Service im Backend gemacht hat. Er hat eine Conda-Umgebungsdatei erstellt. Ein Dockerfile, und dann ist das hier, wie es den Docker-Build-Prozess ausführt. Es hat auch einen Scan durchgeführt, einen Sicherheitsscan, sodass du alle CVEs sehen kannst. Und es sagt dir, wann das erstellt wurde.

Wave und Seqera Containers können viel mehr als das, aber das ist eine Art einfacher Anwendungsfall, der am häufigsten vorkommt. Und ich sollte sagen, dass diese Images für mindestens fünf Jahre gehostet werden. Du kannst also diese URLs in deine Pipelines einbauen und wissen, dass sie nicht so bald verschwinden werden.

Ich habe also meine URL für mein Docker-Image für cowpy.

Ich kann jetzt "docker pull" mit dieser URL machen, und es wird alle verschiedenen Schichten holen und dieses Image herunterladen, sodass es lokal für mich verfügbar ist.

## 1.2. Den Container verwenden, um cowpy als einmaligen Befehl auszuführen

Okay, jetzt versuchen wir, es tatsächlich zu verwenden. Ich werde jetzt einen "docker run"-Befehl anstelle von "docker pull" verwenden, und ich werde das Flag "--rm" verwenden, das Docker einfach sagt, diesen Container herunterzufahren, sobald er fertig ist mit dem, was ich ihn gebeten habe zu tun. Und dann gebe ich den Identifier für den Container ein, der einfach eine URI ist.

Und dann am Ende gebe ich den Befehl an, den Docker innerhalb des Containers ausführen soll, der aus diesem Image generiert wurde. Ich sage einfach cowpy, was der Name des Tools ist, das von Conda Forge installiert wurde und innerhalb des Images verfügbar ist.

Ich drücke Enter und da haben wir es. Wir haben cowpy auf einem System ausgeführt. Wir haben eine kleine Kuh, die uns einige Informationen gibt.

Beachte, dass cowpy nicht auf meinem lokalen System installiert ist. Wenn ich es also einfach ohne all das Docker-Zeug ausführe, sagt es: Befehl nicht gefunden. Das hat also ein Image gepullt. Es hat einen Container mit Docker erstellt, und dann ist es in diesen Container gegangen und hat diesen Befehl für uns ausgeführt und uns die Ausgabe zurück zu unserem Terminal gegeben. Sehr, sehr cool.

## 1.3. Den Container verwenden, um cowpy interaktiv auszuführen

Okay, wir gehen jetzt einen Schritt weiter und führen diesen Container interaktiv aus und schauen uns um, damit wir sehen können, was innerhalb des Containers passiert.

Wenn ich also zurückgehe und meinen run-Befehl nehme und ich werde cowpy am Ende dort loswerden, weil ich eigentlich nicht cowpy ausführen möchte. Ich möchte ein Bash-Terminal ausführen.

Und dann gehe ich zurück hierher und mache "-it", was für Interactive und Terminal oder TTY steht, und ich drücke Enter.

Und jetzt siehst du, dass sich die Eingabeaufforderung, der Teil vor dem, was ich tippe, geändert hat. Das war die Codespaces-Eingabeaufforderung, wo es das Verzeichnis anzeigte, und jetzt steht da base und root und tmp. Ich bin jetzt also innerhalb des Containers, und wenn ich "ls" mache, siehst du, dass die Dateien, die ich in diesem Verzeichnis sehe, sich von den Dateien unterscheiden, die ich in meinem Workspace habe.

Und tatsächlich kann ich keine der Dateien aus meinem lokalen Codespaces-Workspace oder meiner lokalen Festplatte innerhalb des Docker-Containers sehen. Die Docker-Container-Laufzeitumgebung ist vollständig isoliert und kann keine Dateien aus einem Host-Dateisystem außerhalb lesen oder schreiben.

Ich kann jedoch die Software sehen, die innerhalb des Containers installiert ist, und sie ausführen. Ich kann also cowpy ausführen und wir können ein bisschen mehr darüber sehen, wie man cowpy verwendet. Hier kann ich "cowpy 'Hello World'" machen und das sagt ihm, mein Zitat tatsächlich in eine kleine Sprechblase zu setzen. Und du kannst auch verschiedene Arten von Kühen ausführen, es muss also keine Kuh sein. Du kannst ein "-c" machen. Und ich bin in Schweden, also wähle ich einen Elch. Sehr schön. Hab ihm ein paar Geweihe gegeben.

Und es gibt eine ganze Reihe verschiedener, mit denen du herumspielen kannst, die in den Trainingsdokumenten beschrieben sind.

## 1.3.4. Daten in den Container mounten

Okay. Es wäre schön, wenn wir cowpy auf den Dateien in unserem Dateisystem ausführen könnten.

Natürlich ist es nicht sehr nützlich, nur den Container zu haben und keinen Zugriff auf irgendetwas. Es mag sicher und reproduzierbar sein, aber es ist nicht sehr nützlich.

Wie machen wir das also? Ich werde aus diesem Docker-Container herauskommen, indem ich exit tippe, und du kannst sehen, dass die Eingabeaufforderung uns sagt, dass wir jetzt wieder in unseren regulären Codespaces sind.

Und ich werde denselben Befehl noch einmal ausführen. Aber dieses Mal werde ich hier einige zusätzliche Flags hinzufügen. Und das wichtige ist "-v", was für das Mounten eines Volumes steht, was im Grunde wie ein Teil eines Festplattenspeichers ist.

Das "-v" nimmt zwei Teile: Es gibt einen String und dann einen Doppelpunkt und einen String. Und der erste Teil ist das lokale Dateisystem, das in den Container gemountet werden soll. Und dann ist der zweite Teil, wo das innerhalb des Containers landen soll.

Ich möchte jetzt einfach mein gesamtes lokales Dateisystem hier laden. Also "." ist das aktuelle Arbeitsverzeichnis. Ich mache also einfach "." und dann ":", und dann werden wir das in ein neues Verzeichnis innerhalb des Containers namens "my_project" legen. Das könnte wirklich alles heißen.

Und dann führe ich es wieder aus.

Im Arbeitsverzeichnis, wo ich abgelegt werde, das ist /tmp, sind die Dateien nicht da. Aber wenn ich "ls my_project" mache, da haben wir es: Alle dieselben Dateien, die wir lokal in unseren Codespaces hatten, sind jetzt innerhalb des Containers unter diesem Pfad verfügbar.

Das ist Lese- und Schreibzugriff, sodass ich neue Dateien in diesem Verzeichnis erstellen kann und sie werden in meinem Host-Dateisystem angezeigt. Dieses spezielle Verzeichnis verhält sich dann genau so, als wäre ich außerhalb des Containers, sodass ich jetzt lesen und schreiben und Dinge tun kann.

## 1.3.5. Die gemounteten Daten verwenden

Okay, lass uns einfach beweisen, dass wir das tun können. Ich mache "cat /my_project/data/greetings.csv". Wenn du dich erinnerst, sieht dieser Dateiinhalt so aus. Ich kann das jetzt zu cowpy pipen und die Kuh wird die verschiedenen Ausgaben dieser Datei in ihrer kleinen Sprechblase ausdrucken, was irgendwie lustig ist.

Du siehst also, wir können jetzt die Software im Container verwenden, um mit den Dateien auf unserem Host-System zu interagieren.

Okay, lass uns wieder rausgehen und wir machen mit dem Rest des Trainingsmaterials weiter.

## 2. Container in Nextflow verwenden

Das ist also wirklich cool, Container zu verwenden. Hoffentlich macht das Sinn. Und du kannst den Wert dieser Container sehen und warum das nützlich ist, um Analysesoftware auszuführen.

Aber wie machen wir diesen ganzen Prozess innerhalb von Nextflow? Wir wollen nicht selbst eine Menge Docker-Befehle ausführen. Wir wollen einfach Nextflow das alles für uns handhaben lassen.

Also lass uns das durcharbeiten. Wir werden einen neuen Prozess zu unserer Pipeline hinzufügen, um cowpy auszuführen. Okay, also lass uns ein neues Modul für unseren neuen Prozess erstellen. Gehe also in modules, nennen wir es cowPy.nf, und dann kopiere ich den Code aus dem Trainingsmaterial hier.

Aber du siehst, der Prozess ist sehr einfach. Er sieht ähnlich aus wie die, die wir bisher gemacht haben. Wir haben einen input-Block mit einem path, der unsere Eingabedatei ist, und auch einen value hier, sodass dies ein Zeichen sein wird, also könnten wir wieder einen Elch verwenden, wenn wir wollen.

Und dann eine Ausgabe, die eine einzelne Datei hier ist, ein path und dann ein script. Und wir machen dasselbe, was wir interaktiv innerhalb des Containers gemacht haben: Wir machen "cat", um die Eingabedatei zu lesen. Wir pipen diesen Inhalt zu cowpy. Wir wählen ein bestimmtes Zeichen basierend auf dieser Eingabe, wir schreiben in eine Ausgabedatei namens cowpy input file, die dann zur Ausgabe geechot wird.

Großartig. Lass uns das einbinden. Also include \{ cowpy \} from "./modules/cowpy.nf", habe ich es cowpy genannt? Ja.

Und dann rufen wir unseren neuen Prozess hier unten im main-Block des Workflows auf. Also führen wir cowpy aus. Und wir nehmen unseren neuen cowpy-Prozess und wir sagen collectGreetings.out.

Und dann, wenn du dich erinnerst, gab es zwei Ausgaben für dieses Modul. Eine namens outfile und eine namens report. Die VS Code-Erweiterung schlägt diese automatisch für uns vor und wir wollen .outfile.

Du kannst immer in diesen Prozess hier springen. Entweder du hoverst darüber und es sollte dir schnell zeigen, was die Ausgaben waren. Und wir können auch Befehl-Klick darauf machen und es wird die Moduldatei öffnen, wenn du mehr Details sehen möchtest.

Also hier sind wir. Das ist die outfile dort, und das ist der path. Das wird jetzt also die Eingabedatei für unseren cowpy-Prozess sein. Fantastisch.

Wenn du dich erinnerst, hat ein cowpy-Prozess zwei Eingaben. Wir hatten auch den value-Kanal für das Zeichen. Also können wir hier "params.character" hinzufügen. Ich hätte das hart codieren können, wenn ich wollte, aber lass es uns zu einer CLI-Option machen, sodass wir dash, dash character machen können.

Richtig. Ich muss jetzt den Eingabeparameter definieren, den wir gerade aufgerufen haben, und ihm einen Standardwert geben. Also character, String. Und ich mag den Elch, also setze ich ihn standardmäßig auf moose.

Richtig, lass uns versuchen, es auszuführen. Wenn ich also Nextflow run hello containers mache, werden wir sehen, was passiert.

Ich hätte dash resume verwenden können, wenn ich die alten work-Verzeichnisse herumliegen hätte. Und wieder wären diese ersten Prozesse gecacht worden und es wäre etwas schneller gewesen, aber es sollte im Grunde dasselbe sein.

Jetzt können wir sofort sehen, dass es einen Fehler geworfen hat, als es zu unserem neuen Prozess kam. Es sagt uns hier, dass es einen Fehler beim Ausführen des cowpy-Prozesses gab und es mit einem Exit-Status 127 beendet wurde. Das ist der Befehl, den es versucht hat auszuführen. Es sieht richtig aus, es sieht so aus, wie wir es erwartet haben. Es nimmt diesen Ausgabedateinamen, der ungefähr richtig aussieht, es führt es mit einem moose-Zeichen aus und versucht es zu speichern.

Aber du kannst den Befehlsfehler hier sehen, der sagt, cowpy-Befehl nicht gefunden. Und das macht Sinn, weil wir Nextflow noch nicht gesagt haben, einen Container zu verwenden. Wir haben ihm nur den cowpy-Befehl gegeben. Und wie ich vorhin sagte, ist cowpy nicht auf unserem lokalen System installiert. Als es also versuchte, es auszuführen, ist es fehlgeschlagen.

## 2.3.1. Einen Container für cowpy angeben

Wir müssen Nextflow sagen, dass ein Container verfügbar ist und es ihn verwenden kann. Wie machen wir das also?

Wenn wir in unser Modul hier gehen, werden wir eine neue Deklaration oben hinzufügen, die "container" heißt. Und wir werden das dann auf einen String setzen.

Wenn du dich erinnerst, kann ich in Seqera Containers diese URL kopieren und ich lasse sie einfach hier in Anführungszeichen fallen.

Jetzt gehe zurück und versuche es noch einmal auszuführen.

Mal sehen, ob es diesmal funktioniert.

Leider schlägt es auf genau dieselbe Weise fehl, obwohl wir jetzt einen Container für den Prozess definiert haben. Um unser Docker-Image zu verwenden, müssen wir Nextflow sagen, dass es die Docker-Nutzung aktivieren soll, wenn wir den Workflow ausführen.

Und wir werden das tun, indem wir eine neue Konfigurationsdatei erstellen. Ich werde also sagen touch nextflow.config.

Das ist ein spezieller Dateiname, bei dem, wenn er sich im Arbeitsverzeichnis befindet, während ich die Pipeline starte, er automatisch geladen wird. Wenn ich also in diese Nextflow.config-Datei gehe, siehst du, dass sie tatsächlich bereits existiert, was ich vergessen hatte. Und wir haben docker.enabled hier bereits drin, aber es ist auf false gesetzt, was der Standard ist.

Wenn ich das also auf equals True ändere, docker.enabled. Und es gibt Referenzdokumente für all diese Config-Scopes in den Nextflow-Dokumenten. Und auch wenn ich mit einer VS Code-Erweiterung darüber hovere, zieht es die Dokumentation spezifisch dafür ein und sagt mir, was es bedeutet und wie man es setzt.

Wir haben es jetzt also auf true gesetzt, und wenn ich Nextflow wieder ausführe, wird Nextflow jetzt wissen, dass es dieses Docker-Image für uns holen soll, falls wir es nicht bereits lokal haben, und dann diesen Prozess mit dieser Container-Umgebung ausführen.

Und wir können sehen, dass es erfolgreich ausgeführt wurde und wir haben ein kleines Häkchen neben cowpy. Fantastisch. Wenn ich hochgehe und im results-Verzeichnis nachschaue, ist die Datei noch nicht da. Und das liegt daran, dass wir diese Ausgabedatei noch veröffentlichen müssen, genau wie alle anderen.

Also gehen wir zum published-Block innerhalb des Workflows, sagen mycowpy equals cowpy.out.

Und dann hier unten im output-Block, mycowpy, geschweifte Klammern path. Hoppla. Hello containers. Mode, copy.

Wenn ich jetzt noch einmal ausführe, sollte es auf genau dieselbe Weise laufen. Ich hätte dash resume verwenden können und ich vergesse es jedes Mal. Und dann gehe ich hoch und jetzt haben wir eine neue Datei namens cowpy-COLLECTED erstellt, und da ist mein Elch, der BONJOUR, HELLO, HOLÀ sagt. Fantastisch.

Jetzt könnte ich natürlich auch jetzt "--character" übergeben. Was sind die verschiedenen Optionen? Ich glaube, es gibt einen Turkey? Also kann ich character Turkey verwenden. Es wird auf genau dieselbe Weise laufen. Ich habe eine weitere Gelegenheit verpasst, dash resume zu verwenden, und jetzt, wenn wir unsere Datei laden, haben wir jetzt einen Turkey. Fantastisch.

## 2.3.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat

Okay. Letzte kleine Sache. Lass uns diesen Befehl einfach schnell noch einmal ausführen, diesmal resume, und einen kurzen Blick in das work-Verzeichnis werfen, um zu sehen, was Nextflow unter der Haube macht, damit all das für uns funktioniert.

Diesmal ist es super schnell, lass uns in dieses work-Verzeichnis gehen, cd work/. Wenn du dich erinnerst, haben wir hier eine Reihe von Punkt-Dateien und die, an der wir in diesem Fall interessiert sind, ist die, von der ich sagte, dass wir sie fast nie ansehen müssen, genannt .command.run.

Wenn ich code dot command run mache, wird es im Editor geöffnet. Und ich kann in dieser Datei suchen und wenn ich nach unten scrolle, sollte ich Docker run sehen. Und du siehst, Nextflow macht den docker run-Befehl für uns, wenn Docker in einer Config aktiviert ist. Es hat eine ganze Reihe verschiedener Flags und Dinge hier, aber du kannst das "-v"-Flag sehen, das wir selbst verwendet haben, als wir es ausgeführt haben. Und du kannst sehen, dass es das lokale Workspace-Verzeichnis in den Container mountet, sodass der Container auf unsere Eingabedateien zugreifen und die Ausgaben speichern kann. Und dann am Ende führt es auch .command.sh aus, was das generierte Skript ist, das den cowpy-Befehl drin hat.

Und du kannst sehen, dass Nextflow die Workflow-Logik nimmt, die das Zeug ist, das uns tatsächlich interessiert, das spezifisch für unsere Analyse ist, und es macht all das clevere Zeug hinter den Kulissen, damit Docker auf unserem System funktioniert.

Und es macht das auf eine wirklich portable Weise, sodass ein Endbenutzer\*in der Pipeline die Technologie, die er\*sie verwendet, austauschen kann: Docker, Singularity, Apptainer, Conda. Das spielt für die Pipeline-Logik keine Rolle, aber Nextflow wird alle zugrunde liegenden Infrastrukturanforderungen handhaben, sodass es überall läuft.

Und das ist wirklich die Superkraft von Nextflow. Reproduzierbarkeit und Portabilität. Und mit Nextflow kannst du deinen Workflow tatsächlich teilen und andere Leute können ihn auf ihren Systemen ausführen und es wird einfach funktionieren.

Das ist eine wirklich, wirklich schwierige Sache zu tun, und jetzt weißt du auch, wie du das mit deinen Workflows machst.

Okay, das war's für dieses Kapitel. Wenn du zum Ende des Kurses gehst, findest du ein Quiz über Container. Hoffentlich hat das alles Sinn gemacht. Es ist eine wirklich coole Art, mit Analysen zu arbeiten. Und wenn du neu bei Containern bist, hoffe ich, dass ich dich überzeugt habe, dass es der richtige Weg ist, und du wirst nie zurückblicken.

Aber damit mach vielleicht eine kleine Pause, und du kommst in ein paar Minuten zu mir zurück, um den letzten Teil sechs von Hello Nextflow durchzugehen, der ganz um Konfiguration geht.

Vielen Dank.
