# Teil 5: Hallo Container - Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../05_hello_containers.md) zurück.

    Die im Transkript gezeigten Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern in den Materialien.

## Willkommen

Hallo, willkommen zu Teil Fünf des Hello Nextflow Trainingskurses.

Dieses Kapitel heißt Hallo Container. Wir werden darüber sprechen, wie Nextflow mit Tools wie Docker und Singularity integriert wird, um Software-Container zu verwenden und Software für die Nutzenden deiner Pipeline bereitzustellen.

Das bedeutet, dass Leute, die deine Pipeline ausführen, nicht selbst alle verschiedenen Tools installieren müssen. Nextflow erledigt das für sie.

Container sind eine extrem leistungsstarke Technologie und entscheidend für Reproduzierbarkeit und Benutzerfreundlichkeit. Wir beginnen mit einer kurzen Einführung in Container selbst, führen einige Docker-Befehle manuell aus und integrieren dann dieselben Container in unsere Nextflow-Pipeline.

Okay. Lass uns anfangen.

Wie zuvor starten wir mit dem Laden des Trainingsmaterials. Gehe zu training.nextflow.io. Hello Nextflow, Kapitel Fünf, Hallo Container.

Ich wechsle in meine Codespaces-Umgebung und auf der linken Seite sehen wir hello containers dot nf.

Wie zuvor ist dies dasselbe Skript, mit dem wir das vorherige Kapitel vier beendet haben, es sollte also vertraut aussehen.

Wir haben unsere Befehlszeilenparameter zur Angabe der Eingabedatei und des Batch-Namens. Wir binden unsere drei Module ein, und wir haben unseren Workflow, in dem wir die drei Prozesse ausführen.

## 0. Aufwärmen: hello-containers.nf ausführen

Du kannst diesen Workflow gerne erneut ausführen und überprüfen, ob er die Ausgaben erzeugt, die du erwartest. Vorerst werde ich ihn tatsächlich schließen und ins Terminal eintauchen.

## 1. Einen Container 'manuell' verwenden

Zu Beginn dieses Kapitels werden wir eine kleine Wiederholung der Container-Technologie machen. Wenn du sehr vertraut mit Docker oder Singularity oder anderen Container-Technologien bist, dann betrachte dies als Auffrischung oder überspringe es gerne komplett.

Nextflow unterstützt viele verschiedene Arten von Container-Technologien. Dazu gehören Docker, Singularity, Podman, Shifter, Charliecloud und mehr.

In diesem Training konzentrieren wir uns auf Docker. Das ist in den Code Spaces vorinstalliert und ist eine der beliebtesten Container-Technologien, besonders wenn du auf deinem eigenen Computer oder Laptop entwickelst.

Wenn du in einer akademischen Umgebung auf einem gemeinsam genutzten HPC arbeitest, wirst du vielleicht feststellen, dass Singularity verfügbar ist und nicht Docker. Das ist in Ordnung. Alle Konzepte sind genau gleich. Ein paar der manuellen Befehle sind unterschiedlich, aber wenn du Docker verstehst, verstehst du auch Singularity.

Tatsächlich ist Singularity auch in der Code Spaces-Umgebung installiert. Wenn du möchtest, kannst du also versuchen, dieselben Aufgaben mit Singularity anstelle von Docker zu erledigen.

Okay, was ist also Container-Technologie? Die Idee hinter Docker ist, dass es ein Image aus einer entfernten Quelle holen kann. Es auf deine lokale Maschine herunterladen und dann einen Container basierend auf diesem Image erstellen kann.

Dieser laufende Container ist ein bisschen wie eine virtuelle Maschine, die auf deinem Computer läuft. Er ist von deiner Umgebung isoliert und kommt vorpaketiert mit einem Betriebssystem und einer Reihe verfügbarer Software.

## 1.1. Das Container-Image herunterladen

Die Syntax, die wir benötigen, um ein bereits existierendes Image zu holen, ist "docker pull". Ich werde das also in mein Terminal tippen, aber jetzt brauchen wir ein Image zum Experimentieren.

Du kannst Images selbst erstellen. Du kannst sie auf öffentlichen Registries wie Docker Hub oder quay.io finden. Aber eine wirklich gute Möglichkeit, schnell Images zu bekommen, ist die Verwendung von Seqera Containers.

Dies ist ein kostenloser Community-Service, den wir 2024 erstellt haben und den du ohne Login oder ähnliches nutzen kannst.

Wenn du zu seqera.io/containers gehst oder oben auf Containers klickst, wird dir eine Suchoberfläche präsentiert und du kannst den Namen eines beliebigen Tools eingeben, das in Conda oder im Python Package Index verfügbar ist.

Standardmäßig durchsucht es die Bioconda- und Conda Forge-Kanäle, aber du kannst jedem Conda-Kanal ein Präfix voranstellen, wenn du möchtest.

Zum Spaß verwenden wir cowpy. Ich werde cowpy eingeben. Es gibt mir Ergebnisse aus dem Python Package Index und Conda Forge. Ich werde darauf klicken, um es zu meinem Container hinzuzufügen. Ich könnte hier mehrere Pakete hinzufügen, wenn ich wollte. Wähle Docker, wähle linux/amd64 und klicke auf Get Container.

Dies erstellt das Image bei Bedarf für mich, falls es noch nicht erstellt wurde, und gibt mir eine URL, die ich kopieren kann.

Wenn du interessiert bist, kannst du auf View Build Details klicken, und das führt dich zu einer Seite, die die verwendete Conda-Umgebungsdatei und das vollständige Build-Protokoll für den Build zeigt, zusammen mit den Ergebnissen des Sicherheitsscans.

Wenn ich zurück zu meinen Code Spaces gehe, kann ich jetzt diesen Container-Namen einfügen und Enter drücken.

Docker lädt jetzt alle verschiedenen Schichten innerhalb dieses Container-Images herunter und teilt uns dann mit, dass dieses Image zur Verwendung verfügbar ist.

## Ein Singularity-Image herunterladen

Wenn du Singularity verwendest, ist der Prozess im Grunde derselbe. Wir wählen unsere Image-Pakete aus, wählen cowpy. Jetzt wählen wir Singularity anstelle von Docker und klicken auf Get Container. Das gibt uns eine Image-URL mit oras://. Oder wenn du möchtest, kannst du https:// verwenden, indem du dieses Kästchen ankreuzt. Kopiere diese URL. Jetzt gehe zu Code Spaces. Wir haben tatsächlich Apptainer in diesem Space installiert, was dasselbe ist wie Singularity, aber sie sind aneinander aliasiert. Also werde ich apptainer pull machen und dann werde ich es cowpy sif nennen, aber du kannst es nennen, wie du willst. Füge die URL ein. Und das wird dieses Image für mich herunterladen.

Ich könnte ls -lh machen und cowpy.sif sehen

Singularity unterscheidet sich von Docker darin, dass Singularity alle Images in flachen Dateien speichert, während Docker eine Registry hat, in der es alle Schichten separat auf deinem Host-Rechner aufbewahrt, und es hat einen laufenden Daemon, um das alles zu verfolgen.

## 1.2. Den Container verwenden, um cowpy als einmaligen Befehl auszuführen

Okay, lass uns zurück zu Docker gehen. Wir können jetzt versuchen, dieses Image, das wir erstellt haben, auszuführen, indem wir docker run machen.

Ich werde dash dash rm machen, was nur eine einmalige Ausführung des Images durchführt. Und ich werde die Image-URL einfügen. Und dann beendest du dies mit einem Befehl, den du ausführen möchtest.

Das Image, das wir generiert haben, hatte cowpy installiert, also lass uns cowpy versuchen.

Da haben wir es. Es hat unseren Befehl ausgeführt. Ich habe cowpy nicht lokal installiert. Du kannst sehen, wenn ich versuche es auszuführen, existiert es nicht. Jedoch habe ich es in diesem Befehl mit Docker ausgeführt und es hat korrekt diese Ausgabe erzeugt.

## 1.3. Den Container verwenden, um cowpy interaktiv auszuführen

Wir können noch weiter gehen, wenn wir möchten, und einen Container interaktiv starten und darin herumsehen. Wieder mache ich "docker run dash dash rm". Jetzt werde ich dash it machen, was Docker sagt, dass wir ein interaktives Terminal wollen. Ich mache wieder die Image-URL, und dieses Mal, anstatt cowpy zu machen, werde ich bin bash machen, weil der Befehl, den wir ausführen wollen, bash ist.

Das bringt uns in diesen laufenden Container und du kannst sehen, dass sich die Eingabeaufforderung jetzt geändert hat.

Wenn ich LS slash mache, kannst du sehen, dass die Verzeichnisse hier anders sind.

Wenn ich hier rechts ein zweites Terminal öffne, das einfach in GitHub Code Spaces läuft, und LS slash mache, siehst du, dass wir Verzeichnisse wie workspaces und temp haben, während es drüben in Docker anders ist.

Diese Umgebung ist also innerhalb von Docker vollständig getrennt und von meiner Host-Umgebung isoliert. Das ist eine gute Sache, weil das die Ausführung dieses Befehls im Docker-Image isoliert und es zwischen verschiedenen Personen auf verschiedenen Host-Systemen reproduzierbar hält.

Wenn du Daten von deinem Host-System innerhalb des Docker-Images verwenden möchtest, musst du das explizit in den Container mounten.

Wir werden das gleich machen.

## 1.3.2. Die gewünschten Tool-Befehle ausführen

Aber zuerst, lass uns sehen, ob wir cowpy ausführen können. Da nochmal, der Befehl ist jetzt direkt auf der Befehlszeile verfügbar, und wir können anfangen, komplexere Dinge zu tun und Argumente zu übergeben. Hello containers und anstelle der Kuh, lass uns den Tux-Pinguin nehmen. Mal sehen, was wir noch haben.

Lass uns Käse machen. Wunderbar. Wie wäre es mit Drache und Kuh? Ziemlich gut.

## 1.3.3. Den Container verlassen

Okay. Ich kann nicht viel mehr machen, weil ich keine Daten in diesem Container habe. Also lass uns aus diesem laufenden Image aussteigen und sehen, ob wir einige Daten in den Container mounten können. Ich kann das machen, indem ich Strg D mache oder exit tippe. Okay, ich bin jetzt zurück in meinem regulären GitHub Code Space.

## 1.3.4. Daten in den Container mounten

Um einige Daten in den Docker-Container zu mounten, muss ich dash V verwenden. Also werde ich meinen vorherigen Docker-Befehl nehmen, zum Anfang zurückgehen und dash v machen. Ich werde "." für das aktuelle lokale Arbeitsverzeichnis machen, und dann einen Doppelpunkt, um zu sagen, wo das im Host-Verzeichnis gemountet werden soll, und mache slash data. Das mountet also dieses bestimmte Verzeichnis in den Container bei slash data.

Wenn ich jetzt LS slash mache, können wir sehen, dass wir ein neues Verzeichnis namens data haben, und wenn ich LS data mache, kannst du alle Dateien sehen, die wir hier in der Seitenleiste haben. Fantastisch.

## 1.3.5. Die gemounteten Daten verwenden

Jetzt können wir anfangen, einige der Dateien zu verwenden, die auf dem Host-System sind, innerhalb des Docker-Images. Also kann ich cat data greetings csv sagen. Wenn du dich erinnerst, das ist unsere CSV-Datei mit unseren verschiedenen Grüßen von vorhin, und ich kann das zu cowpy pipen. Fantastisch. Jetzt kommen wir irgendwohin.

Okay. Das reicht für die interaktive Ausführung von Docker. Hoffentlich hast du jetzt ein Gefühl dafür, was Docker ungefähr ist und wie man es sowohl verwendet, um einen Befehl einmalig auszuführen, als auch ein Image interaktiv zu verwenden. Wenn du Singularity verwendest, sind die Befehle alle sehr ähnlich, außer dass du Dinge wie apptainer exec oder apptainer run oder singularity exec oder singularity run machst.

## 2. Container in Nextflow verwenden

Als Nächstes gehen wir zurück zu unserem Nextflow-Workflow und sehen, wie man diese Technologie innerhalb der Nextflow-Pipeline verwendet.

Lass uns das Terminal schließen und Hello Containers wieder öffnen.

## 2.1. Ein cowpy-Modul schreiben

Um bei unserem cowpy-Beispiel zu bleiben, lass uns einen neuen Prozess in unserem Workflow erstellen, der cowpy verwendet. Lass uns zu Modulen gehen, eine neue Datei erstellen und sie cowpy nf nennen. Ich werde jetzt ein bisschen schummeln und den Code aus dem Trainingsmaterial kopieren und auf Speichern drücken. Und schauen wir uns das an.

Das ist also ein einfacher Prozess. Hoffentlich verstehst du jetzt, wie die Bausteine eines Prozesses aussehen. Wir haben wieder unser publishDir, das zu results geht. Wir haben zwei Eingaben, eine Eingabedatei und einen String namens character. Wir haben eine Ausgabe cowpy input file, und wir haben ein Skript, das genau so aussieht, wie das, was wir vor einer Sekunde manuell in unserem Docker-Image ausgeführt haben: cat, um eine Datei zu drucken, das auf cowpy zu pipen, zu sagen, welche Art von cowpy-Charakter wir verwenden wollen, und das in die Ausgabedatei auszugeben, die wir hier als Ausgabe übergeben.

## 2.2. cowpy zum Workflow hinzufügen

Okay, lass uns zurück zu unserem Workflow gehen, diesen neuen Prozess importieren. Also cowpy from modules cowpy nf. Lass uns einen neuen Parameter erstellen, damit wir angeben können, welchen Charakter wir wollten. Sagen wir standardmäßig Turkey. Und dann lass uns diesen neuen Prozess am Ende des Workflows aufrufen,

cowpy. Und lass uns die Ausgabe hier von Collect Greetings verwenden. Also collect greetings out, out file hier. Und dann brauchen wir ein zweites Argument, das die neuen params sind, die wir gerade gemacht haben. params dot character.

## 2.2.4. Den Workflow ausführen, um zu überprüfen, dass er funktioniert

Okay, lass uns sehen, ob unser neuer Prozess funktioniert. Nextflow run hello containers. Das sollte diese ersten drei Prozesse ausführen und dann versuchen, cowpy am Ende auszuführen.

Wir haben einen Fehler. Was es hier sagt, cowpy hatte einen Fehler und es hatte einen Exit-Status 127 und tatsächlich, Befehl sh cowpy command not found.

Wir haben Nextflow nicht gesagt, dass wir ein Docker-Image für cowpy verfügbar haben, also hat es versucht, es auf unserem Host-System auszuführen, und wir haben cowpy nicht auf unserem Host-System installiert, also hat es einen Fehler ausgelöst.

## 2.3. Einen Container verwenden, um es auszuführen

Was wir also tun müssen, ist Nextflow zu sagen, dass wir einen Container verfügbar haben. Lass uns zu unserem cowpy-Prozess gehen und wir werden eine neue Direktive oben im Prozess hinzufügen, die container heißt.

Wir finden dann unser Image, kopieren die URL und setzen das in einen String.

Das reicht allein nicht aus, weil eine Nextflow-Pipeline mehrere Möglichkeiten haben kann, Software anzugeben. Ich könnte auch conda conda-forge cowpy machen, zum Beispiel. Und Nextflow muss wissen, welche dieser Technologien du verwenden möchtest.

## 2.3.2. Die Verwendung von Docker über die nextflow.config-Datei aktivieren

Um also mit aktiviertem Docker auszuführen, werden wir uns etwas vorgreifen und die Nextflow-Config-Datei verwenden, was etwas ist, das wir im nächsten Kapitel ausführlicher behandeln werden. Du kannst in diesem Verzeichnis sehen, dass wir eine Datei namens Nextflow Config haben, und hier hast du bereits docker.enabled False.

Wir werden das auf True ändern, um Docker zu aktivieren, und dann können wir versuchen, den Workflow erneut auszuführen.

## 2.3.3. Den Workflow mit aktiviertem Docker ausführen

Nextflow run hello containers nf und dieses Mal wurde cowpy erfolgreich ausgeführt. Lass uns in Results nachsehen. cowpy collected test und da ist unser Truthahn. Wunderbar.

Im Hintergrund dort wusste Nextflow also, dass es einen Container für diesen Prozess verfügbar hatte.

Es hat das Image geholt und die Befehle für uns ausgeführt.

## 2.3.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat

Wenn du neugierig bist, können wir tatsächlich genau sehen, was es gemacht hat, indem wir in das Arbeitsverzeichnis schauen. Wenn ich code work mache, und dann den Hash und dann command run, was, wenn du dich erinnerst, die tatsächliche Datei ist, die für diese Aufgabe ausgeführt wird, können wir hineingehen und nach einer Funktion namens NXF launch suchen. Und hier kannst du den genauen Docker-Befehl sehen, den Nextflow verwendet hat, der viel wie das aussieht, was wir manuell im Terminal vorhin gemacht haben. Docker run. Dieses Host-Verzeichnis in den Container binden, und dann die Container-URL angeben.

Es gibt hier also keine Magie. Es ist nur so, dass Nextflow automatisch die schwere Arbeit für dich erledigt auf eine Weise, die bedeutet, dass du einfach Container in deiner Pipeline angeben kannst, die dann für jeden anderen, der deinen Workflow ausführt, leicht verfügbar sind. Und diese Leute müssen nicht mehr darüber nachdenken, Software zu verwalten, um deine Analyse-Pipeline auszuführen.

Sehr, sehr einfach, sehr bequem und auch wirklich reproduzierbar. Rundum gut.

Okay, großartige Arbeit. Das ist das Ende von Kapitel Fünf. Schließ dich uns im nächsten Video für Teil sechs an, das der letzte Teil dieses Hello Nextflow Trainings ist, wo wir ausführlicher über Nextflow-Konfiguration sprechen werden.

Wir sehen uns im nächsten Video.

[Nächstes Video-Transkript :octicons-arrow-right-24:](06_hello_config.md)
