# Orientierung - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtiger Hinweis"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../00_orientation.md) zurück.

## Willkommen

Hallo und willkommen bei Hello Nextflow. Mein Name ist Phil Ewels. Ich bin Product Manager für Open Source Software bei Seqera, dem Unternehmen hinter Nextflow.

Dieser Kurs ist eine praktische Einführung in die Entwicklung von Workflows mit Nextflow. Er ist für Leute konzipiert, die komplett neu bei Nextflow sind und ihre eigenen Pipelines entwickeln möchten.

Die Beispiele sind alle einfache Textverarbeitung, damit du dich auf die Nextflow-Konzepte konzentrieren kannst, ohne Fachwissen zu benötigen – nur etwas Vertrautheit mit der Kommandozeile.

Wir gehen die Grundlagen von Nextflow durch: Prozesse schreiben, sie zu mehrstufigen Workflows verbinden, Software-Abhängigkeiten mit Containern verwalten und Pipelines für verschiedene Rechenumgebungen konfigurieren. Am Ende hast du eine funktionierende Pipeline von Grund auf erstellt.

Dieser Kurs konzentriert sich auf die _Entwicklung_ von Pipelines. Wenn du nur bestehende Pipelines _ausführen_ möchtest, ohne zu tief in den Code einzutauchen, haben wir einen kürzeren "Nextflow Run"-Kurs, der besser für dich geeignet sein könnte.

Sobald du hier die Grundlagen beherrschst, haben wir auch Folgekurse, die diese Konzepte auf echte wissenschaftliche Analysen anwenden. Wir zeigen dir, wie du die Pipelines und Best Practices der nf-core-Community nutzt.

Wenn du nicht weiterkommst, geh zu community.seqera.io. Dort gibt es ein aktives Community-Forum mit einem Bereich speziell für Trainingsfragen. Du kannst es jederzeit nutzen, aber wir veranstalten auch vierteljährlich Trainingswochen mit Leuten, die speziell zum Helfen da sind. Wenn du also während einer dieser Wochen am Training teilnimmst, sei auf keinen Fall schüchtern und frag nach Hilfe.

Du kannst auch versuchen, Seqera AI um Hilfe zu bitten. Es ist großartig darin, Nextflow-Code zu erklären und dir beim Debuggen zu helfen.

Wenn du bereit bist, Nextflow im großen Maßstab auszuführen, ist Seqera Platform der beste Ort dafür. Es läuft auf deiner Infrastruktur ohne Vendor Lock-in, mit allem von Pipeline-Starts über Echtzeit-Monitoring bis hin zu interaktiven Analyseumgebungen. Aber für jetzt konzentrieren wir uns einfach auf die Grundlagen.

Gut, lass uns anfangen.

## training.nextflow.io

Okay. Das Erste, was zu beachten ist: Alle Trainingskurse auf training.nextflow.io sind sehr interaktiv. Die Idee ist, dass du dem Trainingsmaterial und meinen Anweisungen folgst, und wir gehen das Trainingsmaterial gemeinsam durch. Du brauchst also zwei Dinge: deinen Laptop und diese Website geöffnet. Und das ist so ziemlich alles.

Das ist die Homepage, wie sie heute aussieht, wenn ich das aufnehme. Du kannst sehen, dass es eine Übersicht über die verschiedenen Dinge gibt, den Hintergrund und die verschiedenen Kurse, die wir haben, und die Liste wächst ständig.

Nextflow for newcomers ist, wo wir sind. Es gibt zwei Kurse hier: Nextflow Run, das ist ein anderer Kurs, und Hello Nextflow, um den es uns geht.

Du kannst auch alle verschiedenen Kurse in der Seitenleiste sehen. Ich kann zu Hello Nextflow springen, und wir können alle verschiedenen Kapitel sehen, die wir zusammen durcharbeiten werden.

Es gibt ein paar andere wichtige Dinge zu beachten. Erstens ist das Trainingsmaterial versioniert, du kannst das hier oben sehen. Es steht 3.0 latest, was zum Zeitpunkt meiner Aufnahme die neueste stabile Version ist. Das wird sich im Laufe der Zeit ändern. Wir bringen neue Kurse heraus und aktualisieren das Material im Laufe der Zeit. Wenn es also 3.1 oder 3.2 ist, mach dir nicht zu viele Sorgen. Wenn es 4.0 ist, gibt es wahrscheinlich ein neues Video, und du solltest vielleicht danach suchen, weil es wahrscheinlich bedeutende Updates geben wird.

Ein weiteres Dropdown oben ist dieses für die Sprache. Das ist brandneu für Version 3.0. Wir haben das zuvor übersetzte Material genommen, das von Menschen von Hand gemacht wurde, und haben es in ein LLM eingespeist und diese ganze neue Infrastruktur für die Pflege verschiedener Übersetzungen des Trainingsmaterials mit LLM-Übersetzung eingerichtet.

Jetzt haben wir also all diese fantastischen Übersetzungen hier. Wenn du also auf Koreanisch folgen möchtest, kannst du die ganze Website auf Koreanisch laden und dort mitmachen. Dasselbe gilt für all diese anderen Sprachen, Hindi und Deutsch und so weiter. Ich werde auf Englisch folgen. Das ist die primäre Sprache, in der wir das Material schreiben.

Ein paar andere Buttons: Wenn du den Light Mode magst, anstelle des Dark Mode, kannst du der Website im Light Mode folgen, hier oben.

Und dann ist alles, was wir uns ansehen, in einem einzigen GitHub-Repository, das Open Source ist, genannt nextflow-io/training. Und wenn du diesen Button irgendwann klickst, geht es zum GitHub-Repository. Darauf kommen wir gleich zurück.

## GitHub Codespaces einrichten

Okay, jetzt hast du das also im Browser-Tab geöffnet. Lass uns zu Hello Nextflow gehen und reinklicken. Du kannst auf der Intro-Seite sehen, dass es uns ein paar der Anforderungen, die Übersicht und den Lektionsplan von ungefähr dem, was wir abdecken werden, zeigt, und dann tauchen wir in die ersten Schritte ein.

Es gibt verschiedene Möglichkeiten, wie du dieses interaktive Tutorial machen kannst. Wenn du möchtest, kannst du das gerne lokal auf deinem eigenen Computer mit deiner eigenen Nextflow-Installation machen. Wenn wir zu Environment Options klicken, kannst du sehen, dass es mehr Details gibt, wie du das entweder mit lokalen Devcontainern machen kannst, oder du kannst auch einfach die gesamte Software lokal installieren, mit manueller Installation.

Wir arbeiten daran, das mit Seqera Studios gut zum Laufen zu bringen, das ist also eine weitere Option. Aber die gängigste Option im Moment ist die Verwendung von GitHub Codespaces.

Codespaces richtet eine Sandbox-Umgebung auf einem Remote-Server ein, der von GitHub betrieben wird. Und es ist kostenlos für eine bestimmte Menge an Nutzung, was normalerweise für das Training ausreicht. Und es richtet dir eine VS Code-Instanz ein, eine IDE, wo du auf alle Dateien aus dem Repository zugreifen, Nextflow ausführen und alles machen kannst. Und wir haben Codespaces für dich vorkonfiguriert. Es hat also alles, was du brauchst.

Das Schöne daran ist, dass es nur ein Klick ist, um einen Codespace einzurichten. Es ist für alle gleich, und wir wissen, dass du alle Voraussetzungen bereits installiert hast, also ist es schön und schnell.

Das Erste, was zu tun ist, ist zu "Getting Started" zu gehen. Suche nach diesem Button, der sagt _Open in Codespaces_. Ich werde Cmd + Klick machen, um es in einem neuen Tab zu öffnen, und es bringt uns zu GitHub.

So sieht es aus. Wir können sehen, wir haben alle Optionen hier für dich eingestellt. Wenn du möchtest, kannst du auf "change options" klicken. Einige Dinge, die du hier tun kannst: Du kannst zum Beispiel eine größere Instanz-Maschine wählen, falls du feststellst, dass sie abstürzt, weil ihr der Speicher ausgeht oder so etwas. Oder spezifische Versionen des Trainingsmaterials einstellen. Aber normalerweise kannst du einfach mit dem gehen, was wir hier eingerichtet haben, und du kannst es sehen. In diesem Fall verwendet es das 3.0-Release.

Also werde ich auf "create new Codespace" klicken. Und das bringt mich rein.

Beachte auch, es steht "no Codespace to resume" dort. Wenn ich zuvor einen Codespace erstellt habe, bringt mich das erneute Klicken auf diesen Button im Trainingsmaterial zur gleichen Seite und es listet alle Codespaces auf, die ich bereits laufen habe. Dann kannst du einfach direkt zurück in sie springen und dort weitermachen, wo du aufgehört hast. Es macht also nichts, wenn du deinen Laptop geschlossen hast.

Sie schalten sich automatisch nach ein paar Minuten Inaktivität ab, aber das ist kein Problem. Du kannst sie einfach neu starten.

Sobald du einen neuen Codespace startest, wird er auf dieser Seite so sitzen und für eine ganze Weile laden. Jetzt ist also ein guter Zeitpunkt für eine kurze Pause. Vielleicht hast du vergessen, auf die Toilette zu gehen, oder du möchtest eine Tasse Tee, bevor wir anfangen? Geh jetzt, während du darauf wartest, denn es wird eine Weile dort drehen.

Nur kurz, während wir warten, dass es lädt, werde ich auch zu github.com/codespaces gehen und nur zeigen, das ist die Übersichtsseite, wo du alle verschiedenen Codespaces sehen kannst, die du gerade laufen hast.

Du kannst sehen, ich habe hier einen für nextflow-io/training. Keine Änderungen, weil ich noch nichts darin gemacht habe. Die Menge an Ressourcen, die er verwendet, und du kannst sehen, im Moment richtet er sich ein. Ich kann hier hingehen, dieses kleine Dropdown klicken und auf "delete" klicken. Wenn du also versehentlich mehrere Codespaces eingerichtet hast und einige nicht verwendest, kannst du die alten löschen und aufräumen.

Schließlich noch eine weitere Möglichkeit, hier reinzukommen. Wenn wir zum GitHub-Repository gehen. Und das funktioniert für jedes GitHub-Repository. Klicke auf "code". Du hast Befehle zum lokalen Klonen des Repositories. Und es gibt einen Tab namens Codespaces. Und wieder kannst du einen neuen erstellen, und du kannst die sehen, die bereits laufen.

Wenn du also vergisst, wie du deinen Codespace erstellt hast, kannst du auf diese Weise immer wieder zurückkommen.

## Die VS Code-Oberfläche

Okay, der Build ist fertig und es fängt jetzt an, die GitHub Codespaces zu laden. Es dauert nicht immer so lange, also mach dir keine Sorgen. Es ist nur das erste Mal, wenn du den Codespace erstellst. Wenn du in einen zurückspringst, der bereits existiert, ist es viel schneller.

Sei nicht zu ungeduldig, wenn das das erste Mal ist, es ist noch nicht fertig, obwohl es anfängt, uns eine Oberfläche zu geben.

Aber während wir darauf warten, dass die letzten Dinge eingerichtet werden, werde ich dich einfach durch die Oberfläche führen, falls du mit VS Code etwas unvertraut bist.

Erstens gibt es die Chat-Seitenleiste für KI-Sachen, die wir nicht brauchen. Also werde ich die schließen, das loswerden und etwas Platz freimachen.

Links haben wir den Datei-Explorer, der uns alle Dateien im Git-Repository zeigt, das ist der Workspace, den wir erstellt haben. Beachte, das sind keine lokalen Dateien. Das ist alles auf dem Remote-Server, wo wir arbeiten. Du kannst lokale Dateien per Drag & Drop ziehen und so, aber zum größten Teil werden wir heute nicht daran denken. Wir arbeiten einfach rein remote.

Es gibt andere Tools in dieser Seitenleiste, zum Beispiel Suche. Du kannst also alle Dateien in einem Repository auf einmal durchsuchen. Und wenn wir Entwicklungsarbeit am Training-Repo machen würden, könnten wir Integration mit Quellcodeverwaltung mit Git und Debugging und anderen Dingen machen.

Andere Dinge sind, es gibt ein Haupt-Code-Bearbeitungsfenster hier oben, das gerade eine Vorschau der Readme geladen hat, die für das Trainingsmaterial ist. In diesem Fall betrachtet es also Markdown, aber normalerweise wird das ein Code-Editor sein.

Und darunter haben wir das Terminal, wo wir alle unsere Befehle ausführen und direkt mit Nextflow interagieren werden.

Alles im Codespace ist vorinstalliert, also ist der Nextflow-Befehl bereits da und so weiter und so fort.

Okay. Wenn du so weit kommst, sollte es ungefähr fertig sein. Du kannst jetzt sehen, es hat den Nextflow Language Server heruntergeladen und einige Extensions für uns in VS Code eingerichtet, einschließlich der Nextflow-Extension, die nützlich sein wird. Also kann ich das schließen und ich kann die README.md schließen.

Und jetzt kannst du sehen, ich habe etwas mehr auf der linken Seite. Ich bin hier etwas reingezoomt, aber wenn ich rauszoome, kannst du sehen, dass einer der Buttons Nextflow mit dem Nextflow-Icon sagt. Und das hat einige nette Sachen hier drin, um das Projekt zu erkunden und so, worauf wir später zurückkommen werden.

Okay. Falls du jemals eines dieser Panels verlierst, diese Buttons oben rechts sind wirklich nützlich und diese zeigen und verstecken einfach Dinge. Das zeigt und versteckt den Explorer, zeigt und versteckt das Terminal unten. Und so weiter.

Ich werde diese ziemlich oft verwenden, weil ich sehr reingezoomt bin, also versuche ich, dir zu helfen, den ganzen Text auf meinem Bildschirm zu sehen, und es ist also nützlich, mit dem Terminal Vollbild gehen zu können und es dann zu verstecken, wenn wir uns Code ansehen. Aber die meiste Zeit kannst du einfach all dieses Zeug gleichzeitig offen haben.

Okay, was gibt es noch zu sehen? Nicht viel mehr. Beachte, dass Nextflow, wie ich sage, installiert ist. Also kann ich "nextflow -version" eingeben und es sollte anzeigen, welche Version wir installiert haben.

Es gibt auch andere Sachen, die hier installiert sind. Am Ende jedes Kapitels haben wir zum Beispiel eine Reihe von Quizfragen auf der Website. Und du kannst die auch im Terminal machen, wenn du möchtest, indem du "quiz" eingibst.

Es gibt einige andere Tastaturkürzel, die ich verwenden werde, nur falls du neugierig bist. Zum Beispiel habe ich gerade Cmd+K auf meinem Mac gedrückt und das hat das Terminal geleert, um die gesamte vorherige Ausgabe loszuwerden. Das ist also schön, um die Dinge sauber zu halten. Wenn du siehst, dass ich das mache, so mache ich es.

Und auch, wenn du neu im Terminal bist, denk daran, dass du Tab zur automatischen Vervollständigung verwenden kannst, was ich viel machen werde, um Pfade automatisch zu vervollständigen.

Ich kann also auf der linken Seite hier sehen, es gibt einen Ordner namens Hello Nextflow, den wir durcharbeiten werden. Wenn ich "ls" mache, um Dateien aufzulisten, kann ich "hel" machen, Tab drücken, vervollständigt automatisch. Und das ist also ein sehr schneller Weg, um Pfade zu vervollständigen.

## Nur den Hello Nextflow-Ordner öffnen

Okay. Das ist großartig. Es gibt allerdings viele Sachen in diesem Repository.

Es gibt alle Dateien zum Generieren der Website, und es gibt mehrere verschiedene Kurse hier drin, und du kannst es von dieser Route aus machen und einfach in den "Hello Nextflow"-Ordner klicken. Aber es ist schön, sich wirklich nur darauf zu konzentrieren.

Du kannst das als deinen Workspace mit ein bisschen Herumklicken hier setzen und ein Projektverzeichnis einstellen und so. Aber der einfachste Weg ist, "code" einzugeben, das ist der CLI-Befehl zum Starten von VS Code, und dann "hello-nextflow".

Das öffnet einen neuen Browser-Tab und du kannst den alten schließen. Und es sieht genau gleich aus. Aber jetzt kannst du sehen, wir sind in diesem Unterverzeichnis und alle anderen Dateien sind unsichtbar, und wir haben ein saubereres Setup.

Du kannst hier sehen, dass auch das aktuelle Arbeitsverzeichnis jetzt im Hello Nextflow-Ordner ist. Also schön und sauber. Wir müssen uns keine Sorgen machen, am falschen Ort zu sein. Okay.

## Neue Nextflow-Syntax für 2026

Es gibt eine besondere Sache, die ich an dieser Stelle erwähnen muss. Gerade jetzt, Anfang 2026, fangen wir an, verschiedene Features in Nextflow einzubringen, und eines der großen neuen ist ein neuer Sprachsyntax-Parser innerhalb von Nextflow.

Grundsätzlich die Engine, die deine Nextflow-Dateien liest und das zur Laufzeit versteht. Es gibt einige Änderungen an der Syntax, und es ist wirklich wichtig, dass du Nextflow mit dem richtigen aktivierten Syntax-Parser verwendest.

Zwei Dinge brauchst du dafür. Du brauchst eine aktuelle Version von Nextflow und du musst sicherstellen, dass es aktiviert ist.

Wenn ich "nextflow -version" nochmal mache, siehst du, dass der Codespace mit 25.10.2 läuft und 25.10 ist die Mindestversion, um dieses Zeug verwenden zu können.

Wenn du 26.04 verwendest, was für mich noch nicht rausgekommen ist, aber bald wird. Dann wird das den neuen Syntax-Parser standardmäßig ausführen, und du musst nichts anderes tun.

Aber wenn du 25.10 ausführst, musst du den strict syntax parser aktivieren, wie er genannt wird, oder v2 syntax parser.

Das wird mit einer Umgebungsvariable gemacht. Es ist bereits im Codespace gesetzt, also musst du nichts tun. Aber wenn du lokal ausführst, musst du das setzen, und ich kann das verifizieren, indem ich "echo $NXF_SYNTAX_PARSER" mache, und es sollte auf v2 gesetzt sein.

Wenn du also lokal ausführst, mach einfach "export NXF_SYNTAX_PARSER=v2". So einfach ist das. Aber denk daran, das zu tun, denn sonst wirst du einige seltsame Diskrepanzen und Fehler sehen, während wir weitermachen.

Wenn du dir bei irgendetwas von diesem Zeug rund um Nextflow-Version und Syntax-Parser unsicher bist, erstens, denk daran, du musst dir keine Sorgen machen, wenn du in Codespaces bist. Alles sollte richtig eingerichtet sein. Aber zweitens, wenn du zum Nextflow-Trainingsmaterial gehst, wenn du runtergehst, spricht es über Versionsanforderungen, es gibt einen Link hier, der dich runter zur Hilfeseite rund um Explore Versions bringt, und das geht im Detail durch alles.

Es lohnt sich, das zu lesen, wenn du einen Moment hast, denn es hilft zu klären, was einige der verschiedenen Begriffe sind, die du hören könntest, wenn du anfängst, Nextflow zu verwenden. Dinge wie DSL1, DSL2, Syntax Parser eins, Syntax Parser zwei und so weiter. Es lohnt sich also, einfach mal einen Blick darauf zu werfen, und das wiederholt einiges von dem, was ich gerade gesagt habe.

Es ist auch wirklich nützlich, wenn du zuvor Nextflow-Code geschrieben hast und für eine Auffrischung zurückkommst. Es sagt dir einige der Dinge, die sich ändern, und verlinkt dich zu Teilen der Nextflow-Dokumentation, die dir sagen, wie du deinen Nextflow-Code aktualisierst.

## Kursdateien

Okay. Letzte Sache, um uns vertraut zu machen, ist einfach die Dateien zu sehen, die in diesem Verzeichnis sind. Du kannst entweder in der Seitenleiste schauen oder oft im Trainingsmaterial verwenden wir den tree-Befehl, -L, das ist die Anzahl der Ebenen, in die geschaut werden soll. Wir sagen zwei, und wenn ich das Vollbild mache, siehst du, das spiegelt genau wider, was wir in der Seitenleiste dort sehen, aber es schließt versteckte Dateien aus, die mit einem Punkt beginnen.

Die \*.nf-Dateien, das steht für Nextflow. Das sind also die Nextflow-Skriptdateien, und es gibt eine Starter-Datei hier für jedes der verschiedenen Kapitel des Trainingsmaterials, die wir öffnen und erkunden und dann bearbeiten werden.

Wir werden diese Dateien ändern, während wir weitermachen, und am Ende jedes Kapitels sollten die Dateien also ziemlich gleich aussehen wie am Anfang des Kapitels für das nächste. Aber wir geben dir diese verschiedenen Dateien, damit du immer irgendwie frisch anfangen kannst und dir nicht zu viele Sorgen machen musst, die Syntax zu vermasseln.

Wenn du mit etwas vergleichen musst, das definitiv funktionieren sollte, kannst du im solutions-Ordner nachschauen, und das ist wie ein Endzustand für jedes der Kapitel, also kannst du vergleichen, was du geschrieben hast, mit dem, was dort ist.

Es gibt ein data-Verzeichnis. Das hat nur eine greetings.csv-Datei, die wir als Beispiel-Eingabedaten in einem Teil des Kurses verwenden werden, und Dinge wie eine Config-Datei und einige Parameter, die wir später im Kurs beschreiben werden.

## Abschluss

Okay, jetzt läuft hoffentlich alles. Dein Bildschirm sieht genauso aus wie meiner und du verstehst, wie du an alles rankommst und was all die verschiedenen Dateien sind.

Wenn du zum Ende der Seite bei Getting Started scrollst, kleine Checkbox, solltest du sagen, dass ich verstehe, was ich tue. Meine Umgebung ist eingerichtet und läuft und du hast dein Arbeitsverzeichnis richtig auf den "Hello Nextflow"-Ordner gesetzt.

Wenn du all das abgehakt hast und sie grün aussehen, können wir zum nächsten Video und zum nächsten Kapitel weitermachen, das ist Teil eins. Hello World. Bis gleich.
