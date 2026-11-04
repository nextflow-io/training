# Teil 6: Hello Config - Video-Transkript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../06_hello_config.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur zur Orientierung und enthalten möglicherweise nicht alle Abschnittsnummern aus den Materialien.

## Willkommen

Hi und willkommen zurück zu Teil 6 von Hello Nextflow. Dieser Abschnitt dreht sich komplett um Configs, und es ist der letzte Teil dieses Kurses.

Nextflow ist besonders gut in zwei Dingen: Reproduzierbarkeit und Portabilität. Bei Configs sehen wir, wie die zweite davon richtig glänzt. Die Fähigkeit, eine Nextflow-Pipeline so zu konfigurieren, dass sie auf verschiedene Arten läuft und auf unterschiedlichen Systemen funktioniert, ohne den zugrunde liegenden Pipeline-Code bearbeiten zu müssen.

Diese Superkraft ermöglicht es, dass Nextflow-Pipelines von anderen Leuten an anderen Orten wiederverwendet werden können, oder über verschiedene Infrastrukturen hinweg, auf die du selbst Zugriff hast.

Das bedeutet, du kannst Pipeline-Code auf deinem Laptop entwickeln, ihn in die Cloud pushen, auf deinem HPC ausführen, und es ist derselbe Pipeline-Code, der überall läuft.

In diesem Abschnitt werden wir ein paar Themen durchgehen. Wir beginnen damit, wie Nextflow mit Config-Dateien umgeht, wo er sie lädt, wie du sie schreibst und strukturierst, und diese Trennung zwischen der Pipeline selbst und dem, was in eine Config-Datei gehört.

Dann gehen wir zu einigen häufigen Anwendungsfällen über, wie z.B. das Ändern des Speicherorts für Ausgabedateien, und auch wie man die Pipeline auf verschiedenen Infrastrukturen zum Laufen bringt, sowohl mit verschiedenen Arten von Software-Packaging als auch beim Einreichen von Jobs an verschiedene Infrastrukturen.

## Config-Dateihierarchien

Okay, lass uns anfangen. Wenn es um das Laden von Config-Dateien geht, kann Nextflow aus vielen verschiedenen Orten laden, was eine gute Sache ist und auch eine etwas riskante Sache sein kann, weil es manchmal etwas schwierig sein kann zu wissen, woher es eine Config-Datei bekommt und in welcher Reihenfolge es Dinge lädt.

Daher empfehle ich dir wirklich, auf diesen Link hier zu klicken, der uns zu den Nextflow-Docs führt. Auf dieser Konfigurationsseite werden die wichtigsten Orte aufgelistet, von denen Config geladen wird, und wichtig ist die Reihenfolge, in der diese Dinge geladen werden.

Du kannst sehen, du kannst eine Config-Datei in dein Nextflow-Home-Verzeichnis legen, was normalerweise ".nextflow" in deinem Home-Verzeichnis ist. Und diese Datei wird immer von jedem Nextflow-Lauf auf deinem System geladen.

Der nächste Ort ist eine Datei im Root deines Pipeline-Repositorys oder -Verzeichnisses namens "nextflow.config".

Danach eine weitere Datei namens "nextflow.config", aber diesmal im Verzeichnis, von dem aus du Nextflow startest: dem Launch-Verzeichnis.

Schließlich kannst du Config-Dateipfade auf der Kommandozeile mit einem "-c"-Argument angeben, und du kannst das mehrmals machen. Sie werden in der Reihenfolge angewendet, in der du sie angibst.

Du kannst Config-Dateien an all diesen Orten bereitstellen, wenn du möchtest, und sie werden iterativ geladen, wobei jede die vorherige nur in den Config-Bereichen überschreibt, wo sie kollidieren.

Das ist ein wirklich mächtiges System, weil es bedeutet, dass du sinnvolle Standards setzen und dann schrittweise immer spezifischer werden kannst, während du dich auf diese Config eingrenzst.

## 0. Aufwärmen: hello-config.nf ausführen

Okay, lass uns das schließen und in unsere Codespaces springen und loslegen. Wie zuvor habe ich hier aufgeräumt, ich habe meine vorherigen Results-Verzeichnisse entfernt, meine Nextflow- und Work-Verzeichnisse und so weiter. Mach dir keine Sorgen, wenn du diese Dateien noch herumliegen hast. Ich bin nur sehr stark eingezoomt und daher werden die Dinge sonst sehr schnell unübersichtlich.

Wir werden mit hello-config.nf arbeiten, der letzten Datei in unserem Verzeichnis, und dies sollte dort anknüpfen, wo wir im vorherigen Abschnitt aufgehört haben.

Wir haben also unsere vier verschiedenen Prozesse, die aus Moduldateien eingebunden werden. Wir haben unsere Pipeline-Parameter, unseren Workflow-Block, wo wir die verschiedenen Prozesse aufrufen und die Kanäle zusammenfügen, die Ausgabekanäle veröffentlichen, und dann den Output-Block ganz unten, wo wir definieren, wo diese Dateien gespeichert werden sollen und wie sie kopiert werden sollen.

Wir haben auch bereits eine "nextflow.config"-Datei aus dem letzten Kapitel, wo wir Docker aktivieren, und wir werden heute in diese Datei aufbauen.

Wie zuvor haben wir den Ausgabepfad in diesem Hauptskript auf hello config geändert, damit er nicht mit früheren Ergebnissen kollidiert, die du generiert hast.

Okay, lass uns kurz überprüfen, ob noch alles wie erwartet funktioniert. Ich öffne ein Terminal und mache nextflow run hello-config.nf. Nextflow lädt hoch. Sollte unsere vier verschiedenen Prozesse ausführen. Etwas schöne ASCII-Kunst mit cowpy generieren und dann unsere Ergebnisse in unseren Results-Dateien in diesem Verzeichnis speichern.

Ich kann kurz hier reinschauen, nur um sicherzustellen, dass diese Dateien so aussehen, wie wir erwarten, und sicher genug, da ist unser riesiger Turkey. Super.

## 1.1. Standardwerte in nextflow.config verschieben

Das Erste, was wir jetzt machen werden, ist, einige Dinge aus unserem Skript in unsere Config-Datei zu verschieben.

Und was uns interessiert, sind hauptsächlich die Parameter zu diesem Zeitpunkt. Wir wollen die Standardwerte in die Config-Datei bringen, damit klarer ist, was die Standards sind, und es für Leute einfacher ist, sie zu überschreiben.

Ich werde diesen Params-Block hier aus dem Skript nehmen und ihn in die Config-Datei einfügen. Und wir müssen hier etwas vorsichtig sein, weil die Syntax zwischen Config und Skripten momentan etwas unterschiedlich ist. Die Config-Datei kann keine Typdeklarationen annehmen, weil wir diese Params nicht wirklich definieren, wir referenzieren sie nur. Also werde ich diese loswerden.

Aber ansonsten ist es sehr ähnlich. Wir haben einen Params-Block und dann haben wir unsere verschiedenen Eingabeparameter, Batch-Parameter, Character-Parameter.

Ich kann jetzt zurück zu meinem Skript gehen und ich muss diese Standards nicht mehr definieren, weil diese Werte jetzt in meiner Nextflow-Config-Datei sind.

Allerdings lasse ich die Parameternamen und ihre Typen stehen, damit Nextflow diese Information kennt und weiterhin die gesamte Typsicherheit und alles durchführen kann.

Okay. Wir speichern diese Dateien und überprüfen schnell, ob noch alles genauso funktioniert wie vorher. Es sollte keine Änderungen geben. Wir haben die Werte gleich gelassen. Wir haben nur verschoben, wo sie definiert wurden.

Super.

## 1.2. Eine laufspezifische Konfigurationsdatei verwenden

Bisher haben wir Nextflow aus demselben Verzeichnis gestartet, in dem wir unser Pipeline-Skript haben. Also unser Launch-Verzeichnis und unser Pipeline-Verzeichnis sind quasi dasselbe.

Um zu zeigen, wie wir verschiedene Config-Dateien mit verschiedenen Launch-Verzeichnissen haben können, werden wir jetzt ein neues Unterverzeichnis erstellen.

Ich sage also mkdir, und wir werden es tux-run nennen.

Und dann werde ich cd machen, ins Verzeichnis tux-run wechseln. Und beachte, dass unser Arbeitsverzeichnis jetzt nicht mehr im selben Verzeichnis wie die Pipeline-Skripte ist.

Okay, lass uns eine neue "nextflow.config"-Datei erstellen. Also touch Nextflow config, und öffnen wir sie einfach in VS Code. Du kannst auch in der Seitenleiste hier sehen, dass wir jetzt in diesem Unterverzeichnis sind.

Jetzt können wir denselben Params-Block nehmen, den wir in der obersten nextflow.config hatten, das hier rüberkopieren und jetzt können wir diese Werte ändern.

Erstens sind die Daten jetzt ein anderer relativer Pfad, weil wir in einem Unterverzeichnis sind, also müssen wir das aktualisieren. Und dann werden wir batch in experiment ändern und character von Turkey in tux ändern.

Jetzt speichern wir das und probieren es aus. Genau wie bei data muss ich jetzt ../ sagen, um zum Skript zu kommen. Also Hello config. Und ich drücke Enter.

Der Pipeline-Code hat sich überhaupt nicht geändert, aber jetzt werden wir zwei Config-Sets laden, und die Launch-Dir-Config-Datei sollte die Standards überschreiben, die in der Pipeline nextflow.config gesetzt wurden, und wir sollten verschiedene Ergebnissets bekommen.

Sicher genug, in unserem Verzeichnis hier, in tux-run, kannst du sehen, dass wir ein .nextflow-Verzeichnis und ein Work-Verzeichnis haben, und das liegt daran, dass diese immer in deinem Launch-Verzeichnis erstellt werden. Diese sind also unterschiedlich zu den Work- und Results-Verzeichnissen, die wir von früheren Läufen hatten.

Wenn ich jetzt in Results schaue, können wir unsere Collected sehen und da ist unser kleiner Tux-Charakter. Du kannst sehen, dass diese Parameter richtig interpretiert wurden.

## 1.3. Eine Parameterdatei verwenden

Okay. Als ich vorhin über die verschiedenen Config-Dateien gesprochen habe, die geladen werden könnten, habe ich einen anderen Ort ausgelassen, von dem wir Config bekommen können.

Du kannst sie von der Kommandozeile bekommen, wie wir mit Doppelstrich-Parameternamen gesehen haben, aber wir können auch eine YAML- oder JSON-Datei mit nur Params bereitstellen.

Die Config-Datei kann viele verschiedene Arten von Bereichen haben, aber diese Dateien enthalten nur Parameter, und es ist eine schöne benutzerfreundliche Art, viele Parameter auf einmal zu übergeben, und vielleicht eine etwas reproduzierbarere Art, weil du sie in eine Datei schreibst, also ist es einfach, sie zu einem späteren Zeitpunkt zu bekommen.

Gehen wir also zurück zu unserem Terminal und stellen sicher, dass wir ein Verzeichnis nach oben gehen, damit ich nicht mehr im Unterverzeichnis bin, und ich schaue mir die YAML-Datei an, die wir hier haben, namens test-params.yaml.

Wenn ich also einfach code test-params.yaml mache, kannst du sehen, das ist nur eine reguläre YAML-Datei. Nichts Besonderes daran. Mit den Schlüsseln, die unsere Parameternamen sind, mit der YAML-Formatierung also einem Doppelpunkt hier, und dann einem Wert.

Beachte, dass dies kein Nextflow-Code ist, also können wir hier keine Dinge wie Variablen einfügen. Das sind nur statische Werte.

Da JSON tatsächlich als YAML parst, können wir auch eine test-params.json-Datei haben, die sehr ähnlich aussieht. Es ist nur ein anderes Datenformat.

Wir haben also zwei verschiedene Testdateien hier und haben leicht unterschiedliche Variablen.

Okay, wie geben wir diese an Nextflow? Es ist sehr einfach. Wir machen Nextflow run hello config, wie zuvor. Und anstelle von "-c" für Config-Datei oder dem Laden dieser Standard-Dateinamen, machen wir -params-file. Einzelner Bindestrich, weil es eine Kern-Nextflow-Option ist.

Und dann den Pfad für diese Datei übergeben. Ich werde also "-params-file test-params.yaml" machen, und wir werden sehen, ob diese richtig geladen werden.

Okay. Es ist gelaufen. Erinnern wir uns kurz, was in dieser YAML-Datei war. Also der Batch war auf YAML gesetzt, so sollte er also heißen, und er sollte einen Stegosaurus haben. Gehen wir also hoch und schauen in Results. Und wir haben COLLECTED-yaml. Schauen wir mal, ob wir einen Stegosaurus haben. Fantastisch, ein Stegosaurus mit Hut. Das gefällt uns.

Das hat also wirklich gut funktioniert und ist genau dasselbe mit der JSON-Datei. Wir tauschen einfach die Dateierweiterung hier aus und Nextflow weiß, wie es das lesen soll.

Und in diesem Fall sollten wir einen Batch namens JSON haben und wir sollten eine Schildkröte haben. Schauen wir mal. Wunderbar. Eines meiner Lieblings-CLI-Tools.

## 2.1. Das Ausgabeverzeichnis mit -output-dir anpassen

Okay, das war hauptsächlich über Eingaben in die Pipeline und das Ändern von Parametern. Was ist mit den Ausgaben?

Obwohl wir die Unterverzeichnisse mit Params geändert haben, hast du vielleicht bemerkt, dass alle unsere Dateien immer noch nach Results gehen.

Wir können dieses Basisverzeichnis, in das alle Dateien veröffentlicht werden, mit einem Kommandozeilen-Flag namens -output-dir ändern. Wenn ich also Nextflow run hello config mache und dann -output-dir, und wir nennen es "custom-outdir-cli". Kann nicht tippen. Nur damit wir uns erinnern, woher diese Dateien kamen.

Das ist eine Kern-Nextflow-Option und eine sehr neue. Dies wurde erst kürzlich hinzugefügt, und das ist eine der Sachen, die wir mit dem neuen Language-Parser und allem machen können.

Es ist etwas umständlich zu tippen. Du kannst es auch einfach "-o" nennen, wenn du willst. Wenn ich also zurückgehe. Kann das einfach auf "-o" verkürzen, was etwas einfacher ist.

Okay. Wir führen das aus. Wir haben nichts in unserer Pipeline oder auch in unserer Config an diesem Punkt geändert, und es sollte hoffentlich alle unsere Ergebnisse in ein anderes Top-Level-Verzeichnis speichern. Und du kannst dir vorstellen, du kannst dies auf grundsätzlich jeden Pfad setzen, den du willst.

Es ist gerade oben angekommen. Wir haben ein custom-outdir-cli, und alle Dateien sind dort genau gleich organisiert, mit ihren gleichen Unterverzeichnissen und Dateinamen. Das ist also eine wirklich einfache Möglichkeit, einfach zu ändern, wo die Pipeline ihre Ergebnisse veröffentlicht, ohne zu viel darüber nachzudenken, wie diese Ergebnisse organisiert sind.

## 2.1.2. Fest kodierte Pfade aus dem Output-Block entfernen

Wenn ich in dieses Verzeichnis schaue, können wir sehen, dass wir immer noch ein Unterverzeichnis namens Hello Config haben, was sich jetzt etwas redundant anfühlt.

Laden wir also unser Skript nochmal und können jetzt dieses Unterverzeichnis aus dem Output-Block am Ende entfernen. Weil wir es nicht wirklich mehr brauchen. Also können wir das jetzt einfach machen, das hier löschen. Und dann, wenn es nur das ist, kannst du das entweder komplett löschen oder als leeren String belassen. Ich werde es vorerst als leeren String belassen, weil wir zurückkommen und einige verschiedene Dinge an seiner Stelle einfügen werden. Aber wenn dir Unterverzeichnisse egal sind, ist es am saubersten, die Path-Deklaration dort komplett zu entfernen.

Okay, speichern wir. Probieren wir es schnell nochmal aus. Ich werde tatsächlich mein "custom-outdir-cli"-Verzeichnis entfernen, damit wir nicht durch vorhandene Dateien dort verwirrt werden. Denn denk daran, wenn du Dinge veröffentlichst, werden die Dateien, die bereits da waren, nicht entfernt. Es werden nur neue hinzugefügt. Führen wir diesen Befehl nochmal aus, custom-outdir-cli.

Und wenn du jetzt "ls custom-outdir-cli" machst, gibt es kein Verzeichnis mehr, das Hello Config heißt.

## 2.2.1. outputDir in der Konfigurationsdatei setzen

Okay, das Kommandozeilen-Flag hier, "-o" oder "-output-dir" ist gut. Aber wie setzt man Standards dafür in der Config? Wie machen wir das?

Ich öffne die "nextflow.config"-Datei, schließe alles andere und werde das los. Wir können hier eine neue Config-Option hinzufügen, die ich gerade von der Trainingsmaterial-Website kopiert habe, und sie heißt outputDir.

Sie ist nicht unter irgendeinem Bereich. Nicht unter params oder so. Sie ist auf oberster Ebene, und wir können dies auf einen String setzen. Eine einfache Sache wäre, es einfach zu etwas anderem als Results als fest kodiertem String zu ändern. Aber weil dies in einer Nextflow-Config-Datei ist, können wir hier etwas cleverer sein und auch Variablen einschließen.

Und du kannst hier sehen, dass wir eine Params-Variable eingefügt haben, params.batch, die Teil dieses Strings ist. Das bedeutet, wir können Variablen wiederverwenden, die von anderen Orten kommen. Und in diesem Fall, wenn wir --batch machen, wenn wir die Nextflow-Pipeline ausführen, werden wir ein Unterverzeichnis in unserem benutzerdefinierten Pfad bekommen, basierend darauf, wie der Batch-Name war.

Okay, probieren wir das aus und schauen uns kurz an, wie die Ergebnisse aussehen. Wenn ich also Nextflow run hello config und --batch my_run mache. Erinnern wir uns, wie die Config aussah. Also es ist custom-outdir-config.

Tree custom-outdir-config. Und du kannst sehen, der Batch hieß my_run. Und dann haben wir dieses Unterverzeichnis namens my_run. Also dieser dynamische Dateipfad hat funktioniert.

Und nicht nur das, es ging nicht mehr in ein Standard-Results-Verzeichnis, und ich musste nichts auf der Kommandozeile angeben, um das Basisverzeichnis zu ändern. Wir haben also erfolgreich den Standardwert für das Standard-outputDir zurückgesetzt.

## 2.2.2. Unterverzeichnisse mit Batch- und Prozessnamen

Okay, gehen wir noch etwas weiter. Das ist eine dynamische Variable innerhalb der Config-Datei. Wie sieht es mit dem Skript selbst aus? Bisher hatten wir diese Pfade hier und diese können auch dynamisch sein. Anstatt also einfach etwas fest zu kodieren, können wir einige geschweifte Klammern einfügen und etwas Dynamisches einsetzen.

Zum Beispiel haben wir unsere Prozesse namens sayHello. Wir könnten sayHello.name machen, was ein Attribut des Prozesses ist, was irgendwie langweilig ist, weil es in diesem Fall nur "sayHello" ist. Aber es ist variabel.

Das gibt dir also eine Idee. Wir können das hier einfügen und convertToUpper.name, collectGreetings.name, collectGreetings.name nochmal und cowpy sagen.

Wenn wir jetzt ausführen, wird das Basisverzeichnis immer noch custom-outdir-config sein. Und es wird in einem Unterverzeichnis namens params.batch sein, aber die Unterverzeichnisse darunter sollten nach Prozessnamen organisiert sein.

Probieren wir das einfach aus und sehen, ob es funktioniert. Ich werde also das vorherige Verzeichnis entfernen, damit wir nicht verwirrt werden, und genau denselben Nextflow-Run-Befehl verwenden.

Es sollte auf die gleiche Weise laufen. Ich könnte dash resume bei all diesen verwenden, um es etwas schneller zu machen und die zuvor berechneten Ergebnisse zu verwenden. Wenn ich jetzt tree custom-outdir-config mache, kannst du sehen, es ist nicht in Results, es ist in unserem Basisverzeichnis mit dem Batch-Namen. Und du kannst sehen, alle Ergebnisse sind jetzt in Unterverzeichnissen organisiert, die nach dem Prozess benannt sind. Wir haben also zwei verschiedene Stellen, an denen wir hier dynamische Ausgabepfade definieren.

Okay. Letzte Sache, fügen wir diese Zwischenordner zurück, die wir vorher hatten, weil sie irgendwie nett waren. Intermediates.

Und wir können auch ein bisschen über dieses params.batch nachdenken, vielleicht als Pipeline-Entwickler\*in mochte ich wirklich, das im Unterverzeichnis zu haben, aber wenn Endnutzer\*innen der Pipeline "-o" oder -output-dir auf der CLI setzen, überschreibt es komplett diese gesamte Anweisung, und wir verlieren dieses Unterverzeichnis.

Was wir also machen können, ist diesen dynamischen Pfad aus der outputDir-Config zu nehmen, die überschrieben werden würde, und ihn in den Ausgabepfad zu legen, der nicht überschrieben wird.

Wir können also params.batch slash intermediates slash sayHello.name machen, und all das in einem doppelt zitierten String, damit es von Nextflow interpoliert wird.

Kann das jetzt kopieren, hoppla. Das zu den anderen Prozessen runterkopieren. Denk daran, sie alle in Anführungszeichen zu setzen. Und intermediates von diesen speziellen Ausgaben zu entfernen.

Okay? Es sieht jetzt etwas komplexer aus, aber du kannst sehen, wir bauen wirklich eine schön organisierte Ausgabeverzeichnisstruktur in unserem Code auf.

Und was wirklich schön ist, ist, dass diese zusätzliche Komplexität im Code nicht an die CLI durchgereicht wird. Wir können also unseren Befehl mit -output-dir und welchen Batch-Variablen auch immer ausführen, einfach darüber nachdenken, wie die Pipeline ausgeführt wird, und nicht wirklich zu viel darüber nachdenken, was im Code ist. Und unsere Ausgabedateien werden wirklich schön auf eine sehr gut organisierte Weise konstruiert, was nett für Leute ist, die die Pipeline verwenden.

Super. Während ich das schreibe, merke ich, dass ich einen Fehler gemacht habe. Mal sehen, ob jemand mich hier ertappt hat. Wir haben collectGreetings.name, also ist etwas schiefgegangen. Und ja, sicher genug, ich habe aus Versehen vergessen, diese in geschweifte Klammern zu setzen.

Denk also daran, vorsichtig zu sein, wenn du deinen Code schreibst, und stelle sicher, dass du Nextflow sagst, was eine Variable ist und was nur ein String ist. Denn es wird genau das tun, was du ihm sagst. Und nicht mehr. Wie alle guten Computer. Okay, das sollte es beheben.

## 2.3. Den Publish-Modus auf Workflow-Ebene setzen

Es gibt einen Teil dieses Skripts, den ich immer noch nicht liebe, nämlich die Tatsache, dass wir mode copy immer wieder schreiben, und wenn es eine Sache gibt, die wir nicht mögen, dann ist es uns zu wiederholen.

Wir können das also etwas aufräumen, indem wir das nehmen und in die Config verschieben. Und tatsächlich können wir es für die gesamte Pipeline auf einmal setzen. Also müssen wir es nicht mehrmals sagen.

Wir gehen zu unserer Config-Datei und wir haben einen neuen Bereich hier namens workflow. Und wir können entweder geschweifte Klammern machen oder wir können Punkt-Notation verwenden. Macht keinen Unterschied. Die Trainingsmaterial-Website verwendet Punkt-Notation. Ich kann output sagen und wir können mischen und abgleichen, also mode equals copy. Super.

Und jetzt können wir hier zurückgehen und diese löschen. Jetzt könnten wir sie an Ort und Stelle lassen. Die Config überschreibt im Grunde, was hier geschrieben ist, aber da wir es in der Pipeline-Level-Config haben, und diese beiden Dateien zusammen ausgeliefert werden, gibt es wirklich keinen Grund, es zweimal zu machen.

Okay. Nur eine Überprüfung von uns selbst, weil wir anscheinend Fehler machen. Lass uns das nochmal ausführen und einfach überprüfen, dass wir den Copy-Modus für das Veröffentlichen von Dateien richtig verwenden. Wir werden also das Skript nochmal ausführen und dieses Mal haben wir die Ergebnisse in ein Verzeichnis namens config-output-mode gelegt, schauen, wie die Dateien dort aussehen.

Und dann, wenn ich "ls -l" mache, um batch anzuschauen, und wir können zum Beispiel cowpy anschauen. Und wir sollten sehen, ja, dass dies eine richtige Datei hier ist, die kein Soft-Link ist, also wurde dieses Config-Attribut richtig angewendet.

## 3. Eine Software-Packaging-Technologie auswählen

Okay. Bisher haben wir uns auf die Eingaben und Ausgaben konzentriert, die Dateien, mit denen der Workflow läuft. Aber wie sieht es mit der Infrastruktur aus? Ich habe am Anfang gesagt, dass Nextflow es dir erlaubt, dieselbe Pipeline auf verschiedenen Computing-Setups auszuführen. Wie sieht das also aus?

Um das zu zeigen, werden wir von der Verwendung von Docker zum Ausführen von cowpy wechseln, und stattdessen werden wir Conda verwenden, um dasselbe zu tun.

Ich kann das ganz einfach machen. Wenn ich zu code gehe, "nextflow.config". Wenn du dich erinnerst, oben haben wir docker.enabled früher definiert, im letzten Kapitel, damit wir den Container mit cowpy verwenden konnten.

Ich werde Nextflow sagen, Docker nicht zu verwenden. Setze das auf false. Und ich werde sagen Conda enabled equals true. Also sage Nextflow, bitte verwende Conda.

Jetzt ist es nicht genug, Conda einfach nur zu aktivieren. Genau wie bei Docker müssen wir Nextflow sagen, wo es die Software bekommen kann, die es braucht.

Wenn wir also in die Module hier hüpfen. Und das cowpy-Skript öffnen. Wir können sehen, dass wir oben eine Container-Deklaration haben. Und der Container wird von Docker verwendet, aber auch Singularity, Apptainer und vielen anderen Software-Tools.

Aber er kann nicht für Conda verwendet werden, also haben wir eine separate Deklaration namens "conda", und wir könnten einfach "cowpy" schreiben. Und das wird es der Conda-Paket-Auflösung überlassen, den besten Weg zu finden, das zu lösen, gemäß deiner lokalen Conda-Umgebung.

Oder es ist gute Praxis zu tun, was die Trainingsmaterial-Website sagt zu tun, nämlich einen bestimmten Conda-Channel mit seiner Doppelpunkt-Notation zu definieren, und definitiv eine bestimmte Version der Software zu definieren, damit jede Person, die die Pipeline ausführt, dieselbe Version bekommt.

Beachte, dass Container in dieser Hinsicht etwas überlegen sind, denn wenn du etwas mit Conda installierst, wird es immer noch alle Abhängigkeiten für dieses Paket herausfinden, und die können sich im Laufe der Zeit ändern. Das nennt sich Dependency Drift.

Container hingegen sperren den gesamten Stack aller Software-Abhängigkeiten komplett ein, also kannst du etwas zuversichtlicher sein, dass A, es funktionieren wird, und B, es reproduzierbar sein wird.

Wenn du also in der Lage bist, Docker oder Singularity oder Apptainer zu verwenden, würde ich das definitiv empfehlen.

Was jetzt schön ist, ist, dass die Moduldatei, die von Pipeline-Entwickler\*innen geschrieben wird, jetzt sowohl Container als auch Conda hat, und so sagen wir der Person, die diese Pipeline ausführt, es ist uns egal, welche Software-Packaging-Lösung du verwendest. Es wird sowohl mit Docker als auch mit Conda funktionieren, und hier bekommst du die Software in beiden Fällen.

Wir können das Terminal hochziehen und lass uns das ausprobieren. Also Nextflow run hello config --batch conda. Und das erste Mal, wenn das mit Conda läuft, wird es ein bisschen langsam sein, wenn es zu diesem bestimmten Prozess kommt, weil es "conda install" ausführen muss.

Und es erstellt eine spezielle Conda-Umgebung nur für diesen einen Prozess. Es verwendet also nicht meine globale Conda-Umgebung, die ich auf meinem Terminal habe. Es erstellt eine nur für diesen einen Prozess. Das ist gut, weil es Dinge wie Abhängigkeitskonflikte zwischen verschiedenen Prozessen in deinem Workflow vermeidet. Wenn deine Prozesse Tools haben, die verschiedene Versionen von Python oder so brauchen, ist das okay, weil sie verschiedene Conda-Umgebungen verwenden.

Nextflow cached diese Conda-Umgebungen lokal, du kannst sehen, es sagt dir, wo dieser Pfad ist, er ist im Work-Verzeichnis hier. Und beim nächsten Mal, wenn ich dieses Skript mit Conda ausführe, wird es viel schneller sein, weil es diese existierende Conda-Umgebung finden und einfach wiederverwenden wird. Aber das erste Mal müssen wir es holen, auflösen, alle Abhängigkeiten herunterladen und alles einrichten.

Okay, super, es ist gelaufen. Wir können uns kurz erinnern, wofür die Pipeline derzeit konfiguriert ist. Wenn wir in die Config-Datei schauen, war es "custom-outdir-config" gerade für mich. Schauen wir mal, ob ich zu diesem Basisverzeichnis gehe. Und ich habe --batch conda gemacht. Da ist unser Conda-Unterverzeichnis. Es hat also funktioniert und da ist unsere cowpy-Ausgabe.

Es hat also cowpy geholt, es auf meinem lokalen System mit Conda installiert und den Prozess ausgeführt. Und was toll ist, ist, dass ich als Endnutzer\*in überhaupt nicht über das Software-Management nachdenken musste. Nextflow hat es einfach für mich sortiert. Ich habe gesagt, ich muss Conda auf diesem System verwenden. Die Pipeline-Entwickler\*in hat gesagt, welche Pakete ich brauche. Und Nextflow hat den Rest gemacht. Sehr mächtig.

Beachte, dass du tatsächlich eine Mischung verschiedener Technologien verwenden kannst. Ich kann also Docker für bestimmte Prozesse aktivieren und Conda für andere Prozesse, oder sagen, dass einige Prozesse einfach welche lokale Software ich installiert habe verwenden sollen. Das ist ziemlich ungewöhnlich, aber es ist möglich, und in einigen Fällen, zum Beispiel, wenn du bestimmte Software verwendest, die schwierig in Docker zu packen sein könnte, hast du einen Ausweg.

## 4. Eine Ausführungsplattform auswählen

Das war also Software-Packaging. Der andere Teil der Portabilität auf andere Systeme ist, wo die Jobs tatsächlich laufen. Im Moment laufe ich im Grunde auf meinem_Laptop oder in diesem Codespaces, was ein einzelner Computer ist. Es gibt nichts Ausgefallenes. Nextflow ist ein bisschen clever beim Parallelisieren der Jobs, so gut es kann, aber es ist alles auf einem System.

Wenn du auf einem HPC läufst, hast du wahrscheinlich irgendeine Art von Job-Scheduler wie SLURM oder PBS oder so, und du reichst Jobs an diesen Scheduler ein und er verteilt alle Jobs an verschiedene Compute-Nodes.

Eine andere Art zu laufen ist in der Cloud. Vielleicht verwendest du also AWS Batch oder Azure Cloud oder Google. Und diese funktionieren alle in einem ähnlichen System, wo du einen Scheduler hast und du Jobs einreichst und sie werden an verschiedene Orte zur Berechnung eingereicht.

In der fernen Vergangenheit, als ich anfing, Bioinformatik zu machen, war die Software aller Leute für die Durchführung von Analysen sehr an ihre Computerinfrastruktur gebunden, was es fast unmöglich machte zu replizieren.

Aber mit dieser Config-Trennung in Nextflow und mit Nextflows Fähigkeit, mit sehr vielen verschiedenen Compute-Infrastruktur-Backends zu interagieren, ist es sehr einfach, unsere Pipeline zu nehmen, ohne den Pipeline-Code überhaupt zu modifizieren, und das einfach auszutauschen.

## 4.1. Ein anderes Backend anvisieren

Wenn wir zu unserer "nextflow.config"-Datei gehen, können wir jetzt etwas Config auf Prozess-Ebene einfügen. Wenn ich oben den Process-Bereich einfüge, kann ich den Executor setzen, und hier ist er auf local gesetzt, was der Standard ist.

Beachte, dass wir, weil dies auf Prozess-Ebene ist, Dinge auf verschiedene Prozesse ausrichten können. Und du kannst Executors tatsächlich prozessspezifisch einrichten und eine Hybrid-Ausführung haben, wo einige Jobs lokal laufen könnten, wo auch immer der Nextflow-Job ausgeführt wird. Einige werden an verschiedene HPCs eingereicht und einige könnten an die Cloud eingereicht werden. Du kannst so clever sein, wie du willst.

Es ist jetzt sehr schwierig, das in einer Trainingsumgebung wie dieser zu demonstrieren, weil ich keinen HPC habe, an den ich einreichen kann. Aber was ich machen kann, ist, wenn ich slurm eintippe, können wir ein bisschen schummeln und du kannst ein Gefühl dafür bekommen.

Und das ist wirklich am interessantesten für Leute, die daran gewöhnt sind, auf SLURM zu laufen und wissen, wie die SLURM-Header aussehen. Aber wenn ich Nextflow run mache, hello config. Es wird scheitern, weil es versuchen wird, Jobs an einen Cluster einzureichen, der nicht existiert. Also werden wir irgendeine Fehlermeldung bekommen, dass sbatch nicht verfügbar ist.

Ja, geschrieben. Das ist das Tool. Das ist das CLI-Tool, das du verwendest, um Jobs an einen SLURM-Cluster einzureichen. Aber was wir machen können, ist, wir können in unser Work-Verzeichnis hier gehen, per Befehl klicken, dieses Verzeichnis öffnen und die .command.run anschauen. Und du kannst oben in der .command.run-Datei sehen, wir haben unsere sbatch-Header, die einem theoretischen SLURM-Cluster sagen, wie diese Job-Einreichung zu handhaben ist.

Und du kannst sehen, dass Nextflow clever ist, es macht alle richtigen Dinge. Es ist nur so, dass wir keinen Cluster hatten, an den wir einreichen konnten.

## 5. Compute-Ressourcen-Zuweisungen kontrollieren

Was ist noch unterschiedlich zwischen verschiedenen Computing-Infrastrukturen? Eine andere Sache ist, wie viel verfügbare Ressourcen du hast, und tatsächlich ist es in vielen Compute-Umgebungen eine Anforderung, dass du angeben musst, wie viele CPUs und wie viel Speicher ein Job braucht.

Auch hier abstrahiert Nextflow das für uns, sodass es nicht mehr spezifisch für einen einzelnen Compute-Umgebungstyp ist, und wir können hier in den Process-Level-Bereich eingeben. CPUs equals eins, memory equals zwei Gigabyte. Unsere Pipeline ist nicht sehr anspruchsvoll, also sollte das in Ordnung sein.

Ich habe diese Zahlen hier nur geraten, aber woher weißt du, was eine vernünftige Menge an Ressourcen zu verwenden ist? Es ist eine ziemlich schwierige Aufgabe, rauszugehen und all diese verschiedenen Prozesse einer großen Pipeline mit vielen Proben durchzugehen und zu verstehen, was die Ressourcennutzung war.

Ein guter Ansatz dafür ist also, diese Werte zu Beginn auf hohe Zahlen zu setzen, nur damit deine Pipeline ohne Fehler läuft, und dann Nextflow zu bitten, einen Nutzungsbericht für dich zu generieren.

Das ist super einfach zu machen, ich werde also zurück zu einem Terminal gehen. Oh, ich muss daran denken, das zurück auf local zu setzen, damit meine Pipeline tatsächlich läuft. Und ich sage Nextflow run, und ich werde ein Kommandozeilen-Flag verwenden -with-report.

Und ich kann das leer lassen und es wird einen Standard-Dateinamen geben, aber ich werde ihm einen bestimmten Dateinamen geben, damit das an einem bestimmten Ort gespeichert wird.

Enter drücken, und die Pipeline läuft genau wie normal, aber wenn sie fertig ist, wird sie einen schönen HTML-Bericht für mich generieren.

In der Seitenleiste hier habe ich diese HTML-Datei. Wenn ich das lokal laufen lassen würde, würde ich es einfach öffnen. Weil ich in Codespaces bin, werde ich mit Rechtsklick darauf klicken und Download klicken, was es auf meinen lokalen Computer herunterlädt. Und ich kann es einfach im Webbrowser öffnen.

Nextflow kann so einen Bericht für jede Pipeline generieren und er hat einige wirklich schöne Informationen. Es ist also gute Praxis, diese Dinge immer zu speichern. Es sagt uns, wann wir gelaufen sind, wo wir gelaufen sind, ob es erfolgreich war oder nicht, welche Parameter verwendet wurden, was der CLI-Befehl war, solche Sachen.

Und es gibt auch diese Plots über Ressourcennutzung. Es sagt uns also, welcher Prozentsatz der CPU-Aufrufe für jeden Prozess als Boxplot verwendet wurde, weil es viele Aufgaben für jeden Prozess gibt, also können wir die Verteilung sehen.

Du kannst unsere Prozesse hier sehen, cowpy und collectGreetings hatten nur eine einzige Aufgabe, also ist es nur eine einzelne Linie. Und wir haben sowohl CPU als auch Speicher und Job-Dauer, und sie waren sehr schnell.

Wenn du übrigens Seqera Platform verwendest, bekommst du dieselben Plots in die Platform-Oberfläche eingebaut, ohne etwas tun zu müssen. Du bekommst diese Information also immer griffbereit.

Okay, wir können diesen Bericht verwenden und bei einem echten Lauf, und ein Gefühl dafür bekommen, wie viele CPUs und wie viel Speicher von unserer Pipeline verwendet werden, und zurückkommen und diese Werte zurück in unsere Config-Datei einfügen, sodass wir beim nächsten Mal vielleicht nicht ganz so viel anfordern. Und wir können etwas schlanker sein.

Du kannst jetzt richtig clever werden beim Konfigurieren von Pipeline-Config-Dateien. Und wenn du Seqera Platform verwendest, halte Ausschau nach einem kleinen Knopf, der wie eine Glühbirne aussieht. Denn wenn du darauf klickst, wird es eine hochoptimierte Config-Datei generieren, die speziell auf deine Daten, deinen Lauf und deine Pipeline zugeschnitten ist. Um es auf die effizienteste Weise möglich auszuführen.

Aber für jetzt werde ich sagen, dass eigentlich die Standard-Anzahl von CPUs, die Nextflow gab, in Ordnung war, aber nur ein Gigabyte Speicher brauchen.

## 5.3. Ressourcen-Zuweisungen für einen bestimmten Prozess setzen

Im wirklichen Leben ist es ziemlich ungewöhnlich, dass alle Prozesse in deiner Pipeline dieselben Anforderungen brauchen werden. Du hast vielleicht etwas wie MultiQC als Reporting-Tool, das sehr wenig an Ressourcen braucht und ziemlich schnell läuft.

Und dann hast du vielleicht etwas, das ein Referenzgenom indiziert oder ein Alignment macht oder irgendeinen anderen Job macht. Es ist egal, was es ist, das eine Menge Ressourcen braucht. Und für diese verschiedenen Job-Einreichungen an einen Scheduler willst du also verschiedene Mengen an Ressourcen geben.

Unter diesem Process-Bereich können wir eine Config definieren, die bestimmte Prozesse auf verschiedene Arten anvisiert.

Hier verwenden wir withName, wir können auch Labels verwenden, und diese können ein Muster verwenden, um einen oder mehrere Prozesse anzuvisieren. Hier sagen wir nur, alle Prozesse, die einen Namen cowpy haben, setzen auf zwei Gigabyte Speicher und zwei CPUs, und weil dies ein spezifischerer Selektor ist als der oberste Process-Level, wird dies in diesen Fällen überschrieben, also kannst du hier eine schöne Config-Datei aufbauen, die wirklich alle deine verschiedenen Prozesse in deiner Pipeline anpasst, um sie wirklich effizient zu machen.

## 5.5. Ressourcen-Limits hinzufügen

Als Pipeline-Entwickler\*in kenne ich die Tools wahrscheinlich ziemlich gut, und ich möchte, dass alles so schnell und so gut wie möglich läuft. Es könnte also sein, dass ich ziemlich hohe Zahlen für einige davon eingebe, weil ich weiß, dass es viel schneller laufen wird, wenn ich cowpy 20 CPUs statt zwei gebe.

Das ist in Ordnung, bis du auf deinem Laptop oder auf GitHub Actions Continuous Integration-Test läufst oder auf irgendeinem anderen System, das vielleicht nicht 20 verfügbare CPUs hat.

Wenn du jetzt versuchst, die Pipeline auszuführen, wird sie abstürzen, weil Nextflow sagen wird, ich kann diesen Job nirgendwohin einreichen. Ich habe die verfügbaren Ressourcen nicht.

Um diesen harten Crash zu vermeiden, können wir jetzt etwas mehr Config hinzufügen, die spezifisch für unser System ist, namens Resource Limits. Und das sieht so aus. Es ist wieder unter dem Process-Bereich.

Und Resource Limits, du kannst im Grunde die Obergrenze von dem angeben, was du verfügbar hast. Es ist hier eine Map, und du kannst innerhalb dieser Map den Speicher, die CPUs und die Zeit setzen.

Was jetzt passiert, ist, wenn Nextflow eine Aufgabe von einem Prozess einreicht, schaut es, was angefordert wurde, und macht im Grunde einfach ein Minimum zwischen dem und dem. Wenn wir also 20 CPUs angefordert haben, aber nur vier verfügbar sind, wird es vier anfordern. Die Pipeline stürzt nicht ab und sie verwendet so nah wie möglich an dem, wofür sie von Pipeline-Entwickler\*innen designed wurde.

## 6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln

Okay. Ich habe gesagt, dass die Resource Limits hier systemspezifisch sein könnten, und vielleicht habe ich eine Nextflow-Config-Datei in meiner Pipeline, und ich weiß, dass Leute das an einer Reihe verschiedener Orte verwenden werden. Anstatt jetzt jeden zu zwingen, jedes Mal seine eigene Nextflow-Config-Datei zu erstellen, kann ich verschiedene Voreinstellungen von Konfiguration zusammen in Config-Profile gruppieren.

Ich werde hier ein bisschen runterscrollen und dann direkt nach params, weil die Reihenfolge der Config-Datei hier wichtig ist, die Config-Datei wird sequenziell geladen, also werde ich diese Profile nach allem anderen setzen, damit es die zuvor definierten Params überschreibt. Und ich werde diese Profile vom Trainingsmaterial einfügen.

Es gibt also einen neuen Top-, Top-Level-Bereich namens profiles. Wir können hier beliebige Namen haben. Also haben wir my_laptop und univ_hpc. Und hier können wir sehen, dass wir die anderen gleichen Config-Parameter setzen, die wir vorher hatten. Jetzt innerhalb nur eines Profils. Wir haben also einen lokalen Executor für das Laufen auf my_laptop und ich reiche an einen SLURM-Cluster auf dem HPC ein.

Ich verwende Docker lokal, Conda auf dem HPC, und das HPC-System hat viel höhere Resource Limits.

Ich kann jetzt die Pipeline mit der -profile CLI-Option ausführen, sagen, welches Profil ich verwenden möchte. Ich werde also my_laptop verwenden, und Nextflow wird alle Config innerhalb dieses Profil-Bereichs anwenden. Das kann ich jetzt ausprobieren. Es ist derselbe Befehl wie zuvor. Nextflow run hello config, und ich mache dash profile, einzelner Bindestrich, weil es die Kern-Nextflow-Option ist, dash profile my_laptop.

Es wird jetzt Batch alle Config-Optionen anwenden. Oh, und du kannst sehen, ich habe vorher gesagt, das könnte passieren, dass die Process-Anforderung nach vier CPUs fragte und ich nur zwei auf dieser Codespaces-Instanz habe.

Das ist also eine gute Gelegenheit, die Process Resource Limits einfach auszuprobieren, und zu sagen, dass ich nur zwei CPUs auf my_laptop habe, oder in diesem Codespaces. Wenn wir es jetzt nochmal ausführen, sollte es diese Anforderung auf zwei begrenzen und hoffentlich wird die Pipeline laufen. Super.

## 6.2. Ein Profil für Testparameter erstellen

Beachte, dass diese Profile nicht nur Konfiguration über ihre Infrastruktur haben müssen. Du kannst Gruppierungen von beliebiger Config hier haben, einschließlich Parametern.

Eine andere Sache, die du also sehr oft in Pipelines von Leuten sehen wirst, ist ein Test-Profil, das Parameter enthält, die du normalerweise auf Benutzer\*innenbasis einreichen würdest. Aber hier haben wir im Grunde verschiedene sinnvolle Standards, wenn ich Testfälle ausführen möchte.

Und das ist toll, weil ich nicht notwendigerweise all diese Dinge angeben muss, die erforderliche Parameter sein könnten. Andernfalls kann ich einfach dash profile test sagen und es wird einfach out of the box laufen.

Etwas zu beachten ist, dass Profile auch mehr als eins kombiniert werden können. Ich kann also profile my_laptop hier machen, und dann auch test hinzufügen. Ich mache profile nicht zweimal. Ich mache nur eine kommagetrennte Liste hier ohne Leerzeichen. Und es wird diese Profile in Reihenfolge anwenden. Es wird also die Config vom my_laptop-Profil nehmen und dann die Test-Config obendrauf anwenden.

Wirklich praktisch und du kannst sehen, wie du hier viele sinnvolle Standardgruppen einrichten kannst, um es einfach zu machen, deine Pipeline auszuführen.

## 6.3. nextflow config verwenden, um die aufgelöste Konfiguration zu sehen

Hoffentlich habe ich dich überzeugt, dass Nextflow-Config-Auflösung mächtig ist, aber ich würde dir nicht vorwerfen, wenn du an diesem Punkt ein bisschen Augenverdreher bekommst, nachdem ich etwa 20 verschiedene Wege gesagt habe, Config bereitzustellen und all diese verschiedenen Schichten wie eine Zwiebelschale zu geben.

Wenn du dir also jemals unsicher bist, was die endgültige aufgelöste Config für Nextflow ist, wisse, dass es einen Befehl namens "nextflow config" gibt, und wir können das ausführen und es wird uns sagen, was die aufgelöste Konfiguration an unserem aktuellen Standort ist.

Wenn ich es also hier ausführe, findet es die "nextflow.config"-Datei im aktuellen Arbeitsverzeichnis, und es verarbeitet alle verschiedenen Configs, und es gibt mir die aufgelöste Ausgabe.

Beachte, dass die Nextflow-Config-Datei auch die profile CLI-Option nehmen kann. Wenn ich es also sage, in my_laptop- und Test-Profilen aufzulösen, und du kannst sehen, es hat auch die Resource Limits hier von der my_laptop-Config-Option angewendet und auch die Params gesetzt, die im Test waren.

Das ist also eine nette Art, einfach zu erforschen, wie die Config-Auflösung funktioniert, wenn du dir überhaupt unsicher bist.

## Abschluss

Okay, das war's. Das ist Nextflow-Config in einer Nussschale. Du kannst eine Menge Sachen mit Config machen. Es ist wirklich mächtig. Aber das sind die meisten der häufigen Anwendungsfälle, denen du begegnen wirst, und diese Konzepte gelten für alle verschiedenen Optionen.

Klopfe dir selbst auf die Schulter, denn das ist das Ende des Hello Nextflow-Trainingskurses. Du bist hoffentlich jetzt zuversichtlich sowohl beim Schreiben deiner eigenen Nextflow-Pipeline von Grund auf, beim Konfigurieren und Ausführen, und du kennst alle Ins und Outs und die Dinge, auf die du achten musst.

Es gibt noch ein Quiz, das du auf der Config-Trainingsseite ausprobieren kannst. Geh also runter und probiere das aus und stelle sicher, dass du all diese Teile über die Config verstanden hast.

Und schließe dich uns im letzten Video an, nur für einen kurzen Abschluss über einige der nächsten Schritte, die nach diesem Trainingskurs gut sein könnten.

Danke, dass du bei uns geblieben bist. Gut gemacht und wir sehen uns im nächsten Video.
