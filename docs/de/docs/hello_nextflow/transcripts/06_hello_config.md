# Teil 6: Hello Config - Videotranskript

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Wichtige Hinweise"

    Diese Seite zeigt nur das Transkript. Für vollständige Schritt-für-Schritt-Anleitungen kehre zum [Kursmaterial](../06_hello_config.md) zurück.

    Die im Transkript angegebenen Abschnittsnummern dienen nur zur Orientierung und umfassen möglicherweise nicht alle Abschnittsnummern im Material.

## Willkommen

Hi und willkommen zurück zu Teil Sechs von Hello Nextflow. In diesem Abschnitt geht es um Configs, und es ist der letzte Teil dieses Kurses.

Nextflow ist besonders gut in zwei Dingen: Reproduzierbarkeit und Portabilität. Bei Configs sehen wir, wie die zweite dieser Eigenschaften wirklich glänzt. Die Fähigkeit, eine Nextflow-Pipeline so zu konfigurieren, dass sie auf verschiedene Arten läuft und auf verschiedenen Systemen funktioniert, ohne den zugrunde liegenden Pipeline-Code bearbeiten zu müssen.

Diese Superkraft ermöglicht es, dass Nextflow-Pipelines von anderen Leuten an anderen Orten wiederverwendet werden können, oder über verschiedene Infrastrukturen hinweg, auf die du selbst Zugriff haben könntest.

Das bedeutet, du kannst Pipeline-Code auf deinem Laptop entwickeln, ihn in die Cloud pushen, ihn auf deinem HPC ausführen, und es ist derselbe Pipeline-Code und er läuft überall.

In diesem Abschnitt werden wir einige Themen durchgehen. Wir beginnen damit, wie Nextflow mit Config-Dateien umgeht, woher es sie lädt, wie du sie schreibst und strukturierst, und diese Trennung zwischen der Pipeline selbst und dem, was in eine Config-Datei gehören sollte.

Dann gehen wir auf einige häufige Anwendungsfälle ein, wie zum Beispiel das Ändern, wo Ausgabedateien gespeichert werden, und auch, wie man die Pipeline auf verschiedenen Infrastrukturen zum Laufen bringt, sowohl mit verschiedenen Arten von Software-Packaging als auch beim Einreichen von Jobs an verschiedene Infrastrukturen.

## Config-Dateihierarchien

Okay, legen wir los. Wenn es um das Laden von Config-Dateien geht, kann Nextflow von vielen verschiedenen Orten laden, was eine gute Sache ist und auch eine etwas riskante Sache sein kann, weil es manchmal etwas schwierig sein kann zu wissen, woher es eine Config-Datei bekommt und in welcher Reihenfolge es Dinge lädt.

Ich empfehle dir wirklich, auf diesen Link hier zu klicken, der uns zur Nextflow-Dokumentation führt. Und auf dieser Konfigurationsseite werden die wichtigsten Orte aufgelistet, von denen Config geladen wird, und wichtig, die Reihenfolge, in der diese Dinge geladen werden.

Du kannst sehen, du kannst eine Config-Datei in dein Nextflow-Home-Verzeichnis legen, das typischerweise ".nextflow" in deinem Home-Verzeichnis ist. Und diese Datei wird immer von jedem Nextflow-Lauf auf deinem System geladen.

Der nächste Ort ist eine Datei im Root deines Pipeline-Repositorys oder -Verzeichnisses namens "nextflow.config".

Danach eine weitere Datei namens "nextflow.config", aber diesmal im Verzeichnis, von dem aus du Nextflow startest: dem Launch-Verzeichnis.

Schließlich kannst du Config-Dateipfade auf der Kommandozeile mit einem "-c"-Argument angeben, und du kannst das mehrmals tun. Und sie werden in der Reihenfolge angewendet, in der du sie angibst.

Du kannst Config-Dateien an all diesen Orten bereitstellen, wenn du möchtest, und sie werden iterativ geladen, wobei jede die vorherige nur in den Config-Scopes überschreibt, wo sie sich überschneiden.

Das ist ein wirklich mächtiges System, weil es bedeutet, dass du sinnvolle Standardwerte setzen und dann immer spezifischer werden kannst, je mehr du dich auf diese Config einengst.

## 0. Aufwärmen: hello-config.nf ausführen

Okay, lass uns das schließen und in unsere Codespaces springen und loslegen. Wie zuvor habe ich hier aufgeräumt, ich habe meine vorherigen Results-Verzeichnisse, mein Nextflow- und meine Work-Verzeichnisse usw. entfernt. Mach dir keine Sorgen, wenn du diese Dateien noch herumliegen hast. Es ist nur, weil ich sehr stark gezoomt bin und sonst die Dinge sehr schnell unübersichtlich werden.

Wir werden mit hello-config.nf arbeiten, der letzten Datei in unserem Verzeichnis, und dies sollte dort weitermachen, wo wir im vorherigen Abschnitt aufgehört haben.

Wir haben also unsere vier verschiedenen Prozesse, die aus Moduldateien eingebunden werden. Wir haben unsere Pipeline-Parameter, unseren Workflow-Block, wo wir die verschiedenen Prozesse aufrufen und die Kanäle zusammenfügen, die Ausgabekanäle veröffentlichen, und dann den Output-Block am Ende, wo wir definieren, wo diese Dateien gespeichert werden sollen und wie sie kopiert werden sollen.

Wir haben auch bereits eine "nextflow.config"-Datei aus dem letzten Kapitel, wo wir Docker aktivieren, und wir werden diese Datei heute erweitern.

Wie zuvor haben wir den Ausgabepfad in diesem Hauptskript auf hello config geändert, nur damit er nicht mit vorherigen Ergebnissen kollidiert, die du generiert hast.

Okay, lass uns schnell überprüfen, ob alles noch so funktioniert, wie wir es erwarten. Öffne ein Terminal und führe nextflow run hello-config.nf aus. Nextflow lädt. Sollte unsere vier verschiedenen Prozesse ausführen. Etwas schöne ASCII-Kunst mit cowpy generieren und dann unsere Ergebnisse in unseren Results-Dateien in diesem Verzeichnis speichern.

Ich kann hier schnell reinschauen, nur um sicherzustellen, dass diese Dateien so aussehen, wie wir es erwarten, und tatsächlich, da ist unser riesiger Truthahn. Großartig.

## 1.1. Standardwerte in nextflow.config verschieben

Jetzt werden wir als Erstes anfangen, einige Dinge aus unserem Skript in unsere Config-Datei zu verschieben.

Und was uns hauptsächlich interessiert, sind die Parameter zu diesem Zeitpunkt. Wir wollen die Standardwerte in die Config-Datei übernehmen, damit klarer ist, was die Standardwerte sind, und es für Leute einfacher ist, sie zu überschreiben.

Ich werde einfach diesen params-Block hier aus dem Skript nehmen und ihn in die Config-Datei einfügen. Und wir müssen hier etwas vorsichtig sein, denn im Moment ist die Syntax zwischen Config und Skripten etwas unterschiedlich. Die Config-Datei kann keine Typdeklarationen annehmen, weil wir diese params nicht wirklich definieren, wir referenzieren sie nur. Also werde ich diese loswerden.

Aber ansonsten ist es sehr ähnlich. Wir haben einen params-Block und dann haben wir unsere verschiedenen Eingabeparameter, Batch-Parameter, Character-Parameter.

Ich kann jetzt zurück zu meinem Skript gehen und ich muss diese Standardwerte nicht mehr definieren, weil diese Werte jetzt in meiner Nextflow-Config-Datei sind.

Allerdings lasse ich die Parameternamen und ihre Typen stehen, damit Nextflow diese Informationen kennt und immer noch die gesamte Typsicherheit und alles durchführen kann.

Okay. Wir speichern diese Dateien und überprüfen schnell, ob alles noch genauso funktioniert wie zuvor. Es sollte hier keine Änderungen geben. Wir haben die Werte gleich gehalten. Wir haben nur verschoben, wo sie definiert wurden.

Großartig.

## 1.2. Eine laufspezifische Konfigurationsdatei verwenden

Nun, bisher haben wir Nextflow aus demselben Verzeichnis gestartet, in dem wir unser Pipeline-Skript haben. Also sind unser Launch-Verzeichnis und unser Pipeline-Verzeichnis irgendwie dasselbe.

Um zu zeigen, wie wir verschiedene Config-Dateien mit verschiedenen Launch-Verzeichnissen haben können, werden wir jetzt ein neues Unterverzeichnis erstellen.

Also werde ich mkdir sagen, und wir werden es tux-run nennen.

Und dann werde ich cd, das Verzeichnis wechseln in tux-run. Und beachte, dass unser Arbeitsverzeichnis jetzt nicht mehr im selben Verzeichnis wie die Pipeline-Skripte ist.

Okay, lass uns eine neue "nextflow.config"-Datei erstellen. Also touch nextflow config, und lass uns sie einfach in VS Code öffnen. Du kannst auch in der Seitenleiste hier sehen, dass wir jetzt in diesem Unterverzeichnis sind.

Jetzt können wir denselben params-Block nehmen, den wir in unserer obersten nextflow.config hatten, dies hier kopieren und jetzt können wir diese Werte ändern.

Erstens, die Daten sind jetzt ein anderer relativer Pfad, weil wir in einem Unterverzeichnis sind, also müssen wir das aktualisieren. Und dann werden wir batch zu experiment ändern, und wir werden den character von Turkey zu tux ändern.

Jetzt klicke auf Speichern dort, und lass es uns ausprobieren. Genau wie bei data muss ich jetzt ../ sagen, um zum Skript zu gelangen. Also ist es Hello config. Und ich drücke Enter.

Der Pipeline-Code hat sich überhaupt nicht geändert, aber jetzt werden wir zwei Sätze von Config laden, und die Launch-Dir-Config-Datei sollte die Standardwerte überschreiben, die in der Pipeline nextflow.config gesetzt wurden, und wir sollten verschiedene Ergebnissätze bekommen.

Tatsächlich, innerhalb unseres Verzeichnisses hier, innerhalb von tux-run, kannst du sehen, dass wir ein .nextflow-Verzeichnis und ein work-Verzeichnis haben, und das liegt daran, dass diese immer in deinem Launch-Verzeichnis erstellt werden. Diese sind also anders als die work- und results-Verzeichnisse, die wir von früheren Läufen hatten.

Wenn ich jetzt in results schaue, können wir unsere collected sehen und da ist unser kleiner tux-Character. Du kannst also sehen, dass diese Parameter richtig interpretiert wurden.

## 1.3. Eine Parameterdatei verwenden

Okay. Vorhin, als ich über die verschiedenen Config-Dateien sprach, die geladen werden konnten, habe ich einen anderen Ort ausgelassen, von dem wir Config bekommen können.

Du kannst es von einer Kommandozeile bekommen, wie wir mit Doppel-Minus-Parameternamen gesehen haben, aber wir können auch eine YAML- oder eine JSON-Datei bereitstellen, nur mit params.

Die Config-Datei kann alle verschiedenen Arten von Scopes haben, aber diese Dateien sind nur Parameter, und es ist eine schöne benutzerfreundliche Art, viele Parameter auf einmal bereitzustellen, und vielleicht eine etwas reproduzierbarere Art, weil du sie in eine Datei schreibst, also ist es einfach, sie zu einem späteren Zeitpunkt zu bekommen.

Also lass uns zurück zu unserem Terminal gehen und nur bevor wir es vergessen, stelle sicher, dass wir ein Verzeichnis nach oben gehen, also bin ich nicht mehr im Unterverzeichnis, und ich werde mir die YAML-Datei anschauen, die wir hier haben, genannt test-params.yaml.

Wenn ich also einfach code test-params.yaml mache, kannst du sehen, das ist nur eine reguläre YAML-Datei. Nichts Besonderes daran. Mit den Schlüsseln als unsere Parameternamen, mit der YAML-Formatierung also ein Doppelpunkt hier, und dann ein Wert.

Beachte, dass dies kein Nextflow-Code ist, also können wir hier keine Dinge wie Variablen einfügen. Das sind nur statische Werte.

Auch weil JSON tatsächlich als YAML geparst wird, können wir auch eine test-params.json-Datei haben, die sehr ähnlich aussieht. Es ist nur ein anderes Datenformat.

Wir haben also zwei verschiedene Testdateien hier und wir haben leicht unterschiedliche Variablen.

Okay, wie geben wir diese also an Nextflow? Es ist sehr einfach. Wir machen Nextflow run hello config, wie zuvor. Und anstelle von "-c" für Config-Datei oder dem Laden dieser Standard-Dateinamen machen wir -params-file. Einzelner Bindestrich, weil es eine Kern-Nextflow-Option ist.

Und dann den Pfad für diese Datei übergeben. Also werde ich "-params-file test-params.yaml" machen, und wir werden sehen, ob diese richtig geladen werden.

Okay. Es ist gelaufen. Lass uns uns nur daran erinnern, was in dieser YAML-Datei war. Also war der batch auf YAML gesetzt, so sollte es genannt werden, und es sollte einen Stegosaurus haben. Also lass uns hochgehen und in results schauen. Und wir haben COLLECTED-yaml. Also lass uns sehen, ob wir einen Stegosaurus haben. Fantastisch, ein Stegosaurus mit Hut. Das ist es, was wir mögen.

Das hat also wirklich gut funktioniert, und es ist genau dasselbe mit der JSON-Datei. Wir tauschen einfach die Dateierweiterung hier aus und Nextflow weiß, wie man das liest.

Und in diesem Fall sollten wir einen batch namens JSON haben und wir sollten eine Schildkröte haben. Lass uns mal schauen. Wunderbar. Eines meiner Lieblings-CLI-Tools.

## 2.1. Das Ausgabeverzeichnis mit -output-dir anpassen

Okay, das war hauptsächlich über Eingaben in die Pipeline und das Ändern von Parametern. Was ist mit den Ausgaben?

Nun, obwohl wir die Unterverzeichnisse mit params geändert haben, hast du vielleicht bemerkt, dass alle unsere Dateien immer noch nach results gehen.

Wir können dieses Basisverzeichnis ändern, in das alle Dateien veröffentlicht werden, mit einem Kommandozeilen-Flag namens -output-dir. Wenn ich also Nextflow run hello config mache, und dann mache ich -output-dir, und wir werden es "custom-outdir-cli" nennen. Kann nicht tippen. Nur damit wir uns erinnern, woher diese Dateien kamen.

Das ist eine Kern-Nextflow-Option und es ist eine sehr neue. Diese wurde erst kürzlich hinzugefügt, und das ist eine der Dinge, die wir mit dem neuen Language-Parser und allem machen können.

Es ist etwas umständlich zu tippen. Du kannst es auch einfach "-o" nennen, wenn du möchtest. Wenn ich also einfach zurückgehe. Kann das einfach zu "-o" verkürzen, was etwas einfacher ist.

Okay. Wir führen das aus. Wir haben nichts in unserer Pipeline oder sogar in unserer Config zu diesem Zeitpunkt geändert, und es sollte hoffentlich alle unsere Ergebnisse in ein anderes Top-Level-Verzeichnis speichern. Und du kannst dir vorstellen, du kannst dies auf im Grunde jeden Pfad setzen, den du möchtest.

Es ist gerade oben angekommen. Wir haben ein custom-outdir-cli, und alle Dateien sind dort auf genau die gleiche Weise organisiert, mit ihren gleichen Unterverzeichnissen und Dateinamen. Das ist also eine wirklich einfache Möglichkeit, einfach zu ändern, wohin die Pipeline ihre Ergebnisse veröffentlicht, ohne zu viel darüber nachzudenken, wie diese Ergebnisse organisiert sind.

## 2.1.2. Hartcodierte Pfade aus dem Output-Block entfernen

Wenn ich in dieses Verzeichnis schaue, können wir sehen, dass wir immer noch ein Unterverzeichnis namens Hello Config haben, was sich jetzt etwas redundant anfühlt.

Also lass uns einfach unser Skript wieder laden und wir können jetzt dieses Unterverzeichnis aus dem Output-Block am Ende entfernen. Weil wir es nicht wirklich mehr brauchen. Also können wir das jetzt einfach löschen, das hier löschen. Und wenn es nur das ist, kannst du das entweder komplett löschen oder als leeren String belassen. Ich werde es vorerst als leeren String belassen, weil wir zurückkommen und in Zukunft einige verschiedene Dinge an seine Stelle setzen werden. Aber wenn dir Unterverzeichnisse egal sind, ist es am saubersten, die path-Deklaration dort komplett zu entfernen.

Okay, lass uns speichern. Nur schnell nochmal ausprobieren. Ich werde tatsächlich mein "custom-outdir-cli"-Verzeichnis entfernen, damit wir nicht durch bereits vorhandene Dateien dort verwirrt werden. Denn denk daran, wenn du Dinge veröffentlichst, entfernt es nicht die Dateien, die bereits dort waren. Es fügt nur neue hinzu. Lass uns diesen Befehl nochmal ausführen, custom-outdir-cli.

Und jetzt, wenn du "ls custom-outdir-cli" machst, gibt es dort kein Verzeichnis mehr namens Hello Config.

## 2.2.1. outputDir in der Konfigurationsdatei setzen

Okay, das Kommandozeilen-Flag hier, "-o" oder "-output-dir" ist gut. Aber wie sieht es mit dem Setzen von Standardwerten dafür in der Config aus? Wie machen wir das?

Ich öffne die "nextflow.config"-Datei, schließe alles andere und werde das los. Wir können hier eine neue Config-Option hinzufügen, die ich gerade von der Trainingsmaterial-Website kopiert habe, und sie heißt outputDir.

Sie ist nicht unter irgendeinem Scope. Sie ist nicht unter params oder irgendetwas. Sie ist Top-Level, und wir können dies auf einen String setzen. Nun, eine einfache Sache zu tun ist, es einfach auf irgendetwas anderes als results als hartcodierten String zu ändern. Aber weil dies in einer Nextflow-Config-Datei ist, können wir hier etwas cleverer sein und auch Variablen einschließen.

Und du kannst hier sehen, dass wir eine params-Variable eingeschlossen haben, params.batch, die Teil dieses Strings ist. Das bedeutet, dass wir Variablen wiederverwenden können, die von anderen Orten kommen. Und in diesem Fall, wenn wir --batch machen, wenn wir Nextflow Pipeline ausführen, werden wir ein Unterverzeichnis in unserem benutzerdefinierten Pfad basierend darauf bekommen, wie der Batch-Name war.

Okay, also lass uns das ausprobieren und einfach einen kurzen Blick darauf werfen, wie die Ergebnisse aussehen. Wenn ich also Nextflow run hello config mache und --batch my_run. Erinnern wir uns daran, wie die Config aussah. Also ist es custom-outdir-config.

Tree custom-outdir-config. Und du kannst sehen, der batch hieß my_run. Und dann haben wir dieses Unterverzeichnis namens my_run. Also hat dieser dynamische Dateipfad funktioniert.

Und nicht nur das, es ging nicht mehr in ein Standard-results-Verzeichnis, und ich musste nichts auf der Kommandozeile angeben, um das Basisverzeichnis zu ändern. Wir haben also erfolgreich den Standardwert für das Standard-outputDir zurückgesetzt.

## 2.2.2. Unterverzeichnisse mit Batch- und Prozessnamen

Okay, lass uns das etwas weiter treiben. Das ist eine dynamische Variable innerhalb der Config-Datei. Wie sieht es mit dem Skript selbst aus? Nun, bisher hatten wir diese Pfade hier und diese können auch dynamisch sein. Anstatt also nur etwas hartzucodieren, können wir einige geschweifte Klammern einfügen und etwas Dynamisches einfügen.

Zum Beispiel haben wir unsere Prozesse namens sayHello. Wir könnten sayHello.name machen, was ein Attribut des Prozesses ist, was irgendwie langweilig ist, weil es in diesem Fall nur "sayHello" ist. Aber es ist variabel.

Das gibt dir also eine Idee. Wir können das hier einfügen und convertToUpper.name sagen, collectGreetings.name, collectGreetings.name nochmal, und cowpy.

Wenn wir jetzt ausführen, wird das Basisverzeichnis immer noch custom-outdir-config sein. Und es wird in einem Unterverzeichnis namens params.batch sein, aber die Unterverzeichnisse darunter sollten nach Prozessnamen organisiert sein.

Lass uns das einfach ausprobieren und sehen, ob es funktioniert. Also werde ich das vorherige Verzeichnis entfernen, damit wir nicht verwirrt werden, und genau denselben Nextflow Run-Befehl verwenden.

Es sollte auf die gleiche Weise laufen. Ich könnte dash resume bei all diesen verwenden, um es etwas schneller zu machen und die zuvor berechneten Ergebnisse zu verwenden. Wenn ich jetzt tree custom-outdir-config mache, kannst du sehen, es ist nicht in results, es ist in unserem Basisverzeichnis mit dem Batch-Namen. Und du kannst sehen, alle Ergebnisse sind jetzt innerhalb von Unterverzeichnissen organisiert, die nach dem Prozess benannt sind. Wir haben also zwei verschiedene Orte, an denen wir hier dynamische Ausgabepfade definieren.

Okay. Letzte Sache, lass uns diese Zwischenordner zurückfügen, die wir vorher hatten, weil sie irgendwie nett waren. Intermediates.

Und wir können auch ein bisschen über dieses params.batch nachdenken, vielleicht als Pipeline-Entwickler\*in mochte ich es wirklich, das im Unterverzeichnis zu haben, aber wenn Endnutzer\*innen der Pipeline "-o" oder -output-dir auf der CLI setzen, überschreibt es komplett diese gesamte Anweisung, und wir verlieren dieses Unterverzeichnis.

Was wir also tun können, ist, wir können diesen dynamischen Pfad aus der outputDir-Config nehmen, die überschrieben würde, und ihn in den Output-Pfad einfügen, der nicht überschrieben wird.

Wir können also params.batch slash intermediates slash sayHello.name machen, und all das in einem doppelt zitierten String, damit es von Nextflow interpoliert wird.

Kann jetzt kopieren, hoppla. Diese zu den anderen Prozessen kopieren. Denk daran, sie alle in Anführungszeichen zu setzen. Und intermediates von diesen speziellen Ausgaben zu entfernen.

Okay? Es sieht jetzt etwas komplexer aus, aber du kannst sehen, wir fangen wirklich an, eine schön organisierte Ausgabeverzeichnisstruktur in unserem Code aufzubauen.

Und was wirklich schön ist, ist, dass diese zusätzliche Komplexität im Code nicht zur CLI durchgeht. Wir können also unseren Befehl mit -output-dir und welchen Batch-Variablen auch immer ausführen, nur darüber nachdenken, wie man die Pipeline ausführt und nicht wirklich zu viel darüber nachdenken, was im Code ist. Und unsere Ausgabedateien werden wirklich schön auf eine sehr gut organisierte Weise konstruiert, was für Leute, die die Pipeline verwenden, im Grunde nett ist.

Großartig. Während ich das schreibe, merke ich, dass ich einen Fehler gemacht habe. Mal sehen, ob mich hier jemand erwischt hat. Wir haben collectGreetings.name, also ist etwas schiefgelaufen. Und ja, tatsächlich, ich habe versehentlich vergessen, diese in geschweifte Klammern zu setzen.

Also denk daran, sei vorsichtig, wenn du deinen Code schreibst, und stelle sicher, dass du Nextflow sagst, was eine Variable ist und was nur ein String ist. Denn es wird genau das tun, was du ihm sagst. Und nichts mehr. Wie alle guten Computer. Okay, das sollte es beheben.

## 2.3. Den Publish-Modus auf Workflow-Ebene setzen

Es gibt einen Teil dieses Skripts, den ich immer noch nicht liebe, und das ist die Tatsache, dass wir mode copy immer und immer wieder schreiben, und wenn es eine Sache gibt, die wir nicht mögen, dann ist es, uns zu wiederholen.

Wir können das also etwas aufräumen, indem wir das nehmen und in die Config verschieben. Und tatsächlich können wir es für die gesamte Pipeline auf einmal setzen. Wir müssen es also nicht mehrmals sagen.

Wir gehen zu unserer Config-Datei und wir haben hier einen neuen Scope namens workflow. Und wir können entweder geschweifte Klammern machen oder wir können Punkt-Notation verwenden. Macht keinen Unterschied. Die Trainingsmaterial-Website verwendet Punkt-Notation. Ich kann output sagen und wir können mischen und anpassen, also mode equals copy. Großartig.

Und jetzt können wir hier zurückgehen und diese löschen. Wir könnten sie an Ort und Stelle lassen. Die Config überschreibt im Grunde, was hier geschrieben ist, aber da wir es in der Pipeline-Level-Config haben, und diese beiden Dateien zusammen ausgeliefert werden, gibt es keinen Grund, es wirklich zweimal zu tun.

Okay. Nur eine Plausibilitätsprüfung für uns selbst, weil wir anscheinend Fehler machen. Lass uns das nochmal ausführen und einfach überprüfen, dass wir den Copy-Modus korrekt zum Veröffentlichen von Dateien verwenden. Wir werden also das Skript nochmal ausführen und diesmal haben wir die Ergebnisse in ein Verzeichnis namens config-output-mode gesteckt, schauen, wie die Dateien dort aussehen.

Und dann, wenn ich "ls -l" mache, um batch anzuschauen, und wir können zum Beispiel cowpy anschauen. Und wir sollten sehen, ja, dass dies hier eine richtige Datei ist, die kein Soft-Link ist, also wurde dieses Config-Attribut richtig angewendet.

## 3. Eine Software-Packaging-Technologie auswählen

Okay. Bisher haben wir uns auf die Eingaben und die Ausgaben konzentriert, die Dateien, mit denen der Workflow läuft. Aber wie sieht es mit der Infrastruktur aus? Ich sagte am Anfang, dass Nextflow es dir ermöglicht, dieselbe Pipeline auf verschiedenen Computing-Setups auszuführen. Wie sieht das also aus?

Um das zu zeigen, werden wir von der Verwendung von Docker zum Ausführen von cowpy wechseln, und stattdessen werden wir Conda verwenden, um dasselbe zu tun.

Ich kann das sehr einfach machen. Wenn ich zu code gehe, "nextflow.config". Wenn du dich erinnerst, haben wir oben docker.enabled früher definiert, und im letzten Kapitel, damit wir den Container mit cowpy darin verwenden konnten.

Ich werde Nextflow sagen, Docker nicht zu verwenden. Setze das auf false. Und ich werde sagen Conda enabled equals true. Also sage Nextflow, bitte verwende Conda.

Nun, nur Conda zu aktivieren ist nicht genug für sich. Genau wie wir es mit Docker gemacht haben, müssen wir Nextflow sagen, wo es die Software bekommen kann, die es braucht.

Wenn wir also in die Module hier hüpfen. Und das cowpy-Skript öffnen. Wir können sehen, wir haben oben eine Container-Deklaration. Und der Container wird von Docker verwendet, aber auch von Singularity, Apptainer und vielen der anderen Software-Tools.

Aber er kann nicht für Conda verwendet werden, also haben wir eine separate Deklaration namens "conda", und wir könnten einfach "cowpy" schreiben. Und das wird es der Conda-Paketauflösung überlassen, den besten Weg herauszufinden, das gemäß deiner lokalen Conda-Umgebung zu lösen.

Oder es ist gute Praxis zu tun, was die Trainingsmaterial-Website sagt, nämlich einen spezifischen Conda-Kanal mit seiner Doppelpunkt-Notation zu definieren, und definitiv eine spezifische Version der Software zu definieren, damit jede Person, die die Pipeline ausführt, dieselbe Version bekommt.

Beachte, dass Container in dieser Hinsicht etwas überlegen sind, denn wenn du etwas mit Conda installierst, wird es immer noch alle Abhängigkeiten für dieses Paket herausfinden, und die können sich im Laufe der Zeit ändern. Genannt Dependency Drift.

Container jedoch sperren den gesamten Stack aller Software-Abhängigkeiten ganz nach unten, sodass du etwas zuversichtlicher sein kannst, dass A, es funktionieren wird, und B, es reproduzierbar sein wird.

Wenn du also in der Lage bist, Docker oder Singularity oder Apptainer zu verwenden, würde ich das definitiv empfehlen.

Was jetzt schön ist, ist, dass die Moduldatei, die von der Pipeline-Entwickler\*in geschrieben wurde, jetzt sowohl Container als auch Conda hat, und so sagen wir der Person, die diese Pipeline ausführt, es ist uns egal, welche Software-Packaging-Lösung du verwendest. Es wird sowohl mit Docker als auch mit Conda funktionieren, und hier bekommst du die Software in beiden Fällen.

Wir können das Terminal hochziehen und lass uns das ausprobieren. Also Nextflow run hello config --batch conda. Und das erste Mal, wenn dies mit conda läuft, wird es bei diesem speziellen Prozess etwas langsam sein, weil es "conda install" ausführen muss.

Und es erstellt eine spezielle Conda-Umgebung nur für diesen einen Prozess. Es verwendet also nicht meine globale Conda-Umgebung, die ich auf meinem Terminal habe. Es erstellt eine nur für diesen einen Prozess. Das ist gut, weil es Dinge wie Abhängigkeitskonflikte zwischen verschiedenen Prozessen in deinem Workflow vermeidet. Wenn deine Prozesse Tools haben, die verschiedene Versionen von Python oder solche Dinge benötigen, ist das okay, weil sie verschiedene Conda-Umgebungen verwenden.

Nextflow cached diese Conda-Umgebungen lokal, du kannst sehen, es sagt dir, wo dieser Pfad ist, er ist im work-Verzeichnis hier. Und wenn ich dieses Skript das nächste Mal mit Conda ausführe, wird es viel schneller sein, weil es diese bestehende Conda-Umgebung finden und einfach wiederverwenden wird. Aber das erste Mal müssen wir es holen, auflösen, alle Abhängigkeiten herunterladen und alles einrichten.

Okay, großartig, es ist gelaufen. Wir können uns nur daran erinnern, wofür die Pipeline derzeit konfiguriert ist. Wenn wir in die Config-Datei schauen, war es "custom-outdir-config" gerade für mich. Schauen wir, ob ich zu diesem Basisverzeichnis gehe. Und ich habe --batch conda gemacht. Da ist unser conda-Unterverzeichnis. Es hat also funktioniert und da ist unsere cowpy-Ausgabe.

Es hat also cowpy geholt, es auf meinem lokalen System mit conda installiert und den Prozess ausgeführt. Und was großartig ist, ist, dass ich als Endnutzer\*in überhaupt nicht über irgendein Software-Management nachdenken musste. Nextflow hat es einfach für mich sortiert. Ich sagte, ich muss Conda auf diesem System verwenden. Die Pipeline-Entwickler\*in sagte, welche Pakete ich brauchte. Und Nextflow hat den Rest gemacht. Sehr mächtig.

Beachte, dass du tatsächlich eine Mischung verschiedener Technologien verwenden kannst. Ich kann also Docker für bestimmte Prozesse aktivieren und conda für andere Prozesse, oder sagen, dass einige Prozesse einfach welche lokale Software verwenden sollten, die ich installiert hatte. Das ist ziemlich ungewöhnlich, aber es ist möglich, und in einigen Fällen, zum Beispiel, wenn du bestimmte Software verwendest, die schwierig in Docker zu packen sein könnte, hast du einen Ausweg.

## 4. Eine Ausführungsplattform auswählen

Das ist also Software-Packaging. Der andere Teil der Portabilität auf andere Systeme ist, wo die Jobs tatsächlich laufen. Im Moment laufe ich im Grunde auf meinem Laptop oder in diesem Codespaces, was ein einzelner Computer ist. Es gibt nichts Ausgefallenes. Nextflow ist ein bisschen clever beim Parallelisieren der Jobs, so gut es kann, aber es ist alles auf einem System.

Wenn du nun auf einem HPC läufst, hast du wahrscheinlich eine Art Job-Scheduler wie SLURM oder PBS oder so etwas, und du reichst Jobs an diesen Scheduler ein und er verteilt alle Jobs an verschiedene Compute-Nodes.

Eine andere Art zu laufen ist in der Cloud. Vielleicht verwendest du AWS Batch oder Azure Cloud oder Google. Und diese funktionieren alle in einem ähnlichen System, wo du einen Scheduler hast und du Jobs einreichst und sie werden an verschiedene Orte zur Berechnung eingereicht.

Nun, in der fernen Vergangenheit, als ich anfing, Bioinformatik zu machen, war die Software aller Leute zum Ausführen von Analysen sehr an ihre Recheninfrastruktur gebunden, was es fast unmöglich machte zu replizieren.

Aber mit dieser Config-Trennung in Nextflow und mit Nextflows Fähigkeit, mit sehr vielen verschiedenen Compute-Infrastruktur-Backends zu interagieren, ist es sehr einfach, unsere Pipeline zu nehmen, ohne den Pipeline-Code überhaupt zu modifizieren, und das einfach auszutauschen.

## 4.1. Ein anderes Backend anvisieren

Wenn wir also zu unserer "nextflow.config"-Datei gehen, und wir können jetzt etwas Prozess-Level-Config einfügen. Wenn ich also oben den Prozess-Scope einfüge und ich den Executor setzen kann, und hier ist er auf local gesetzt, was der Standard ist.

Beachte, weil dies Prozess-Level ist, können wir Dinge auf verschiedene Prozesse ausrichten. Und so kannst du tatsächlich Executors so einrichten, dass sie prozessspezifisch sind und eine hybride Ausführung haben, wo einige Jobs lokal laufen könnten, wo auch immer der Nextflow-Job ausgeführt wird. Einige werden an verschiedene HPC eingereicht und einige könnten an die Cloud eingereicht werden. Du kannst so clever sein, wie du möchtest.

Nun, es ist sehr schwierig, dies in einer Trainingsumgebung wie dieser zu demonstrieren, weil ich keinen HPC habe, an den ich einreichen kann. Aber was ich tun kann, ist, wenn ich slurm eintippe, können wir ein bisschen schummeln und du kannst ein Gefühl dafür bekommen.

Und das ist wirklich am interessantesten für Leute, die es gewohnt sind, auf SLURM zu laufen und wissen, wie die SLURM-Header aussehen. Aber wenn ich Nextflow run mache, hello config. Es wird fehlschlagen, weil es versuchen wird, Jobs an einen Cluster einzureichen, der nicht existiert. Wir werden also irgendeine Fehlermeldung über sbatch bekommen, das nicht verfügbar ist.

Ja, geschrieben. Das ist das Tool. Das ist das CLI-Tool, das du verwendest, um Jobs an einen slurm-Cluster einzureichen. Aber was wir tun können, ist, wir können in unser work-Verzeichnis hier gehen, per Befehl klicken, dieses Verzeichnis öffnen und die .command.run anschauen. Und du kannst oben in der .command.run-Datei sehen, wir haben unsere sbatch-Header, die einem theoretischen SLURM-Cluster sagen, wie diese Job-Einreichung zu handhaben ist.

Und du kannst also sehen, dass Nextflow clever ist, es macht all die richtigen Dinge. Es ist nur so, dass wir keinen Cluster hatten, an den wir einreichen konnten.

## 5. Compute-Ressourcenzuweisungen steuern

Was ist noch anders zwischen verschiedenen Computing-Infrastrukturen? Eine andere Sache ist, wie viele verfügbare Ressourcen du hast, und tatsächlich ist es in vielen Compute-Umgebungen eine Anforderung, dass du angeben musst, wie viele CPUs und wie viel Speicher ein Job benötigt.

Wieder abstrahiert Nextflow das für uns, sodass es nicht mehr spezifisch für einen einzelnen Compute-Umgebungstyp ist, und wir können hier im Prozess-Level-Scope eingeben. CPUs equals eins, memory equals zwei Gigabyte. Unsere Pipeline ist nicht sehr anspruchsvoll, also sollte das in Ordnung sein.

Nun, ich habe diese Zahlen hier nur geraten, aber woher weißt du, was eine sinnvolle Menge an Ressourcen ist? Es ist eine ziemlich schwierige Aufgabe, all diese verschiedenen Prozesse einer großen Pipeline mit vielen Proben durchzugehen und zu verstehen, wie die Ressourcennutzung war.

Ein guter Ansatz dafür ist also, diese Werte zu Beginn auf hohe Zahlen zu setzen, nur damit deine Pipeline ohne Fehler läuft, und dann Nextflow zu bitten, einen Nutzungsbericht für dich zu generieren.

Das ist super einfach zu machen, also werde ich zurück zu einem Terminal gehen. Oh, ich muss mich daran erinnern, das zurück auf local zu setzen, damit meine Pipeline tatsächlich läuft. Und ich werde sagen Nextflow run, und ich werde ein Kommandozeilen-Flag -with-report verwenden.

Und ich kann das leer lassen und es wird einen Standard-Dateinamen geben, aber ich werde ihm einen spezifischen Dateinamen geben, damit das an einem spezifischen Ort gespeichert wird.

Enter drücken, und die Pipeline läuft genau wie normal, aber wenn sie fertig ist, wird sie einen schönen HTML-Bericht für mich generieren.

Also in der Seitenleiste hier habe ich diese HTML-Datei. Wenn ich das lokal laufen lassen würde, würde ich es einfach öffnen. Ich bin, weil ich in Codespaces bin, werde ich mit der rechten Maustaste darauf klicken und auf Download klicken, was es auf meinen lokalen Computer herunterladen wird. Und ich kann es einfach leicht im Webbrowser öffnen.

Nextflow kann einen Bericht wie diesen für jede Pipeline generieren und er hat einige wirklich schöne Informationen. Es ist also gute Praxis, diese Dinge immer zu speichern. Er sagt uns, wann wir gelaufen sind, wo wir gelaufen sind, ob es erfolgreich war oder nicht, welche Parameter verwendet wurden, was der CLI-Befehl war, solche Dinge.

Und es gibt auch diese Plots über Ressourcennutzung. Er sagt uns also, welcher Prozentsatz der CPU-Aufrufe für jeden Prozess als Box-Plot hier verwendet wurde, weil es viele Aufgaben für jeden Prozess gibt, sodass wir die Verteilung sehen können.

Du kannst unsere Prozesse hier sehen, cowpy und collectGreetings hatten nur eine einzige Aufgabe, also ist es nur eine einzelne Linie. Und wir haben sowohl CPU als auch Speicher und Job-Dauer, und sie waren sehr schnell.

Wenn du übrigens Seqera Platform verwendest, bekommst du dieselben Plots in die Platform-Oberfläche eingebaut, ohne etwas tun zu müssen. Du bekommst also immer diese Informationen zur Hand.

Okay, wir können also diesen Bericht verwenden und bei einem echten Lauf ein Gefühl dafür bekommen, wie viele CPUs und wie viel Speicher von unserer Pipeline verwendet werden, und zurückkommen und diese Werte zurück in unsere Config-Datei einfügen, sodass wir das nächste Mal vielleicht nicht ganz so viel anfordern. Und wir können etwas schlanker sein.

Nun, du kannst wirklich clever werden beim Konfigurieren von Pipeline-Config-Dateien. Und wieder, wenn du Seqera Platform verwendest, halte Ausschau nach einem kleinen Button, der wie eine Glühbirne aussieht. Denn wenn du darauf klickst, wird es eine hochoptimierte Config-Datei generieren, die speziell auf deine Daten, deinen Lauf und deine Pipeline zugeschnitten ist. Um sie auf die effizienteste Weise möglich auszuführen.

Aber für jetzt werde ich sagen, dass tatsächlich die Standard-Anzahl von CPUs, die Nextflow gab, in Ordnung war, und wir nur ein Gigabyte Speicher brauchen.

## 5.3. Ressourcenzuweisungen für einen bestimmten Prozess setzen

Nun, im wirklichen Leben ist es ziemlich ungewöhnlich, dass alle Prozesse in deiner Pipeline dieselben Anforderungen haben werden. Du könntest etwas wie MultiQC als Reporting-Tool haben, das sehr wenig in Bezug auf Ressourcen benötigt und ziemlich schnell läuft.

Und dann hast du vielleicht etwas, das ein Referenzgenom indiziert oder ein Alignment macht oder einen anderen Job macht. Es ist egal, was es ist, das viele Ressourcen benötigt. Und für diese verschiedenen Job-Einreichungen an einen Scheduler möchtest du also verschiedene Mengen an Ressourcen geben.

Unter diesem Prozess-Scope können wir eine Config definieren, die spezifische Prozesse auf verschiedene Arten anvisiert.

Hier verwenden wir withName, wir können auch Labels verwenden, und diese können ein Muster verwenden, um einen oder mehrere Prozesse anzuvisieren. Hier sagen wir nur, alle Prozesse, die einen Namen cowpy haben, setzen auf zwei Gigabyte Speicher und zwei CPUs, und weil dies ein spezifischerer Selektor ist als der Top-Level-Prozess, wird dies in diesen Fällen überschrieben, sodass du hier eine schöne Config-Datei aufbauen kannst, die wirklich alle deine verschiedenen Prozesse in deiner Pipeline anpasst, um sie wirklich effizient zu machen.

## 5.5. Ressourcenlimits hinzufügen

Nun, als Pipeline-Entwickler\*in kenne ich die Tools wahrscheinlich ziemlich gut, und ich möchte, dass alles so schnell und so gut wie möglich läuft. Es könnte also sein, dass ich ziemlich hohe Zahlen für einige davon eingebe, weil ich weiß, dass es viel schneller laufen wird, wenn ich cowpy 20 CPUs anstelle von zwei gebe.

Das ist in Ordnung, bis du versuchst, auf deinem Laptop oder auf GitHub Actions Continuous Integration Test oder einem anderen System zu laufen, das vielleicht keine 20 CPUs verfügbar hat.

Wenn du jetzt versuchst, die Pipeline auszuführen, wird sie abstürzen, weil Nextflow sagen wird, ich kann diesen Job nirgendwo einreichen. Ich habe die verfügbaren Ressourcen nicht.

Um diesen harten Absturz zu vermeiden, können wir jetzt etwas mehr Config hinzufügen, die spezifisch für unser System ist, genannt Ressourcenlimits. Und das sieht so aus. Es ist wieder unter dem Prozess-Scope.

Und Ressourcenlimits, du kannst im Grunde die Obergrenze dessen angeben, was du verfügbar hast. Es ist hier eine Map, und du kannst innerhalb dieser Map den Speicher, die CPUs und die Zeit setzen.

Was jetzt passiert, ist, wenn Nextflow eine Aufgabe von einem Prozess einreicht, schaut es, was angefordert wird, und es macht im Grunde nur ein Minimum zwischen dem und dem. Wenn wir also 20 CPUs angefordert haben, aber nur vier verfügbar sind, wird es vier anfordern. Die Pipeline stürzt nicht ab und sie verwendet so nah wie möglich an dem, was von der Pipeline-Entwickler\*in entworfen wurde.

## 6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln

Okay. Ich sagte, dass die Ressourcenlimits hier systemspezifisch sein könnten, und vielleicht habe ich eine Nextflow-Config-Datei in meiner Pipeline, und ich weiß, dass Leute dies an einer Reihe verschiedener Orte verwenden werden. Anstatt jetzt alle zu zwingen, jedes Mal ihre eigene Nextflow-Config-Datei zu erstellen, kann ich verschiedene Voreinstellungen von Konfiguration zusammen in Config-Profile gruppieren.

Ich werde hier ein bisschen nach unten scrollen und dann einfach nach params, weil die Reihenfolge der Config-Datei hier wichtig ist, die Config-Datei wird sequenziell geladen, also werde ich diese Profile nach allem anderen einfügen, damit es die zuvor definierten params überschreibt. Und ich werde diese Profile vom Trainingsmaterial einfügen.

Es gibt also einen neuen Top-Level-Scope namens profiles. Wir können hier beliebige Namen haben. Wir haben also my_laptop und univ_hpc. Und hier können wir sehen, wir setzen die anderen gleichen Config-Parameter, die wir vorher hatten. Jetzt nur innerhalb eines Profils. Wir haben also einen lokalen Executor zum Laufen auf my_laptop und ich reiche an einen SLURM-Cluster auf dem HPC ein.

Ich verwende Docker lokal, conda auf dem HPC, und das HPC-System hat viel höhere Ressourcenlimits.

Jetzt kann ich die Pipeline mit der -profile CLI-Option ausführen, sagen, welches Profil ich verwenden möchte. Ich werde also my_laptop verwenden, und Nextflow wird alle Config innerhalb dieses Profil-Scopes anwenden. Ich kann das jetzt ausprobieren. Es ist derselbe Befehl wie zuvor. Nextflow run hello config, und ich mache dash profile, einzelner Bindestrich, weil es die Kern-Nextflow-Option ist, dash profile my_laptop.

Es wird jetzt diese Config-Option batch-weise anwenden. Oh, und du kannst sehen, ich sagte vorher, das könnte passieren, dass die Prozessanforderung, es nach vier CPUs fragte und ich nur zwei auf dieser Codespaces-Instanz habe.

Das ist also eine gute Gelegenheit, einfach die Prozess-Ressourcenlimits auszuprobieren und zu sagen, dass ich nur zwei CPUs auf my_laptop habe, oder in diesem Codespaces. Wenn wir es jetzt nochmal ausführen, sollte es diese Anforderung auf zwei begrenzen und hoffentlich wird die Pipeline laufen. Großartig.

## 6.2. Ein Profil von Testparametern erstellen

Beachte, dass diese Profile nicht nur Konfiguration über ihre Infrastruktur haben müssen. Du kannst hier Gruppierungen jeder Config haben, einschließlich Parameter.

Eine andere Sache, die du also sehr oft in den Pipelines der Leute sehen wirst, ist ein Test-Profil, das Parameter enthält, die du normalerweise pro Nutzer\*in einreichen würdest. Aber hier haben wir im Grunde verschiedene sinnvolle Standardwerte für wenn ich Testfälle ausführen möchte.

Und das ist großartig, weil ich nicht unbedingt all diese Dinge angeben muss, die möglicherweise erforderliche Parameter sind. Ansonsten kann ich einfach dash profile test sagen und es wird einfach out of the box laufen.

Etwas zu beachten ist, dass Profile auch mehr als eins kombiniert werden können. Ich kann also hier profile my_laptop machen, und dann auch test hinzufügen. Ich mache nicht zweimal profile. Ich mache nur eine kommagetrennte Liste hier ohne Leerzeichen. Und es wird diese Profile in Reihenfolge anwenden. Es wird also die Config vom my_laptop-Profil nehmen, und dann wird es die Test-Config oben drauf anwenden.

Wirklich praktisch und du kannst sehen, wie du hier viele sinnvolle Standardgruppen einrichten kannst, um es einfach zu machen, deine Pipeline auszuführen.

## 6.3. nextflow config verwenden, um die aufgelöste Konfiguration zu sehen

Hoffentlich habe ich dich überzeugt, dass Nextflow-Config-Auflösung mächtig ist, aber ich würde es dir nicht verübeln, wenn du an diesem Punkt ein bisschen schielst, nachdem ich etwa 20 verschiedene Wege gesagt habe, Config bereitzustellen und all diese verschiedenen Schichten wie eine Zwiebelschale zu geben.

Wenn du dir also jemals unsicher bist, was die endgültige aufgelöste Config für Nextflow ist, wisse, dass es einen Befehl namens "nextflow config" gibt, und wir können das ausführen und es wird uns sagen, was die aufgelöste Konfiguration an unserem aktuellen Standort ist.

Wenn ich es also hier ausführe, findet es die "nextflow.config"-Datei im aktuellen Arbeitsverzeichnis, und es verarbeitet alle verschiedenen Config, und es gibt mir die aufgelöste Ausgabe.

Beachte, dass die Nextflow-Config-Datei auch die profile CLI-Option nehmen kann. Wenn ich ihr also sage, in my_laptop- und test-Profilen aufzulösen, und du kannst sehen, es hat auch die Ressourcenlimits hier von der my_laptop-Config-Option angewendet und auch die params gesetzt, die im Test waren.

Das ist also eine schöne Art, einfach zu erkunden, wie die Config-Auflösung funktioniert, wenn du dir überhaupt unsicher bist.

## Abschluss

Okay, das war's. Das ist Nextflow-Config in einer Nussschale. Du kannst viel mit Config machen. Es ist wirklich mächtig. Aber das sind die meisten der häufigen Anwendungsfälle, die du vorfinden wirst, und diese Konzepte gelten für alle verschiedenen Optionen.

Klopf dir selbst auf die Schulter, denn das ist das Ende des Hello Nextflow Trainingskurses. Du bist hoffentlich jetzt zuversichtlich sowohl beim Schreiben deiner eigenen Nextflow-Pipeline von Grund auf, beim Konfigurieren und Ausführen, und du kennst all die Ins und Outs und die Dinge, auf die du achten musst.

Es gibt noch ein Quiz, das du auf der Config-Trainingsseite ausprobieren kannst. Geh also runter und probiere das aus und stelle sicher, dass du all diese Teile über die Config verstanden hast.

Und schließ dich uns im letzten Video an, nur für einen kurzen Abschluss über einige der nächsten Schritte, die nach diesem Trainingskurs gut zu tun sein könnten.

Danke, dass du bei uns geblieben bist. Gut gemacht und ich sehe dich im nächsten Video.
