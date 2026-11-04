# Part 3: Hello Workflow - Transcripció del vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importants"

    Aquesta pàgina mostra només la transcripció. Per a instruccions completes pas a pas, torneu al [material del curs](../03_hello_workflow.md).

    Els números de secció mostrats a la transcripció es proporcionen només amb finalitats indicatives i poden no incloure tots els números de secció dels materials.

## Benvinguda i resum

Hola, i benvinguts de nou a la part tres de Hello Nextflow. Aquesta part s'anomena Hello Workflow, i és en aquesta part del curs on realment comencem a justificar el nom de pipeline o workflow.

Agafarem el nostre script de pipeline simple fins ara amb el seu únic procés, i començarem a afegir processos addicionals i veurem com Nextflow gestiona aquesta orquestració i el flux de dades a través del pipeline.

Tornem als nostres code spaces. Veureu que he esborrat tots els meus directoris .nextflow\* i els directoris de treball i tot per intentar mantenir-ho net. No us preocupeu si encara teniu aquests fitxers de les parts anteriors del curs.

Treballarem a partir d'un fitxer anomenat hello-workflow.nf. Com abans, això bàsicament representa l'script que hem construït fins aquest punt, i ens dóna un punt de partida net. I de nou, a la sortida podem veure que el camí ara és hello_workflow. Així que els fitxers publicats haurien d'anar a un subdirectori diferent a la vostra carpeta de resultats.

Per recapitular on som fins ara, tenim un únic procés aquí, amb una entrada greeting, una sortida greeting file. I després el simple script Bash, que només fa una comanda echo a un fitxer.

Tenim una única entrada de workflow, el bloc params aquí, on diem que espera un path, i el valor per defecte és data/greetings.csv, que és aquest fitxer d'aquí dalt.

Després al workflow mateix, tenim un bloc main. Estem creant un canal. Estem analitzant el CSV en files i després prenent el primer element de cada array, i estem passant aquest canal a aquest procés, que després genera tres tasques, i estem publicant des del workflow, les sortides d'aquest procés.

I finalment, al bloc output, estem dient a Nextflow que publiqui aquests fitxers d'aquest canal al directori anomenat hello_workflow. I que copiï aquests fitxers en lloc de fer enllaços simbòlics.

## 1. Afegir un segon pas al workflow

D'acord, en aquesta part afegirem un segon procés al nostre workflow. Agafarem les sortides del procés sayHello, i les processarem en un segon pas, que convertirà totes les lletres dins d'aquests fitxers convertToUppercase.

Això és només un exemple senzill, és només un processament de cadenes simple de nou, però us mostra com podem prendre la lògica, dins del workflow.

Farem servir una comanda bash anomenada "tr" per això, que és l'abreviatura de translate. És una comanda Unix que existeix des de sempre. Si no la coneixeu, no us culpo. No crec que l'hagués utilitzat mai abans de la formació, però podeu provar-la molt ràpidament al terminal. Si faig "echo 'hello world'" i després pipe a 'tr' i després entre cometes dieu rang de caràcters, així que A a Z, minúscules, i després voleu fer A a Z majúscules. I només diu, tradueix aquestes lletres a aquestes lletres.

I quan premo enter, podeu veure que ara ha posat tot en majúscules. Molt bé si us agrada cridar a la gent.

Així que aquest és un estil molt simple de comanda bash que utilitzarem al nostre segon procés.

## 1.2. Escriure el pas de conversió a majúscules com a procés Nextflow

Així que si torno al meu script, faré una mica de trampa i simplement copiaré el codi de la documentació de la formació. Però podeu veure exactament què està passant.

Tenim un nou procés aquí. Aquest l'hem anomenat convertToUpper, però podríem anomenar-lo com vulguem.

Tenim una única entrada path, com vam fer abans. No és un canal de valor, és un canal de path. I després una única sortida.

Al bloc script fem "cat" al fitxer d'entrada. I podem posar això entre claus si volem. i que pren aquesta variable. I executem la mateixa comanda bash al pipe i escrivim els resultats a un fitxer amb aquest nom de fitxer, i això és recollit pel path de sortida.

Ara necessitem fer alguna cosa amb aquest nou procés. Així que anirem al workflow on construïm la lògica diferent d'un workflow, i després d'aquest primer procés, executarem el nostre segon procés. Així que convertToUpper és el nom del procés aquí.

Pren una entrada així que no podem simplement cridar-lo per si mateix. Volem processar la sortida del primer procés. Així que igual que vam fer amb això, sayHello out on estem publicant aquests resultats. Volem utilitzar aquests mateixos resultats aquí com a entrada, així que podem copiar-los i posar-los allà.

Volem el procés sayHello ".out", i Nextflow sap que això significa un simple registre de sortida única aquí, que és aquest fitxer. Així que això després es passarà com a entrada a un segon procés.

## 1.5. Configurar la publicació de sortida del workflow

D'acord. I finalment, perquè realment desem els resultats d'aquest segon procés, també necessitem publicar-los des del workflow, i després definir-los al bloc output, mateixa sintaxi que abans. Així que podem copiar això i dir second outputs, o com vulgueu anomenar-ho.

Preneu el nom del procés que ens interessa, convertToUpper out, i després aquí baix al bloc output. Afegiu això i podríem fer els mateixos atributs aquí. Així que també volem aquests fitxers al subdirectori Hello Workflow, i també volem copiar-los.

Genial. Provem d'executar-ho. Així que si obro el terminal i faig "nextflow run hello-workflow.nf", i veurem què fa. Veurem si es veu diferent de les parts anteriors.

Així que llança Nextflow. A la documentació, diu que ho feu amb "-resume", però vaig esborrar tot el meu directori de treball, així que no hauria fet cap diferència aquí. Però si ho vau fer, llavors això també funcionarà.

I es veu gairebé exactament igual. Però podeu veure ara que hi ha una segona línia de sortida aquí, on podeu veure el nom del segon procés que acabem d'afegir. I efectivament, podeu veure que s'ha executat tres vegades amb èxit.

Brillant. Si tingués els meus directoris de treball anteriors al voltant i hagués fet això amb "-resume", aquests haurien estat, emmagatzemats en memòria cau només el primer pas del pipeline. Perquè aquestes sortides eren exactament les mateixes, així que Nextflow hauria sabut reutilitzar-les de nou.

I així podeu veure com podeu utilitzar -resume per construir iterativament el vostre workflow, pas a pas, si ho necessiteu.

D'acord, fem una ullada al directori de resultats aquí dalt i vegem si ha funcionat. Podem veure que tenim alguns fitxers més aquí dalt. Tenim els nostres fitxers originals com vam fer abans del primer procés. I efectivament, tenim els nostres fitxers upper i les lletres estan totes en majúscules, així que ha funcionat. És realment bonic de veure.

També és interessant només comprovar dins d'aquests directoris de treball. Com abans el hash aquí correspon als directoris de treball. Així que si miro a "ls work", i després expandeixo això, veurem els diferents fitxers aquí.

Veiem el fitxer de sortida del primer procés, que s'ha introduït aquí com a entrada. I podem veure el nou fitxer de sortida que es va generar.

Ara si faig això amb "-la" per llistar i mostrar tots els fitxers, veurem algunes coses més. Primer, veureu que aquest fitxer és en realitat un enllaç simbòlic al primer procés. Això és bàsicament sempre un enllaç simbòlic si pot ser, per estalviar espai de fitxer. No estem publicant els fitxers aquí i només fa referència a aquest fitxer d'una primera tasca a una segona tasca perquè tot estigui encapsulat dins d'aquest directori de treball, i segur i aïllat de tot el demés.

I això ha d'estar allà perquè si mirem el fitxer .command.sh, així que si faig "cat work/b8/56\*", podeu veure que les parts de fitxer aquí són relatives, així que està fent cat d'aquest fitxer d'entrada, que s'ha enllaçat simbòlicament al mateix directori de treball.

Així que així és com es veurà cada directori de treball. Quan el mireu a Nextflow, tindreu tots els fitxers d'entrada allà preparats a aquest directori de treball. I després també tindreu qualsevol fitxer de sortida que es va crear. Així que això és genial. Això es veu com esperàvem.

## 2.1. Definir la comanda de recollida i provar-la al terminal

D'acord, tornem al nostre workflow. Quin és el següent pas que volem fer?

Tenim dos processos ara i estan prenent aquest fitxer CSV, analitzant-lo i dividint-lo. I després tenim tres tasques per a cadascun d'aquests processos i Nextflow gestiona la paral·lelització de tot això, així que tot s'executa costat a costat quan és possible.

Aquesta manera de dividir el treball per executar coses en paral·lel és molt comuna. I l'invers d'això és després recollir tot de nou. Així que això és el que farem amb el nostre procés final al workflow és que en tindrem un tercer aquí, que pren aquestes tres sortides diferents i les combina totes en un únic fitxer.

Podem fer això força simplement en un terminal, només per tenir una idea de com es veurà això.

Si vaig a la carpeta de resultats. Així, "cd results/hello_workflow/", i tenim tots els fitxers UPPER aquí. Puc simplement utilitzar "cat", que utilitzem per imprimir el contingut d'aquest fitxer, i podeu donar múltiples fitxers a "cat" i llegirà un després de l'altre.

Així que puc dir "UPPER-\*", que em dóna la mateixa llista de tres noms de fitxer amb expansió Bash. I puc dir combined.txt. Crec que a la documentació, llista els noms de fitxer exactes, però està fent el mateix.

Ara, si utilitzo "cat combined.txt", podem veure que tenim el contingut del fitxer dels tres fitxers.

Així que això és bàsicament tot el que aquest procés farà és que intentarem donar-li tots els diferents fitxers de sortida d'un procés anterior en una única tasca de procés, i després farem "cat" junts i desarem el fitxer de sortida.

## 2.2. Crear un nou procés per fer el pas de recollida

D'acord, així que afegim el nostre nou procés. Enganxaré això dels materials de formació, i podeu veure que ens ha deixat una mica d'exercici per al lector aquí amb aquests signes d'interrogació. Però podeu veure l'esquema general del procés és bàsicament el que acabem de fer al terminal, on estem fent "cat" d'un munt de fitxers d'entrada i escrivint-lo a un fitxer de sortida aquí anomenat collected, i després la sortida espera aquest únic path de nou.

Així que necessitem algun tipus d'entrada aquí i seran un conjunt de paths. Així que de nou, definim un canal d'entrada path i anomenem-lo input_files. Ara, això prèviament ens ha donat un únic path aquí, però un path també pot tenir múltiples fitxers aquí, encara que sigui una única declaració.

Copiaré això aquí baix perquè volem fer "cat" d'aquests fitxers. I podríeu pensar que tenim alguns problemes aquí amb imprimir un array o coses així, però Nextflow és generalment força sensible quan es tracta d'això. I si se li dóna un canal amb múltiples fitxers com aquest, els posarà tots junts amb separadors d'espai. Així que això ens donarà la sintaxi correcta.

Això és genial. Així que ara connectem el nostre nou procés. Vaig al workflow. Diré combine the outputs, el nom del nou procés, i igual que abans. Prendré aquest procés anterior, convertToUpper i faré ".out".

Genial. Provem-ho i vegem si funciona al terminal. Si només torno un parell de directoris enrere i després torno a executar la comanda Nextflow, i veurem què passa.

Així que el workflow s'ha llançat i ara podeu veure que tenim tres noms de procés diferents, que és genial. Els dos primers es veuen igual que abans, i el tercer nou s'executa, que és bo.

No obstant això, hi ha alguna cosa una mica estranya aquí. Volíem combinar aquests fitxers de sortida en un únic fitxer, i tanmateix aquest procés podem veure que s'ha executat tres vegades, no una.

Efectivament, si anem a un d'aquests directoris de treball. I fem "cat work/" "collected", llavors veurem. Només hi ha una única paraula aquí, no tres.

I així el que ha passat és que Nextflow ha continuat aquesta paral·lelització igual que va fer als passos anteriors. I aquest procés ens va donar un canal amb tres elements, i aquests tres elements de canal es van passar al nostre procés posterior, que va generar tres tasques de procés.

Bàsicament va intentar recollir tres vegades separades i cada vegada només tenia un únic fitxer, així que només va fer cat d'un únic fitxer a una sortida, i de fet, podem veure això al fitxer .command.sh també.

Si faig .command.sh, podem veure que només té un únic nom de fitxer aquí i només un únic fitxer es va preparar a aquest directori de treball.

## 2.3. Afegir el pas de recollida al workflow

Així que d'alguna manera necessitem dir a Nextflow que reuneixi totes aquestes sortides d'un procés anterior i les doni a aquest procés posterior com un únic element de canal, en lloc de tres.

Ho fem amb un operador de canal anomenat _collect_.

Aquest és un operador súper útil, que veureu als pipelines Nextflow tot el temps. Aquest és un canal aquí, aquest canal de sortida, igual que el que vam crear a la part superior. I així podem afegir operadors de canal igual que vam fer abans. Podem simplement fer punt, i després en aquest cas, collect, parèntesis.

I això és tot el que necessitem. Això després manipularà aquest canal abans que es passi a aquest procés.

Si voleu veure què li està passant, també podem veure-ho aquí. Així que aquí, això no està relacionat amb executar aquest procés en absolut, així que podria posar-ho en qualsevol punt després d'executar aquest procés. Però prenem el mateix canal de sortida, i l'estem mirant amb .view, i després l'estem mirant de nou amb .collect.view.

I quan executem això, ens mostrarà les dues estructures diferents d'aquest canal, abans i després de collect. Així que provem-ho ara. D'acord, acabo d'allunyar una mica perquè algunes de les sortides són força llargues, però si executo el pipeline, veurem si funciona.

Espero que un tercer procés s'executi només una vegada, perquè està recollint les sortides i efectivament, podeu veure collectGreetings com un d'un. Així que això ha executat només una tasca.

I després si mirem les declaracions view, tenim tres declaracions view per als tres elements d'abans, amb un camí de fitxer a cadascun.

I després d'aquesta declaració collect, això només s'ha activat una vegada perquè hi ha un únic element en aquest canal. I ara tenim aquesta llista de tres camins de fitxer diferents.

Això és exactament el que esperàvem. I podeu veure esperançadament, això és bàsicament l'invers d'aquest operador "map" que vam fer per anar dels arrays CSV a elements de canal separats. Ara estem prenent elements de canal separats i tornant-los a posar en un únic array.

Genial, podem netejar aquestes declaracions view. Ja no les necessitem. Podem passar al següent pas.

Abans d'anar més lluny, i abans que m'oblidi, afegiré una nova declaració publish aquí. Third output. Podeu anomenar això alguna cosa més semàntica i descriptiva al vostre workflow. I després afegiré això al bloc output de nou i diré path 'hello_workflow' mode 'copy'. Només perquè el fitxer de sortida generat per aquest procés es desi a la nostra carpeta de resultats aquí dalt.

Només per comprovar ràpidament que funciona. Hauria de ser una mica més net ara perquè no tenim aquestes declaracions view. I veurem si obtenim el nostre nou fitxer de sortida aquí dalt. Una d'una tasca executada, tenim un nou fitxer anomenat collected, i ara tenim totes tres paraules. Fantàstic. Què segueix?

## 3. Passar paràmetres addicionals a un procés

D'acord. A continuació mirarem de gestionar múltiples entrades en un únic procés. Fins ara podeu veure que tots els nostres processos només estan prenent una cosa com a entrada. Tots tenen una única línia sota la seva entrada.

Ho demostrarem permetent que Nextflow especifiqui un identificador de lot diferent perquè potser executeu aquest workflow múltiples vegades i podeu donar-li un ID de lot diferent cada vegada.

Simplement afegiré una segona línia a l'entrada aquí per a collectGreetings. I l'anomenaré "val", perquè això és una cadena. Ara és un valor, no un path, i l'anomenaré "batch_name".

Després editaré l'script aquí baix per utilitzar aquesta variable, i intentaré posar-la al mateix lloc que el material de formació. Així que la poso al mig d'aquest camí de fitxer COLLECTED-$\{batch_name\}-output.

Encara no hem acabat. Recordeu que hem de dir a Nextflow quins seran els noms dels fitxers de sortida. Així que també hem de fer el mateix aquí dalt: COLLECTED-$\{batch_name\}-output.txt".

Fantàstic. Nextflow ara està rebent una segona entrada de variable i l'està interpolant a l'script i la sortida.

Una última cosa, ara hem de trobar on s'està cridant això, i hem de passar la segona entrada al procés. Això és igual que qualsevol altra entrada a una funció en qualsevol altre llenguatge.

Igual que vam fer abans a la formació, utilitzaré l'especial "params" aquí, i l'anomenarem "params.batch" perquè puguem tenir una opció CLI --batch. I ara podeu veure que el nostre procés aquí té dues entrades separades només separades per comes, que s'estan passant.

És realment important aconseguir l'ordre correcte, així que l'ordre dels arguments aquí per a channel i després el param ha de coincidir. El channel i el batch name allà. Això és només coincidència posicional.

D'acord. Puc executar aquest pipeline ara directament amb --batch, però primer fem el correcte i definim-ho a l'entrada aquí a Params. Així que l'afegiré a batch i després direm que és una cadena i donem-li un valor per defecte. Així que simplement anomenem-lo batch. D'acord? Ara provem d'executar el workflow.

--batch Trio. Crec que diu al material de formació, però podríem utilitzar qualsevol cadena que vulguem allà. I esperançadament veurem que el fitxer de sortida de resultats apareix aquí.

I efectivament, COLLECTED-trio-output - això ha funcionat correctament. Ha reanomenat el nostre fitxer. I podeu imaginar ara que això és útil perquè si executo això de nou amb un nom de lot diferent, com replicate_two, llavors ens donarà un nom de lot diferent aquí dalt.

I i llavors no sobreescriurà els fitxers de sortida en aquest cas. Així que això és bonic.

## 4. Afegir una sortida al pas de recollida

D'acord, així que ara tenim múltiples entrades al nostre procés aquí. Però què passa si volem crear múltiples sortides? El nostre exemple aquí és que crearem un informe per a aquest procés, només dient quants fitxers es van recollir.

I ho farem amb una comanda echo aquí. Així que podem dir echo. There were, copiaré això del material de formació, perquè no hàgiu de veure'm escriure-ho.

There were $\{count_greetings\} greetings in this batch, i desar això a un nou fitxer ara anomenat $\{batch_name\}, així que mateixa variable, podem reutilitzar-la tantes vegades com vulguem, report.txt.

## 4.1.1. Comptar el nombre de salutacions recollides

Necessitem calcular això d'alguna manera. Podríem fer aquesta lògica a l'script Bash si volguéssim, utilitzant lògica Bash. No obstant això, també podem simplement fer scripting directament dins del codi Nextflow, sempre que estigui dins del bloc script al procés i per sobre de la secció entre cometes.

Qualsevol cosa aquí no s'inclourà a l'script renderitzat final, i només serà executat per Nextflow quan renderitzi una tasca.

Així que aquí només estem fent una mica de lògica. Estem creant una nova variable anomenada count_greetings. Prenem el canal input files aquí, i estem cridant .size() sobre ell.

D'acord, aquesta funció em donarà un número aquí a aquesta variable, i ara el nostre avís ha desaparegut perquè aquesta variable s'està definint.

D'acord, així que estem creant aquest segon fitxer al directori de treball, però necessitem dir a Nextflow que l'esperi com a sortida publicada d'aquest procés. Així que ho fem amb exactament la mateixa sintaxi que vam fer per al primer fitxer.

Diem path perquè és, de nou, podríem estar publicant una variable aquí si volguéssim amb "val", però direm "path". I després el nom de fitxer esperat. Noteu que no està ressaltat aquí. Això és perquè vaig utilitzar cometes simples. He d'utilitzar cometes dobles.

## 4.1.2. Emetre el fitxer d'informe i nomenar sortides

D'acord, això és genial. I ara podríem començar a accedir a aquestes sortides aquí baix igual que vaig fer aquí. Però ara és un array d'objectes diferents, així que podria fer collectGreetings.out[0] per obtenir el primer, o un per obtenir el segon, que és el nostre nou informe.

Però realment no m'agrada fer això gaire perquè és força fàcil equivocar-se amb el recompte d'índex. I et sents allà comptant línies molt i afegeixes una nova sortida i de sobte tot es trenca. Així que

és molt més bonic referenciar tot per nom en lloc d'això. I podem fer això amb una clau especial aquí anomenada "emit".

Així que podem anomenar això com vulguem. Anomenem-ho emit outfile, i emit reports. Si definiu aquests i podeu fer-ho en un o molts, depèn de vosaltres. Ara puc anar aquí baix i en lloc d'això puc anar punt out punt reports i només cridar-ho per nom, que és molt més fàcil d'entendre el vostre codi quan el llegiu, i és més segur davant de canvis al codi.

He afegit el .out.report aquí, però en realitat necessito tenir dues sortides diferents que s'estan publicant. Així que reanomenaré com alguna cosa més interessant com collected i report i és això el que vaig anomenar? Vaig anomenar-ho out file, perdó. Així que aquest nom emit aquí outfile i report. perquè estem publicant dos canals de sortida diferents i així necessitem referenciar-los tots dos al bloc publish.

Després també necessitem definir-los al bloc output. Així que vaig reanomenar això collected, i de nou, per a reports, una mica verbós aquí, però és realment útil quan entres a llegir un nou workflow, veure totes les diferents sortides aquí, tots els diferents canals llistats costat a costat, i hi ha maneres de fer això menys verbós, que tocarem més endavant.

D'acord, provem-ho i executem el nostre workflow i vegem què passa.

Esperançadament ara hauria d'executar-se bàsicament igual que abans. I obtindrem un nou fitxer de sortida aquí dalt anomenat replicate_two, report. I aquí està. S'ha obert i diu que hi ha tres salutacions al lot, que és el que esperàvem, així que és perfecte.

Si vaig al directori de treball aquí només per demostrar-vos que es va executar al codi Nextflow, en lloc de l'script bash, puc anar a cat work/ command.sh, i veureu aquí que només està fent echo d'aquesta cadena directament. There were three greetings in this batch, i així aquesta variable va ser interpolada per Nextflow. Es va calcular al bloc script abans que escrivís el fitxer .command.sh. Així que el càlcul de la variable resultant està bàsicament codificat en dur en això abans que s'executi al vostre entorn de càlcul en aquest cas.

I així podeu veure aquesta separació entre l'script. Bloc aquí i qualsevol cosa per sobre d'ell. Espero que tingui sentit.

## Conclusió i qüestionari

D'acord, això és el final d'aquesta part de Hello Nextflow. Així que com abans, aneu i feu una ullada al qüestionari. Feu-ho a la pàgina web o al CLI, repasseu algunes de les preguntes i només comproveu que heu entès part del material que hem cobert. Vegeu si hi ha alguna cosa allà que destaqui alguna cosa que potser no hàgiu entès. No gaires preguntes. Bonic i fàcil de fer. O podeu fer-ho a la pàgina web aquí baix també.

I feu una petita pausa, una petita volta i torneu i uniu-vos a nosaltres a la part quatre de Hello Nextflow, on parlarem de mòduls. Moltes gràcies.
