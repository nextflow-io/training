# Part 2: Hello Channels - Transcripció del vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importants"

    Aquesta pàgina mostra només la transcripció. Per a instruccions completes pas a pas, torneu al [material del curs](../02_hello_channels.md).

    Els números de secció mostrats a la transcripció es proporcionen només amb finalitats indicatives i poden no incloure tots els números de secció dels materials.

## Benvinguda

Hola i benvinguts de nou a la Part 2 de Hello Nextflow. Aquest capítol s'anomena Hello Channels.

Els canals són com l'adhesiu del vostre pipeline de Nextflow. Són les parts que mantenen units tots els diferents processos, que Nextflow utilitza per passar tota la informació i orquestrar el vostre workflow.

Hi ha una altra part dels canals que són els operadors. Aquests són bàsicament funcions que podem utilitzar sobre els canals per modificar-ne els continguts. Submergim-nos en VS Code i vegem on som.

Estic molt ampliat en aquest VS Code, així que per mantenir les coses netes i ordenades, he eliminat tots els fitxers _.nextflow\*_ i el directori _work/_ i el _results/_ i tot del Capítol U. I estic començant de nou aquí. Però no us preocupeu massa per això. Si no voleu, podeu deixar aquests fitxers. No causaran cap problema.

Començarem treballant amb _hello-channels.nf_ per a aquest capítol, i si obro això, hauria de semblar molt similar al fitxer amb el qual estàvem treballant anteriorment. Pot ser que diferents parts estiguin en diferents parts de l'script, però tot hauria de ser bàsicament el mateix.

Una cosa que és diferent és que el path al bloc output aquí ara és _hello_channels_ per a aquesta part, el que significa que els fitxers de resultats s'emmagatzemaran en un subdirectori diferent als vostres resultats si encara ho teniu allà. Així que hauria de ser un lloc net i agradable per començar sense confondre's amb les sortides.

D'acord, així que recordem ràpidament què fa aquest script quan executem aquest workflow. Fem _"nextflow run hello-channels.nf"_. Podem fer _"--input myinput"_, i quan executem això, utilitzarà aquest paràmetre, params.input, que es va passar com a variable per al procés sayHello aquí dalt, que va a greeting i es desa a output.txt. I podem veure això al fitxer de resultats. Genial.

## 1. Proporcionar entrades variables mitjançant un canal explícitament

Això és bonic. Però és, és força simplista. Tenim una variable en aquest paràmetre, que va a un procés que s'executa una vegada, i no escala realment. I no podem donar-li molts fitxers diferents per crear aquí. No podem donar-li moltes salutacions diferents. Només en tenim una.

En realitat, Nextflow tracta d'escalar la vostra anàlisi. Així que probablement voleu que faci més d'una cosa. I ho fem amb _canals_.

Els canals són un concepte una mica únic per a moltes persones que comencen amb Nextflow. Ve d'aquest tipus de conceptes de programació funcional, i pot trigar una mica de temps a entendre'l, però una vegada ho cliqueu, realment desbloquegen el poder de Nextflow i és clau per a com escriviu els vostres workflows.

## 1.1. Crear un canal d'entrada

Comencem prenent aquest script i fent que utilitzi un _canal_ en lloc de només un _param_.

Anem al workflow, que és on està tota la nostra lògica de workflow sobre encadenar coses. I entraré aquí i crearé un nou canal.

Crear un nou canal.

I l'anomenaré "_greeting_ch"_. Aquesta és la convenció de fer "_\_ch"_ així, només perquè pugueu recordar que aquesta variable és un canal. Però podeu anomenar-lo com vulgueu.

I després diré equals, i faré _"channel.of"._

Channel és com l'espai de noms per a tot el que té a veure amb canals. "c" minúscula si heu estat utilitzant Nextflow abans. I el _".of"_ és quelcom anomenat Channel factory, que és bàsicament una manera de crear un canal.

Hi ha moltes Channel factories diferents. Si faig només "." aquí, podeu veure que VS Code està suggerint moltes d'elles, però _".of"_ és la més simple i només pren una entrada aquí.

Així que puc fer uns parèntesis i diré _"Hello Channels!"_.

Genial. Tinc un canal. Fantàstic. Puc prémer desar, podria executar-lo de nou, però no passarà res interessant. VS Code m'ha donat una línia d'advertència taronja aquí sota i m'ha dit que això està configurat: heu creat això, però mai l'heu utilitzat realment per a res. Aquest canal no s'està consumint.

D'acord, així que com l'utilitzem? Molt simple. Prendré això, ho copiaré, i eliminaré _params.input_ i posaré _"greeting_ch"_ aquí en el seu lloc. Així que passarem aquest canal com a entrada a sayHello.

Tingueu en compte que he codificat aquesta cadena de moment. Això és una mica un pas enrere després del nostre bonic paràmetre que vam utilitzar al final del capítol anterior, però només manté les coses simples per començar perquè pugueu veure la lògica.

D'acord, aniré al meu terminal i executaré el workflow de nou. Sense cap _"--input"_ aquesta vegada, i s'executarà i utilitzarà aquest canal que hem creat i esperem que hauríem de tenir un fitxer aquí dalt a _results/hello_channels/_ i ara diu "Hello Channels!". Fantàstic. Així que això és el que esperàvem del nostre canal aquí. Genial.

## 1.4. Utilitzar view() per inspeccionar els continguts del canal

Una cosa més a afegir aquí, només una introducció ràpida a una altra funció que podem utilitzar en canals anomenada "_.view"_.

Això és anàleg a la comanda _print_ en Python o altres llenguatges als quals podríeu estar acostumats, i només bolca els continguts d'aquest canal al terminal quan l'executem.

Així que feu "_.view"_, i després si torno a executar el workflow, hauria d'imprimir al terminal quins són els continguts d'aquest canal, en el moment en què el vam crear.

Efectivament, podeu veure que s'ha imprès al terminal aquí. _"Hello Channels!"_.

Tingueu en compte que podeu trencar aquestes coses a través de línies si voleu, i de fet, el formatador automàtic de Nextflow intentarà fer-ho per vosaltres. L'espai en blanc no és realment important aquí, així que podeu encadenar aquestes coses una darrere l'altra.

## 2. Modificar el workflow per executar-se amb múltiples valors d'entrada

D'acord, així que el nostre canal té una cosa que és bonica, però és bàsicament el mateix que abans. Així que fem-ho una mica més complicat. Afegim algunes coses més al nostre canal.

La Channel factory "_.of()"_ pot prendre múltiples elements, així que escrivim alguns més. Farem _Hello, Bonjour, Hej_. I després podem executar aquest workflow de nou i veurem què passa.

Hauria d'executar-se de nou. I hem imprès ara. _"Hello", "Bonjour"_ i _"Hej"_ al terminal amb la nostra declaració view. Fantàstic.

## 2.1.2. Executar la comanda i mirar la sortida del registre

Podríeu pensar que hem acabat en aquest punt. Però en realitat hi ha una mica de trampa aquí, que ens farà ensopegar. Si mirem el nostre fitxer de sortida aquí. Podeu veure que té _"Hello"_ dins, però no té cap de les altres sortides. De fet, només és aquesta.

Si executem aquest workflow múltiples vegades, fins i tot podríem veure que de vegades té _"Bonjour"_, de vegades té _"Hej"_. És una mica aleatori.

Si mirem el terminal, podem veure que s'ha executat tres vegades i podem veure les diferents sortides de view. Però si vaig al directori work, puc fer _"cat work"_. Posar aquest hash i expandir-lo i _output.txt_. Podeu veure que aquest fitxer al directori work és diferent del directori results, i aquest és _"Hej"._ Així que hi ha alguna cosa que no funciona del tot bé aquí.

I la clau és que tenim tres tasques que s'han executat. La sortida de Nextflow intenta resumir això a mesura que avança el processament, de manera que no ocupi completament tot el vostre terminal, i aquest registre ANSI utilitza codis d'escapament ANSI, bàsicament ha sobreescrit les altres tasques. Així que només us mostra l'última que va passar a ser actualitzada.

## 2.1.3. Executar la comanda de nou amb l'opció -ansi-log false

Hi ha algunes coses que podem fer per entendre això una mica millor. Podem mirar al directori work mateix i podeu veure tots els diferents directoris work allà, però això és una mica confús perquè estarà barrejat amb diferents execucions de Nextflow.

O podem dir a Nextflow que no utilitzi els codis d'escapament ANSI.

Així que si executo la comanda de nou, però aquesta vegada dic _"-ansi-log false"_ per desactivar-lo, també podria utilitzar les variables d'entorn _$NO_COLOR_ o _"$NXF_ANSI_LOG=false"_. Llavors utilitza el tipus de registre de Nextflow més antic sense cap d'aquests codis d'escapament. Només imprimeix directament a un terminal sense actualitzacions intel·ligents.

I ara podem veure tots tres d'aquests processos que s'han executat. I cadascun d'ells el seu propi hash de tasca. I si anem a aquests directoris work, veurem les tres salutacions diferents que vam especificar.

Així que això té una mica més de sentit ara. Esperem que entengueu que Nextflow estava fent això, només estava sent una mica intel·ligent amb el que us mostrava al terminal amb aquests directoris work.

No obstant això, això està solucionat per a un problema amb els directoris work, però no ha solucionat un problema amb el fitxer de sortida. Encara només tenim un fitxer de sortida que diu _"Hello"_.

## 2.2. Assegurar que els noms dels fitxers de sortida seran únics

Ara per entendre això, hem de tornar al nostre script de workflow. Estem generant el nostre canal aquí, el passem al nostre procés, i si mirem el procés, estem escrivint la salutació a un fitxer anomenat _"output.txt"_ i passant aquest fitxer de sortida de tornada al bloc output aquí baix, publicant-lo.

No obstant això, cada tres vegades que aquest procés s'executa aquestes tres tasques diferents. Totes generen un fitxer anomenat _"output.txt"_, tots aquests fitxers de sortida es publiquen al directori results, i tots se sobreescriuen mútuament. Així que qualsevol fitxer de resultats que obtingueu allà és només l'últim que es va generar, però va esborrar tots els altres. Això no és realment el que volem.

## 2.2.1. Construir un nom de fitxer de sortida dinàmic

Hi ha diferents maneres de gestionar això, però la més simple per ara és només crear diferents noms de fitxer únics. Així que cada vegada que la tasca s'executi amb una salutació diferent, generarà un fitxer de sortida diferent, que ja no xocarà quan es publiqui. I després obtindrem tres fitxers de sortida únics.

Ho fem exactament de la mateixa manera. Podem utilitzar aquesta variable en qualsevol lloc dins del bloc script i podem utilitzar-la múltiples vegades.

Així que puc enganxar-la aquí, _"$\{greeting\}\_output.txt"_, i després també necessito enganxar-la aquí dalt perquè ja no estem creant un fitxer anomenat _output.txt_. Així que si no actualitzo això, Nextflow fallarà amb un error dient que esperava un fitxer, que mai es va generar.

Així que necessito fer el mateix allà i necessito utilitzar cometes dobles, no cometes simples, perquè aquesta variable s'entengui.

D'acord, provem-ho i vegem si ha funcionat. Executarem el workflow de nou. Esperem que ens mostri les tres tasques diferents dins dels tres directoris work diferents. I efectivament, podeu veure a la carpeta results aquí dalt a l'esquerra. Ara tenim tres fitxers diferents amb tres noms de fitxer diferents i cadascun amb els continguts diferents que esperem. Així que els fitxers ja no se sobreescriuen mútuament, i tot està allà com esperem.

Aquesta és una mica una configuració trivial per la qual hem passat aquí, però subratlla alguns dels conceptes clau que necessiteu entendre sobre com funciona la publicació de fitxers, i algunes de les coses en què podríeu caure com a trampes. Així que esperem que pugueu evitar això als vostres propis workflows.

També val la pena assenyalar que el que hem fet aquí és una mica poc pràctic en situacions de la vida real. Hem pres algunes dades d'entrada i estem utilitzant aquestes dades, però també estem anomenant el fitxer després d'aquestes dades, cosa que normalment no podeu fer.

Així que en pipelines de Nextflow més madurs i reals, sovint passareu un objecte meta amb totes les metadades associades amb una mostra donada. Després podeu crear noms de fitxer dinàmics basats en això, que és molt més pràctic.

Si esteu interessats en com fer això amb les millors pràctiques, hi ha una missió secundària a _training.nextflow.io_, que tracta específicament sobre metadades i mapes meta, així que podeu aprofundir allà per a més detalls.

## 3. Proporcionar múltiples entrades mitjançant un array

D'acord. A continuació explorarem una mica sobre com s'estructuren els canals i com difereixen d'altres tipus d'estructures de dades en el llenguatge de codificació. I pensaré una mica sobre com podria potencialment utilitzar un array, que podria ser un concepte familiar si veniu d'altres llenguatges.

Puc utilitzar un array en un canal? Provem-ho. Crearé un array, i he copiat això de la documentació, _"greetings_array"_ i _"Hello", "Bonjour"_ i _"Holà"_. I després posaré això aquí en lloc de les meves cadenes codificades. Així que diré "channel.of" _"greetings_array"_, passant aquest array a un canal. Provem-ho.

Obrir el terminal, i executar el pipeline.

D'acord. Podeu veure que la declaració view aquí va imprimir el nostre array com s'esperava, però després tot aquest text vermell, o no serà vermell si encara teniu _"-ansi-log"_ desactivat, però tot aquest text vermell ens està dient que alguna cosa va malament.

Ja no tenim una marca de verificació verda aquí. Tenim una creu vermella, i si només faig això una mica més ample perquè sigui més fàcil de llegir, Nextflow ens està dient què va malament.

Així que desglossem això secció per secció. Diu que l'error va ser causat per, i després la raó de l'error, que són fitxers de sortida que falten. Així que bàsicament aquest bloc output va dir que aquest fitxer hauria de ser creat i no ho va ser. A continuació diu que aquesta és la comanda que es va executar. Així que això és bàsicament els continguts d'aquest fitxer _.command.sh_. Així és com es veia després que totes aquestes variables s'haguessin posat.

I podeu veure aquí que la nostra comanda echo en realitat només s'ha executat una vegada i ha utilitzat tot l'array, però en una representació de cadena, que no és realment el que volíem.

I després la comanda va sortir així, i aquest era el directori work on podem anar i veure els fitxers per entendre una mica més.

D'acord. Així que el que va passar llavors va ser. Nextflow només va passar tot aquest array com un únic element de canal al procés, el que va significar que el procés només s'ha executat una vegada. Tenia una tasca i no va utilitzar les dades en una estructura que esperàvem.

## 3.2. Utilitzar un operador per transformar els continguts del canal

Així que necessitem fer alguna cosa a aquest canal primer, abans que pugui ser utilitzat. I això està preparant l'escenari per utilitzar operadors, que són funcions especials que podem utilitzar en canals per manipular els continguts del canal.

En aquest cas, utilitzarem quelcom anomenat _flatten_. Que passem al final del canal aquí. Així que creem el canal i després executem _flatten_. I de nou, si hi passem per sobre, ens mostra la documentació per a aquesta comanda directament a VS Code, que és molt útil. També podeu trobar tota aquesta documentació al lloc web de Nextflow, la documentació.

Podria simplement executar aquest codi ara i veure si funciona, però també és una bona oportunitat per introduir com fer codi dinàmic dins d'operadors i dins del codi de Nextflow, que s'anomenen closures.

Així que afegiré de nou una comanda view aquí abans d'executar _flatten_. I aquí aquesta té aquestes claus ondulades, que és la closure dinàmica. I només hi ha algun codi arbitrari dins d'aquí que s'executarà, dins del context d'un operador view.

Aquí, això està dient pren la salutació, que és l'entrada de l'operador view, i això és aquí. Podria anomenar això com volgués, podria anomenar això _"foo"_ i només necessito referir-m'hi com a _"foo"_ més tard. I després dic amb això, retorna això.

I després estableix retornant una cadena que diu abans del flatten per a una variable. Molt simple.

Ara afegiré un altre d'aquests exactament igual, però diré després de _flatten_.

Així que el que això fa, perquè això s'executa en seqüència, veureu com es veu el canal abans d'executar _flatten_, i després de nou després d'executar _flatten_.

I després aquest canal greeting encara es crea, així que encara es passarà al procés. I esperem que ara el workflow s'executi. Provem-ho.

Genial. Així que primer de tot és que el pipeline no va fallar aquesta vegada. Teníem tres processos que s'han executat correctament i tenim una petita marca de verificació. I després podem veure que les nostres declaracions view van funcionar.

Tenim abans de _flatten_, que és aquest array que vam veure abans de la fallada, i després tenim tres vegades el després de _flatten_ es va cridar on tenim _"Hello", "Bonjour",_ i tots aquests altres tres elements separats a l'array, que ara són com esperàvem, tres elements separats al canal.

I podeu veure que l'operador _view_ es va executar tres vegades. I això és perquè aquest canal després de _flatten_ ara té tres elements. I així l'operador es crida tres vegades.

Molt ràpidament, només mencionaria que quan estava creant Channel factories abans, vaig fer _"."_, i després vam veure que hi havia moltes maneres diferents de crear canals, i una d'elles s'anomena "_fromList"_. I això està realment dissenyat específicament per fer aquesta mateixa operació. Així que podríem haver fet simplement from list greetings away, i això funcionarà. És una sintaxi lleugerament neta i més agradable. Però per als propòsits d'aquesta demostració, volíem fer-ho una mica més pas a pas perquè poguéssiu veure com s'està manipulant el canal i com diferents operadors poden canviar el que hi ha al contingut d'un canal.

## 4. Llegir valors d'entrada des d'un fitxer CSV

D'acord, com podem fer això una mica més realista? Probablement no voldreu estar creant molts codis al vostre pipeline de Nextflow amb arrays codificats. Probablement voldreu prendre les dades de fora quan llanceu, i aquestes dades gairebé certament estaran en fitxers.

Així que la següent cosa que farem és que replicarem això, però en lloc de prendre les dades d'un únic paràmetre CLI o d'una cadena o array codificat, les prendrem d'un fitxer.

Així que desfer-nos del nostre greetings away. I ara canviarem aquesta Channel factory de nou. Acabo de dir que hi havia un munt per triar i n'hi ha una anomenada _".fromPath"_. I li diré que, en aquest cas, prengui _params.input_, que està tornant a la nostra entrada que estàvem utilitzant abans.

Ara aquest paràmetre no està realment preparat per ser utilitzat encara. Encara estem dient que és una cadena i està codificat aquí amb un valor per defecte, però podríem sobreescriure aquesta cadena. Ara volem que això sigui un fitxer en el seu lloc. Així que el tipus és diferent. Ja no és un _String_. És un _Path_.

I després podem establir el valor per defecte si volem, de nou, a un Path. I si miro a explorar a l'esquerra, podeu veure en aquest repositori, en aquest directori de treball, tinc un directori anomenat data. Tinc un fitxer allà anomenat _"greetings.csv"._

Així que puc simplement establir el valor per defecte aquí a _"data/greetings.csv"_. Ara, quan executi aquest pipeline de nou sense cap opció de línia de comandes, utilitzarà aquest valor per defecte. Sap que és un path, així que sap que hauria de gestionar-lo com un path i no una cadena.

I després passarà això a una Channel factory des d'aquest _params.input_ i crearà el nostre canal, que després s'utilitzarà en aquest procés anomenat _sayHello_. Provem-ho.

D'acord. Ha fallat. No us preocupeu. Això era esperat. I si esteu seguint el material de formació, veureu que era esperat allà també. Vegem què està passant aquí.

Ha intentat executar el pipeline. Ha intentat executar el procés, i ha obtingut un error força similar al que vam veure abans.

Aquí diu: hem intentat executar _echo_, però en lloc de fer echo dels continguts d'aquest fitxer CSV, només va fer echo del path. I podeu veure que és el path absolut complet aquí a aquest fitxer CSV.

I després efectivament, perquè va intentar escriure això a aquest path realment complicat, no sabia realment què fer. I estava fora de l'abast del directori work del procés.

Vaig esmentar al principi que Nextflow encapsula cada tasca executada dins d'un directori work especial. I si intenteu escriure a dades, que estan fora d'aquest directori work, Nextflow us aturarà com a precaució de seguretat. I això és el que ha passat aquí. Intentem escriure a un path absolut i Nextflow va fallar i ens va impedir.

## 4.2. Utilitzar l'operador splitCsv() per analitzar el fitxer

D'acord, donem una ullada a aquest canal i vegem com es veu. Podem fer _".view",_ i he copiat això del lloc web. Així que _.view_, i tenim una closure dinàmica aquí i diem un nom de variable "_csv"_ com a entrada. Així que això són els continguts del canal, i diem abans de splitCsv, i així és com es veu.

Si l'executo de nou, encara fallarà, però ens mostrarà què hi ha dins d'aquest canal. No és particularment emocionant. És aquesta variable _path_. Així que podeu veure que només és una cadena aquí perquè s'està imprimint a un terminal, però és un objecte _path_, que conté la informació i les metadades sobre aquest fitxer.

No volem passar les metadades del fitxer a l'entrada. Volem passar els continguts d'aquest fitxer. Si mirem el fitxer _greetings.csv_, podeu veure aquí que té aquestes variables diferents aquí. _Hello, Bonjour, Holà_ de nou. I aquestes són les coses que realment volem passar al nostre procés, no només el fitxer mateix com un únic objecte.

Així que necessitem analitzar aquest fitxer CSV. Necessitem desempaquetar-lo, arribar als continguts del fitxer CSV, i després passar els continguts dins del canal al procés.

Com probablement podeu dir pel missatge de registre, volem utilitzar el _splitCsv_, que és un altre operador, un altre operador de canal. Així que si faig "_dot" "s"_, i després podeu veure que està auto-suggerit. Ups, _splitCsv_ i uns parèntesis.

I després després de _splitCsv_, posaré una altra declaració _view_ només perquè puguem veure com es veu després. Executem el pipeline i vegem què tenim.

D'acord. Encara ha fallat, però d'una manera nova i emocionant, que és progrés.

Aquesta vegada de nou, tenim algun problema amb el nostre script, que s'ha renderitzat. Ara. Ja no tenim el path final, però tenim un array de variables, que sembla molt l'error que teníem abans quan estàvem passant un array com a entrada fixa.

Amb el nostre registre de l'operador view, podem veure abans que _splitCsv_ era el path. I efectivament, després de _splitCsv_, tenim tres sortides diferents i cadascuna d'aquestes sortides sembla molt cada una de les files del fitxer _greetings.csv_, que té sentit.

Així que el que ha passat aquí és que Nextflow ha analitzat aquest fitxer CSV donant-nos tres objectes, un array per a cada línia del fitxer CSV. Així que després tres vegades hem passat un array de variables al canal en lloc d'un únic valor de cadena.

D'acord, així que l'última vegada que vam tenir aquest problema, vam utilitzar _flatten_. Només molt ràpidament. Provem flatten i vegem què passa.

Puc anomenar aquestes variables, el que sigui. Així que l'anomenaré _myarray_ perquè ja no és realment un CSV. Provem d'executar-lo de nou i vegem què passa amb _flatten_.

Així que aquesta vegada executarem, vam analitzar el CSV en tres objectes array, i després el vam aplanar. I aquesta vegada va passar. I el pipeline de Nextflow s'ha executat. No obstant això podeu veure que _flatten_ realment va a fons i aplana tot. I així obtenim tres entrades d'array independents per a cada fila. I així va executar el procés tres vegades cada fila d'un CSV. I ara tenim un munt de fitxers de resultats, i 123, 456, i tot tipus de coses, no només aquesta primera columna del CSV, que és el que realment volíem.

## 4.3. Utilitzar l'operador map() per extreure les salutacions

Així que com arribem només a la primera columna? Si flatten és massa simplista aquí, necessitem un operador més complex on puguem realment personalitzar i dir-li què volem del CSV.

Per fer això, utilitzarem _map_. Bàsicament _map_ només diu, executa algun codi, alguna funció sobre cada element que em donin i fes algun tipus de transformació sobre ell. I perquè és tan flexible, el veureu aparèixer en el codi de Nextflow tot el temps.

Per si mateix, no fa res. Així que no volem parèntesis regulars, volem una closure aquí i necessitem dir-li què fer. Així que diré _"row"_, perquè se li estan donant files del CSV, així que és un nom de variable lògic. És l'entrada. I vull retornar només el primer element d'aquest array.

Els arrays a Nextflow es basen en zero, així que direm només el primer element, que és la fila zero. Si volguéssim la segona columna, podria ser un o la tercera columna ser dos, i així successivament. Podem retornar el que vulguem aquí, però retornaré només el primer valor.

I ara, podem executar el pipeline de nou i veure si fa el que esperem.

Efectivament, després de _splitCsv_ tenim els nostres arrays, i després després del _map,_ tenim les nostres cadenes netes i agradables, només _"Hello", "Bonjour"_ i _"Holà"_. I el pipeline ara està fent el que volem. Fantàstic.

Així que podem desfer-nos de totes aquestes comandes view ara. Ja no les necessitem.

## Resum

Hem acabat el nostre tipus de depuració i aquest és el codi amb què acabem. Prenent el nostre paràmetre CLI anomenat _input_, que està classificat com un _Path_. Nextflow troba el path, el carrega, i entén el fitxer CSV. Retorna totes les diferents files. I després mapegem només el primer element d'aquesta fila al canal que tipus de dóna els continguts del canal, que es passa al procés.

I el procés s'executa sobre cada element del canal, que són tres. I executa el procés tres vegades, donant-li tres tasques. I aquests resultats es publiquen després des del workflow, recollits per la sortida del procés. Publicats des d'un workflow i desats al bloc output a un subdirectori anomenat _"hello_channels"_.

Molt guai. Ara estem arribant a alguna cosa que s'assembla més a un pipeline de Nextflow de la vida real que podríeu executar per a alguna anàlisi real.

## Conclusió

D'acord. Esperem que ara estigueu aconseguint una sensació del que són els canals i operadors de Nextflow i com els operadors treballen sobre els canals i com podeu crear-los.

Els canals, com vaig dir al principi d'aquest vídeo, són l'adhesiu de Nextflow. I podeu veure aquí que podem prendre diferents entrades i manipular-les i prendre aquestes dades i després passar-les a la lògica de workflow posterior.

I aquest bloc workflow aquí és realment on construïu tota aquesta paral·lelització i tota la lògica intel·ligent, i expliqueu a Nextflow com construir el vostre DAG de workflow, i com orquestrar el vostre pipeline.

Els canals no són el concepte més fàcil d'entendre. Així que feu una pausa, penseu una mica sobre això, potser llegiu el material de nou, i assegureu-vos realment que teniu aquests conceptes clars perquè això és clau per a la vostra comprensió de Nextflow i com millor entengueu els canals i els diferents operadors de canal i les diferents Channel factories. Més diversió tindreu escrivint Nextflow i més poderosos seran els vostres pipelines.

Això no és el mateix que la programació regular en Python o altres llenguatges. No estem utilitzant declaracions _if_ aquí, això és programació de flux funcional utilitzant canals i operadors. Així que és una mica diferent, però també és súper poderós.

Aquest és el final d'aquest capítol. Aneu i feu una pausa ràpida i us veuré al següent vídeo per a la part tres on passarem per Hello Workflow, i parlarem una mica més sobre els workflows.

Igual que el capítol anterior, hi ha algunes preguntes de qüestionari a la part inferior de la pàgina web aquí, així que podeu fer una ràpida passada per aquestes i assegurar-vos que enteneu totes les diferents parts del material que acabem de fer. I a part d'això, us veuré al següent vídeo. Moltes gràcies.

D'acord.
