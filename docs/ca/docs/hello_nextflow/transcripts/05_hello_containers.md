# Part 5: Hello Containers - Transcripció del vídeo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Notes importants"

    Aquesta pàgina mostra només la transcripció. Per a instruccions completes pas a pas, torneu al [material del curs](../05_hello_containers.md).

    Els números de secció mostrats a la transcripció es proporcionen només amb finalitats indicatives i poden no incloure tots els números de secció dels materials.

## Benvinguda i context

Hola, i benvinguts de nou a Hello Nextflow. Aquesta és la part cinc anomenada Hello Containers. I en aquesta part del curs, parlarem sobre com encapsular els requisits de programari per a un pipeline perquè les persones que executin el pipeline no hagin de pensar en instal·lar el programari.

Si heu estat treballant en bioinformàtica tant de temps com jo, potser recordeu el que sovint anomeno els mals temps antics, on quan volíeu executar el pipeline d'algú altre o replicar el seu treball, passàveu hores o dies intentant instal·lar totes les diferents eines de programari que utilitzaven, a les mateixes versions, intentant compilar-les a la vostra màquina, i era un malson. Era realment difícil.

Si estàveu executant en un HPC, potser heu utilitzat mòduls d'entorn on els administradors de sistemes intentaven instal·lar programari per a vosaltres, cosa que estava bé, però encara imperfecta.

Però ara tenim millors maneres de fer-ho. Nextflow té suport integrat per a diferents tecnologies de contenidors de programari. Docker és la més comuna. És la que utilitzarem avui. Funciona bé a Codespaces. Funciona bé al vostre ordinador local i funciona bé al núvol.

Però també Singularity o Apptainer, que són molt comuns en sistemes HPC i efectivament funcionen exactament de la mateixa manera. O Podman, Shifter, hi ha un munt d'altres que són tots molt similars.

L'únic extra, que és una mica similar però no del tot, que Nextflow suporta és Conda. I Nextflow pot gestionar entorns Conda per a vosaltres per procés, cosa que és molt millor que fer els vostres propis entorns Conda. I de nou, pot venir amb un pipeline.

Començarem aquest capítol parlant una mica sobre tecnologies de contenidors i Docker i com funcionen. I farem la primera meitat només manualment a Docker perquè entengueu què està passant sota el capó i com funciona això. Perquè això és realment important per entendre què està fent Nextflow i com entendre què està fent el vostre workflow quan s'està executant.

Així que. Saltem al nostre Codespaces. Ara he netejat tot de nou, però si saltem a Hello Containers, hauríeu de veure que tots els nostres scripts i tot són allà igual que al final del capítol de mòduls. Així que tenim els nostres diferents mòduls aquí, que vaig crear al directori de mòduls.

Encara hi són. Han d'estar-hi perquè pugui executar-se. i el workflow i la sortida són tots iguals excepte que hem canviat el camí de publicació de sortida a Hello Containers, perquè els vostres fitxers acabin en aquell directori.

Podem executar això ara per comprovar que funciona si voleu, o podem continuar amb el terminal.

## 1. Utilitzar un contenidor 'manualment'

Utilitzarem Docker per gestionar els nostres contenidors, i puc comprovar que està instal·lat al meu Codespaces fent "docker -v", que em mostra la versió que està instal·lada i tot, i que està funcionant correctament.

Ara els contenidors i Docker tenen dos conceptes que són realment importants. Un s'anomena imatge, i un s'anomena contenidor. La imatge és la instantània, si voleu, de tot el sistema de fitxers que utilitzareu, i el contenidor és l'entorn en execució. Així que creeu un contenidor utilitzant una imatge.

Un cop esteu dins d'aquell contenidor, normalment funciona com un sistema operatiu complet. Està tallat del món exterior. Està separat de tot la resta, i això és bo. Així és com obtenim una reproducibilitat tan bona amb Nextflow.

Perquè per a tasques executades dins d'un contenidor, no estan contaminades per cap fitxer de configuració al vostre sistema local. Cap altra influència externa, s'executen al seu propi petit entorn aïllat. Els fitxers es produeixen llavors d'una manera molt, molt reproducible perquè esteu utilitzant les mateixes biblioteques subjacents, totes les mateixes dependències, exactament el mateix programari per a cada persona executant en cada entorn informàtic diferent. El que francament crec que és fantàstic i increïble que funcioni. I fins i tot, fins i tot avui dia encara em sorprèn que això sigui possible.

## 1.1. Descarregar la imatge del contenidor

Així que provarem d'utilitzar algunes imatges Docker i Docker, quan l'executeu al vostre sistema, té un registre docker al vostre ordinador, o en aquest cas, al code space, que fa un seguiment de totes les diferents imatges que s'han descarregat i utilitzat en el passat, i les diferents capes de les quals estan construïdes.

Podem veure quines imatges tenim localment amb Docker fent "docker image ls". I en aquest cas podeu veure que hi ha un munt d'imatges Docker aquí, que totes tenen a veure amb la configuració d'aquest Codespaces. Tot a veure amb contenidors de desenvolupament i coses. Així que no us heu de preocupar massa per elles, però a mesura que afegim més imatges i les descarreguem, a mesura que aquest curs avança, podeu comprovar aquesta llista i veureu que el registre local està fent un seguiment de totes aquestes coses que hem descarregat.

Però agafarem una de nova fent "docker pull". I això li diu a Docker que obtingui una nova imatge del web.

Llavors posem l'URI per a aquell contenidor. Ara això podria ser una imatge docker que heu construït localment i després heu pujat a internet. podria ser una imatge que algú altre ha fet. Hi ha moltes, moltes, moltes maneres diferents de fer imatges Docker, però possiblement una de les maneres més simples és externalitzar això, i aconseguir que algú altre ho faci per a vosaltres.

I el que utilitzarem en aquest tutorial és un servei de Seqera anomenat Seqera Containers.

Ara, Seqera Containers és totalment gratuït, i utilitza una peça de programari de codi obert que desenvolupem anomenada Wave, que es va construir per gestionar contenidors d'una manera complementària a Nextflow. I gestiona molts dels casos d'ús comuns amb els quals ens trobem tractant amb Nextflow.

És molt comú que el programari que necessitem estigui empaquetat a Conda, a canals Bioconda, o conda-forge o altres canals més específics de domini. I Wave i Seqera Containers és realment bo construint imatges a partir d'això.

Així que puc anar a aquesta interfície web i jugarem amb el paquet anomenat "cowpy". Així que escric el nom del paquet que vull. Cerca, l'ha trobat a l'índex de paquets Python, així que puc utilitzar-lo. O si espero una mica més, està cercant bioconda i conda-forge. I podeu veure, puc especificar qualsevol canal conda aquí. Així que si voleu trobar un canal Nvidia o qualsevol altra cosa, això també hauria de funcionar.

I llavors puc especificar si vull que construeixi una imatge docker per a mi o una imatge singularity i també quina arquitectura de CPU vull utilitzar. Així que amd64 o arm64.

I un cop es llisten els resultats de bioconda, ara puc veure totes les diferents versions que estan disponibles també. Ho posaré. I ara podria continuar cercant i obtenir més paquets de Conda si vull i compondre aquest contenidor com vulgui, però només vull aquest. Així que faré clic a Get Container.

Ara, algú altre ja ha sol·licitat el mateix contenidor abans i es retorna des d'un registre, així que simplement l'obtenim immediatament. Però si ningú més mai hagués demanat aquest paquet de programari o aquesta combinació de paquets de programari, Wave i Seqera Containers el construirien sobre la marxa per a nosaltres.

Podem copiar aquesta URL i també podem veure els detalls de construcció. I això ens mostra el que el servei va fer al backend. Va crear un fitxer d'entorn conda. Un fitxer docker, i llavors això és, executant el procés de construcció docker. També va executar un escaneig, un escaneig de seguretat, així que podeu veure qualsevol CVE. I us diu quan es va crear això.

Wave i Seqera Containers poden fer molt més que això, però aquest és un cas d'ús simple, que és el més comú. I hauria de dir que aquestes imatges s'allotgen durant almenys cinc anys. Així que podeu construir aquestes URL als vostres pipelines i saber que no desapareixeran aviat.

Així que tinc la meva URL per a la meva imatge docker per a cowpy.

Ara puc fer "docker pull" d'aquella URL, i obtindrà totes les diferents capes i descarregarà aquesta imatge perquè estigui disponible per a mi localment.

## 1.2. Utilitzar el contenidor per executar cowpy com una comanda única

D'acord, ara provem d'utilitzar-lo realment. Així que ara utilitzaré una comanda "docker run" en lloc de "docker pull", i utilitzaré la bandera "--rm", que simplement li diu a Docker que tanqui aquest contenidor un cop hagi acabat de fer el que li he demanat. I llavors poso l'identificador per al contenidor, que és només un URI.

I llavors al final, especifico la comanda que vull que Docker executi dins del contenidor generat a partir d'aquesta imatge. Només diré cowpy, que és el nom de l'eina que està instal·lada des de Conda Forge, que està disponible dins de la imatge.

Premeré enter i aquí ho teniu. Hem executat cowpy en un sistema. Tenim una petita vaca donant-nos alguna informació.

Ara noteu que cowpy no està instal·lat al meu sistema local. Així que si l'executo sense totes les coses de Docker, diu, comanda no trobada. Així que això ha descarregat una imatge. Ha creat un contenidor utilitzant Docker, i llavors ha entrat en aquell contenidor i ha executat aquesta comanda per a nosaltres i ens ha donat la sortida de tornada al nostre terminal. Molt, molt genial.

## 1.3. Utilitzar el contenidor per executar cowpy interactivament

D'acord, anirem un pas més enllà ara i executarem aquest contenidor interactivament i furgarem una mica, perquè puguem veure què està passant dins del contenidor.

Així que si torno enrere i agafo la meva comanda run i eliminaré cowpy al final allà, perquè realment no vull executar cowpy. Vull executar un terminal Bash.

I llavors tornaré aquí i faré "-it", que significa Interactiu i Terminal o TTY, i premeré enter.

I ara podeu veure el prompt, la part abans que escrigui, ha canviat. Aquest era el prompt de Codespaces on deia el directori, i ara diu base i roots i tmp. Així que ara estic dins del contenidor, i si faig "ls", veureu que els fitxers que veig en aquest directori són diferents dels fitxers que tinc al meu espai de treball.

I de fet, no puc veure cap dels fitxers del meu espai de treball de codespaces local o del meu disc local dins del contenidor Docker. El temps d'execució del contenidor docker, està completament aïllat i no pot escriure ni llegir cap fitxer d'un sistema de fitxers amfitrió fora.

Puc, però, veure el programari que està instal·lat dins del contenidor i executar-lo. Així que puc executar cowpy i podem veure una mica més sobre com utilitzar cowpy. Aquí puc fer "cowpy 'Hello World'" i això li diu, li diu que posi realment la meva cita dins d'una petita bombolla de parla. I també podeu executar diferents tipus de vaques, així que no ha de ser una vaca. Podeu fer un "-c". I estic a Suècia, així que triaré un ant. Molt bonic. Li he donat unes banyes.

I hi ha un munt de diferents que podeu provar, que podeu veure descrits als documents de formació.

## 1.3.4. Muntar dades al contenidor

D'acord. Seria bonic si poguéssim executar cowpy als fitxers del nostre sistema de fitxers.

Per descomptat, no és súper útil només tenir el contenidor i cap accés a res en absolut. Pot ser segur i reproducible, però no és gaire útil.

Així que com ho fem? Sortiré d'aquest contenidor Docker escrivint exit, i podeu veure que el prompt ens diu que ara estem de tornada al nostre Codespaces regular.

I executaré la mateixa comanda de nou. Però aquesta vegada afegiré algunes banderes addicionals aquí enrere. I la important és "-v", que significa muntar un volum, que és com bàsicament una part, part d'un espai de disc.

El "-v" pren dues parts: hi ha com una cadena i llavors dos punts i una cadena. I la primera part és el sistema de fitxers local, que hauria de ser muntat al contenidor. I llavors la segona part és on això hauria d'acabar dins del contenidor.

Ara només vull carregar tot el meu sistema de fitxers local aquí. Així que "." és el directori de treball actual. Així que només faré "." i llavors ":", i llavors posarem això en un nou directori dins del contenidor anomenat "my_project". Això realment podria anomenar-se qualsevol cosa.

I llavors executaré de nou.

Al directori de treball on estic deixat, que és /tmp, els fitxers no hi són. Però si faig "ls my_project", aquí ho tenim: tots els mateixos fitxers que teníem localment al nostre Codespaces ara estan disponibles dins del contenidor en aquell camí.

Aquest és accés de lectura i escriptura així que puc crear nous fitxers en aquest directori i apareixeran al meu sistema de fitxers amfitrió. Aquest directori particular, llavors es comporta exactament com si estigués fora del contenidor així que ara puc llegir i escriure i fer coses.

## 1.3.5. Utilitzar les dades muntades

D'acord, només provem que podem fer això. Faig "cat /my_project/data/greetings.csv". Si recordeu aquest contingut de fitxer es veu així. Ara puc canalitzar això a cowpy i la vaca imprimirà les diferents sortides d'aquell fitxer a la seva petita bombolla de parla, cosa que és una mica divertida.

Així que podeu veure, ara podem utilitzar el programari al contenidor per interactuar amb els fitxers al nostre sistema amfitrió.

D'acord, tornem a sortir i continuarem amb la resta del material de formació.

## 2. Utilitzar contenidors a Nextflow

Així que això és realment genial utilitzant contenidors. Esperem que tingui sentit. I podeu veure el valor d'aquests contenidors i per què això és útil per executar programari d'anàlisi.

Però com fem tot aquest mateix procés dins de Nextflow? No volem estar executant un munt de comandes Docker nosaltres mateixos. Volem només deixar que Nextflow gestioni tot això per a nosaltres.

Així que treballem en això. Afegirem un nou procés al nostre pipeline, per executar cowpy. D'acord, així que creem un nou mòdul per al nostre nou procés. Així que anem a mòduls, anomenem-lo cowPy.nf, i llavors copiaré el codi del material de formació aquí.

Però podeu veure que el procés és molt simple. Es veu molt com els que hem fet fins ara, tenim un bloc d'entrada amb un camí, que és el nostre fitxer d'entrada, i també un valor aquí perquè això serà un caràcter, així que podríem utilitzar un ant de nou si volem.

I llavors una sortida, que és un únic fitxer aquí, un camí i llavors un script. I estem fent el mateix que vam fer interactivament dins del contenidor: estem fent "cat" per llegir el fitxer d'entrada. Estem canalitzant aquell contingut a cowpy. Estem triant un caràcter específic basat en aquella entrada, estem escrivint a un fitxer de sortida anomenat cowpy input file, que llavors s'ecoitza a la sortida.

Genial. Incloem això. Així que include \{ cowpy \} from "./modules/cowpy.nf", l'he anomenat cowpy? Sí.

I llavors cridem el nostre nou procés aquí baix al bloc principal del workflow. Així que executem cowpy. I agafarem el nostre nou procés cowpy i direm collectGreetings.out.

I llavors si recordeu, hi havia dues sortides per a aquest mòdul. Una anomenada outfile i una anomenada report. L'extensió VS Code està auto-suggerint aquestes per a nosaltres i volem el .outfile.

Sempre podeu saltar a aquest procés aquí. O bé passeu el cursor per sobre i hauria de mostrar-vos ràpidament quines eren les sortides. I també podem fer command click en ell i obrirà el fitxer del mòdul si voleu veure amb més detall.

Així que aquí anem. Aquest és l'outfile allà, i aquest és el camí. Així que ara això serà el fitxer d'entrada per al nostre procés cowpy. Fantàstic.

Ara si recordeu, un procés cowpy té dues entrades. També teníem el canal de valor per al caràcter. Així que podem afegir "params.character" aquí. Podria haver codificat això si hagués volgut, però fem-ho una opció CLI perquè puguem fer dash, dash character.

Bé. Ara necessito definir el paràmetre d'entrada que acabem de cridar i donar-li un valor per defecte. Així que character, String. I m'agrada l'ant, així que el definiré a moose per defecte.

Bé, provem d'executar-ho. Així que si faig Nextflow run hello containers, veurem què passa.

Podria haver utilitzat dash resume si hagués tingut els vells directoris de treball per allà. I de nou, aquests primers processos haurien estat, emmagatzemats en memòria cau i hauria estat una mica més ràpid, però hauria de ser bàsicament el mateix.

Ara podem veure directament que ha llançat un error quan va arribar al nostre nou procés, ens està dient aquí que hi va haver un error executant el procés cowpy i va sortir amb un estat de sortida 127. Aquesta és la comanda que va intentar executar. Es veu bé, es veu com esperàvem. Està agafant aquell nom de fitxer de sortida, que es veu correcte, l'està executant amb un caràcter moose i intentant desar-lo.

Però podeu veure que l'error de comanda aquí està dient que la comanda cowpy no s'ha trobat. I això té sentit perquè encara no hem dit realment a Nextflow que utilitzi un contenidor. Només li hem donat la comanda cowpy. I com he dit abans, cowpy no està instal·lat al nostre sistema local. Així que quan va intentar executar-lo, va fallar.

## 2.3.1. Especificar un contenidor per a cowpy

Necessitem dir a Nextflow que hi ha un contenidor disponible i que pot utilitzar-lo. Així que com ho fem?

Si saltem al nostre mòdul aquí, afegirem una nova declaració a la part superior anomenada "container". I llavors definirem això a una cadena.

Ara, si recordeu, a Seqera Containers, puc copiar aquella URL i només la deixo caure entre cometes aquí.

Ara tornem i provem d'executar-ho de nou.

Deixeu-me veure si funciona aquesta vegada.

Malauradament, falla exactament de la mateixa manera, encara que ara hem definit un contenidor perquè el procés s'executi. Així que per utilitzar la nostra imatge docker, necessitem dir a Nextflow que habiliti l'ús de Docker quan executem el workflow.

I ho farem creant un nou fitxer de configuració. Així que diré touch nextflow.config.

Aquest és un nom de fitxer especial on si està al directori de treball mentre inicio el pipeline, es carregarà automàticament. Així que si vaig a aquest fitxer Nextflow dot config, podeu veure que realment ja existeix, cosa que havia oblidat. I tenim docker.enabled aquí ja, però està definit a false, que és el valor per defecte.

Així que si canvio això a equals True en lloc d'això, docker.enabled. I hi ha documents de referència per a tots aquests àmbits de configuració als documents de Nextflow. I també podeu veure que quan passo el cursor per sobre amb una extensió VS Code, extreu els documents específics d'això i em diu què significa i com definir-lo.

Així que ara l'hem definit a true, i si executo Nextflow de nou, Nextflow ara sabrà obtenir aquella imatge docker per a nosaltres si encara no la tenim localment, i llavors executar aquell procés amb aquell entorn de contenidor.

I així podem veure que s'ha executat amb èxit i tenim un petit tic al costat d'un cowpy. Fantàstic. Si pujo i miro al directori de resultats, el fitxer encara no hi és. I això és perquè encara necessitem, publicar aquest fitxer de sortida igual que tots els altres.

Així que anem al bloc published dins del workflow, diguem mycowpy equals cowpy.out.

I llavors aquí baix al bloc de sortida, mycowpy, claus ondulades path. Ups. Hello containers. Mode, copy.

Si executo de nou ara, hauria d'executar-se exactament de la mateixa manera. Podria haver utilitzat dash resume i oblido cada vegada. I llavors pujo i ara tenim un nou fitxer creat anomenat cowpy-COLLECTED, i allà està el meu ant dient BONJOUR, HELLO, HOLA Fantàstic.

Ara per descomptat també podria passar ara "--character". Quines són les diferents opcions? Crec que hi ha un Turkey? Així que puc utilitzar character Turkey. S'executarà exactament de la mateixa manera. He perdut una altra oportunitat d'utilitzar dash resume, i ara si carreguem el nostre fitxer i ara tenim un Turkey. Fantàstic.

## 2.3.4. Inspeccionar com Nextflow va llançar la tasca contenidoritzada

D'acord. Última petita cosa. Executem aquesta comanda de nou ràpidament, reprenem aquesta vegada, i donem una ullada ràpida al directori de treball per veure què és el que Nextflow està fent sota el capó per fer que tot això funcioni per a nosaltres.

Aquesta vegada és súper ràpid, anem a aquest directori de treball, cd work/. Ara si recordeu tenim un munt de fitxers punt aquí i el que ens interessa en aquest cas és el que vaig dir que gairebé mai necessitem mirar, anomenat .command.run.

Si faig code dot command run, l'obrirà a l'editor. I puc cercar en aquest fitxer i si faig scroll cap avall hauria de veure Docker run. I podeu veure que Nextflow està fent la comanda docker run per a nosaltres, quan Docker està habilitat a la configuració. Té un munt de diferents banderes i coses aquí, però podeu veure la bandera "-v" que vam utilitzar nosaltres mateixos quan estàvem executant. I podeu veure que està muntant el directori de l'espai de treball local al contenidor, perquè el contenidor pugui accedir als nostres fitxers d'entrada i desar les sortides. I llavors al final, també està executant .command.sh, que és l'script generat, que té la comanda cowpy dins.

I així podeu veure que Nextflow està agafant la lògica del workflow, que és el que realment ens importa, que és específic de la nostra anàlisi, i està fent totes les coses intel·ligents entre bastidors per fer que Docker funcioni al nostre sistema.

I ho està fent d'una manera realment portable perquè un usuari final del pipeline pugui canviar la tecnologia que està utilitzant: Docker, Singularity, Apptainer, Conda. Això realment no importa a la lògica del pipeline, però Nextflow gestionarà totes les necessitats d'infraestructura subjacents, perquè s'executi a qualsevol lloc.

I això és realment el superpoder de Nextflow. És reproducibilitat i portabilitat. I amb Nextflow podeu realment compartir el vostre workflow i altres persones poden executar-lo als seus sistemes i simplement funcionarà.

Això és una cosa realment, realment difícil de fer, i ara sabeu com fer-ho també amb els vostres workflows.

D'acord, això és tot per a aquest capítol. si baixeu al final d'un curs, trobareu un qüestionari de nou sobre alguns contenidors. Esperem que tot tingués sentit. És una manera realment genial de treballar amb anàlisi. I si sou nous als contenidors, espero haver-vos convençut que és el camí a seguir, i mai mirareu enrere.

Però amb això, feu una petita pausa potser, i us uniu a mi en un parell de minuts per repassar la part final sis del Hello Nextflow, que tracta tot sobre configuració.

Moltes gràcies.
