# Bölüm 3: İş akışı çalıştırmalarını yönetme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Pipeline'ları çalıştırdıkça ve yeniden çalıştırdıkça, çalıştırma geçmişi ve eski `work/` dizinleri birikmeye başlar.
[Bölüm 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work)'de zaten tamamlanmış işleri atlamak için `-resume` seçeneğini kullandınız.
Burada bir çalıştırma hakkında rapor oluşturmayı, [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) ile geçmiş çalıştırmaların geçmişini incelemeyi ve artık ihtiyaç duymadığınız eski work dizinlerini [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) ile silmeyi öğreneceksiniz.

---

## 1. Pipeline raporları oluşturma

Nextflow, bir çalıştırma hakkında birkaç farklı türde rapor oluşturabilir; her biri kendi `-with-*` bayrağıyla eklenir: bir çalıştırma raporu (`-with-report`), bir çalıştırma zaman çizelgesi (`-with-timeline`), bir görev izleme dosyası (`-with-trace`) ve bir iş akışı diyagramı (`-with-dag`).
Burada ilk ikisini oluşturacağız; diğerleri için Nextflow referansındaki [Çalıştırma raporları](https://nextflow.io/docs/latest/reports.html) bölümüne bakın.

### 1.1. Çalıştırma raporu oluşturma

Pipeline tamamlandıktan sonra HTML raporu oluşturmak için herhangi bir `nextflow run` komutuna `-with-report` ekleyin:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow, raporu çalışma dizininde `report-<timestamp>.html` adlı bir dosyaya yazar.
Bir tarayıcıda açarak çalıştırma özetini, her görevin durumu ve çalışma süresiyle birlikte yer aldığı bir tabloyu ve süreç bazında ayrıştırılmış kaynak kullanım grafiklerini görebilirsiniz.

**Tasks** sekmesi, pipeline'ın çalıştırdığı her görevi süreç adı, durumu ve kaynak kullanımıyla birlikte listeler:

![Çalıştırma raporu görev tablosu](img/execution_report_tasks.png)

Rapor, özellikle bir pipeline beklenenden uzun sürdüğünde veya bir görev başarısız olduğunda çok işe yarar: görev tablosu, zamanın tam olarak nerede harcandığını ve hangi görevlerin başarılı ya da başarısız olduğunu gösterir.

### 1.2. Çalıştırma zaman çizelgesi oluşturma

Her görevin ne zaman çalıştığını Gantt grafiği tarzında görüntülemek için bir çalıştırmaya `-with-timeline` ekleyin:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow, zaman çizelgesini `timeline-<timestamp>.html` adlı bir dosyaya yazar.
Bir tarayıcıda açarak her görev için, ne zaman çalıştığına ve ne kadar sürdüğüne göre konumlandırılmış ve boyutlandırılmış çubukları görebilirsiniz:

![Çalıştırma zaman çizelgesi](img/execution_timeline.png)

Zaman çizelgesi, [Bölüm 1](./01_run_nextflow.md#31-run-the-workflow)'deki yelpaze açılıp kapanma şeklini bir bakışta görünür kılar: üç `sayHello` görevi paralel olarak çalışır, ardından üç `convertToUpper` görevi, sonra `collectGreetings` ve `cowpy` birbiri ardına çalışır; çünkü her biri kendinden önceki her şeye bağımlıdır.

### Özetle

`-with-report` ile HTML çalıştırma raporu ve `-with-timeline` ile çalıştırma zaman çizelgesi oluşturmayı, ayrıca Nextflow'un desteklediği diğer rapor türlerini nerede bulacağınızı öğrendiniz.

### Sırada ne var?

Geçmiş çalıştırmaların geçmişini nasıl inceleyeceğinizi öğrenin.

---

## 2. Geçmiş çalıştırmaların günlüğünü inceleme

İster bir pipeline geliştiriyor olun ister onu üretim ortamında çalıştırıyor olun, bir noktada geçmiş çalıştırmalar hakkında bilgi aramanız gerekecektir.

### 2.1. Geçmiş dosyası

Bir Nextflow iş akışını her başlattığınızda, geçerli çalışma dizinindeki `.nextflow` adlı gizli bir dizin altında bulunan `history` adlı bir günlük dosyasına bir satır yazılır.

??? abstract "Dosya içeriği"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Her satır, bu dizinden başlatılan bir çalıştırmanın zaman damgasını, süresini, çalıştırma adını, durumunu, revizyon kimliğini, oturum kimliğini ve tam komut satırını gösterir.

Son iki satıra bakın: bunlar aynı komutun iki ayrı çağrısıdır (biri düz, biri `-resume` ile) ve aynı oturum kimliğini paylaşırlar.
Oturum kimliği yalnızca gerçekten yeni bir çalıştırma başlattığınızda değişir; `-resume` kullanmak onu korur ve Nextflow hangi cache'i yeniden kullanacağını bu sayede bilir.

### 2.2. Daha kullanışlı bir görünüm için `nextflow log` kullanma

Ham geçmiş dosyasını okumak işe yarar, ancak `nextflow log` aynı bilgileri bir başlıkla birlikte biçimlendirir:

```bash
nextflow log
```

??? success "Komut çıktısı"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow, `-resume` için kullandığı önbellekleme bilgilerini oturum kimliğine göre anahtarlanmış şekilde `.nextflow/cache` altında gruplar.
Bu nedenle, geçmiş bir çalıştırmayı araştırmanız veya temizlemeniz gerektiğinde doğru çalıştırma adını ya da oturum kimliğini buradan bulmak her zaman ilk adımdır.

### Özetle

Nextflow'un geçmiş çalıştırmaların kaydını nerede tuttuğunu ve bunu `nextflow log` ile nasıl inceleyeceğinizi öğrendiniz.

### Sırada ne var?

Artık ihtiyaç duymadığınız eski work dizinlerini nasıl kaldıracağınızı öğrenin.

---

## 3. Eski work dizinlerini silme

Her çalıştırma, çıktılarını `results/` dizinine kopyaladıktan sonra bile görev dizinlerini `work/` altında bırakır.
Geliştirme sürecinde yeterince pipeline çalıştırıldığında bu alt dizinler birikmeye başlar; bu nedenle Nextflow, artık ihtiyaç duyulmayan dizinleri kaldırmak için `nextflow clean` komutunu sağlar.

### 3.1. Silme ölçütlerini belirleme

`nextflow clean`, neyin kaldırılacağını seçmek için birkaç farklı yöntemi destekler; tam liste için [referans belgelerine](https://www.nextflow.io/docs/latest/reference/cli.html#clean) bakın.
Burada, çalıştırma adını kullanarak belirli bir çalıştırmadan önceki her şeyi sileceksiniz.

`nextflow log` kullanarak saklamak istediğiniz en son çalıştırmayı bulun; [2.2'deki örnekte](#22-use-nextflow-log-for-a-friendlier-view) bu, `-resume` çalıştırmasından önceki son düz çalıştırma olan `elegant_panini`'dir.
Çalıştırma adı, konsolda `Launching (...)` satırında veya `nextflow log` çıktısının `RUN NAME` sütununda gösterilen, makine tarafından oluşturulmuş iki parçalı dizedir.

### 3.2. Deneme çalıştırması yapma

Gerçekten silmeden önce belirli bir komutun ne sileceğini kontrol etmek için önce `-n` ekleyin:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Komut çıktısı"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

Bu 16 görev dizinidir: `turkey` çalıştırmasından 8 görev ve `tux` çalıştırmasından 8 görev; bu dört süreçli pipeline'ın iki tam çalıştırması için tam olarak beklenen sayıdır.
`elegant_panini` çalıştırmasının kendisi ve `-resume` çalıştırmasının ondan yeniden kullandığı önbelleğe alınmış görevler olduğu gibi bırakılır.

Çıktınız farklı dizin adları listeleyecektir; kaç satır gördüğünüz ise kaç çalıştırma yaptığınıza bağlıdır. Hiç satır görmüyorsanız ya çalıştırma adı günlüğünüzdeki bir adla eşleşmiyordur ya da ondan önce silinecek bir şey yoktur.

### 3.3. Silme işlemini gerçekleştirme

Deneme çalıştırması doğru göründüğünde, aynı komutu `-n` yerine `-f` ile yeniden çalıştırın:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Komut çıktısı"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean`, görev dizinlerini boşaltır ancak iki karakterli üst dizinleri (örneğin `e5/`) yerinde bırakır.

!!! warning "Uyarı"

    Geçmiş çalıştırmalara ait work dizinlerini silmek, onları Nextflow'un cache'inden kaldırır ve yalnızca orada depolanan çıktıları siler.
    Bu durum, Nextflow'un ilgili süreçleri yeniden çalıştırmadan çalıştırmaya devam etme yeteneğini ortadan kaldırır; bu nedenle yalnızca devam ettirmeye ihtiyaç duymayacağınızdan emin olduğunuz çalıştırmaları temizleyin.
    Bu aynı zamanda, `work/` dizinine veya `symlink` yayımlama moduna güvenmek yerine, önem verdiğiniz her şeyi `mode 'copy'` ile `results/` dizinine yayımlamanın neden değerli olduğunu da açıklar.

### Özetle

`nextflow clean` ile eski work dizinlerini nasıl kaldıracağınızı ve bunu yapmanın söz konusu çalıştırmalardan devam etme yeteneğini nasıl ortadan kaldırdığını öğrendiniz.

### Sırada ne var?

[Bölüm 4](./04_remote_repositories.md)'te GitHub gibi uzak depolardan doğrudan pipeline çalıştırmayı öğrenin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- `-with-report` ile HTML çalıştırma raporu ve `-with-timeline` ile çalıştırma zaman çizelgesi oluşturma
- `nextflow log` ile geçmiş çalıştırmaların geçmişini inceleme
- `nextflow clean` ile eski work dizinlerini kaldırma ve buna eşlik eden devam ettirme ödünleşimini anlama
