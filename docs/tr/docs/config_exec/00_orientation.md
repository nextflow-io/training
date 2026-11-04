# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim Ortamını Başlatın

GitHub Codespaces üzerinde sunduğumuz hazır ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Ortam yüklenirken okumaya devam edebilmek için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (kullandığınız cihaza göre sağ tıklama, ctrl+tıklama veya cmd+tıklama kullanabilirsiniz).
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam Temelleri

Bu eğitim ortamı, kursu tamamlamak için gereken tüm yazılım, kod ve verileri içermektedir; dolayısıyla herhangi bir şey kurmanıza gerek yoktur.

Codespace, VSCode arayüzüyle yapılandırılmıştır. Bu arayüz; bir dosya sistemi gezgini, bir kod düzenleyici ve bir terminal kabuğu içermektedir.
Kurs boyunca verilen tüm talimatlar (örneğin "dosyayı açın", "kodu düzenleyin" veya "bu komutu çalıştırın"), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunmaktadır.

Bu kursu kendi başınıza tamamlıyorsanız, daha fazla ayrıntı için [ortam temelleri](../envsetup/01_setup.md) sayfasını incelemenizi öneririz.

### Sürüm Gereksinimleri

Bu kurs, v2 sözdizimi ayrıştırıcısının etkin olduğu (25.10+ sürümünde varsayılan) Nextflow 25.10.2 veya daha yeni bir sürümünü gerektirmektedir.
Yerel veya özel bir ortam kullanıyorsanız, [burada](../info/nxf_versions.md) belgelenen doğru ayarları kullandığınızdan emin olun.

## Çalışmaya Hazırlanın

Codespace'iniz çalışmaya başladıktan sonra, başlamadan önce yapmanız gereken iki şey vardır: çalışma dizininizi ayarlamak ve sağlanan materyallere göz atmak.

### Çalışma Dizinini Ayarlayın

Codespace varsayılan olarak tüm eğitim kurslarının kök dizininde açılır.
Bu kurs için `config-exec/` dizinine geçin:

```bash
cd config-exec/
```

Ardından VSCode'u bu dizine odaklanacak şekilde ayarlayın; böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak tam yolu kullanarak her zaman geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/config-exec
    ```

### Sağlanan Materyalleri İnceleyin

Kurs materyallerini sol taraftaki dosya gezginini kullanarak veya `tree` komutuyla inceleyebilirsiniz.
Tam yapıyı görmek için terminalden aşağıdaki komutu çalıştırın:

```bash
tree . -L 2
```

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

**`main.nf`** ve **`modules/`** dosyaları, [Nextflow Run](../nextflow_run/index.md) kursundaki çok adımlı pipeline'ın aynısıdır; **`nextflow.config`** dosyası ise orada gördüğünüz yapılandırmanın aynısıdır.
Bu alıştırmalar boyunca her ikisini de genişleteceksiniz.

**`data/`** dizini, pipeline'ın okuduğu CSV girdi dosyasını içermektedir.

## Hazırlık Kontrol Listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım

Tüm kutuları işaretleyebildiyseniz, başlayabilirsiniz.

**[Bölüm 1: Hesaplama ortamınıza uyarlama](./01_packaging_and_execution.md) sayfasına geçmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
