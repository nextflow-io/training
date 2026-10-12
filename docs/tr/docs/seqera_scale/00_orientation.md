# Başlarken

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Eğitim Ortamını Başlatın

GitHub Codespaces üzerinde sunduğumuz hazır ortamı kullanmak için aşağıdaki "Open in GitHub Codespaces" düğmesine tıklayın. Diğer seçenekler için [Ortam seçenekleri](../envsetup/index.md) sayfasına bakın.

Ortam yüklenirken okumaya devam edebilmek için eğitim ortamını yeni bir tarayıcı sekmesinde veya penceresinde açmanızı öneririz (kullandığınız cihaza göre sağ tıklama, ctrl+tıklama veya cmd+tıklama kullanabilirsiniz).
Kursu tamamlamak için bu talimatları paralel olarak açık tutmanız gerekecektir.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Ortam Temelleri

Bu eğitim ortamı, kursu tamamlamak için gerekli tüm yazılım, kod ve verileri içermektedir; bu nedenle kendiniz herhangi bir şey yüklemenize gerek yoktur.

Codespace, bir dosya sistemi gezgini, kod düzenleyici ve terminal kabuğu içeren VSCode arayüzüyle yapılandırılmıştır.
Kurs boyunca verilen tüm talimatlar (örneğin "dosyayı açın", "kodu düzenleyin" veya "bu komutu çalıştırın"), aksi belirtilmedikçe VSCode arayüzünün bu üç bölümüne atıfta bulunmaktadır.

Bu kursu kendi başınıza tamamlıyorsanız, daha fazla ayrıntı için lütfen [ortam temelleri](../envsetup/01_setup.md) sayfasını inceleyin.

## Çalışmaya Hazırlanın

Codespace'iniz çalışmaya başladıktan sonra, başlamadan önce yapmanız gereken iki şey vardır: çalışma dizininizi ayarlamak ve sağlanan materyallere göz atmak.

### Çalışma Dizinini Ayarlayın

Codespace varsayılan olarak tüm eğitim kurslarının kök dizininde açılır.
Bu kurs için `seqera-scale/` dizinine geçin:

```bash
cd seqera-scale/
```

Ardından VSCode'u bu dizine odaklanacak şekilde ayarlayın; böylece dosya gezgini kenar çubuğunda yalnızca ilgili dosyalar görünür:

```bash
code .
```

!!! tip "İpucu"

    Herhangi bir nedenle bu dizinden çıkarsanız (örneğin codespace'iniz uyku moduna geçerse), Github Codespaces eğitim ortamında çalıştığınızı varsayarak her zaman tam yolu kullanarak geri dönebilirsiniz:

    ```bash
    cd /workspaces/training/seqera-scale
    ```

### Sağlanan Materyalleri İnceleyin

Kurs materyallerini soldaki dosya gezginini kullanarak veya `tree` komutuyla inceleyebilirsiniz.
Tam yapıyı görmek için terminalden aşağıdaki komutu çalıştırın:

```bash
tree -a .
```

??? abstract "Dizin içeriği"

    ```console
    .
    └── .seqera_config
    ```

**`.seqera_config`** dosyası, `tw` CLI'ını Seqera erişim tokenınız ve çalışma alanınızla yapılandırmak için 3. bölümde dolduracağınız bir taslaktır.

## Hazırlık Kontrol Listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Ortamım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım

Tüm kutuları işaretleyebildiyseniz, başlamaya hazırsınız.

**[Bölüm 1: Pipeline'ları web arayüzünden başlatma](./01_run_with_seqera.md) sayfasına geçmek için bu sayfanın sağ alt köşesindeki oka tıklayın.**
