# Yerel Devcontainer'lar

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Yerel bir Docker kurulumunuz varsa veya kurmaya hazırsanız, bu materyallerle yerel olarak çalışmanın en kolay yolu Visual Studio Code'un devcontainer özelliğini kullanmaktır. Bu yaklaşım, manuel kurulum gerektirmeden gerekli tüm araçları ve bağımlılıkları sağlar.

## Gereksinimler

Yerel devcontainer kurulumunu kullanmak için şunlara ihtiyacınız olacak:

- [Visual Studio Code](https://code.visualstudio.com/)
- Yerel bir Docker kurulumu, örneğin:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (Windows/macOS için)
  - [Docker Engine](https://docs.docker.com/engine/install/) (Linux için)
  - [Colima](https://github.com/abiosoft/colima) (macOS için alternatif)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (Docker Desktop'a dahildir, ancak diğer Docker kurulumlarında ayrı kurulum gerekebilir)
- VS Code için [Dev Containers eklentisi](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

Devcontainer'ı açmaya çalışmadan önce Docker kurulumunuzun çalışıyor olması gerekir.

Docker buildx'in kullanılabilir olduğunu doğrulamak için şunu çalıştırın:

```bash
docker buildx version
```

Bu komut başarısız olursa, devam etmeden önce buildx eklentisini kurmanız gerekecektir.

## Kurulum Talimatları

VS Code devcontainer'larını kullanarak yerel ortamınızı kurmak için bu adımları izleyin:

### VS Code'da "Dev Containers" eklentisini kurun

- VS Code'u açın
- Extensions'a gidin (Ctrl+Shift+X veya macOS'ta Cmd+Shift+X)
- "Dev Containers" arayın
- "Install"a tıklayın

![Installing Dev Containers extension in VS Code](img/install_extension.png)

### Depoyu klonlayın:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Depoyu VS Code'da açın:

- VS Code'u başlatın
- Menüden **File -> Open Folder** seçin
- Az önce klonladığınız training deposu klasörüne gidin ve seçin
- **Open**'a tıklayın

### Konteyner'da Yeniden Aç

VS Code tarafından "Reopen in Container" sorulursa, tıklayın. Alternatif olarak:

- F1'e basın (veya Ctrl+Shift+P / macOS'ta Cmd+Shift+P)
- "Dev Containers: Reopen in Container" yazın
- **Önemli**: Yapılandırma seçmeniz istendiğinde, **local-dev** devcontainer yapılandırmasını seçin

![Reopen in Container prompt](img/reopen_prompt.png)

![Selecting local configuration](img/select_local_config.png)

Konteyner'ın oluşturulmasını bekleyin. İlk seferde tüm gerekli bileşenleri indirip kurduğu için birkaç dakika sürebilir.

Konteyner oluşturulup çalıştırıldığında, gerekli tüm araçların yüklü olduğu tam yapılandırılmış bir ortama sahip olacaksınız:

- Java
- Nextflow
- Docker
- Git
- Ve eğitim için gereken diğer tüm bağımlılıklar

![VS Code with devcontainer running](img/running_container.png)

## Devcontainer'ları Kullanmanın Avantajları

Devcontainer yaklaşımını kullanmak birçok avantaj sunar:

- **Tutarlılık**: Farklı makinelerde tutarlı bir geliştirme ortamı sağlar
- **Basitlik**: Tüm bağımlılıklar önceden yüklenmiş ve yapılandırılmıştır
- **İzolasyon**: Geliştirme ortamı yerel sisteminizden izole edilmiştir
- **Tekrarlanabilirlik**: Devcontainer'ı kullanan herkes aynı kurulumu alır
- **Manuel kurulum yok**: Java, Nextflow ve diğer araçları manuel olarak kurmanıza gerek yok

## Ortamınızı Kontrol Etme

Devcontainer'ınız çalıştıktan sonra, her şeyin doğru şekilde kurulduğunu doğrulamak için şunu çalıştırabilirsiniz:

```bash
nextflow info
```

Bu, ortamınızın düzgün yapılandırıldığını doğrulayan Nextflow sürümünü ve çalışma zamanı bilgilerini görüntülemelidir.

## Sorun Giderme

Devcontainer kurulumunda sorunlarla karşılaşırsanız:

1. Devcontainer'ı açmadan önce Docker kurulumunuzun (Docker Desktop, Colima, Docker Engine vb.) çalıştığından emin olun
2. İstendiğinde **local-dev** yapılandırmasını seçtiğinizi doğrulayın
3. `docker buildx version` çalıştırarak Docker buildx'in kurulu ve çalışır durumda olduğunu doğrulayın
4. Konteyner oluşturma başarısız olursa, "Dev Containers: Rebuild Container" komutunu çalıştırarak yeniden oluşturmayı deneyin
5. Kalıcı sorunlar için [VS Code Dev Containers sorun giderme rehberine](https://code.visualstudio.com/docs/devcontainers/troubleshooting) bakın
