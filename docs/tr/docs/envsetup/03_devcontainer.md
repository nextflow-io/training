# Yerel Devcontainer'lar

Yerel bir Docker kurulumunuz varsa veya kurulum yapmaktan çekinmiyorsanız, bu materyallerle yerel olarak çalışmanın en kolay yolu Visual Studio Code'un devcontainer özelliğini kullanmaktır. Bu yaklaşım, manuel kurulum gerektirmeden gerekli tüm araçları ve bağımlılıkları sağlar.

## Gereksinimler

Yerel devcontainer kurulumunu kullanmak için şunlara ihtiyacınız olacak:

- [Visual Studio Code](https://code.visualstudio.com/)
- Yerel bir Docker kurulumu, örneğin:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (Windows/macOS için)
  - [Docker Engine](https://docs.docker.com/engine/install/) (Linux için)
  - [Colima](https://github.com/abiosoft/colima) (macOS için alternatif)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (Docker Desktop'ta dahildir, ancak diğer Docker kurulumlarında ayrı kurulum gerekebilir)
- VS Code için [Dev Containers eklentisi](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

Devcontainer'ı açmayı denemeden önce Docker kurulumunuzun çalışıyor olması gerekir.

Docker buildx'in mevcut olduğunu doğrulamak için şunu çalıştırın:

```bash
docker buildx version
```

Bu komut başarısız olursa, devam etmeden önce buildx eklentisini kurmanız gerekecektir.

## Kurulum Talimatları

VS Code devcontainer'larını kullanarak yerel ortamınızı kurmak için şu adımları izleyin:

### VS Code'da "Dev Containers" eklentisini kurun

- VS Code'u açın
- Eklentiler'e gidin (Ctrl+Shift+X veya macOS'ta Cmd+Shift+X)
- "Dev Containers" araması yapın
- "Install" (Kur) düğmesine tıklayın

![VS Code'da Dev Containers eklentisini kurma](img/install_extension.png)

### Depoyu klonlayın:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Depoyu VS Code'da açın:

- VS Code'u başlatın
- Menüden **File -> Open Folder** (Dosya -> Klasör Aç) seçeneğini seçin
- Az önce klonladığınız training deposu klasörüne gidin ve seçin
- **Open** (Aç) düğmesine tıklayın

### Container'da Yeniden Açın

VS Code tarafından "Reopen in Container" (Container'da Yeniden Aç) istemi görürseniz, üzerine tıklayın. Alternatif olarak:

- F1'e basın (veya Ctrl+Shift+P / macOS'ta Cmd+Shift+P)
- "Dev Containers: Reopen in Container" yazın
- **Önemli**: Bir yapılandırma seçmeniz istendiğinde, **local-dev** devcontainer yapılandırmasını seçin

![Container'da Yeniden Aç istemi](img/reopen_prompt.png)

![Yerel yapılandırma seçimi](img/select_local_config.png)

Container'ın oluşturulmasını bekleyin. İlk seferinde gerekli tüm bileşenleri indirip kurduğu için bu birkaç dakika sürebilir.

Container oluşturulup çalışmaya başladığında, aşağıdakiler dahil olmak üzere gerekli tüm araçların kurulu olduğu tam yapılandırılmış bir ortama sahip olacaksınız:

- Java
- Nextflow
- Docker
- Git
- Ve eğitim için gerekli diğer tüm bağımlılıklar

![Devcontainer çalışırken VS Code](img/running_container.png)

## Devcontainer Kullanmanın Avantajları

Devcontainer yaklaşımını kullanmak çeşitli avantajlar sunar:

- **Tutarlılık**: Farklı makinelerde tutarlı bir geliştirme ortamı sağlar
- **Basitlik**: Tüm bağımlılıklar önceden kurulu ve yapılandırılmıştır
- **İzolasyon**: Geliştirme ortamı yerel sisteminizden izole edilmiştir
- **Tekrarlanabilirlik**: Devcontainer'ı kullanan herkes aynı kurulumu elde eder
- **Manuel kurulum gerektirmez**: Java, Nextflow ve diğer araçları manuel olarak kurmanıza gerek yoktur

## Ortamınızı Kontrol Etme

Devcontainer'ınız çalışmaya başladığında, her şeyin doğru şekilde kurulduğunu şunu çalıştırarak doğrulayabilirsiniz:

```bash
nextflow info
```

Bu, Nextflow sürümünü ve çalışma zamanı bilgilerini görüntülemeli ve ortamınızın düzgün şekilde yapılandırıldığını doğrulamalıdır.

## Sorun Giderme

Devcontainer kurulumunda sorunlarla karşılaşırsanız:

1. Devcontainer'ı açmadan önce Docker kurulumunuzun (Docker Desktop, Colima, Docker Engine vb.) çalıştığından emin olun
2. İstendiğinde **local-dev** yapılandırmasını seçtiğinizi kontrol edin
3. `docker buildx version` komutunu çalıştırarak Docker buildx'in kurulu ve çalıştığını doğrulayın
4. Container oluşturulamazsa, "Dev Containers: Rebuild Container" komutunu çalıştırarak yeniden oluşturmayı deneyin
5. Kalıcı sorunlar için [VS Code Dev Containers sorun giderme kılavuzuna](https://code.visualstudio.com/docs/devcontainers/troubleshooting) başvurun
