SHELL := /bin/bash

install:
	curl -s "https://get.sdkman.io" | bash
	source ${HOME}/.sdkman/bin/sdkman-init.sh && sdk install asciidoctorj
	curl  https://repo1.maven.org/maven2/com/puravida-software/asciidoctor/asciidoctor-quizzes/0.1.0/asciidoctor-quizzes-0.1.0.jar \
		-o ${HOME}/.sdkman/candidates/asciidoctorj/current/lib/asciidoctor-quizzes-0.1.0.jar


training:
	mkdir -p build
	cd asciidocs && asciidoctorj -cp ${HOME}/.sdkman/candidates/asciidoctorj/current/lib/asciidoctor-quizzes-0.1.0.jar -D ../build index.adoc && cp -r img ../build/img
	

handson:
	mkdir -p build/hands-on
	pushd handson-docs
	asciidoctorj -D ../build/hands-on index.adoc
	cp -r img ../build/hands-on/img
	popd
