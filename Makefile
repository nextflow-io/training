SHELL := /bin/bash
CWD := $(shell cd -P -- '$(shell dirname -- "$0")' && pwd -P)

make-antora:
	docker build -t custom-antora docker/

clean-nf-training:
	docker run --rm -it -v $(CWD):/antora alpine /bin/ash -c "rm -rf ./build/nf-training"

nf-training: clean-nf-training
	docker run --rm -it -v $(CWD):/antora --env DOCSEARCH_ENABLED=true --env DOCSEARCH_ENGINE=lunr --env FORCE_SHOW_EDIT_PAGE_LINK=true custom-antora nf-training-playbook.yml --stacktrace

clean: clean-nf-training

build: nf-training

all: make-antora clean nf-training

netlify:
	yarn global add asciidoctor-plantuml asciidoctor-kroki @antora/cli@3.0.0 @antora/site-generator@3.0.0 @antora/lunr-extension
	DOCSEARCH_ENABLED=true DOCSEARCH_ENGINE=lunr FORCE_SHOW_EDIT_PAGE_LINK=true  antora nf-training-playbook.yml --stacktrace