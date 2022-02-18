PWD := $(shell pwd)

clean:
	rm -rf asciidocs/index.html
	rm -rf handson-docs/index.html
 
docs: asciidocs/index.html handson-docs/index.html
asciidocs/index.html: asciidocs/*.adoc
	cd asciidocs; asciidoctor index.adoc
	mkdir -p ../seqera-website-ver2/static/training
	cp -r asciidocs/img ../seqera-website-ver2/static/training
	cp asciidocs/index.html ../seqera-website-ver2/static/training

handson-docs/index.html: handson-docs/*.adoc
	cd handson-docs; asciidoctor index.adoc
	cd handson-docs/solutions; asciidoctor *.adoc
	mkdir -p ../seqera-website-ver2/static/training/handson
	mkdir -p ../seqera-website-ver2/static/training/handson/solutions
	cp -r handson-docs/img ../seqera-website-ver2/static/training/handson/
	cp handson-docs/*.{html,css} ../seqera-website-ver2/static/training/handson
	cp handson-docs/solutions/*.{html,css} ../seqera-website-ver2/static/training/handson/solutions

docker-clean:
	docker run --rm -it -v ${PWD}/asciidocs/:/documents/ asciidoctor/docker-asciidoctor rm -rf ./build
docker:
	docker run --rm -it -v ${PWD}/asciidocs/:/documents/ asciidoctor/docker-asciidoctor asciidoctor -D build index.adoc
	docker run --rm -it -v ${PWD}/asciidocs/:/documents/ asciidoctor/docker-asciidoctor cp -r ./img ./build