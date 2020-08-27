clean:
	rm -rf asciidocs/index.html
	rm -rf handson-docs/index.html
 
docs: asciidocs/index.html handson-docs/index.html
asciidocs/index.html: asciidocs/*.adoc
	cd asciidocs; asciidoctor index.adoc
	mkdir -p ../seqera-website/public/training
	cp -r asciidocs/img ../seqera-website/public/training
	cp asciidocs/index.html ../seqera-website/public/training

handson-docs/index.html: handson-docs/*.adoc
	cd handson-docs; asciidoctor index.adoc
	cd handson-docs/solutions; asciidoctor *.adoc
	mkdir -p ../seqera-website/public/training/handson
	mkdir -p ../seqera-website/public/training/handson/solutions
	cp -r handson-docs/img ../seqera-website/public/training/handson/
	cp handson-docs/*.{html,css} ../seqera-website/public/training/handson
	cp handson-docs/solutions/*.{html,css} ../seqera-website/public/training/handson/solutions

upload:
	rm -rf nf-training/{.nextflow*,work}
	aws s3 sync --acl public-read --delete --include "nf-training/{data,0?_*}" nf-training s3://seqeralabs.com/public/nf-training
