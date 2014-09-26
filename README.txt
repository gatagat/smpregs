smpregs
=======

build kmers.so:
	python setup_kmers.py build
	mv build/lib*/*.so ./
	rm -r build

encode annotations (file with genomic annotations for dm3 is included):
	./encode_annotations.py data/genomic-annotations-dm3.txt data/genomic-annotations-dm3.enc data/genomic-annotations-dm3.fa

get random regions of the same length:
	./smpregs.py -r data/S2-spec.bed > out

get random regions with approx. the same GC content (at most 5bps diffrent):
	./smpregs.py -r data/S2-spec.bed GC:threshold=5 > out

get random regions with the same genomic annotation at position 201:
	./smpregs.py -r data/S2-spec.bed GAPos:pos=201 > out

get random regions with approx. the same genomic annotation histogram:
	./smpregs.py -r data/S2-spec.bed GAHist:threshold=5 > out

combining multiple filters:
	./smpregs.py -r data/S2-spec.bed -n data/genomic-annotations-dm3.fa -g dm3 GAPos:pos=201,GC:threshold=5 > out
